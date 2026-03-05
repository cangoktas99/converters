/*----------------------------------------------------------------------------*/
/*  gmf2of.cpp                                                                */
/*  Convert a libMeshb binary mesh (.meshb) + optional solution (.solb) to   */
/*  an OpenFOAM polyMesh case directory with fields.                         */
/*                                                                            */
/*  Usage:                                                                    */
/*    gmf2of <input_base> <output_case_dir> [time_dir=0]                     */
/*           [-r <refCase> <refTime>] [-c <createPatchDict>]                 */
/*                                                                            */
/*  -r  copy field names, dimensions and BC types from a reference OF case   */
/*  -c  merge boundary refs into named patches using a createPatchDict file  */
/*----------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "libmeshb7.h"

/*----------------------------------------------------------------------------*/
/* Utilities                                                                  */
/*----------------------------------------------------------------------------*/

static void mkdirp(const std::string &path)
{
#if defined(_WIN32) || defined(_WIN64)
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

static std::string trim(const std::string &s)
{
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

/* FoamFile header writer */
static void writeFoamHeader(std::ofstream &out,
                            const std::string &cls,
                            const std::string &location,
                            const std::string &object,
                            int nCells = -1,
                            int nFaces = -1,
                            int nPoints = -1,
                            int nInternalFaces = -1)
{
    out << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        << "| =========                 |                                                 |\n"
        << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        << "|  \\\\    /   O peration     |                                                 |\n"
        << "|   \\\\  /    A nd           |                                                 |\n"
        << "|    \\\\/     M anipulation  |                                                 |\n"
        << "\\*---------------------------------------------------------------------------*/\n"
        << "FoamFile\n{\n"
        << "    version     2.0;\n"
        << "    format      ascii;\n"
        << "    class       " << cls << ";\n"
        << "    location    \"" << location << "\";\n"
        << "    object      " << object << ";\n";
    if (nCells >= 0) {
        out << "    note        \"nPoints:" << nPoints
            << "  nCells:"   << nCells
            << "  nFaces:"   << nFaces
            << "  nInternalFaces:" << nInternalFaces << "\";\n";
    }
    out << "}\n"
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

/*----------------------------------------------------------------------------*/
/* Face key for deduplication: sorted vertex list                             */
/*----------------------------------------------------------------------------*/

using FaceKey = std::vector<int>; /* sorted vertex indices */

static FaceKey makeFaceKey(const std::vector<int> &verts)
{
    FaceKey k = verts;
    std::sort(k.begin(), k.end());
    return k;
}

struct FaceRecord {
    std::vector<int> verts;  /* in owner-CCW order */
    int owner;
    int neighbour;           /* -1 = boundary */
    int patchRef;            /* 0 = internal, else patch ref from boundary elements */
};

/*----------------------------------------------------------------------------*/
/* Element face templates (0-based local node indices, CCW from owner)        */
/*----------------------------------------------------------------------------*/

static std::vector<std::vector<int>> tetFaceTemplates()
{
    return {
        {0,2,1},   /* opposite n3, outward */
        {0,1,3},   /* opposite n2, outward */
        {1,2,3},   /* opposite n0, outward */
        {0,3,2},   /* opposite n1, outward */
    };
}

static std::vector<std::vector<int>> hexFaceTemplates()
{
    return {
        {0,3,2,1}, /* bottom */
        {4,5,6,7}, /* top    */
        {0,1,5,4}, /* front  */
        {1,2,6,5}, /* right  */
        {2,3,7,6}, /* back   */
        {3,0,4,7}, /* left   */
    };
}

static std::vector<std::vector<int>> prismFaceTemplates()
{
    return {
        {0,2,1},   /* bottom */
        {3,4,5},   /* top    */
        {0,1,4,3}, /* side 0-1 */
        {1,2,5,4}, /* side 1-2 */
        {2,0,3,5}, /* side 2-0 */
    };
}

static std::vector<std::vector<int>> pyramidFaceTemplates()
{
    return {
        {0,3,2,1}, /* base */
        {0,1,4},   /* front tri  */
        {1,2,4},   /* right tri  */
        {2,3,4},   /* back  tri  */
        {3,0,4},   /* left  tri  */
    };
}

/*----------------------------------------------------------------------------*/
/* createPatchDict reader                                                     */
/*----------------------------------------------------------------------------*/

struct PatchDef {
    std::string name;
    std::string type;          /* "patch", "wall", "symmetry" */
    std::vector<int> refs;     /* e.g. {1,2,3,4,6} */
};

/* Parse an OpenFOAM createPatchDict file and extract
   patch name, type, and list of source Patch_N names → ref numbers. */
static std::vector<PatchDef> readCreatePatchDict(const std::string &path)
{
    std::vector<PatchDef> defs;
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "Warning: cannot open createPatchDict: %s\n",
                     path.c_str());
        return defs;
    }
    /* Simple state machine:
       Look for "name <patchName>;", "type <patchType>;",
       "patches ( Patch_N ... );" inside each { } block. */
    std::string line;
    int depth = 0;
    bool inPatches = false;
    PatchDef cur;
    bool haveCur = false;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty() || t[0] == '/' || t.rfind("//", 0) == 0) continue;

        /* Track brace depth inside the patches( ... ) list */
        for (char c : t) {
            if (c == '(') depth++;
            if (c == ')') depth--;
        }

        /* "name  patchName;" */
        if (t.find("name") == 0 && t.find("name ") != std::string::npos) {
            auto val = t.substr(t.find("name") + 4);
            val = trim(val);
            if (!val.empty() && val.back() == ';') val.pop_back();
            val = trim(val);
            cur.name = val;
            haveCur = true;
        }
        /* "type  wall;" inside patchInfo { } */
        if (haveCur && t.find("type") == 0 &&
            t.find("type ") != std::string::npos &&
            t.find("type") < 8) {
            auto val = t.substr(t.find("type") + 4);
            val = trim(val);
            if (!val.empty() && val.back() == ';') val.pop_back();
            val = trim(val);
            cur.type = val;
        }
        /* "patches ( Patch_1 Patch_2 ... );" — may span multiple lines */
        auto pp = t.find("patches");
        if (pp != std::string::npos && pp < 12) {
            inPatches = true;
        }
        if (inPatches) {
            /* Extract Patch_N tokens */
            std::istringstream ss(t);
            std::string tok;
            while (ss >> tok) {
                if (tok.rfind("Patch_", 0) == 0) {
                    /* Strip trailing ; or ) */
                    while (!tok.empty() && (tok.back()==';'||tok.back()==')'||tok.back()=='('))
                        tok.pop_back();
                    if (tok.rfind("Patch_", 0) == 0) {
                        int ref = std::atoi(tok.c_str() + 6);
                        if (ref > 0) cur.refs.push_back(ref);
                    }
                }
            }
            if (t.find(';') != std::string::npos ||
                (t.find(')') != std::string::npos && depth <= 1)) {
                inPatches = false;
            }
        }
        /* End of a patch entry: we see a closing } at depth 1 and have data */
        if (depth <= 1 && haveCur && !cur.name.empty() && !cur.refs.empty()) {
            if (cur.type.empty()) cur.type = "patch";
            defs.push_back(cur);
            cur = PatchDef();
            haveCur = false;
        }
    }
    return defs;
}

/*----------------------------------------------------------------------------*/
/* Reference case field metadata                                              */
/*----------------------------------------------------------------------------*/

struct RefFieldMeta {
    std::string name;
    bool        isVector;
    std::string dims;
    /* Per-patch BC content (patchName → lines inside { },
       with nonuniform Lists stripped, "value ..." lines stripped) */
    std::map<std::string, std::string> patchBCs;
};

static std::string readFoamClass(const std::string &path)
{
    std::ifstream in(path);
    if (!in) return "";
    std::string line;
    bool inHeader = false;
    while (std::getline(in, line)) {
        if (line.find("FoamFile") != std::string::npos) { inHeader = true; continue; }
        if (inHeader && line.find("class") != std::string::npos) {
            auto pos = line.find("class");
            auto val = line.substr(pos + 5);
            size_t start = val.find_first_not_of(" \t");
            if (start == std::string::npos) return "";
            val = val.substr(start);
            if (!val.empty() && val.back() == ';') val.pop_back();
            while (!val.empty() && (val.back() == ' ' || val.back() == '\t'))
                val.pop_back();
            return val;
        }
        if (line.rfind("// *", 0) == 0) break;
    }
    return "";
}

static std::string readDimensions(const std::string &path)
{
    std::ifstream in(path);
    if (!in) return "[0 0 0 0 0 0 0]";
    std::string line;
    while (std::getline(in, line)) {
        auto pos = line.find("dimensions");
        if (pos != std::string::npos) {
            auto bra = line.find('[', pos);
            auto ket = line.find(']', bra);
            if (bra != std::string::npos && ket != std::string::npos)
                return line.substr(bra, ket - bra + 1);
        }
    }
    return "[0 0 0 0 0 0 0]";
}

/* Read the boundaryField block from a field file and return per-patch
   BC entries.  Nonuniform List blocks and bare "value ..." lines
   (that OpenFOAM caches) are stripped — we'll replace them with
   interpolated data from the adapted solution.                               */
static std::map<std::string, std::string>
readPatchBCs(const std::string &path)
{
    std::map<std::string, std::string> result;
    std::ifstream in(path);
    if (!in) return result;
    std::string line;

    /* 1. Seek to "boundaryField" */
    bool found = false;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.rfind("boundaryField", 0) == 0) { found = true; break; }
    }
    if (!found) return result;
    /* Skip opening brace */
    while (std::getline(in, line))
        if (line.find('{') != std::string::npos) break;

    /* 2. Parse patch entries at brace depth 1 */
    int depth = 0;
    std::string curPatch;
    std::string curContent;
    bool skipping = false;
    int skipParenDepth = 0;
    bool skipSeenOpen = false;
    bool skipWaitSemicolon = false;

    while (std::getline(in, line)) {
        std::string t = trim(line);

        /* End of boundaryField */
        if (depth == 0 && t == "}") break;

        if (depth == 0 && !t.empty() && t != "{" && t != "}") {
            /* New patch name (may be quoted like ".*") */
            curPatch = t;
            /* Remove trailing semicolons if any */
            while (!curPatch.empty() && curPatch.back() == ';')
                curPatch.pop_back();
            curPatch = trim(curPatch);
            curContent.clear();
            skipping = false;
            skipWaitSemicolon = false;
            continue;
        }

        /* Track braces */
        for (char c : t) {
            if (c == '{') depth++;
            if (c == '}') depth--;
        }

        if (depth == 0 && t.find('}') != std::string::npos) {
            /* Closing brace of current patch */
            if (!curPatch.empty())
                result[curPatch] = curContent;
            curPatch.clear();
            curContent.clear();
            skipping = false;
            skipWaitSemicolon = false;
            continue;
        }

        if (curPatch.empty() || depth < 1) continue;

        /* --- Skip nonuniform List blocks --- */
        if (skipping) {
            for (char c : t) {
                if (c == '(') { skipParenDepth++; skipSeenOpen = true; }
                if (c == ')') skipParenDepth--;
            }
            if (skipSeenOpen && skipParenDepth <= 0) {
                skipWaitSemicolon = true;
                skipping = false;
            }
            continue;
        }
        if (skipWaitSemicolon) {
            skipWaitSemicolon = false;
            if (trim(line) == ";") continue;
            /* Not a semicolon — fall through */
        }

        /* Detect nonuniform List */
        if (t.find("nonuniform") != std::string::npos &&
            t.find("List<") != std::string::npos) {
            skipping = true;
            skipParenDepth = 0;
            skipSeenOpen = false;
            for (char c : t) {
                if (c == '(') { skipParenDepth++; skipSeenOpen = true; }
                if (c == ')') skipParenDepth--;
            }
            if (skipSeenOpen && skipParenDepth <= 0) {
                skipping = false;
                skipWaitSemicolon = true;
            }
            continue;
        }

        /* Strip "value uniform ..." lines — we'll write our own */
        {
            std::string lt = trim(line);
            /* Match "value" but not "valueFraction", "freestreamValue", etc.
               "value" must be followed by space/tab. */
            if (lt.find("value") == 0 &&
                lt.size() > 5 && (lt[5] == ' ' || lt[5] == '\t')) {
                continue; /* skip this line */
            }
        }

        curContent += line + "\n";
    }

    return result;
}

static std::vector<RefFieldMeta> readRefFieldMeta(const std::string &refDir)
{
    std::vector<RefFieldMeta> metas;
    DIR *dir = opendir(refDir.c_str());
    if (!dir) return metas;
    struct dirent *ent;
    while ((ent = readdir(dir)) != nullptr) {
        std::string name = ent->d_name;
        if (name == "." || name == ".." || name[0] == '.') continue;
        if (name == "uniform") continue;
        std::string fpath = refDir + "/" + name;
        std::string cls = readFoamClass(fpath);
        if (cls != "volScalarField" && cls != "volVectorField") continue;
        RefFieldMeta m;
        m.name     = name;
        m.isVector = (cls == "volVectorField");
        m.dims     = readDimensions(fpath);
        m.patchBCs = readPatchBCs(fpath);
        metas.push_back(m);
    }
    closedir(dir);
    std::sort(metas.begin(), metas.end(),
              [](const RefFieldMeta &a, const RefFieldMeta &b) {
                  return a.name < b.name;
              });
    return metas;
}

/*----------------------------------------------------------------------------*/
/* Merged patch info                                                          */
/*----------------------------------------------------------------------------*/

struct PatchInfo {
    std::string name;
    std::string type;      /* "patch", "wall", "symmetry" */
    int startFace;         /* global face index in allFaces */
    int nFaces;
    int bndIdx;            /* start index into boundaryFaces vector */
};

/*----------------------------------------------------------------------------*/
/* Main                                                                       */
/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::fprintf(stderr,
            "Usage: gmf2of <input_base> <output_case_dir> [time_dir=0] "
            "[-r <refCase> <refTime>] [-c <createPatchDict>]\n"
            "  -r  copy field names, dimensions, BCs from a reference case\n"
            "  -c  merge refs into named patches via createPatchDict\n");
        return 1;
    }

    std::string inBase, caseDir, timeDir = "0";
    std::string refCase, refTime, createPatchDictFile;
    {
        std::vector<std::string> posArgs;
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-r" && i + 2 < argc) {
                refCase = argv[++i];
                refTime = argv[++i];
            } else if (arg == "-c" && i + 1 < argc) {
                createPatchDictFile = argv[++i];
            } else {
                posArgs.push_back(arg);
            }
        }
        if (posArgs.size() < 2) {
            std::fprintf(stderr, "Error: need at least <input_base> and <output_case_dir>\n");
            return 1;
        }
        inBase  = posArgs[0];
        caseDir = posArgs[1];
        if (posArgs.size() >= 3) timeDir = posArgs[2];
    }

    std::string meshFile = inBase + ".meshb";
    std::string solFile  = inBase + ".solb";
    std::string meshDir  = caseDir + "/constant/polyMesh";
    std::string fieldDir = caseDir + "/" + timeDir;

    /* ---- Read createPatchDict if provided ---- */
    std::vector<PatchDef> patchDefs;
    std::map<int, std::string> refToPatch;   /* ref → merged patch name */
    std::map<std::string, std::string> patchTypeMap; /* name → OF type */
    if (!createPatchDictFile.empty()) {
        patchDefs = readCreatePatchDict(createPatchDictFile);
        for (auto &pd : patchDefs) {
            for (int r : pd.refs)
                refToPatch[r] = pd.name;
            patchTypeMap[pd.name] = pd.type;
            std::printf("  Patch '%s' (type=%s): refs",
                        pd.name.c_str(), pd.type.c_str());
            for (int r : pd.refs) std::printf(" %d", r);
            std::printf("\n");
        }
    }

    /* ---- Open mesh file ---- */
    int ver, dim;
    int64_t msh = GmfOpenMesh(meshFile.c_str(), GmfRead, &ver, &dim);
    if (!msh) {
        std::fprintf(stderr, "Error: cannot open %s\n", meshFile.c_str());
        return 1;
    }
    std::printf("Reading %s (ver=%d, dim=%d) ...\n", meshFile.c_str(), ver, dim);
    if (dim != 3) {
        std::fprintf(stderr, "Error: only 3D meshes supported\n");
        GmfCloseMesh(msh);
        return 1;
    }

    /* ---- Query keyword counts ---- */
    int64_t nVer = GmfStatKwd(msh, GmfVertices);
    int64_t nTet = GmfStatKwd(msh, GmfTetrahedra);
    int64_t nHex = GmfStatKwd(msh, GmfHexahedra);
    int64_t nPri = GmfStatKwd(msh, GmfPrisms);
    int64_t nPyr = GmfStatKwd(msh, GmfPyramids);
    int64_t nBTri = GmfStatKwd(msh, GmfTriangles);
    int64_t nBQad = GmfStatKwd(msh, GmfQuadrilaterals);

    int64_t nCells = nTet + nHex + nPri + nPyr;
    std::printf("  nVer=%lld  nTet=%lld  nHex=%lld  nPri=%lld  nPyr=%lld\n",
                (long long)nVer, (long long)nTet, (long long)nHex,
                (long long)nPri, (long long)nPyr);
    std::printf("  nBndTri=%lld  nBndQad=%lld\n",
                (long long)nBTri, (long long)nBQad);

    /* ---- Read vertices ---- */
    std::vector<std::array<double,3>> pts(nVer);
    {
        GmfGotoKwd(msh, GmfVertices);
        for (int64_t i = 0; i < nVer; i++) {
            int ref;
            GmfGetLin(msh, GmfVertices, &pts[i][0], &pts[i][1], &pts[i][2], &ref);
        }
    }

    /* ---- Read element connectivity (0-based) ---- */
    struct ElemData { int nodes[8]; int ref; };

    auto readElems = [&](int kwd, int nNodes, int64_t n,
                         std::vector<ElemData> &out) {
        out.resize(n);
        if (n == 0) return;
        GmfGotoKwd(msh, kwd);
        for (int64_t i = 0; i < n; i++) {
            int &r = out[i].ref;
            int *nd = out[i].nodes;
            switch (nNodes) {
                case 4: GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&r); break;
                case 5: GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&nd[4],&r); break;
                case 6: GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&nd[4],&nd[5],&r); break;
                case 8: GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],
                                            &nd[4],&nd[5],&nd[6],&nd[7],&r); break;
            }
            for (int k = 0; k < nNodes; k++) nd[k]--;
        }
    };

    std::vector<ElemData> tets, hexes, prisms, pyrs;
    readElems(GmfTetrahedra, 4, nTet, tets);
    readElems(GmfHexahedra,  8, nHex, hexes);
    readElems(GmfPrisms,     6, nPri, prisms);
    readElems(GmfPyramids,   5, nPyr, pyrs);

    /* Boundary surface elements */
    struct BndSurf { std::vector<int> verts; int ref; };
    std::vector<BndSurf> bndSurfs;
    {
        auto readBndElems = [&](int kwd, int nv, int64_t n) {
            if (n == 0) return;
            GmfGotoKwd(msh, kwd);
            for (int64_t i = 0; i < n; i++) {
                BndSurf bs;
                bs.verts.resize(nv);
                int r;
                if (nv == 3)
                    GmfGetLin(msh, kwd, &bs.verts[0],&bs.verts[1],&bs.verts[2],&r);
                else
                    GmfGetLin(msh, kwd, &bs.verts[0],&bs.verts[1],
                                        &bs.verts[2],&bs.verts[3],&r);
                for (int &v : bs.verts) v--;
                bs.ref = r;
                bndSurfs.push_back(bs);
            }
        };
        readBndElems(GmfTriangles,      3, nBTri);
        readBndElems(GmfQuadrilaterals, 4, nBQad);
    }

    std::map<FaceKey,int> bndFaceRef;
    for (auto &bs : bndSurfs)
        bndFaceRef[makeFaceKey(bs.verts)] = bs.ref;

    GmfCloseMesh(msh);

    /* ---- Cell ordering ---- */
    int64_t tetOff = 0;
    int64_t hexOff = nTet;
    int64_t priOff = nTet + nHex;
    int64_t pyrOff = nTet + nHex + nPri;

    /* ---- Build face list via deduplication ---- */
    std::map<FaceKey, FaceRecord> faceMap;

    auto addElemFaces = [&](int64_t cellId, const int *nodes,
                            const std::vector<std::vector<int>> &tmpl) {
        for (auto &localFace : tmpl) {
            std::vector<int> fv;
            fv.reserve(localFace.size());
            for (int li : localFace) fv.push_back(nodes[li]);
            FaceKey key = makeFaceKey(fv);
            auto it = faceMap.find(key);
            if (it == faceMap.end()) {
                FaceRecord fr;
                fr.verts     = fv;
                fr.owner     = (int)cellId;
                fr.neighbour = -1;
                fr.patchRef  = 0;
                auto bit = bndFaceRef.find(key);
                if (bit != bndFaceRef.end())
                    fr.patchRef = bit->second;
                faceMap[key] = fr;
            } else {
                it->second.neighbour = (int)cellId;
            }
        }
    };

    auto ttmpl = tetFaceTemplates();
    auto htmpl = hexFaceTemplates();
    auto ptmpl = prismFaceTemplates();
    auto ytmpl = pyramidFaceTemplates();

    for (int64_t i = 0; i < nTet; i++)
        addElemFaces(tetOff+i, tets[i].nodes,  ttmpl);
    for (int64_t i = 0; i < nHex; i++)
        addElemFaces(hexOff+i, hexes[i].nodes, htmpl);
    for (int64_t i = 0; i < nPri; i++)
        addElemFaces(priOff+i, prisms[i].nodes,ptmpl);
    for (int64_t i = 0; i < nPyr; i++)
        addElemFaces(pyrOff+i, pyrs[i].nodes,  ytmpl);

    /* ---- Separate internal / boundary faces ---- */
    std::vector<FaceRecord> internalFaces, boundaryFaces;
    for (auto &kv : faceMap) {
        if (kv.second.neighbour >= 0)
            internalFaces.push_back(kv.second);
        else
            boundaryFaces.push_back(kv.second);
    }

    std::sort(internalFaces.begin(), internalFaces.end(),
        [](const FaceRecord &a, const FaceRecord &b) {
            if (a.owner != b.owner) return a.owner < b.owner;
            return a.neighbour < b.neighbour;
        });

    /* ---- Build merged patch info ---- */
    /* Assign each boundary face a patch name using the createPatchDict
       mapping (if provided), otherwise each ref gets its own Patch_<ref>. */
    auto patchNameFor = [&](int ref) -> std::string {
        auto it = refToPatch.find(ref);
        if (it != refToPatch.end()) return it->second;
        return "Patch_" + std::to_string(ref);
    };
    auto patchTypeFor = [&](const std::string &name) -> std::string {
        auto it = patchTypeMap.find(name);
        if (it != patchTypeMap.end()) return it->second;
        return "patch";
    };

    /* Determine the desired patch order: use patchDefs order if available,
       then any remaining refs in numeric order. */
    std::vector<std::string> patchOrder;
    std::set<std::string> patchOrderSet;
    for (auto &pd : patchDefs) {
        patchOrder.push_back(pd.name);
        patchOrderSet.insert(pd.name);
    }
    /* Collect all refs from boundary faces */
    std::set<int> allRefs;
    for (auto &bf : boundaryFaces)
        allRefs.insert(bf.patchRef);
    for (int r : allRefs) {
        std::string nm = patchNameFor(r);
        if (!patchOrderSet.count(nm)) {
            patchOrder.push_back(nm);
            patchOrderSet.insert(nm);
        }
    }

    /* Assign each boundary face its patch name, then sort by
       (patch order index, owner) */
    std::map<std::string, int> patchIdx;
    for (int i = 0; i < (int)patchOrder.size(); i++)
        patchIdx[patchOrder[i]] = i;

    std::sort(boundaryFaces.begin(), boundaryFaces.end(),
        [&](const FaceRecord &a, const FaceRecord &b) {
            int ia = patchIdx[patchNameFor(a.patchRef)];
            int ib = patchIdx[patchNameFor(b.patchRef)];
            if (ia != ib) return ia < ib;
            return a.owner < b.owner;
        });

    /* Build PatchInfo array */
    std::vector<PatchInfo> patches;
    {
        int bndIdx = 0;
        int faceStart = (int)internalFaces.size();
        for (auto &pname : patchOrder) {
            PatchInfo pi;
            pi.name = pname;
            pi.type = patchTypeFor(pname);
            pi.bndIdx = bndIdx;
            pi.startFace = faceStart;
            pi.nFaces = 0;
            /* Count faces belonging to this patch */
            while (bndIdx < (int)boundaryFaces.size() &&
                   patchNameFor(boundaryFaces[bndIdx].patchRef) == pname) {
                pi.nFaces++;
                bndIdx++;
            }
            if (pi.nFaces > 0) {
                faceStart += pi.nFaces;
                patches.push_back(pi);
            }
        }
    }

    /* Full ordered face list */
    std::vector<FaceRecord> allFaces;
    allFaces.insert(allFaces.end(), internalFaces.begin(), internalFaces.end());
    allFaces.insert(allFaces.end(), boundaryFaces.begin(), boundaryFaces.end());

    int nAllFaces      = (int)allFaces.size();
    int nInternalFaces = (int)internalFaces.size();

    std::printf("  nAllFaces=%d  nInternalFaces=%d  nBndFaces=%d  nPatches=%d\n",
                nAllFaces, nInternalFaces, (int)boundaryFaces.size(),
                (int)patches.size());
    for (auto &pi : patches)
        std::printf("    %s (%s): %d faces\n",
                    pi.name.c_str(), pi.type.c_str(), pi.nFaces);

    /* ---- Create output directories ---- */
    mkdirp(caseDir);
    mkdirp(caseDir + "/constant");
    mkdirp(meshDir);
    mkdirp(fieldDir);

    /* ---- Write points ---- */
    {
        std::ofstream out(meshDir + "/points");
        writeFoamHeader(out, "vectorField", "constant/polyMesh", "points");
        out << (long long)nVer << "\n(\n";
        out << std::scientific;
        out.precision(10);
        for (auto &p : pts)
            out << "(" << p[0] << " " << p[1] << " " << p[2] << ")\n";
        out << ")\n";
    }
    std::printf("  Wrote points (%lld)\n", (long long)nVer);

    /* ---- Write faces ---- */
    {
        std::ofstream out(meshDir + "/faces");
        writeFoamHeader(out, "faceList", "constant/polyMesh", "faces");
        out << nAllFaces << "\n(\n";
        for (auto &fr : allFaces) {
            out << fr.verts.size() << "(";
            for (int k = 0; k < (int)fr.verts.size(); k++) {
                if (k) out << " ";
                out << fr.verts[k];
            }
            out << ")\n";
        }
        out << ")\n";
    }
    std::printf("  Wrote faces (%d)\n", nAllFaces);

    /* ---- Write owner ---- */
    {
        std::ofstream out(meshDir + "/owner");
        writeFoamHeader(out, "labelList", "constant/polyMesh", "owner",
                        (int)nCells, nAllFaces, (int)nVer, nInternalFaces);
        out << nAllFaces << "\n(\n";
        for (auto &fr : allFaces) out << fr.owner << "\n";
        out << ")\n";
    }

    /* ---- Write neighbour ---- */
    {
        std::ofstream out(meshDir + "/neighbour");
        writeFoamHeader(out, "labelList", "constant/polyMesh", "neighbour",
                        (int)nCells, nAllFaces, (int)nVer, nInternalFaces);
        out << nInternalFaces << "\n(\n";
        for (int i = 0; i < nInternalFaces; i++)
            out << allFaces[i].neighbour << "\n";
        out << ")\n";
    }
    std::printf("  Wrote owner/neighbour\n");

    /* ---- Write boundary ---- */
    {
        std::ofstream out(meshDir + "/boundary");
        writeFoamHeader(out, "polyBoundaryMesh", "constant/polyMesh", "boundary");
        out << (int)patches.size() << "\n(\n";
        for (auto &pi : patches) {
            out << "    " << pi.name << "\n    {\n"
                << "        type            " << pi.type << ";\n"
                << "        nFaces          " << pi.nFaces << ";\n"
                << "        startFace       " << pi.startFace << ";\n"
                << "    }\n";
        }
        out << ")\n";
    }
    std::printf("  Wrote boundary (%d patches)\n", (int)patches.size());

    /* ---- Read solution (.solb) ---- */
    int64_t sol = GmfOpenMesh(solFile.c_str(), GmfRead, &ver, &dim);
    if (!sol) {
        std::printf("No .solb found at %s, skipping field output.\n", solFile.c_str());
        std::printf("Done.\n");
        return 0;
    }
    std::printf("Reading %s ...\n", solFile.c_str());

    int NmbTyp = 0, SolSiz = 0;
    int SolTypArr[GmfMaxTyp] = {};
    bool haveSolAtElems = false;

    for (int kwd : {GmfSolAtTetrahedra, GmfSolAtHexahedra,
                    GmfSolAtPrisms, GmfSolAtPyramids}) {
        int nt=0, ss=0, st[GmfMaxTyp]={};
        int64_t n = GmfStatKwd(sol, kwd, &nt, &ss, st);
        if (n > 0 && !haveSolAtElems) {
            NmbTyp = nt; SolSiz = ss;
            memcpy(SolTypArr, st, nt * sizeof(int));
            haveSolAtElems = true;
        }
    }

    bool haveSolAtVerts = false;
    int64_t nSolVerts = 0;
    if (!haveSolAtElems) {
        int nt=0, ss=0, st[GmfMaxTyp]={};
        nSolVerts = GmfStatKwd(sol, GmfSolAtVertices, &nt, &ss, st);
        if (nSolVerts > 0) {
            NmbTyp = nt; SolSiz = ss;
            memcpy(SolTypArr, st, nt * sizeof(int));
            haveSolAtVerts = true;
        }
    }

    if (NmbTyp == 0) {
        std::printf("No recognisable solution fields in .solb, skipping.\n");
        GmfCloseMesh(sol);
        std::printf("Done.\n");
        return 0;
    }

    /* ---- Read field names from ReferenceStrings ---- */
    std::vector<std::string> refNames(NmbTyp);
    {
        int64_t nRef = GmfStatKwd(sol, GmfReferenceStrings);
        if (nRef > 0) {
            GmfGotoKwd(sol, GmfReferenceStrings);
            for (int64_t r = 0; r < nRef; r++) {
                int kwdCode, subIdx;
                char nameBuf[256] = {};
                GmfGetLin(sol, GmfReferenceStrings,
                          &kwdCode, &subIdx, nameBuf);
                std::string s(nameBuf);
                auto pos = s.rfind(' ');
                if (pos != std::string::npos && pos + 1 < s.size()) {
                    int idx = std::atoi(s.substr(pos + 1).c_str());
                    std::string name = s.substr(0, pos);
                    if (idx >= 1 && idx <= NmbTyp)
                        refNames[idx - 1] = name;
                }
            }
        }
    }
    bool haveRefNames = false;
    for (auto &n : refNames)
        if (!n.empty()) { haveRefNames = true; break; }

    /* ---- Build per-field descriptors ---- */
    struct FieldDesc {
        std::string name;
        bool        isVector;
        int         solOffset;
        std::string dims;
        /* Per-patch BC content from reference case (patchName → BC lines) */
        std::map<std::string, std::string> patchBCs;
    };
    std::vector<FieldDesc> fieldDescs;

    static const std::map<std::string, std::string> knownDims = {
        {"p",       "[1 -1 -2 0 0 0 0]"},
        {"T",       "[0 0 0 1 0 0 0]"},
        {"U",       "[0 1 -1 0 0 0 0]"},
        {"nuTilda", "[0 2 -1 0 0 0 0]"},
        {"nut",     "[0 2 -1 0 0 0 0]"},
        {"k",       "[0 2 -2 0 0 0 0]"},
        {"omega",   "[0 0 -1 0 0 0 0]"},
        {"e",       "[0 2 -2 0 0 0 0]"},
        {"alphat",  "[1 -1 -1 0 0 0 0]"},
    };

    /* ---- Read reference case metadata ---- */
    std::vector<RefFieldMeta> refMetas;
    if (!refCase.empty()) {
        std::string refDir = refCase + "/" + refTime;
        refMetas = readRefFieldMeta(refDir);
        if (refMetas.empty()) {
            std::fprintf(stderr, "Warning: no fields found in %s, "
                         "ignoring -r\n", refDir.c_str());
        } else {
            std::printf("  Reference case: %zu fields from %s\n",
                        refMetas.size(), refDir.c_str());
            for (auto &m : refMetas)
                std::printf("    %s (%s)\n", m.name.c_str(),
                            m.isVector ? "vector" : "scalar");
        }
    }

    if (!refMetas.empty()) {
        /* Path C: use reference case metadata */
        int off = 0;
        int solIdx = 0;
        for (auto &rm : refMetas) {
            FieldDesc fd;
            fd.name     = rm.name;
            fd.isVector = rm.isVector;
            fd.dims     = rm.dims;
            fd.patchBCs = rm.patchBCs;
            fd.solOffset = off;
            if (rm.isVector) {
                if (solIdx < NmbTyp && SolTypArr[solIdx] == GmfVec) {
                    off += 3; solIdx++;
                } else if (solIdx + 2 < NmbTyp &&
                           SolTypArr[solIdx] == GmfSca &&
                           SolTypArr[solIdx+1] == GmfSca &&
                           SolTypArr[solIdx+2] == GmfSca) {
                    off += 3; solIdx += 3;
                } else {
                    std::fprintf(stderr, "Warning: cannot map vector field "
                                 "'%s' at solIdx=%d\n", rm.name.c_str(), solIdx);
                    off += 3; solIdx++;
                }
            } else {
                off += 1; solIdx++;
            }
            fieldDescs.push_back(fd);
        }
        std::printf("  Fields from reference case:");
        for (auto &fd : fieldDescs)
            std::printf(" %s(%s)", fd.name.c_str(),
                        fd.isVector ? "vec" : "sca");
        std::printf("\n");
        if (off != SolSiz) {
            std::fprintf(stderr, "Warning: reference fields consume %d values "
                         "but solb has %d per entry\n", off, SolSiz);
        }
    } else if (haveRefNames) {
        int off = 0;
        for (int t = 0; t < NmbTyp; t++) {
            FieldDesc fd;
            fd.name = refNames[t].empty()
                ? ("field_" + std::to_string(t))
                : refNames[t];
            fd.solOffset = off;
            fd.isVector  = (SolTypArr[t] == GmfVec);
            auto it = knownDims.find(fd.name);
            fd.dims = (it != knownDims.end())
                ? it->second : "[0 0 0 0 0 0 0]";
            fieldDescs.push_back(fd);
            off += fd.isVector ? 3 : 1;
        }
        std::printf("  Fields from ReferenceStrings:");
        for (auto &fd : fieldDescs)
            std::printf(" %s(%s)", fd.name.c_str(),
                        fd.isVector ? "vec" : "sca");
        std::printf("\n");
    } else {
        int nScaFound = 0;
        int off = 0;
        for (int t = 0; t < NmbTyp; t++) {
            if (SolTypArr[t] == GmfSca) {
                FieldDesc fd;
                fd.isVector  = false;
                fd.solOffset = off;
                if (nScaFound == 0) {
                    fd.name = "p"; fd.dims = "[1 -1 -2 0 0 0 0]";
                } else if (nScaFound == 1) {
                    fd.name = "T"; fd.dims = "[0 0 0 1 0 0 0]";
                } else {
                    fd.name = "scalar_" + std::to_string(nScaFound);
                    fd.dims = "[0 0 0 0 0 0 0]";
                }
                fieldDescs.push_back(fd);
                nScaFound++;
                off++;
            } else if (SolTypArr[t] == GmfVec) {
                FieldDesc fd;
                fd.name = "U"; fd.dims = "[0 1 -1 0 0 0 0]";
                fd.isVector  = true;
                fd.solOffset = off;
                fieldDescs.push_back(fd);
                off += 3;
            } else {
                off += SolSiz;
            }
        }
        std::printf("  Fields (legacy convention):");
        for (auto &fd : fieldDescs)
            std::printf(" %s(%s)", fd.name.c_str(),
                        fd.isVector ? "vec" : "sca");
        std::printf("\n");
    }

    /* ---- Read solution data ---- */
    auto readSolBuf = [&](int kwd, int64_t n) -> std::vector<double> {
        std::vector<double> buf(n * SolSiz, 0.0);
        if (n == 0) return buf;
        GmfGotoKwd(sol, kwd);
        for (int64_t i = 0; i < n; i++)
            GmfGetLin(sol, kwd, &buf[i * SolSiz]);
        return buf;
    };

    std::vector<double> flatSol(nCells * SolSiz, 0.0);
    std::vector<double> vertSol; /* kept for boundary face interp */

    if (haveSolAtElems) {
        std::printf("  Reading SolAtElements ...\n");
        auto tetSolBuf = readSolBuf(GmfSolAtTetrahedra, nTet);
        auto hexSolBuf = readSolBuf(GmfSolAtHexahedra,  nHex);
        auto priSolBuf = readSolBuf(GmfSolAtPrisms,     nPri);
        auto pyrSolBuf = readSolBuf(GmfSolAtPyramids,   nPyr);

        auto copySol = [&](const std::vector<double> &src, int64_t offset, int64_t n) {
            for (int64_t i = 0; i < n; i++)
                for (int k = 0; k < SolSiz; k++)
                    flatSol[(offset+i)*SolSiz + k] = src[i*SolSiz + k];
        };
        copySol(tetSolBuf, tetOff, nTet);
        copySol(hexSolBuf, hexOff, nHex);
        copySol(priSolBuf, priOff, nPri);
        copySol(pyrSolBuf, pyrOff, nPyr);
    } else {
        std::printf("  Reading SolAtVertices (%lld), interpolating to cells ...\n",
                    (long long)nSolVerts);
        vertSol = readSolBuf(GmfSolAtVertices, nSolVerts);

        auto interpCells = [&](const std::vector<ElemData> &elems,
                               int nNodes, int64_t cellOff) {
            for (int64_t i = 0; i < (int64_t)elems.size(); i++) {
                double *dst = &flatSol[(cellOff + i) * SolSiz];
                for (int k = 0; k < nNodes; k++) {
                    int vi = elems[i].nodes[k];
                    const double *src = &vertSol[vi * SolSiz];
                    for (int s = 0; s < SolSiz; s++)
                        dst[s] += src[s];
                }
                for (int s = 0; s < SolSiz; s++)
                    dst[s] /= nNodes;
            }
        };
        interpCells(tets,   4, tetOff);
        interpCells(hexes,  8, hexOff);
        interpCells(prisms, 6, priOff);
        interpCells(pyrs,   5, pyrOff);
    }

    GmfCloseMesh(sol);

    /* ---- Compute boundary face values ---- */
    std::vector<std::vector<double>> bndFaceSol(boundaryFaces.size());
    if (!vertSol.empty()) {
        /* Interpolate from vertex data */
        std::printf("  Interpolating boundary face values from vertices ...\n");
        for (size_t i = 0; i < boundaryFaces.size(); i++) {
            bndFaceSol[i].resize(SolSiz, 0.0);
            int nv = (int)boundaryFaces[i].verts.size();
            for (int v : boundaryFaces[i].verts) {
                const double *src = &vertSol[v * SolSiz];
                for (int s = 0; s < SolSiz; s++)
                    bndFaceSol[i][s] += src[s];
            }
            for (int s = 0; s < SolSiz; s++)
                bndFaceSol[i][s] /= nv;
        }
    } else {
        /* Fallback: use owner cell value */
        std::printf("  Using owner cell values for boundary faces ...\n");
        for (size_t i = 0; i < boundaryFaces.size(); i++) {
            bndFaceSol[i].resize(SolSiz, 0.0);
            int oc = boundaryFaces[i].owner;
            const double *src = &flatSol[oc * SolSiz];
            for (int s = 0; s < SolSiz; s++)
                bndFaceSol[i][s] = src[s];
        }
    }

    /* ---- Write field files ---- */
    /* Helper: find BC content for a given patch name in a field's patchBCs.
       Falls back to ".*" regex entry, then to default zeroGradient. */
    auto findPatchBC = [](const std::map<std::string, std::string> &patchBCs,
                          const std::string &patchName) -> std::string {
        auto it = patchBCs.find(patchName);
        if (it != patchBCs.end()) return it->second;
        /* Try regex ".*" */
        it = patchBCs.find("\".*\"");
        if (it != patchBCs.end()) return it->second;
        return "        type            zeroGradient;\n";
    };

    auto writeBoundaryField = [&](std::ofstream &out, const FieldDesc &fd,
                                  bool isVector) {
        out << "boundaryField\n{\n";
        for (auto &pi : patches) {
            out << "    " << pi.name << "\n    {\n";
            std::string bc = findPatchBC(fd.patchBCs, pi.name);
            out << bc;
            /* For symmetry/empty patches, don't write value */
            if (pi.type == "symmetry" || pi.type == "empty") {
                out << "    }\n";
                continue;
            }
            /* Write interpolated value */
            if (isVector) {
                out << "        value           nonuniform List<vector>\n"
                    << pi.nFaces << "\n(\n";
                out << std::scientific;
                out.precision(10);
                for (int i = 0; i < pi.nFaces; i++) {
                    const auto &row = bndFaceSol[pi.bndIdx + i];
                    out << "(" << row[fd.solOffset+0] << " "
                               << row[fd.solOffset+1] << " "
                               << row[fd.solOffset+2] << ")\n";
                }
            } else {
                out << "        value           nonuniform List<scalar>\n"
                    << pi.nFaces << "\n(\n";
                out << std::scientific;
                out.precision(10);
                for (int i = 0; i < pi.nFaces; i++)
                    out << bndFaceSol[pi.bndIdx + i][fd.solOffset] << "\n";
            }
            out << ")\n;\n";
            out << "    }\n";
        }
        out << "}\n";
    };

    auto writeScalarField = [&](const FieldDesc &fd) {
        std::ofstream out(fieldDir + "/" + fd.name);
        writeFoamHeader(out, "volScalarField", timeDir, fd.name);
        out << "dimensions      " << fd.dims << ";\n\n"
            << "internalField   nonuniform List<scalar>\n"
            << nCells << "\n(\n";
        out << std::scientific; out.precision(10);
        for (int64_t c = 0; c < nCells; c++)
            out << flatSol[c*SolSiz + fd.solOffset] << "\n";
        out << ");\n\n";
        writeBoundaryField(out, fd, false);
    };

    auto writeVectorField = [&](const FieldDesc &fd) {
        std::ofstream out(fieldDir + "/" + fd.name);
        writeFoamHeader(out, "volVectorField", timeDir, fd.name);
        out << "dimensions      " << fd.dims << ";\n\n"
            << "internalField   nonuniform List<vector>\n"
            << nCells << "\n(\n";
        out << std::scientific; out.precision(10);
        for (int64_t c = 0; c < nCells; c++) {
            const double *row = &flatSol[c*SolSiz];
            out << "(" << row[fd.solOffset+0] << " "
                       << row[fd.solOffset+1] << " "
                       << row[fd.solOffset+2] << ")\n";
        }
        out << ");\n\n";
        writeBoundaryField(out, fd, true);
    };

    for (auto &fd : fieldDescs) {
        if (fd.isVector)
            writeVectorField(fd);
        else
            writeScalarField(fd);
        std::printf("  Wrote %s\n", fd.name.c_str());
    }

    std::printf("Done.\n");
    return 0;
}
