/*----------------------------------------------------------------------------*/
/*  of2gmf.cpp                                                                */
/*  Convert an OpenFOAM polyMesh + solution fields (p, T, U) to libMeshb     */
/*  binary format (.meshb + .solb).                                           */
/*                                                                            */
/*  Usage:                                                                    */
/*    of2gmf <case_dir> <time_dir> <output_base>                              */
/*                                                                            */
/*  Example:                                                                  */
/*    of2gmf ./myCase 1000 result                                             */
/*    -> writes result.meshb, result.solb                                     */
/*----------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
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
#include <stdexcept>
#include <dirent.h>

#include "libmeshb7.h"

/*----------------------------------------------------------------------------*/
/* Helpers: OpenFOAM file parser                                              */
/*----------------------------------------------------------------------------*/

/* Skip FoamFile header block and the "// * * *" separator line.
   Positions the stream just after that line, ready to read data. */
static void skipFoamHeader(std::istream &in)
{
    std::string line;
    /* Scan until the separator line (starts with "// *") */
    while (std::getline(in, line)) {
        if (line.rfind("// *", 0) == 0)
            return;
    }
    throw std::runtime_error("Malformed OpenFOAM file: separator not found");
}

/* Read the FoamFile header and return the value of the 'class' field.
   Returns "" if the file cannot be opened or does not have a FoamFile header. */
static std::string detectFoamClass(const std::string &path)
{
    std::ifstream in(path);
    if (!in) return "";
    std::string line;
    bool inFoamFile = false;
    while (std::getline(in, line)) {
        if (line.find("FoamFile") != std::string::npos) { inFoamFile = true; continue; }
        if (inFoamFile && line.find("class") != std::string::npos) {
            std::istringstream ss(line);
            std::string key, val;
            ss >> key >> val;
            if (!val.empty() && val.back() == ';') val.pop_back();
            return val;
        }
        if (line.rfind("// *", 0) == 0) break;
    }
    return "";
}

/* Scan <fieldDir> for files whose FoamFile class is volScalarField or
   volVectorField and return their names sorted alphabetically. */
static std::vector<std::string> discoverFields(const std::string &fieldDir)
{
    std::vector<std::string> names;
    DIR *dir = opendir(fieldDir.c_str());
    if (!dir) return names;
    struct dirent *ent;
    while ((ent = readdir(dir)) != nullptr) {
        std::string name = ent->d_name;
        if (name == "." || name == ".." || name[0] == '.') continue;
        std::string cls = detectFoamClass(fieldDir + "/" + name);
        if (cls == "volScalarField" || cls == "volVectorField")
            names.push_back(name);
    }
    closedir(dir);
    std::sort(names.begin(), names.end());
    return names;
}

/* Skip whitespace and comments ('//' to end-of-line) */
static void skipWS(std::istream &in)
{
    while (true) {
        in >> std::ws;
        if (in.peek() == '/') {
            char c1, c2;
            in.get(c1); in.get(c2);
            if (c2 == '/') {
                std::string dummy;
                std::getline(in, dummy);
            } else {
                in.unget(); in.unget();
                break;
            }
        } else {
            break;
        }
    }
}

/* Read the integer count that appears before the opening '(' */
static int readCount(std::istream &in)
{
    int n;
    skipWS(in);
    if (!(in >> n))
        throw std::runtime_error("Expected integer count");
    return n;
}

/* Expect and consume a single character */
static void expect(std::istream &in, char c)
{
    skipWS(in);
    char got;
    if (!(in >> got) || got != c) {
        std::string msg = "Expected '";
        msg += c; msg += "'";
        throw std::runtime_error(msg);
    }
}

/*----------------------------------------------------------------------------*/
/* Read points: returns vector of (x,y,z)                                    */
/*----------------------------------------------------------------------------*/
static std::vector<std::array<double,3>> readPoints(const std::string &path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);
    skipFoamHeader(in);

    int n = readCount(in);
    expect(in, '(');

    std::vector<std::array<double,3>> pts(n);
    for (int i = 0; i < n; i++) {
        expect(in, '(');
        in >> pts[i][0] >> pts[i][1] >> pts[i][2];
        expect(in, ')');
    }
    expect(in, ')');
    return pts;
}

/*----------------------------------------------------------------------------*/
/* Read faces: returns vector of face vertex lists (0-based)                  */
/*----------------------------------------------------------------------------*/
static std::vector<std::vector<int>> readFaces(const std::string &path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);
    skipFoamHeader(in);

    int n = readCount(in);
    expect(in, '(');

    std::vector<std::vector<int>> faces(n);
    for (int i = 0; i < n; i++) {
        skipWS(in);
        int nv;
        if (!(in >> nv))
            throw std::runtime_error("Expected vertex count for face " + std::to_string(i));
        expect(in, '(');
        faces[i].resize(nv);
        for (int j = 0; j < nv; j++)
            in >> faces[i][j];
        expect(in, ')');
    }
    expect(in, ')');
    return faces;
}

/*----------------------------------------------------------------------------*/
/* Read label list (owner or neighbour)                                       */
/* Also extracts nCells and nInternalFaces from the 'note' field if present   */
/*----------------------------------------------------------------------------*/
static std::vector<int> readLabelList(const std::string &path,
                                      int *nCells_out, int *nIntFaces_out)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    if (nCells_out)    *nCells_out    = -1;
    if (nIntFaces_out) *nIntFaces_out = -1;

    /* Search for 'note' line to extract mesh stats */
    std::string line;
    while (std::getline(in, line)) {
        /* Note format: note "nPoints:N nCells:N nFaces:N nInternalFaces:N" */
        auto pos = line.find("note");
        if (pos != std::string::npos) {
            /* Parse nCells and nInternalFaces */
            auto parseVal = [&](const char *key) -> int {
                auto kp = line.find(key);
                if (kp == std::string::npos) return -1;
                kp += strlen(key);
                /* skip ':' */
                while (kp < line.size() && !isdigit((unsigned char)line[kp])) kp++;
                return std::stoi(line.substr(kp));
            };
            if (nCells_out)    *nCells_out    = parseVal("nCells");
            if (nIntFaces_out) *nIntFaces_out = parseVal("nInternalFaces");
        }
        if (line.rfind("// *", 0) == 0) break;
    }

    int n = readCount(in);
    expect(in, '(');

    std::vector<int> data(n);
    for (int i = 0; i < n; i++)
        in >> data[i];
    expect(in, ')');
    return data;
}

/*----------------------------------------------------------------------------*/
/* Boundary patch descriptor                                                  */
/*----------------------------------------------------------------------------*/
struct Patch {
    std::string name;
    std::string type;
    int startFace;
    int nFaces;
};

static std::vector<Patch> readBoundary(const std::string &path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);
    skipFoamHeader(in);

    int np = readCount(in);
    expect(in, '(');

    std::vector<Patch> patches(np);
    for (int i = 0; i < np; i++) {
        skipWS(in);
        in >> patches[i].name;
        expect(in, '{');

        /* Read key-value pairs until '}' */
        std::string key;
        while (skipWS(in), in >> key) {
            if (key == "}") break;
            if (key == "type") {
                in >> patches[i].type;
                /* consume trailing ';' */
                std::string dummy; std::getline(in, dummy);
            } else if (key == "nFaces") {
                in >> patches[i].nFaces;
                std::string dummy; std::getline(in, dummy);
            } else if (key == "startFace") {
                in >> patches[i].startFace;
                std::string dummy; std::getline(in, dummy);
            } else {
                /* skip value */
                std::string dummy; std::getline(in, dummy);
            }
        }
    }
    expect(in, ')');
    return patches;
}

/*----------------------------------------------------------------------------*/
/* Field reading (p, T, U)                                                    */
/* Returns: scalar values per cell (for p,T) or 3-component (for U)          */
/*----------------------------------------------------------------------------*/
struct FieldData {
    bool    present;
    bool    isVector;
    int     nCells;
    std::vector<double> values; /* scalar: size nCells; vector: size nCells*3 */
};

static FieldData readField(const std::string &path, int nCells)
{
    FieldData fd;
    fd.present  = false;
    fd.nCells   = nCells;

    std::string cls = detectFoamClass(path);
    if (cls != "volScalarField" && cls != "volVectorField") return fd;

    std::ifstream in(path);
    if (!in) return fd; /* file not present */

    bool expectVector = (cls == "volVectorField");
    fd.present  = true;
    fd.isVector = expectVector;
    skipFoamHeader(in);

    /* Scan for 'internalField' */
    std::string tok;
    while (in >> tok) {
        if (tok == "internalField") break;
    }
    if (!in) throw std::runtime_error("No internalField in " + path);

    in >> tok; /* "uniform" or "nonuniform" */

    if (tok == "uniform") {
        if (!expectVector) {
            double val;
            in >> val;
            fd.values.assign(nCells, val);
        } else {
            expect(in, '(');
            double ux, uy, uz;
            in >> ux >> uy >> uz;
            expect(in, ')');
            fd.values.resize(nCells * 3);
            for (int c = 0; c < nCells; c++) {
                fd.values[c*3+0] = ux;
                fd.values[c*3+1] = uy;
                fd.values[c*3+2] = uz;
            }
        }
    } else if (tok == "nonuniform") {
        /* skip "List<scalar>" or "List<vector>" */
        std::string dummy; std::getline(in, dummy);

        int n = readCount(in);
        if (n != nCells) {
            std::fprintf(stderr, "Warning: field %s has %d values, expected %d\n",
                         path.c_str(), n, nCells);
        }
        expect(in, '(');
        if (!expectVector) {
            fd.values.resize(n);
            for (int i = 0; i < n; i++) in >> fd.values[i];
        } else {
            fd.values.resize(n * 3);
            for (int i = 0; i < n; i++) {
                expect(in, '(');
                in >> fd.values[i*3+0] >> fd.values[i*3+1] >> fd.values[i*3+2];
                expect(in, ')');
            }
        }
        expect(in, ')');
    } else {
        throw std::runtime_error("Unknown internalField type '" + tok + "' in " + path);
    }

    return fd;
}

/*----------------------------------------------------------------------------*/
/* Cell type identification and node extraction                                */
/*----------------------------------------------------------------------------*/

enum CellType { TET=0, HEX=1, PRISM=2, PYRAMID=3, POLY=4 };

static const int N_NODES[] = {4, 8, 6, 5, 0};

/* For a tet cell: extract 4 nodes in INRIA GMF orientation.
   INRIA convention: face (n0,n1,n2) viewed from n3 appears CCW,
   i.e. the right-hand-rule normal of (n0,n1,n2) points TOWARD n3.

   face[f0] = base tri [a,b,c], apex = 4th vertex.
   If the cell OWNS f0, the OF face normal (a->b->c) points outward
   (away from apex) -> flip b,c so normal points toward apex (INRIA).
   If the cell is NEIGHBOUR of f0, the OF face normal already points
   toward the cell (toward apex) -> keep as-is (already INRIA). */
static void extractTet(const std::vector<std::vector<int>> &faces,
                       const std::vector<int> &owner,
                       const std::vector<int> &cellFaceList,
                       int cellId,
                       int (&out)[8])
{
    /* Find first tri face */
    int f0 = -1;
    for (int f : cellFaceList) {
        if ((int)faces[f].size() == 3) { f0 = f; break; }
    }
    assert(f0 >= 0);

    /* Collect all unique vertices */
    std::set<int> allVerts;
    for (int f : cellFaceList)
        for (int v : faces[f]) allVerts.insert(v);
    assert(allVerts.size() == 4);

    /* Base from f0, apex = remaining vertex */
    int a = faces[f0][0], b = faces[f0][1], c = faces[f0][2];
    int apex = -1;
    for (int v : allVerts)
        if (v != a && v != b && v != c) { apex = v; break; }

    if (owner[f0] == cellId) {
        /* f0 outward: normal (a->b->c) points away from apex; flip for INRIA */
        out[0]=a+1; out[1]=c+1; out[2]=b+1; out[3]=apex+1;
    } else {
        /* f0 inward: normal (a->b->c) points toward apex; already INRIA */
        out[0]=a+1; out[1]=b+1; out[2]=c+1; out[3]=apex+1;
    }
}

/* For a prism: two tri faces (bottom/top) + 3 quad side faces.
   Match bottom vertices to top vertices via side quads. */
static void extractPrism(const std::vector<std::vector<int>> &faces,
                         const std::vector<int> &owner,
                         const std::vector<int> &cellFaceList,
                         int cellId,
                         int (&out)[8])
{
    /* Separate tri and quad faces */
    std::vector<int> triFaces, quadFaces;
    for (int f : cellFaceList) {
        if ((int)faces[f].size() == 3) triFaces.push_back(f);
        else                            quadFaces.push_back(f);
    }
    assert(triFaces.size() == 2 && quadFaces.size() == 3);

    int fBot = triFaces[0];
    /* Bottom vertices */
    int b0 = faces[fBot][0], b1 = faces[fBot][1], b2 = faces[fBot][2];

    /* If cell is neighbour of fBot, flip winding */
    bool flipBot = (owner[fBot] != cellId);
    if (flipBot) std::swap(b1, b2);

    /* For each bottom vertex find matching top vertex via a shared quad face */
    auto matchTop = [&](int bv) -> int {
        for (int fq : quadFaces) {
            const auto &fv = faces[fq];
            /* Check if bv is in this quad */
            int idx = -1;
            for (int k = 0; k < 4; k++)
                if (fv[k] == bv) { idx = k; break; }
            if (idx < 0) continue;
            /* The top vertex is opposite in the quad face.
               For a prism side quad [bvA, bvB, tvB, tvA],
               both top vertices are the ones not in the bottom tri. */
            for (int k = 0; k < 4; k++) {
                int v = fv[k];
                if (v != b0 && v != b1 && v != b2)
                    return v; /* first top vertex found */
            }
        }
        return -1; /* should not happen */
    };

    int t0 = matchTop(b0), t1 = matchTop(b1), t2 = matchTop(b2);

    out[0]=b0+1; out[1]=b1+1; out[2]=b2+1;
    out[3]=t0+1; out[4]=t1+1; out[5]=t2+1;
}

/* For a hex: 6 quad faces. Find bottom, top (no shared verts), match via sides. */
static void extractHex(const std::vector<std::vector<int>> &faces,
                       const std::vector<int> &owner,
                       const std::vector<int> &cellFaceList,
                       int cellId,
                       int (&out)[8])
{
    /* Collect 6 quad faces */
    std::vector<int> qf;
    for (int f : cellFaceList) qf.push_back(f);
    assert(qf.size() == 6);

    int fBot = qf[0];

    /* Find top face: no shared vertices with fBot */
    std::set<int> botVerts(faces[fBot].begin(), faces[fBot].end());
    int fTop = -1;
    for (int i = 1; i < 6; i++) {
        bool share = false;
        for (int v : faces[qf[i]])
            if (botVerts.count(v)) { share = true; break; }
        if (!share) { fTop = qf[i]; break; }
    }
    assert(fTop >= 0);

    /* Collect top vertex set */
    std::set<int> topVerts(faces[fTop].begin(), faces[fTop].end());

    /* Side faces */
    std::vector<int> sideFaces;
    for (int f : qf) if (f != fBot && f != fTop) sideFaces.push_back(f);
    assert(sideFaces.size() == 4);

    /* Bottom vertices in CCW order from fBot (respect orientation) */
    int b0,b1,b2,b3;
    if (owner[fBot] == cellId) {
        b0=faces[fBot][0]; b1=faces[fBot][1];
        b2=faces[fBot][2]; b3=faces[fBot][3];
    } else {
        b0=faces[fBot][0]; b1=faces[fBot][3];
        b2=faces[fBot][2]; b3=faces[fBot][1];
    }

    /* Match each bottom vertex to a top vertex via a shared side face */
    auto matchTop = [&](int bv) -> int {
        for (int fs : sideFaces) {
            const auto &fv = faces[fs];
            bool hasBv = false;
            for (int v : fv) if (v == bv) { hasBv = true; break; }
            if (!hasBv) continue;
            /* Return the top vertex in this side face adjacent to bv */
            int idx = -1;
            for (int k = 0; k < 4; k++) if (fv[k] == bv) idx = k;
            /* Adjacent vertices in the face (CCW): fv[idx-1] and fv[idx+1] */
            int prev = fv[(idx+3)%4], next = fv[(idx+1)%4];
            if (topVerts.count(prev)) return prev;
            if (topVerts.count(next)) return next;
        }
        return -1;
    };

    int t0=matchTop(b0), t1=matchTop(b1), t2=matchTop(b2), t3=matchTop(b3);

    out[0]=b0+1; out[1]=b1+1; out[2]=b2+1; out[3]=b3+1;
    out[4]=t0+1; out[5]=t1+1; out[6]=t2+1; out[7]=t3+1;
}

/* For a pyramid: 1 quad base + 4 tri faces. */
static void extractPyramid(const std::vector<std::vector<int>> &faces,
                           const std::vector<int> &owner,
                           const std::vector<int> &cellFaceList,
                           int cellId,
                           int (&out)[8])
{
    int fBase = -1;
    for (int f : cellFaceList)
        if ((int)faces[f].size() == 4) { fBase = f; break; }
    assert(fBase >= 0);

    std::set<int> allVerts;
    for (int f : cellFaceList)
        for (int v : faces[f]) allVerts.insert(v);
    assert(allVerts.size() == 5);

    std::set<int> baseVerts(faces[fBase].begin(), faces[fBase].end());
    int apex = -1;
    for (int v : allVerts) if (!baseVerts.count(v)) { apex = v; break; }

    int q0,q1,q2,q3;
    if (owner[fBase] == cellId) {
        q0=faces[fBase][0]; q1=faces[fBase][1];
        q2=faces[fBase][2]; q3=faces[fBase][3];
    } else {
        q0=faces[fBase][0]; q1=faces[fBase][3];
        q2=faces[fBase][2]; q3=faces[fBase][1];
    }

    out[0]=q0+1; out[1]=q1+1; out[2]=q2+1; out[3]=q3+1; out[4]=apex+1;
}

/*----------------------------------------------------------------------------*/
/* Main conversion                                                             */
/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::fprintf(stderr,
            "Usage: of2gmf <case_dir> <time_dir> <output_base> [options] [field1 ...]\n"
            "  --all       auto-discover all volScalarField/volVectorField files\n"
            "  -b <ref.meshb>  copy boundary face refs from a reference .meshb\n"
            "  field1 ...  explicit field names (e.g. k omega nut)\n"
            "  (default: p T U)\n"
            "  Reads polyMesh + fields, writes <output_base>.meshb + .solb\n");
        return 1;
    }

    std::string caseDir  = argv[1];
    std::string timeDir  = argv[2];
    std::string outBase  = argv[3];

    std::string meshDir  = caseDir + "/constant/polyMesh";
    std::string fieldDir = caseDir + "/" + timeDir;

    /* ---- Parse optional arguments ---- */
    std::vector<std::string> fieldNames;
    std::string refMeshFile;
    bool doAll = false;
    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--all") {
            doAll = true;
        } else if (arg == "-b" && i + 1 < argc) {
            refMeshFile = argv[++i];
        } else {
            fieldNames.push_back(arg);
        }
    }
    if (doAll) fieldNames = discoverFields(fieldDir);
    if (fieldNames.empty()) fieldNames = {"p", "T", "U"};

    /* ---- Load reference geometry data if provided ---- */
    /* Coordinate-based key for matching vertices across meshes */
    struct CoordKey {
        int64_t ix, iy, iz;
        bool operator<(const CoordKey &o) const {
            if (ix != o.ix) return ix < o.ix;
            if (iy != o.iy) return iy < o.iy;
            return iz < o.iz;
        }
    };
    /* Quantise coordinates for matching.
       Use 1e6 resolution (~1 micron at meter scale) to tolerate the
       ASCII precision loss when coordinates pass through OpenFOAM text
       format.  This is safe since mesh vertices are always well-separated
       relative to this tolerance. */
    auto makeCoordKey = [](double x, double y, double z) -> CoordKey {
        const double S = 1e6;
        return { (int64_t)std::llround(x*S),
                 (int64_t)std::llround(y*S),
                 (int64_t)std::llround(z*S) };
    };

    using FaceKey = std::vector<int>;
    auto makeFaceKey = [](const std::vector<int> &verts) -> FaceKey {
        FaceKey k = verts;
        std::sort(k.begin(), k.end());
        return k;
    };
    using EdgeKey = std::pair<int,int>;
    auto makeEdgeKey = [](int a, int b) -> EdgeKey {
        return a < b ? EdgeKey{a,b} : EdgeKey{b,a};
    };

    /* Data from reference mesh */
    std::map<FaceKey, int> refBndFaceRef;  /* sorted 1-based verts → ref */
    struct RefEdge { int v0, v1, ref; };   /* 1-based ref-mesh vertex IDs */
    std::vector<RefEdge> refEdges;
    std::vector<int> refCorners;           /* 1-based ref-mesh vertex IDs */
    std::vector<int> refRidges;            /* 1-based ref-mesh edge indices (into refEdges) */
    /* Geometry associations: VerticesOnGeometric{Vertices,Edges,Triangles} */
    struct GeomNode { int node, id; };                   /* kw 40 */
    struct GeomEdge { int node, id; double t, gref; };   /* kw 41 */
    struct GeomFace { int node, id; double u, v, gref; };/* kw 42 */
    std::vector<GeomNode> refGeomNodes;
    std::vector<GeomEdge> refGeomEdges;
    std::vector<GeomFace> refGeomFaces;
    /* cad_data (GmfByteFlow) */
    std::vector<unsigned char> refCadData;
    /* Coordinate map: ref-mesh vertex → CoordKey */
    std::map<CoordKey, int> refCoordToIdx; /* coord → 1-based ref vertex */
    std::vector<std::array<double,3>> refPts; /* ref vertex coords (0-based) */
    std::vector<int> refVertexRefs;              /* ref vertex refs (0-based) */

    if (!refMeshFile.empty()) {
        int rver, rdim;
        int64_t rmsh = GmfOpenMesh(refMeshFile.c_str(), GmfRead, &rver, &rdim);
        if (!rmsh) {
            std::fprintf(stderr, "Warning: cannot open ref mesh %s, ignoring -b\n",
                         refMeshFile.c_str());
            refMeshFile.clear();
        } else {
            std::printf("Reading geometry data from %s ...\n", refMeshFile.c_str());

            /* Read reference vertices (for coordinate matching) */
            int64_t nRefVerts = GmfStatKwd(rmsh, GmfVertices);
            if (nRefVerts > 0) {
                refPts.resize(nRefVerts);
                refVertexRefs.resize(nRefVerts);
                GmfGotoKwd(rmsh, GmfVertices);
                for (int64_t i = 0; i < nRefVerts; i++) {
                    GmfGetLin(rmsh, GmfVertices,
                              &refPts[i][0], &refPts[i][1], &refPts[i][2],
                              &refVertexRefs[i]);
                    CoordKey ck = makeCoordKey(refPts[i][0], refPts[i][1], refPts[i][2]);
                    refCoordToIdx[ck] = (int)(i + 1); /* 1-based */
                }
                std::printf("  Read %lld reference vertices\n", (long long)nRefVerts);
            }

            /* Read boundary face refs (triangles + quads) */
            auto readRefBnd = [&](int kwd, int nv) {
                int64_t n = GmfStatKwd(rmsh, kwd);
                if (n == 0) return;
                GmfGotoKwd(rmsh, kwd);
                for (int64_t i = 0; i < n; i++) {
                    int nd[4], ref;
                    if (nv == 3)
                        GmfGetLin(rmsh, kwd, &nd[0],&nd[1],&nd[2],&ref);
                    else
                        GmfGetLin(rmsh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&ref);
                    std::vector<int> verts(nd, nd + nv);
                    std::sort(verts.begin(), verts.end());
                    refBndFaceRef[verts] = ref;
                }
                std::printf("  Read %lld %s from reference mesh\n",
                            (long long)n, nv == 3 ? "triangles" : "quads");
            };
            readRefBnd(GmfTriangles, 3);
            readRefBnd(GmfQuadrilaterals, 4);

            /* Read boundary edges */
            {
                int64_t n = GmfStatKwd(rmsh, GmfEdges);
                if (n > 0) {
                    GmfGotoKwd(rmsh, GmfEdges);
                    refEdges.resize(n);
                    for (int64_t i = 0; i < n; i++) {
                        GmfGetLin(rmsh, GmfEdges,
                                  &refEdges[i].v0, &refEdges[i].v1, &refEdges[i].ref);
                    }
                    std::printf("  Read %lld edges from reference mesh\n", (long long)n);
                }
            }

            /* Read corners */
            {
                int64_t n = GmfStatKwd(rmsh, GmfCorners);
                if (n > 0) {
                    GmfGotoKwd(rmsh, GmfCorners);
                    refCorners.resize(n);
                    for (int64_t i = 0; i < n; i++)
                        GmfGetLin(rmsh, GmfCorners, &refCorners[i]);
                    std::printf("  Read %lld corners from reference mesh\n", (long long)n);
                }
            }

            /* Read ridges */
            {
                int64_t n = GmfStatKwd(rmsh, GmfRidges);
                if (n > 0) {
                    GmfGotoKwd(rmsh, GmfRidges);
                    refRidges.resize(n);
                    for (int64_t i = 0; i < n; i++)
                        GmfGetLin(rmsh, GmfRidges, &refRidges[i]);
                    std::printf("  Read %lld ridges from reference mesh\n", (long long)n);
                }
            }

            /* Read geometry associations (VerticesOnGeometric*) */
            {
                int64_t n = GmfStatKwd(rmsh, GmfVerticesOnGeometricVertices);
                if (n > 0) {
                    refGeomNodes.resize(n);
                    GmfGotoKwd(rmsh, GmfVerticesOnGeometricVertices);
                    for (int64_t i = 0; i < n; i++)
                        GmfGetLin(rmsh, GmfVerticesOnGeometricVertices,
                                  &refGeomNodes[i].node, &refGeomNodes[i].id);
                    std::printf("  Read %lld VerticesOnGeometricVertices\n", (long long)n);
                }
            }
            {
                int64_t n = GmfStatKwd(rmsh, GmfVerticesOnGeometricEdges);
                if (n > 0) {
                    refGeomEdges.resize(n);
                    GmfGotoKwd(rmsh, GmfVerticesOnGeometricEdges);
                    for (int64_t i = 0; i < n; i++)
                        GmfGetLin(rmsh, GmfVerticesOnGeometricEdges,
                                  &refGeomEdges[i].node, &refGeomEdges[i].id,
                                  &refGeomEdges[i].t, &refGeomEdges[i].gref);
                    std::printf("  Read %lld VerticesOnGeometricEdges\n", (long long)n);
                }
            }
            {
                int64_t n = GmfStatKwd(rmsh, GmfVerticesOnGeometricTriangles);
                if (n > 0) {
                    refGeomFaces.resize(n);
                    GmfGotoKwd(rmsh, GmfVerticesOnGeometricTriangles);
                    for (int64_t i = 0; i < n; i++)
                        GmfGetLin(rmsh, GmfVerticesOnGeometricTriangles,
                                  &refGeomFaces[i].node, &refGeomFaces[i].id,
                                  &refGeomFaces[i].u, &refGeomFaces[i].v,
                                  &refGeomFaces[i].gref);
                    std::printf("  Read %lld VerticesOnGeometricTriangles\n", (long long)n);
                }
            }

            /* Read cad_data using dedicated helper (before closing) */
            {
                int nmbByt = 0;
                char *buf = GmfReadByteFlow(rmsh, &nmbByt);
                if (buf && nmbByt > 0) {
                    refCadData.assign(buf, buf + nmbByt);
                    free(buf);
                    std::printf("  Read %d bytes of cad_data from reference mesh\n",
                                nmbByt);
                }
            }

            GmfCloseMesh(rmsh);
        }
    }

    /* ---- Read polyMesh ---- */
    std::printf("Reading %s ...\n", meshDir.c_str());

    auto pts = readPoints(meshDir + "/points");
    auto faces = readFaces(meshDir + "/faces");

    int nCells = -1, nInternalFaces = -1;
    auto owner     = readLabelList(meshDir + "/owner",     &nCells, &nInternalFaces);
    auto neighbour = readLabelList(meshDir + "/neighbour", nullptr, nullptr);
    auto patches   = readBoundary(meshDir + "/boundary");

    int nFaces  = (int)faces.size();
    int nPoints = (int)pts.size();

    if (nCells < 0) {
        /* Fall back: compute from owner array */
        nCells = *std::max_element(owner.begin(), owner.end()) + 1;
    }
    if (nInternalFaces < 0) {
        nInternalFaces = (int)neighbour.size();
    }

    std::printf("  nPoints=%d  nCells=%d  nFaces=%d  nInternalFaces=%d\n",
                nPoints, nCells, nFaces, nInternalFaces);

    /* ---- Build cell→faces map ---- */
    std::vector<std::vector<int>> cellFaces(nCells);
    for (int f = 0; f < nFaces; f++) {
        cellFaces[owner[f]].push_back(f);
        if (f < nInternalFaces)
            cellFaces[neighbour[f]].push_back(f);
    }

    /* ---- Classify cells and extract connectivity ---- */
    std::vector<CellType>  cellType(nCells, POLY);
    /* local index within each type's list */
    std::vector<int>       cellLocalIdx(nCells, -1);

    struct ElemRecord {
        int nodes[8];
        int ref;
    };
    std::vector<ElemRecord> tets, hexes, prisms, pyramids;
    int nPoly = 0;

    for (int c = 0; c < nCells; c++) {
        const auto &cf = cellFaces[c];
        int nFacesC = (int)cf.size();
        int nTri = 0, nQad = 0;
        for (int f : cf) {
            int nv = (int)faces[f].size();
            if (nv == 3) nTri++;
            else if (nv == 4) nQad++;
        }

        ElemRecord er;
        er.ref = 1;
        memset(er.nodes, 0, sizeof(er.nodes));

        if (nFacesC == 4 && nTri == 4) {
            cellType[c] = TET;
            extractTet(faces, owner, cf, c, er.nodes);
            cellLocalIdx[c] = (int)tets.size();
            tets.push_back(er);
        } else if (nFacesC == 6 && nQad == 6) {
            cellType[c] = HEX;
            extractHex(faces, owner, cf, c, er.nodes);
            cellLocalIdx[c] = (int)hexes.size();
            hexes.push_back(er);
        } else if (nFacesC == 5 && nTri == 2 && nQad == 3) {
            cellType[c] = PRISM;
            extractPrism(faces, owner, cf, c, er.nodes);
            cellLocalIdx[c] = (int)prisms.size();
            prisms.push_back(er);
        } else if (nFacesC == 5 && nTri == 4 && nQad == 1) {
            cellType[c] = PYRAMID;
            extractPyramid(faces, owner, cf, c, er.nodes);
            cellLocalIdx[c] = (int)pyramids.size();
            pyramids.push_back(er);
        } else {
            cellType[c] = POLY;
            nPoly++;
        }
    }

    if (nPoly > 0)
        std::printf("Warning: %d polyhedral cell(s) skipped\n", nPoly);

    std::printf("  tets=%zu  hexes=%zu  prisms=%zu  pyramids=%zu\n",
                tets.size(), hexes.size(), prisms.size(), pyramids.size());

    /* ---- Collect boundary faces ---- */
    struct BndFace {
        std::vector<int> verts; /* 1-based */
        int ref;                /* patch index or EGADS face ID */
    };
    std::vector<BndFace> bndTri, bndQad;
    int nRefMatched = 0, nRefFallback = 0;

    for (int pi = 0; pi < (int)patches.size(); pi++) {
        const Patch &p = patches[pi];
        int patchRef = pi + 1; /* fallback: 1-based patch index */
        for (int fi = p.startFace; fi < p.startFace + p.nFaces; fi++) {
            const auto &fv = faces[fi];
            BndFace bf;
            /* OF stores boundary face CCW from owner (outward), keep as-is */
            bf.verts.reserve(fv.size());
            for (int v : fv) bf.verts.push_back(v + 1); /* 0→1-based */

            /* Try to recover original ref from reference mesh */
            if (!refBndFaceRef.empty()) {
                FaceKey key = makeFaceKey(bf.verts);
                auto it = refBndFaceRef.find(key);
                if (it != refBndFaceRef.end()) {
                    bf.ref = it->second;
                    nRefMatched++;
                } else {
                    bf.ref = patchRef;
                    nRefFallback++;
                }
            } else {
                bf.ref = patchRef;
            }

            if ((int)fv.size() == 3) bndTri.push_back(bf);
            else if ((int)fv.size() == 4) bndQad.push_back(bf);
            /* larger polygon faces: skip */
        }
    }
    if (!refBndFaceRef.empty()) {
        std::printf("  Boundary refs: %d matched from reference, %d fallback\n",
                    nRefMatched, nRefFallback);
    }

    /* ---- Write .meshb ---- */
    std::string meshFile = outBase + ".meshb";
    std::printf("Writing %s ...\n", meshFile.c_str());

    int64_t msh = GmfOpenMesh(meshFile.c_str(), GmfWrite, 2, 3);
    if (!msh) {
        std::fprintf(stderr, "Error: cannot create %s\n", meshFile.c_str());
        return 1;
    }

    /* Vertices — use refs from reference mesh if available.
       Build ref-vertex → output-vertex mapping for geometry data. */
    std::vector<int> vertexRef(nPoints, 0);
    /* ref2outVert[r] = 1-based output vertex index for ref vertex r (1-based) */
    std::vector<int> ref2outVert;
    if (!refMeshFile.empty() && !refPts.empty()) {
        ref2outVert.resize((int)refPts.size() + 1, 0);

        /* If vertex counts match, assume direct 1:1 correspondence
           (gmf2of preserves vertex order, createPatch preserves indices) */
        if ((int)refPts.size() == nPoints) {
            std::printf("  Vertex counts match (%d) — using direct index mapping\n",
                        nPoints);
            for (int i = 0; i < nPoints; i++) {
                ref2outVert[i + 1] = i + 1;
                vertexRef[i] = refVertexRefs[i];
            }
        } else {
            /* Fall back to coordinate matching */
            std::map<CoordKey, int> coordToRef;
            for (int i = 0; i < (int)refPts.size(); i++)
                coordToRef[makeCoordKey(refPts[i][0], refPts[i][1], refPts[i][2])]
                    = i + 1;
            int nMatched = 0;
            for (int i = 0; i < nPoints; i++) {
                auto it = coordToRef.find(
                    makeCoordKey(pts[i][0], pts[i][1], pts[i][2]));
                if (it != coordToRef.end()) {
                    int refIdx = it->second;
                    ref2outVert[refIdx] = i + 1;
                    vertexRef[i] = refVertexRefs[refIdx - 1];
                    nMatched++;
                }
            }
            std::printf("  Vertex refs: %d / %d matched by coordinate\n",
                        nMatched, nPoints);
        }
    }
    GmfSetKwd(msh, GmfVertices, nPoints);
    for (int i = 0; i < nPoints; i++)
        GmfSetLin(msh, GmfVertices, pts[i][0], pts[i][1], pts[i][2], vertexRef[i]);

    /* Volume elements */
    if (!tets.empty()) {
        GmfSetKwd(msh, GmfTetrahedra, (int64_t)tets.size());
        for (auto &e : tets)
            GmfSetLin(msh, GmfTetrahedra,
                      e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3], e.ref);
    }
    if (!hexes.empty()) {
        GmfSetKwd(msh, GmfHexahedra, (int64_t)hexes.size());
        for (auto &e : hexes)
            GmfSetLin(msh, GmfHexahedra,
                      e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3],
                      e.nodes[4], e.nodes[5], e.nodes[6], e.nodes[7], e.ref);
    }
    if (!prisms.empty()) {
        GmfSetKwd(msh, GmfPrisms, (int64_t)prisms.size());
        for (auto &e : prisms)
            GmfSetLin(msh, GmfPrisms,
                      e.nodes[0], e.nodes[1], e.nodes[2],
                      e.nodes[3], e.nodes[4], e.nodes[5], e.ref);
    }
    if (!pyramids.empty()) {
        GmfSetKwd(msh, GmfPyramids, (int64_t)pyramids.size());
        for (auto &e : pyramids)
            GmfSetLin(msh, GmfPyramids,
                      e.nodes[0], e.nodes[1], e.nodes[2],
                      e.nodes[3], e.nodes[4], e.ref);
    }

    /* Boundary elements */
    if (!bndTri.empty()) {
        GmfSetKwd(msh, GmfTriangles, (int64_t)bndTri.size());
        for (auto &bf : bndTri)
            GmfSetLin(msh, GmfTriangles,
                      bf.verts[0], bf.verts[1], bf.verts[2], bf.ref);
    }
    if (!bndQad.empty()) {
        GmfSetKwd(msh, GmfQuadrilaterals, (int64_t)bndQad.size());
        for (auto &bf : bndQad)
            GmfSetLin(msh, GmfQuadrilaterals,
                      bf.verts[0], bf.verts[1], bf.verts[2], bf.verts[3], bf.ref);
    }

    /* ---- Write geometry data from reference mesh ---- */
    if (!refMeshFile.empty() && !ref2outVert.empty()) {
        auto &ref2out = ref2outVert; /* alias for brevity */

        /* Write GmfEdges */
        if (!refEdges.empty()) {
            /* Map edges and filter out unmapped ones */
            struct OutEdge { int v0, v1, ref; };
            std::vector<OutEdge> outEdges;
            for (auto &re : refEdges) {
                int ov0 = ref2out[re.v0], ov1 = ref2out[re.v1];
                if (ov0 > 0 && ov1 > 0)
                    outEdges.push_back({ov0, ov1, re.ref});
            }
            if (!outEdges.empty()) {
                GmfSetKwd(msh, GmfEdges, (int64_t)outEdges.size());
                for (auto &e : outEdges)
                    GmfSetLin(msh, GmfEdges, e.v0, e.v1, e.ref);
                std::printf("  Wrote %zu boundary edges\n", outEdges.size());
            }

            /* Write GmfRidges (ridge = edge index; remap to new edge list) */
            if (!refRidges.empty()) {
                /* Build ref-edge-index → out-edge-index map */
                /* refRidges[i] is 1-based index into the original edge list */
                /* We need to find which output edge corresponds to each ref edge */
                std::map<EdgeKey, int> outEdgeIdx; /* sorted pair → 1-based out index */
                for (int i = 0; i < (int)outEdges.size(); i++)
                    outEdgeIdx[makeEdgeKey(outEdges[i].v0, outEdges[i].v1)] = i + 1;

                std::vector<int> outRidges;
                for (int ri : refRidges) {
                    if (ri < 1 || ri > (int)refEdges.size()) continue;
                    int rv0 = ref2out[refEdges[ri-1].v0];
                    int rv1 = ref2out[refEdges[ri-1].v1];
                    if (rv0 > 0 && rv1 > 0) {
                        auto it = outEdgeIdx.find(makeEdgeKey(rv0, rv1));
                        if (it != outEdgeIdx.end())
                            outRidges.push_back(it->second);
                    }
                }
                if (!outRidges.empty()) {
                    GmfSetKwd(msh, GmfRidges, (int64_t)outRidges.size());
                    for (int r : outRidges)
                        GmfSetLin(msh, GmfRidges, r);
                    std::printf("  Wrote %zu ridges\n", outRidges.size());
                }
            }
        }

        /* Write GmfCorners */
        if (!refCorners.empty()) {
            std::vector<int> outCorners;
            for (int rc : refCorners) {
                int ov = ref2out[rc];
                if (ov > 0) outCorners.push_back(ov);
            }
            if (!outCorners.empty()) {
                GmfSetKwd(msh, GmfCorners, (int64_t)outCorners.size());
                for (int c : outCorners)
                    GmfSetLin(msh, GmfCorners, c);
                std::printf("  Wrote %zu corners\n", outCorners.size());
            }
        }
    }

    /* ---- Write geometry associations from reference mesh ---- */
    if (!refMeshFile.empty() && !ref2outVert.empty()) {
        auto mapRefVert = [&](int refVert) -> int {
            if (refVert < 1 || refVert >= (int)ref2outVert.size()) return 0;
            return ref2outVert[refVert];
        };

        /* VerticesOnGeometricVertices */
        if (!refGeomNodes.empty()) {
            std::vector<GeomNode> out;
            for (auto &g : refGeomNodes) {
                int ov = mapRefVert(g.node);
                if (ov > 0) out.push_back({ov, g.id});
            }
            if (!out.empty()) {
                GmfSetKwd(msh, GmfVerticesOnGeometricVertices, (int64_t)out.size());
                for (auto &g : out)
                    GmfSetLin(msh, GmfVerticesOnGeometricVertices, g.node, g.id);
                std::printf("  Wrote %zu VerticesOnGeometricVertices\n", out.size());
            }
        }
        /* VerticesOnGeometricEdges */
        if (!refGeomEdges.empty()) {
            std::vector<GeomEdge> out;
            for (auto &g : refGeomEdges) {
                int ov = mapRefVert(g.node);
                if (ov > 0) out.push_back({ov, g.id, g.t, g.gref});
            }
            if (!out.empty()) {
                GmfSetKwd(msh, GmfVerticesOnGeometricEdges, (int64_t)out.size());
                for (auto &g : out)
                    GmfSetLin(msh, GmfVerticesOnGeometricEdges,
                              g.node, g.id, g.t, g.gref);
                std::printf("  Wrote %zu VerticesOnGeometricEdges\n", out.size());
            }
        }
        /* VerticesOnGeometricTriangles */
        if (!refGeomFaces.empty()) {
            std::vector<GeomFace> out;
            for (auto &g : refGeomFaces) {
                int ov = mapRefVert(g.node);
                if (ov > 0) out.push_back({ov, g.id, g.u, g.v, g.gref});
            }
            if (!out.empty()) {
                GmfSetKwd(msh, GmfVerticesOnGeometricTriangles, (int64_t)out.size());
                for (auto &g : out)
                    GmfSetLin(msh, GmfVerticesOnGeometricTriangles,
                              g.node, g.id, g.u, g.v, g.gref);
                std::printf("  Wrote %zu VerticesOnGeometricTriangles\n", out.size());
            }
        }
    }

    /* Write cad_data (GmfByteFlow) before closing the mesh.
       Note: cad_data embedding is optional — refine can load the .egads
       file directly via the -g flag instead. */
    if (0 && !refCadData.empty()) {
        if (GmfWriteByteFlow(msh, (char *)refCadData.data(),
                             (int)refCadData.size()))
            std::printf("  Wrote %zu bytes of cad_data\n", refCadData.size());
        else
            std::fprintf(stderr, "  Warning: failed to write cad_data\n");
    }

    GmfCloseMesh(msh);
    std::printf("  Wrote %d vertices, %zu tets, %zu hexes, %zu prisms, %zu pyramids\n",
                nPoints, tets.size(), hexes.size(), prisms.size(), pyramids.size());
    std::printf("  Wrote %zu boundary triangles, %zu boundary quads\n",
                bndTri.size(), bndQad.size());

    /* ---- Read solution fields ---- */
    std::vector<FieldData> fields;
    for (auto &name : fieldNames) {
        auto fd = readField(fieldDir + "/" + name, nCells);
        if (fd.present) {
            std::printf("  Field '%s' (%s)\n", name.c_str(),
                        fd.isVector ? "vector" : "scalar");
            fields.push_back(fd);
        } else {
            std::fprintf(stderr, "Warning: field '%s' not found or not a vol field, skipping\n",
                         name.c_str());
        }
    }

    if (fields.empty()) {
        std::printf("No solution fields found in %s, skipping .solb\n",
                    fieldDir.c_str());
        return 0;
    }

    /* Build solution type descriptor */
    std::vector<int> solTypVec;
    for (auto &fd : fields)
        solTypVec.push_back(fd.isVector ? GmfVec : GmfSca);
    int nSolTyp = (int)solTypVec.size();

    /* Compute total doubles per cell */
    int solSiz = 0;
    for (int t : solTypVec)
        solSiz += (t == GmfSca ? 1 : 3);

    /* Build flat per-cell solution array (indexed by OF cell ID) */
    std::vector<double> solAll(nCells * solSiz, 0.0);
    for (int c = 0; c < nCells; c++) {
        int off = c * solSiz, pos = 0;
        for (auto &fd : fields) {
            if (!fd.isVector) {
                solAll[off + pos++] = fd.values[c];
            } else {
                solAll[off + pos+0] = fd.values[c*3+0];
                solAll[off + pos+1] = fd.values[c*3+1];
                solAll[off + pos+2] = fd.values[c*3+2];
                pos += 3;
            }
        }
    }

    /* Build per-type packed buffers using cellLocalIdx mapping */
    auto buildBuf = [&](const std::vector<ElemRecord> &elems,
                        const std::vector<int> &localIdx,
                        CellType type,
                        std::vector<double> &buf) {
        int n = (int)elems.size();
        if (n == 0) return;
        buf.resize(n * solSiz);
        for (int c = 0; c < nCells; c++) {
            if (cellType[c] != type) continue;
            int li = localIdx[c];
            for (int k = 0; k < solSiz; k++)
                buf[li * solSiz + k] = solAll[c * solSiz + k];
        }
    };

    std::vector<double> tetSol, hexSol, priSol, pyrSol;
    buildBuf(tets,    cellLocalIdx, TET,     tetSol);
    buildBuf(hexes,   cellLocalIdx, HEX,     hexSol);
    buildBuf(prisms,  cellLocalIdx, PRISM,   priSol);
    buildBuf(pyramids,cellLocalIdx, PYRAMID, pyrSol);

    /* ---- Write .solb ---- */
    std::string solFile = outBase + ".solb";
    std::printf("Writing %s ...\n", solFile.c_str());

    int64_t sol = GmfOpenMesh(solFile.c_str(), GmfWrite, 2, 3);
    if (!sol) {
        std::fprintf(stderr, "Error: cannot create %s\n", solFile.c_str());
        return 1;
    }

    int *solTypArr = solTypVec.data();

    auto writeSolKwd = [&](int64_t solMsh, int kwd,
                           std::vector<double> &buf, int64_t nElems) {
        if (nElems == 0) return;
        GmfSetKwd(solMsh, kwd, nElems, nSolTyp, solTypArr);
        /* For SolKwd, GmfSetLin takes a single double* pointing to the
           contiguous per-line array of SolSiz doubles. */
        for (int64_t i = 0; i < nElems; i++)
            GmfSetLin(solMsh, kwd, &buf[i * solSiz]);
    };

    writeSolKwd(sol, GmfSolAtTetrahedra,    tetSol, (int64_t)tets.size());
    writeSolKwd(sol, GmfSolAtHexahedra,     hexSol, (int64_t)hexes.size());
    writeSolKwd(sol, GmfSolAtPrisms,        priSol, (int64_t)prisms.size());
    writeSolKwd(sol, GmfSolAtPyramids,      pyrSol, (int64_t)pyramids.size());

    /* ---- Interpolate cell-centered solution to vertices ---- */
    std::vector<double> vertSol(nPoints * solSiz, 0.0);
    std::vector<int> vertCount(nPoints, 0);

    auto addCellToVerts = [&](const std::vector<ElemRecord> &elems,
                              int nNodes, int cellOff) {
        for (int i = 0; i < (int)elems.size(); i++) {
            int ci = cellOff + i;
            for (int k = 0; k < nNodes; k++) {
                int vi = elems[i].nodes[k] - 1; /* 1-based → 0-based */
                for (int s = 0; s < solSiz; s++)
                    vertSol[vi * solSiz + s] += solAll[ci * solSiz + s];
                vertCount[vi]++;
            }
        }
    };

    int tetOff = 0;
    int hexOff = (int)tets.size();
    int priOff = hexOff + (int)hexes.size();
    int pyrOff = priOff + (int)prisms.size();
    addCellToVerts(tets,     4, tetOff);
    addCellToVerts(hexes,    8, hexOff);
    addCellToVerts(prisms,   6, priOff);
    addCellToVerts(pyramids, 5, pyrOff);

    for (int v = 0; v < nPoints; v++) {
        if (vertCount[v] > 0) {
            for (int s = 0; s < solSiz; s++)
                vertSol[v * solSiz + s] /= vertCount[v];
        }
    }

    GmfSetKwd(sol, GmfSolAtVertices, (int64_t)nPoints, nSolTyp, solTypArr);
    for (int64_t i = 0; i < nPoints; i++)
        GmfSetLin(sol, GmfSolAtVertices, &vertSol[i * solSiz]);

    /* ---- Write field names via ReferenceStrings ---- */
    GmfSetKwd(sol, GmfReferenceStrings, (int64_t)fieldNames.size());
    for (int f = 0; f < (int)fieldNames.size(); f++) {
        char nameBuf[256] = {};
        std::snprintf(nameBuf, sizeof(nameBuf), "%s %d",
                      fieldNames[f].c_str(), f + 1);
        GmfSetLin(sol, GmfReferenceStrings,
                  GmfSolAtVertices, 1, nameBuf);
    }

    GmfCloseMesh(sol);
    std::printf("  Solution: %d field(s), %d doubles/cell\n", nSolTyp, solSiz);
    std::printf("  SolAtVertices: %d vertices (cell-averaged)\n", nPoints);
    std::printf("Done.\n");
    return 0;
}
