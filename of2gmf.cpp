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

static FieldData readField(const std::string &path, bool expectVector, int nCells)
{
    FieldData fd;
    fd.present   = false;
    fd.isVector  = expectVector;
    fd.nCells    = nCells;

    std::ifstream in(path);
    if (!in) return fd; /* file not present */

    fd.present = true;
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

/* For a tet cell: extract 4 nodes in correct GMF orientation.
   face[f0] = base tri [a,b,c], apex = 4th vertex.
   If the cell OWNS f0, the face normal (a->b->c CCW) points outward -> correct.
   If the cell is NEIGHBOUR of f0, face is inverted -> flip b and c. */
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
        /* f0 outward from this cell: normal a->b->c points away from apex */
        out[0]=a+1; out[1]=b+1; out[2]=c+1; out[3]=apex+1;
    } else {
        /* f0 inward: flip */
        out[0]=a+1; out[1]=c+1; out[2]=b+1; out[3]=apex+1;
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
            "Usage: of2gmf <case_dir> <time_dir> <output_base>\n"
            "  Reads polyMesh + fields, writes <output_base>.meshb + .solb\n");
        return 1;
    }

    std::string caseDir  = argv[1];
    std::string timeDir  = argv[2];
    std::string outBase  = argv[3];

    std::string meshDir  = caseDir + "/constant/polyMesh";
    std::string fieldDir = caseDir + "/" + timeDir;

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
        int ref;                /* patch index (1-based) */
    };
    std::vector<BndFace> bndTri, bndQad;

    for (int pi = 0; pi < (int)patches.size(); pi++) {
        const Patch &p = patches[pi];
        int ref = pi + 1; /* 1-based patch index */
        for (int fi = p.startFace; fi < p.startFace + p.nFaces; fi++) {
            const auto &fv = faces[fi];
            BndFace bf;
            bf.ref = ref;
            /* OF stores boundary face CCW from owner (outward), keep as-is */
            bf.verts.reserve(fv.size());
            for (int v : fv) bf.verts.push_back(v + 1); /* 0→1-based */
            if ((int)fv.size() == 3) bndTri.push_back(bf);
            else if ((int)fv.size() == 4) bndQad.push_back(bf);
            /* larger polygon faces: skip */
        }
    }

    /* ---- Write .meshb ---- */
    std::string meshFile = outBase + ".meshb";
    std::printf("Writing %s ...\n", meshFile.c_str());

    int64_t msh = GmfOpenMesh(meshFile.c_str(), GmfWrite, 2, 3);
    if (!msh) {
        std::fprintf(stderr, "Error: cannot create %s\n", meshFile.c_str());
        return 1;
    }

    /* Vertices */
    GmfSetKwd(msh, GmfVertices, nPoints);
    for (int i = 0; i < nPoints; i++)
        GmfSetLin(msh, GmfVertices, pts[i][0], pts[i][1], pts[i][2], 0);

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

    GmfCloseMesh(msh);
    std::printf("  Wrote %d vertices, %zu tets, %zu hexes, %zu prisms, %zu pyramids\n",
                nPoints, tets.size(), hexes.size(), prisms.size(), pyramids.size());
    std::printf("  Wrote %zu boundary triangles, %zu boundary quads\n",
                bndTri.size(), bndQad.size());

    /* ---- Read solution fields ---- */
    auto fp   = readField(fieldDir + "/p", false, nCells);
    auto fT   = readField(fieldDir + "/T", false, nCells);
    auto fU   = readField(fieldDir + "/U", true,  nCells);

    bool hasSol = fp.present || fT.present || fU.present;
    if (!hasSol) {
        std::printf("No solution fields found in %s, skipping .solb\n",
                    fieldDir.c_str());
        return 0;
    }

    /* Build solution type descriptor */
    std::vector<int> solTypVec;
    if (fp.present) solTypVec.push_back(GmfSca);
    if (fT.present) solTypVec.push_back(GmfSca);
    if (fU.present) solTypVec.push_back(GmfVec);
    int nSolTyp = (int)solTypVec.size();

    /* Compute total doubles per cell */
    int solSiz = 0;
    for (int t : solTypVec)
        solSiz += (t == GmfSca ? 1 : 3);

    /* Build flat per-cell solution array (indexed by OF cell ID) */
    std::vector<double> solAll(nCells * solSiz, 0.0);
    for (int c = 0; c < nCells; c++) {
        int off = c * solSiz;
        int pos = 0;
        if (fp.present) { solAll[off + pos] = fp.values[c]; pos++; }
        if (fT.present) { solAll[off + pos] = fT.values[c]; pos++; }
        if (fU.present) {
            solAll[off + pos+0] = fU.values[c*3+0];
            solAll[off + pos+1] = fU.values[c*3+1];
            solAll[off + pos+2] = fU.values[c*3+2];
            pos += 3;
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

    GmfCloseMesh(sol);
    std::printf("  Solution: %d field(s), %d doubles/cell\n", nSolTyp, solSiz);
    std::printf("Done.\n");
    return 0;
}
