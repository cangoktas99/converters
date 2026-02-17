/*----------------------------------------------------------------------------*/
/*  gmf2of.cpp                                                                */
/*  Convert a libMeshb binary mesh (.meshb) + optional solution (.solb) to   */
/*  an OpenFOAM polyMesh case directory with fields p, T, U.                 */
/*                                                                            */
/*  Usage:                                                                    */
/*    gmf2of <input_base> <output_case_dir> [time_dir=0]                     */
/*                                                                            */
/*  Example:                                                                  */
/*    gmf2of result ./myCase 0                                                */
/*    -> reads  result.meshb (+ result.solb if present)                      */
/*    -> writes ./myCase/constant/polyMesh/{points,faces,owner,neighbour,     */
/*               boundary}                                                    */
/*       ./myCase/0/{p,T,U}  (only the fields present in .solb)              */
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
#include <stdexcept>
#include <functional>
#include <sys/stat.h>
#include <sys/types.h>

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
/*                                                                            */
/* Tet    [n0,n1,n2,n3]:  face n0,n1,n2 has outward normal (away from n3)   */
/* Hex    [n0..n7]:        bottom=n0..n3, top=n4..n7 (n4 above n0, etc.)     */
/* Prism  [n0..n5]:        bottom tri=n0,n1,n2; top tri=n3,n4,n5             */
/* Pyramid[n0..n4]:        quad base=n0,n1,n2,n3; apex=n4                    */
/*----------------------------------------------------------------------------*/

/* Returns list of faces for the element; each face = {local_node_indices}
   in the CCW-from-owner orientation for the FIRST cell that owns that face.
   Owner is assigned at first-seen; second-seen → neighbour (face is then CW). */

static std::vector<std::vector<int>> tetFaceTemplates()
{
    /* For tet [n0,n1,n2,n3]:
       Face opposite n3 (base): n0,n1,n2 (CCW when viewed from outside, i.e. from away-from-n3)
       Other faces:              each includes n3 */
    return {
        {0,1,2},   /* base, normal away from n3 */
        {0,3,1},   /* left,  CCW from outside */
        {1,3,2},   /* right, CCW from outside */
        {2,3,0},   /* back,  CCW from outside */
    };
}

static std::vector<std::vector<int>> hexFaceTemplates()
{
    /* Hex: bottom=0,1,2,3  top=4,5,6,7  (n4 above n0) */
    return {
        {0,3,2,1}, /* bottom face, CCW from below (outward = downward) */
        {4,5,6,7}, /* top face,    CCW from above (outward = upward)   */
        {0,1,5,4}, /* front  */
        {1,2,6,5}, /* right  */
        {2,3,7,6}, /* back   */
        {3,0,4,7}, /* left   */
    };
}

static std::vector<std::vector<int>> prismFaceTemplates()
{
    /* Prism: bottom tri=0,1,2  top tri=3,4,5 */
    return {
        {0,2,1},   /* bottom, CCW from below */
        {3,4,5},   /* top,    CCW from above */
        {0,1,4,3}, /* side 0-1 */
        {1,2,5,4}, /* side 1-2 */
        {2,0,3,5}, /* side 2-0 */
    };
}

static std::vector<std::vector<int>> pyramidFaceTemplates()
{
    /* Pyramid: quad base=0,1,2,3  apex=4 */
    return {
        {0,3,2,1}, /* base, CCW from below */
        {0,1,4},   /* front tri  */
        {1,2,4},   /* right tri  */
        {2,3,4},   /* back  tri  */
        {3,0,4},   /* left  tri  */
    };
}

/*----------------------------------------------------------------------------*/
/* Main                                                                       */
/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::fprintf(stderr,
            "Usage: gmf2of <input_base> <output_case_dir> [time_dir=0]\n");
        return 1;
    }

    std::string inBase   = argv[1];
    std::string caseDir  = argv[2];
    std::string timeDir  = (argc >= 4) ? argv[3] : "0";

    std::string meshFile = inBase + ".meshb";
    std::string solFile  = inBase + ".solb";
    std::string meshDir  = caseDir + "/constant/polyMesh";
    std::string fieldDir = caseDir + "/" + timeDir;

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

    /* ---- Read vertices (0-based: subtract 1 from GMF 1-based) ---- */
    std::vector<std::array<double,3>> pts(nVer);
    {
        GmfGotoKwd(msh, GmfVertices);
        for (int64_t i = 0; i < nVer; i++) {
            int ref;
            GmfGetLin(msh, GmfVertices, &pts[i][0], &pts[i][1], &pts[i][2], &ref);
        }
    }

    /* ---- Read element connectivity (store as 0-based) ---- */
    struct ElemData {
        int nodes[8];
        int ref;
    };

    auto readElems = [&](int kwd, int nNodes, int64_t n,
                         std::vector<ElemData> &out) {
        out.resize(n);
        if (n == 0) return;
        GmfGotoKwd(msh, kwd);
        for (int64_t i = 0; i < n; i++) {
            int &r = out[i].ref;
            int *nd = out[i].nodes;
            switch (nNodes) {
                case 4:
                    GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&r);
                    break;
                case 5:
                    GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&nd[4],&r);
                    break;
                case 6:
                    GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],&nd[4],&nd[5],&r);
                    break;
                case 8:
                    GmfGetLin(msh, kwd, &nd[0],&nd[1],&nd[2],&nd[3],
                                        &nd[4],&nd[5],&nd[6],&nd[7],&r);
                    break;
            }
            /* Convert to 0-based */
            for (int k = 0; k < nNodes; k++) nd[k]--;
        }
    };

    std::vector<ElemData> tets, hexes, prisms, pyrs;
    readElems(GmfTetrahedra, 4, nTet, tets);
    readElems(GmfHexahedra,  8, nHex, hexes);
    readElems(GmfPrisms,     6, nPri, prisms);
    readElems(GmfPyramids,   5, nPyr, pyrs);

    /* Boundary surface elements (for patch assignment via ref) */
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
                for (int &v : bs.verts) v--;  /* 0-based */
                bs.ref = r;
                bndSurfs.push_back(bs);
            }
        };
        readBndElems(GmfTriangles,      3, nBTri);
        readBndElems(GmfQuadrilaterals, 4, nBQad);
    }

    /* Build lookup: sorted face key → patch ref */
    std::map<FaceKey,int> bndFaceRef;
    for (auto &bs : bndSurfs)
        bndFaceRef[makeFaceKey(bs.verts)] = bs.ref;

    GmfCloseMesh(msh);

    /* ---- Cell ordering: tets first, then hexes, prisms, pyrs ---- */
    /* OF cell ID = tet_offset + local_idx, etc. */
    int64_t tetOff = 0;
    int64_t hexOff = nTet;
    int64_t priOff = nTet + nHex;
    int64_t pyrOff = nTet + nHex + nPri;

    /* ---- Build face list via deduplication ---- */
    std::map<FaceKey, FaceRecord> faceMap;

    /* Helper: add one element's faces.
       cellId (OF 0-based), nodes (0-based), templates */
    auto addElemFaces = [&](int64_t cellId,
                            const int *nodes,
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
                /* Check if this is a boundary surface element */
                auto bit = bndFaceRef.find(key);
                if (bit != bndFaceRef.end())
                    fr.patchRef = bit->second;
                faceMap[key] = fr;
            } else {
                /* Second cell → neighbour */
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

    /* ---- Separate internal from boundary faces ---- */
    std::vector<FaceRecord> internalFaces, boundaryFaces;
    for (auto &kv : faceMap) {
        if (kv.second.neighbour >= 0)
            internalFaces.push_back(kv.second);
        else
            boundaryFaces.push_back(kv.second);
    }

    /* Sort internal faces: by owner asc, then neighbour asc */
    std::sort(internalFaces.begin(), internalFaces.end(),
        [](const FaceRecord &a, const FaceRecord &b) {
            if (a.owner != b.owner) return a.owner < b.owner;
            return a.neighbour < b.neighbour;
        });

    /* Sort boundary faces: by patchRef, then owner */
    std::sort(boundaryFaces.begin(), boundaryFaces.end(),
        [](const FaceRecord &a, const FaceRecord &b) {
            if (a.patchRef != b.patchRef) return a.patchRef < b.patchRef;
            return a.owner < b.owner;
        });

    /* Collect unique patch refs from boundary faces */
    std::map<int,int> refStart; /* patchRef → startFace in boundary section */
    std::map<int,int> refCount;
    int bndBase = (int)internalFaces.size();
    {
        int idx = bndBase;
        for (auto &bf : boundaryFaces) {
            if (!refCount.count(bf.patchRef)) {
                refStart[bf.patchRef] = idx;
                refCount[bf.patchRef] = 0;
            }
            refCount[bf.patchRef]++;
            idx++;
        }
    }

    /* Full ordered face list */
    std::vector<FaceRecord> allFaces;
    allFaces.insert(allFaces.end(), internalFaces.begin(), internalFaces.end());
    allFaces.insert(allFaces.end(), boundaryFaces.begin(), boundaryFaces.end());

    int nAllFaces        = (int)allFaces.size();
    int nInternalFaces   = (int)internalFaces.size();

    std::printf("  nAllFaces=%d  nInternalFaces=%d  nBndFaces=%d\n",
                nAllFaces, nInternalFaces, (int)boundaryFaces.size());

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
        out << (int)refStart.size() << "\n(\n";
        for (auto &kv : refStart) {
            int ref = kv.first;
            out << "    Patch_" << ref << "\n    {\n"
                << "        type            patch;\n"
                << "        nFaces          " << refCount[ref] << ";\n"
                << "        startFace       " << kv.second << ";\n"
                << "    }\n";
        }
        out << ")\n";
    }
    std::printf("  Wrote boundary (%d patches)\n", (int)refStart.size());

    /* ---- Read solution (.solb) ---- */
    int64_t sol = GmfOpenMesh(solFile.c_str(), GmfRead, &ver, &dim);
    if (!sol) {
        std::printf("No .solb found at %s, skipping field output.\n", solFile.c_str());
        std::printf("Done.\n");
        return 0;
    }
    std::printf("Reading %s ...\n", solFile.c_str());

    /* Detect which cell-type solution keywords are present */
    int NmbTyp = 0, SolSiz = 0;
    int SolTypArr[GmfMaxTyp] = {};

    /* Try to get solution metadata from whichever keyword is present */
    int64_t nSolCheck = 0;
    for (int kwd : {GmfSolAtTetrahedra, GmfSolAtHexahedra,
                    GmfSolAtPrisms, GmfSolAtPyramids}) {
        int nt=0, ss=0, st[GmfMaxTyp]={};
        int64_t n = GmfStatKwd(sol, kwd, &nt, &ss, st);
        if (n > 0 && nSolCheck == 0) {
            NmbTyp = nt; SolSiz = ss;
            memcpy(SolTypArr, st, nt * sizeof(int));
            nSolCheck = n;
        }
    }

    if (NmbTyp == 0) {
        std::printf("No recognisable solution fields in .solb, skipping.\n");
        GmfCloseMesh(sol);
        std::printf("Done.\n");
        return 0;
    }

    /* Determine field layout: p, T, U by type order */
    /* We assume the convention from of2gmf:
       SolTypArr = {GmfSca, GmfSca, GmfVec} → offsets 0, 1, 2
       or any subset thereof.                                     */
    int nScaFound = 0, vecOffset = -1;
    std::vector<int> scaOffsets;
    {
        int off = 0;
        for (int t = 0; t < NmbTyp; t++) {
            if (SolTypArr[t] == GmfSca) {
                scaOffsets.push_back(off);
                nScaFound++;
                off++;
            } else if (SolTypArr[t] == GmfVec) {
                vecOffset = off;
                off += 3;
            } else {
                off += SolSiz; /* skip unknown */
            }
        }
    }

    bool hasP = nScaFound >= 1;
    bool hasT = nScaFound >= 2;
    bool hasU = vecOffset >= 0;

    /* Build per-type solution buffers */
    auto readSolBuf = [&](int kwd, int64_t n) -> std::vector<double> {
        std::vector<double> buf(n * SolSiz, 0.0);
        if (n == 0) return buf;
        GmfGotoKwd(sol, kwd);
        /* For SolKwd, GmfGetLin takes a single double* pointing to the
           contiguous per-line array of SolSiz doubles. */
        for (int64_t i = 0; i < n; i++)
            GmfGetLin(sol, kwd, &buf[i * SolSiz]);
        return buf;
    };

    auto tetSol = readSolBuf(GmfSolAtTetrahedra, nTet);
    auto hexSol = readSolBuf(GmfSolAtHexahedra,  nHex);
    auto priSol = readSolBuf(GmfSolAtPrisms,     nPri);
    auto pyrSol = readSolBuf(GmfSolAtPyramids,   nPyr);
    GmfCloseMesh(sol);

    /* Reconstruct flat cell array [0..nCells-1] */
    std::vector<double> flatSol(nCells * SolSiz);
    {
        auto copy = [&](const std::vector<double> &src, int64_t offset, int64_t n) {
            for (int64_t i = 0; i < n; i++)
                for (int k = 0; k < SolSiz; k++)
                    flatSol[(offset+i)*SolSiz + k] = src[i*SolSiz + k];
        };
        copy(tetSol, tetOff, nTet);
        copy(hexSol, hexOff, nHex);
        copy(priSol, priOff, nPri);
        copy(pyrSol, pyrOff, nPyr);
    }

    /* ---- Write field files ---- */
    auto writeScalarField = [&](const std::string &name,
                                const std::string &dims,
                                int solOffset) {
        std::ofstream out(fieldDir + "/" + name);
        writeFoamHeader(out, "volScalarField", timeDir, name);
        out << "dimensions      " << dims << ";\n\n"
            << "internalField   nonuniform List<scalar>\n"
            << nCells << "\n(\n";
        out << std::scientific; out.precision(10);
        for (int64_t c = 0; c < nCells; c++)
            out << flatSol[c*SolSiz + solOffset] << "\n";
        out << ");\n\n"
            << "boundaryField\n{\n"
            << "    \".*\"\n    {\n"
            << "        type            zeroGradient;\n"
            << "    }\n}\n";
    };

    auto writeVectorField = [&](const std::string &name,
                                const std::string &dims,
                                int solOffset) {
        std::ofstream out(fieldDir + "/" + name);
        writeFoamHeader(out, "volVectorField", timeDir, name);
        out << "dimensions      " << dims << ";\n\n"
            << "internalField   nonuniform List<vector>\n"
            << nCells << "\n(\n";
        out << std::scientific; out.precision(10);
        for (int64_t c = 0; c < nCells; c++) {
            const double *row = &flatSol[c*SolSiz];
            out << "(" << row[solOffset+0] << " "
                       << row[solOffset+1] << " "
                       << row[solOffset+2] << ")\n";
        }
        out << ");\n\n"
            << "boundaryField\n{\n"
            << "    \".*\"\n    {\n"
            << "        type            zeroGradient;\n"
            << "    }\n}\n";
    };

    if (hasP) {
        writeScalarField("p", "[1 -1 -2 0 0 0 0]", scaOffsets[0]);
        std::printf("  Wrote p\n");
    }
    if (hasT) {
        writeScalarField("T", "[0 0 0 1 0 0 0]", scaOffsets[1]);
        std::printf("  Wrote T\n");
    }
    if (hasU) {
        writeVectorField("U", "[0 1 -1 0 0 0 0]", vecOffset);
        std::printf("  Wrote U\n");
    }

    std::printf("Done.\n");
    return 0;
}
