# OpenFOAM ↔ GMF Mesh & Solution Converters

## Overview

This project provides two C++ command-line converters between the
[OpenFOAM](https://www.openfoam.com/) polyMesh format and the
[Gamma Mesh Format (GMF)](https://github.com/LoicMarechal/libMeshb)
handled by the **libMeshb** library.  Both directions are supported:

| Program   | Direction                  |
|-----------|----------------------------|
| `of2gmf`  | OpenFOAM → `.meshb`/`.solb` |
| `gmf2of`  | `.meshb`/`.solb` → OpenFOAM |

Converted quantities: mesh topology (vertices + standard cell types) and
solution fields **p** (pressure), **T** (temperature), **U** (velocity).

---

## Part I — Understanding libMeshb and the GMF Format

### What is libMeshb?

**libMeshb** (version 7.80 used here) is a lightweight, single-file C library
written by Loïc Maréchal (INRIA) for reading and writing mesh and solution files
in the *Gamma Mesh Format*.  The library is used throughout the CFD and mesh
generation community — notably by tools such as *Feflo.a*, *Wolf*, and *Gmsh*.

The format comes in two flavours for each file type:

| Extension | Content  | Encoding       |
|-----------|----------|----------------|
| `.mesh`   | Mesh     | ASCII text     |
| `.meshb`  | Mesh     | Binary         |
| `.sol`    | Solution | ASCII text     |
| `.solb`   | Solution | Binary         |

Detection is purely by filename extension
([`libmeshb7.c:583–592`](../libMeshb-7.80/sources/libmeshb7.c)).
The library transparently handles little-endian and big-endian binary files and
supports optional asynchronous I/O for very large meshes.

---

### File Versions

Four versions (1–4) are supported.  The version number controls numeric
precision and pointer width:

| Version | Float size | Integer size | File offset size |
|---------|-----------|-------------|-----------------|
| 1       | 32-bit    | 32-bit      | 32-bit          |
| 2       | 64-bit    | 32-bit      | 32-bit          |
| 3       | 64-bit    | 32-bit      | 64-bit          |
| 4       | 64-bit    | 64-bit      | 64-bit          |

The converters in this project always write **version 2** (double-precision
coordinates, 32-bit indices) — the most widely supported choice.

---

### ASCII File Layout (`.mesh` / `.sol`)

An ASCII mesh file is a plain text file containing named keyword blocks
separated by blank lines:

```
MeshVersionFormatted 2

Dimension 3

Vertices
1982
0.491814 0.357873 0.793755 0
-0.185247 0.576968 0.795482 0
...

Triangles
3960
1 2 3 1
2 4 3 1
...

End
```

Lines beginning with `#` are comments.  Keyword names are case-sensitive
strings matched against the internal keyword table.

---

### Binary File Layout (`.meshb` / `.solb`)

The binary format is a *forward-linked list* of keyword blocks, which allows
O(1) random access to any keyword:

```
[4 bytes]   Endian tag  (1 = native, 16777216 = byte-swapped)
[4 bytes]   Version     (1–4)
[4 bytes]   Keyword code for Dimension
[4/8 bytes] Offset to next keyword block
[4 bytes]   Dimension value (2 or 3)

For each subsequent keyword block:
  [4 bytes]   Keyword code  (enum GmfKwdCod)
  [4/8 bytes] Offset to the next block   ← linked-list pointer
  [4/8 bytes] Number of data lines
  ... data ...
```

The offset pointers mean the library can seek directly to any keyword at
open time without scanning the entire file.

---

### The Keyword Catalogue

The library defines **282 keywords** via the `GmfKwdCod` enum in
[`libmeshb7.h`](../libMeshb-7.80/sources/libmeshb7.h).  Every keyword has
three associated strings in the `GmfKwdFmt` table
([`libmeshb7.c:236–466`](../libMeshb-7.80/sources/libmeshb7.c)):

| Column | Meaning                          | Example         |
|--------|----------------------------------|-----------------|
| `[0]`  | ASCII keyword name               | `"Vertices"`    |
| `[1]`  | Count-field specifier            | `"i"` or `""`  |
| `[2]`  | Per-line data format string      | `"dri"`         |

#### Format string characters

| Char | Meaning                                                  |
|------|----------------------------------------------------------|
| `i`  | Integer (32-bit for ver ≤ 3, 64-bit for ver = 4)         |
| `r`  | Real (32-bit float for ver = 1, 64-bit double otherwise) |
| `d`  | Dimension-expanded real: `dr` → `dim` reals              |
| `c`  | 64-byte string                                           |
| `n`  | Variable-length integer list                             |
| `s`  | Solution-expanded real: `sr` → total solution size reals |
| `h`  | High-order solution real: `hr` → `NmbNod × SolSiz` reals |

A keyword with no count field (`[1] == ""`) is an *information keyword*
(`InfKwd`) storing a single scalar value, e.g. `Dimension` or `Time`.
A keyword whose format is `"sr"` or `"hr"` is a *solution keyword*
(`SolKwd`) with an additional header describing the field composition.

---

### Standard Mesh Elements

The most important mesh keywords, with their per-line format strings:

| Keyword            | Format   | Description                  |
|--------------------|----------|------------------------------|
| `Vertices`         | `dri`    | (x, y[, z], ref)             |
| `Edges`            | `iii`    | 2 vertex indices + ref        |
| `Triangles`        | `iiii`   | 3 indices + ref               |
| `Quadrilaterals`   | `iiiii`  | 4 indices + ref               |
| `Tetrahedra`       | `iiiii`  | 4 indices + ref               |
| `Prisms`           | `iiiiiii`| 6 indices + ref               |
| `Hexahedra`        | `iiiiiiiii`| 8 indices + ref             |
| `Pyramids`         | `iiiiii` | 5 indices + ref               |

> **Important:** GMF uses **1-based** vertex indices throughout.

The `ref` (reference) integer at the end of each line carries domain or
boundary-region tags.  For volume elements it is a subdomain ID; for surface
elements it identifies the geometric boundary.

High-order variants (`TrianglesP2`, `TetrahedraP2`, `HexahedraQ2`, ...) follow
the same pattern with additional mid-edge and mid-face node indices appended.

---

### Solution Files and Field Types

A solution keyword — such as `SolAtVertices` or `SolAtTetrahedra` — has an
*extended header* that describes how many field components are stored per
entity, and what type each component is:

| Code       | Name        | Size per entity (dim = 3) |
|------------|-------------|--------------------------|
| `GmfSca`   | Scalar      | 1                        |
| `GmfVec`   | Vector      | 3                        |
| `GmfSymMat`| Sym. matrix | 6                        |
| `GmfMat`   | Full matrix | 9                        |

**ASCII example** (`out.sol` — 5 vertices, 3 fields: scalar + vector + scalar):

```
MeshVersionFormatted 2
Dimension 3

SolAtVertices
5
3 1 2 1        ← 3 field types: GmfSca(1), GmfVec(2), GmfSca(1)

5 1 6 8 4      ← vertex 1: p=5, U=(1,6,8), T=4
8 2 4 8 2
...

End
```

Solution keywords are available for every element type:
`SolAtVertices`, `SolAtEdges`, `SolAtTriangles`, `SolAtQuadrilaterals`,
`SolAtTetrahedra`, `SolAtPrisms`, `SolAtHexahedra`, `SolAtPyramids`, and
their integer variants (`ISolAt*`) and high-order variants (`HOSolAt*`).

---

### Core API

```c
/* Open / close */
int64_t GmfOpenMesh(const char *file, GmfRead/GmfWrite, &ver, &dim);
int     GmfCloseMesh(int64_t idx);

/* Inspect a keyword */
int64_t GmfStatKwd(int64_t idx, int kwd, ...);
  // For SolKwd: extra args &NmbTyp, &TotSolSiz, SolTypTab[]

/* Navigate to a keyword's data */
int GmfGotoKwd(int64_t idx, int kwd);

/* Declare a keyword block for writing */
int GmfSetKwd(int64_t idx, int kwd, int64_t NmbLin, ...);
  // For SolKwd: extra args NmbTyp, int *SolTypTab

/* Read / write one line at a time */
int GmfGetLin(int64_t idx, int kwd, ...);   // variadic pointers per field
int GmfSetLin(int64_t idx, int kwd, ...);   // variadic values per field
  // EXCEPTION for SolKwd: single double* to contiguous SolSiz-element array

/* Bulk block I/O (fastest path) */
int GmfGetBlock(int64_t idx, int kwd, int64_t beg, int64_t end,
                int nproc, void *cbk_usr, void *cbk_fn, ...);
int GmfSetBlock(int64_t idx, int kwd, int64_t beg, int64_t end,
                int nproc, void *cbk_usr, void *cbk_fn, ...);
```

> **Pitfall discovered during development:** for `SolKwd` keywords, both
> `GmfGetLin` and `GmfSetLin` consume only **one variadic argument** — a
> `double *` pointing to a contiguous array of `SolSiz` doubles.  Passing
> individual `double` values as separate variadic arguments silently
> misinterprets their bit patterns as pointer values, producing corrupted
> output.  Always pass a pointer:
> ```cpp
> // Correct
> GmfSetLin(sol, GmfSolAtTetrahedra, &buf[i * solSiz]);
> GmfGetLin(sol, GmfSolAtTetrahedra, &buf[i * solSiz]);
>
> // Wrong — passes doubles as if they were pointers
> GmfSetLin(sol, GmfSolAtTetrahedra, buf[0], buf[1], buf[2], ...);
> ```

---

## Part II — The Converters

### Project Structure

```
converters/
├── CMakeLists.txt       Build system
├── of2gmf.cpp           OpenFOAM → GMF converter
├── gmf2of.cpp           GMF → OpenFOAM converter
└── README.md            This file
```

The library source lives at `../libMeshb-7.80/sources/` and is compiled
directly into both executables — no pre-built libMeshb installation needed.

---

### Build

```bash
cd converters
cmake -B build
cmake --build build
# Binaries: build/of2gmf  build/gmf2of
```

---

### of2gmf — OpenFOAM → GMF

#### Usage

```
of2gmf <case_dir> <time_dir> <output_base>
```

**Example:**
```bash
of2gmf ./myCase 1000 result
# Writes: result.meshb  result.solb
```

#### What it reads

| Path                                      | Content             |
|-------------------------------------------|---------------------|
| `<case_dir>/constant/polyMesh/points`     | Vertex coordinates  |
| `<case_dir>/constant/polyMesh/faces`      | Face connectivity   |
| `<case_dir>/constant/polyMesh/owner`      | Face owner cells    |
| `<case_dir>/constant/polyMesh/neighbour`  | Internal face neighbours |
| `<case_dir>/constant/polyMesh/boundary`   | Boundary patches    |
| `<case_dir>/<time_dir>/p`                 | Pressure (scalar)   |
| `<case_dir>/<time_dir>/T`                 | Temperature (scalar)|
| `<case_dir>/<time_dir>/U`                 | Velocity (vector)   |

Missing field files are silently skipped.  If no fields are found at all,
no `.solb` is written.

#### Cell type detection

OpenFOAM stores cells implicitly as a set of faces.  The converter
reconstructs cell types by counting the triangular and quadrilateral faces
belonging to each cell:

| Type    | Total faces | Triangles | Quads |
|---------|-------------|-----------|-------|
| Tet     | 4           | 4         | 0     |
| Hex     | 6           | 0         | 6     |
| Prism   | 5           | 2         | 3     |
| Pyramid | 5           | 4         | 1     |

Cells that do not match any of these patterns are polyhedral and are
**skipped with a warning**.

#### Node ordering and orientation

GMF requires element nodes in a specific orientation (outward face normals
using the right-hand rule).  OpenFOAM stores face vertices CCW when viewed
from the *owner* cell — i.e., the face normal points away from the owner.

For each standard type the converter:

- **Tet:** takes one triangular face as the base; if the cell *owns* that face
  the winding is already outward, otherwise it is flipped (swap two vertices).
- **Hex:** locates the two opposite quad faces (no shared vertices), builds the
  bottom face from the owned quad, then matches each bottom vertex to its top
  counterpart via shared side faces.
- **Prism:** finds the two triangular end-caps and matches them vertex-by-vertex
  through the three quad side faces.
- **Pyramid:** takes the single quad as the base with the same ownership-based
  flip logic, and the unique fifth vertex as the apex.

#### Solution layout

Fields are written as `SolAtTetrahedra`, `SolAtHexahedra`, `SolAtPrisms`, and
`SolAtPyramids` — one keyword per element type present in the mesh, each
sharing the same field descriptor.

The descriptor for the full `{p, T, U}` set is:
```c
int SolTyp[] = { GmfSca, GmfSca, GmfVec };  // 1 + 1 + 3 = 5 doubles/cell
```
If only a subset of the fields is present, the descriptor is shortened
accordingly.

---

### gmf2of — GMF → OpenFOAM

#### Usage

```
gmf2of <input_base> <output_case_dir> [time_dir=0]
```

**Example:**
```bash
gmf2of result ./myCase 0
# Writes: ./myCase/constant/polyMesh/{points,faces,owner,neighbour,boundary}
#         ./myCase/0/{p,T,U}
```

#### Mesh reconstruction

The converter reads all volume element keywords (`GmfTetrahedra`,
`GmfHexahedra`, `GmfPrisms`, `GmfPyramids`) and decomposes each element
into its faces using fixed orientation templates:

| Element | Face templates (0-based local node indices)                         |
|---------|---------------------------------------------------------------------|
| Tet     | `{0,1,2}` `{0,3,1}` `{1,3,2}` `{2,3,0}`                          |
| Hex     | `{0,3,2,1}` `{4,5,6,7}` `{0,1,5,4}` `{1,2,6,5}` `{2,3,7,6}` `{3,0,4,7}` |
| Prism   | `{0,2,1}` `{3,4,5}` `{0,1,4,3}` `{1,2,5,4}` `{2,0,3,5}`         |
| Pyramid | `{0,3,2,1}` `{0,1,4}` `{1,2,4}` `{2,3,4}` `{3,0,4}`             |

Faces are deduplicated by their sorted vertex key: the first element to
generate a face becomes its *owner*; the second becomes its *neighbour*.
Faces seen only once are boundary faces.

Cell ordering in the output is fixed: all tets first, then hexes, prisms,
pyramids.  This determines how the flat OF cell-index array maps back to
per-type solution buffers.

#### Boundary patch naming

GMF boundary surface elements (`GmfTriangles`, `GmfQuadrilaterals`) carry a
`ref` integer tag.  Each distinct ref value becomes one OpenFOAM boundary
patch named `Patch_<ref>` with type `patch`.

#### Solution output

Each field file is written as a `nonuniform` list with `zeroGradient`
boundary conditions on all patches.  The solution layout assumed when
reading the `.solb` mirrors what `of2gmf` writes:
- first scalar → `p`
- second scalar → `T`
- first vector → `U`

---

### Supported / Unsupported features

| Feature                         | Status           |
|---------------------------------|------------------|
| ASCII `.mesh` / `.sol` input    | Not supported (binary `.meshb`/`.solb` only) |
| Tet / Hex / Prism / Pyramid     | ✓ Full support   |
| General polyhedra (OF)          | Skipped with warning |
| OpenFOAM `uniform` fields       | ✓ Expanded to constant array |
| OpenFOAM `nonuniform` fields    | ✓ Full read      |
| OpenFOAM `$var` substitution    | Not supported (reads as 0) |
| Binary OpenFOAM format          | Not supported    |
| Parallel decomposed cases       | Not supported    |

---

### Verified Test Case

The converters were tested on the **LS59 turbine blade subsonic case** located at
`/home/ubuntu/Templates/LS59_Turbine_Blade_Subsonic_NS/`:

| Quantity              | Value   |
|-----------------------|---------|
| Vertices              | 93,546  |
| Hexahedral cells      | 46,206  |
| Faces (total)         | 185,391 |
| Internal faces        | 91,845  |
| Boundary patches      | 6       |
| U field (nonuniform)  | 46,206 vectors |

Round-trip result (`of2gmf` → `gmf2of`):
- All vertex coordinates preserved to **full double precision**
- Face counts, owner/neighbour arrays and boundary patch structure
  **exactly reconstructed**
- U-field values preserved to **full double precision**
- Total runtime: **< 0.3 seconds** for the full conversion
