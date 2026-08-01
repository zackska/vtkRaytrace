"""Convert an OpenFOAM AMR field to the uniform grid gpuShadow needs, by splatting cells.

WHY A CONVERSION IS NEEDED AT ALL. `gpuShadow` is a CUDA texture ray-marcher: the field is
bound to a `cudaTextureObject_t` over a `cudaMalloc3DArray` and sampled with `tex3D<float>`,
which is where its throughput comes from (hardware trilinear fetch, one instruction per
sample). CUDA textures can only bind a dense uniform lattice, so an unstructured AMR mesh has
no texture path. The conversion is unavoidable; the cost of it is not.

WHY NOT `vtkResampleToImage`. It probes each sample point into the mesh, so it costs
O(samples x cell-locate) and the locate term grows as AMR refines. Measured on a 5.8 cm
coaxial-atomizer case (`We25_S0_L2`), converting one instant:

    t = 0.0000   barely refined            51 s
    t = 0.0250   jet intact               96 s
    t = 0.1015   2.8 M cells, spray dispersed   686 s

i.e. the cost tracks refinement, not time, and it dominates everything else -- rendering the
resulting grid takes ~4 s.

WHAT THIS DOES INSTEAD. An OpenFOAM `refiner` mesh is a Cartesian octree: a uniform hex block
with `maxRefinement N`, so every cell is an axis-aligned box whose edge is the base size
halved N times, aligned to the base lattice. Each cell therefore maps to a known index box in
the target grid by arithmetic, with no search. That makes the conversion O(cells) and roughly
*flat* in refinement level. Same case, same instant: **3.2 s**, a factor of 214.

TWO THINGS THAT BITE, both learned the hard way:

1. The mesh is **not** all-hex as far as VTK is concerned. AMR leaves hanging nodes on the
   faces of a coarse cell bordering refined neighbours, so cells report 8-26 points. The
   cells are still axis-aligned boxes, so the per-cell box comes from a segmented min/max
   over each cell's point list (`np.minimum.reduceat`), not from assuming 8 corners.
2. `vtkCellCenters` and `vtkCellSizeFilter` are not a shortcut -- they walk cells one at a
   time through the generic cell API and took over 8 minutes on 2.8 M cells. The same numbers
   come out of the raw point/connectivity arrays in ~2 s.

RESOLUTION. The target spacing defaults to the finest cell size, so the lattice aligns with
the octree and nothing is sampled below the resolution the data carries. Sampling finer than
the finest cell does not add information; it interpolates detail into existence.

THIS IS NOT EQUIVALENT TO PROBING, deliberately. The probe path runs `cellDataToPointData`
(averaging cell values onto vertices) and then interpolates trilinearly, smoothing the
interface twice. On a VoF field that dilates the liquid: on the case above it reported a
dark-pixel fraction of 0.0445 against 0.0392 for the splat, about 12 % more liquid than the
solver actually holds. Splatting writes each cell value into the voxels that cell occupies,
so a PLIC interface stays as sharp as the solver left it. Every structure survives either way
(mean |difference| over the rendered image was 0.0089, confined to structure outlines).
"""
import struct

import numpy as np
from vtk.util.numpy_support import vtk_to_numpy


def cell_boxes(ug):
    """Per-cell centre and edge length, from the raw arrays.

    Returns (centres (n,3), edge (n,)). Handles the variable point-count that hanging nodes
    produce, and never touches VTK's per-cell API.
    """
    ncell = ug.GetNumberOfCells()
    pts = vtk_to_numpy(ug.GetPoints().GetData())
    ca = ug.GetCells()
    conn = vtk_to_numpy(ca.GetConnectivityArray())
    seg = vtk_to_numpy(ca.GetOffsetsArray())[:-1]
    lo = np.empty((ncell, 3))
    hi = np.empty((ncell, 3))
    for ax in range(3):
        A = pts[conn, ax]
        lo[:, ax] = np.minimum.reduceat(A, seg)
        hi[:, ax] = np.maximum.reduceat(A, seg)
    return 0.5 * (lo + hi), (hi[:, 0] - lo[:, 0]).astype(np.float64)


def splat(ug, field, bounds, dx=None):
    """Rasterise a cell field onto a uniform lattice.

    ug      : vtkUnstructuredGrid carrying the cell array
    field   : name of the cell array, e.g. 'alpha.water'
    bounds  : (LX, LY, LZ) domain extent, origin assumed at 0
    dx      : target spacing; defaults to the finest cell edge, which aligns the lattice with
              the octree. Passing anything smaller oversamples the data.

    Returns (grid (nz,ny,nx) float32, dx).
    """
    ctr, h = cell_boxes(ug)
    vals = np.nan_to_num(vtk_to_numpy(ug.GetCellData().GetArray(field))).astype(np.float32)
    if dx is None:
        dx = h.min()

    LX, LY, LZ = bounds
    NX = int(round(LX / dx))
    NY = int(round(LY / dx))
    NZ = int(round(LZ / dx))
    grid = np.zeros((NZ, NY, NX), np.float32)

    # cells come in a few discrete sizes (the base and its refinements); each size class is
    # handled with pure broadcasting -- no per-cell Python loop, no search.
    lev = np.round(np.log2(h / dx)).astype(int)
    for L in np.unique(lev):
        m = lev == L
        n = int(2 ** max(L, 0))                       # this cell spans n^3 voxels
        c = ctr[m]
        v = vals[m]
        i0 = np.floor((c[:, 0] - 0.5 * n * dx) / dx + 0.5).astype(np.int32)
        j0 = np.floor((c[:, 1] - 0.5 * n * dx) / dx + 0.5).astype(np.int32)
        k0 = np.floor((c[:, 2] - 0.5 * n * dx) / dx + 0.5).astype(np.int32)
        for di in range(n):
            for dj in range(n):
                for dk in range(n):
                    grid[np.clip(k0 + dk, 0, NZ - 1),
                         np.clip(j0 + dj, 0, NY - 1),
                         np.clip(i0 + di, 0, NX - 1)] = v
    return grid, dx


def write_grid(path, grid, bounds):
    """Write gpuShadow's grid format: 3i dims, 6f bounds, float32 payload in (z,y,x)."""
    LX, LY, LZ = bounds
    NZ, NY, NX = grid.shape
    with open(path, 'wb') as f:
        f.write(struct.pack('3i', NX, NY, NZ))
        f.write(struct.pack('6f', 0, LX, 0, LY, 0, LZ))
        grid.tofile(f)
