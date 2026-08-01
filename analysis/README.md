# analysis/

Study harnesses, not library code. These consume the renderers in `gpu/` (or plain image
data) to answer one specific question, and each carries assumptions from the study it was
written for. Read the header before reusing.

| Script | What it does | Study-specific assumptions |
|---|---|---|
| `ladder_strip.py` | Renders a 5-panel shadowgraph comparison strip across a Weber-number ladder, cropped to a *common* box so panels are like-for-like | `CASES` is hard-coded to `We{16,18,20,22,25}_S0_L2`; expects each case staged as a directory containing `system/blockMeshDict`, `case.foam` and reconstructed time dirs holding `alpha.water` |
| `density_plane.py` | Centre-plane density and species fraction from a **decomposed** OpenFOAM case — reads `processor*/` directly via `vtkPOpenFOAMReader` CaseType 0, so no `reconstructPar` and no serial `polyMesh` needed (worth ~40% of the download when fetching one time directory from object storage) |
| `sigma_gap.py` | Instantaneous singular-value diagnostics on a shadowgraph image set — the relative gap δ = (σ₁−σ₂)/σ₁, the rank for 99% spectral energy, and normalised spectral entropy | Expects the `frames_real.npz` layout (`X`, `meta`, `H`, `W`) built from the Liang (2022) coaxial-swirl dataset; the `AIR2WE` map is calibrated for **S=0 only** |
| `amr_to_grid.py` | Library, not a harness. Converts an OpenFOAM AMR cell field to the uniform grid `gpuShadow` requires, by splatting cells into their voxels instead of probing sample points into the mesh | Assumes a Cartesian octree — a uniform hex block under `topoChanger/refiner`. Not valid for a body-fitted or graded mesh |

```bash
python3 analysis/ladder_strip.py  /path/to/staged/cases     # writes ladder_strip.png there
python3 analysis/sigma_gap.py     /path/to/frames_real.npz  # writes *_sigma_gap.npy alongside
python3 analysis/density_plane.py /path/to/case 0.15        # writes density_plane.png/.npz there
```

Both resolve the `gpuShadow` binary from `$GPUSHADOW_DIR`, falling back to `../gpu`.

## A caveat on `sigma_gap.py`

The metric comes from Zamani Ashtiani & Fukami, *AIAA J* 2026 ([10.2514/1.J067033]),
where δ ranks flow coherence for extreme gust–airfoil interactions. Two things did **not**
carry over when applied to shadowgraphy, and the script reports both variants so you can see
it:

- **δ on a raw frame is useless** — σ₁ is dominated by the static background and the nozzle,
  so δ saturates near 0.95 regardless of regime. Subtract a per-condition mean frame first.
- **δ alone is a weak discriminator even then.** On the Liang S=0 series it is non-monotonic
  in Weber number. What separates regimes is the *whole spectrum*: `r99` rose 13.4 → 14.3 →
  23.9 → 32.3 → 32.5 across We 9/16/25/64/100, and the SWI→Bag boundary at We 16→25 separated
  at AUC 0.977 per frame where the FWI→SWI step gave only 0.756.

So this is the paper's *diagnostic idea*, not its algorithm — full time-dependent bases solve
evolution equations to avoid storing history of large 3-D fields, which is unnecessary at
image sizes where an instantaneous SVD costs microseconds. Frames used were downsampled to
100×80, so `r99` is resolution-limited and the apparent saturation above We≈64 may be an
artefact.

## Getting the field onto a uniform grid

`gpuShadow` is a CUDA texture ray-marcher, so it can only read a dense uniform lattice; an
AMR mesh has no texture path. The conversion is unavoidable, but it is easy to make it the
dominant cost. `vtkResampleToImage` probes every sample point into the mesh, so it scales as
O(samples × cell-locate) and the locate term grows with refinement. On one coaxial-atomizer
case the conversion went **51 s → 686 s** between `t = 0` and `t = 0.1015 s` purely because
the spray dispersed and the octree refined, while rendering the resulting grid stayed at ~4 s.

`amr_to_grid.py` splats instead: every cell of a Cartesian octree is an axis-aligned box on a
known lattice, so it maps to its voxels by arithmetic with no search. Same case, same instant:
**3.2 s** — and roughly flat in refinement rather than growing with it.

Two practical notes. The mesh is *not* all-hex from VTK's point of view — AMR leaves hanging
nodes, so cells report 8–26 points and per-cell boxes must come from a segmented min/max over
each cell's point list. And `vtkCellCenters`/`vtkCellSizeFilter` are not a shortcut; they walk
cells through the generic API and took over 8 minutes on 2.8 M cells where numpy on the raw
arrays takes ~2 s.

The two paths are not equivalent, deliberately. Probing averages cell data onto vertices and
then interpolates trilinearly, smoothing a VoF interface twice and dilating the liquid — on
the case above it reported ~12 % more dark pixels than the solver actually holds. Splatting
keeps a PLIC interface as sharp as the solver left it. Every structure survives either way
(mean |difference| over the rendered image 0.0089, confined to structure outlines).

Default target spacing is the finest cell size, so the lattice aligns with the octree.
Sampling finer does not add information — it interpolates detail into existence.
