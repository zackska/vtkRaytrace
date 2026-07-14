# gpuShadow — GPU synthetic shadowgraph renderer

CUDA port of the acceptance-angle shadowgraph tracer. One kernel, one thread per pixel;
rays march through a 3D texture of the field, using the GPU's hardware trilinear
interpolation for both the field value and its gradient. No BVH, no triangles, no cell
locator. **~300 Mrays/s** — a frame in 3–5 ms vs ~15–30 s for the CPU (VTK) tracer.

## Modes
- **sharp**   — threshold field (alpha) at 0.5 → Snell/TIR refraction + Beer–Lambert absorption (resolved two-phase interface).
- **eikonal** — continuous RK4 (Sharma) ray march through `n = 1 + K·field` (graded-index gas/thermal field; schlieren/BOS).
- **hybrid**  — gas eikonal bend (field1) + Snell interface (field0 = alpha) + absorption (evaporating jets).
- Finite depth-of-field via thin-lens aperture sampling (`naper` rays/pixel, focal plane `zf`, radius `apR`; `naper=1` = infinite DoF).

## Build
```
nvcc -O3 -arch=sm_86 -cudart static -o gpuShadow gpuShadow.cu
```
`-cudart static` → the binary runs on any machine with an NVIDIA driver (no CUDA toolkit install needed).

## I/O
- `grid.bin`: `int32 NX,NY,NZ` ; `float X0,X1,Y0,Y1,Z0,Z1` ; `NX*NY*NZ float32` (k,j,i order).
- `img.bin` : `int32 RESX,RESY` ; `RESX*RESY float32`.
- Args: `grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1]`
  - `outmode` 0 = transmittance (shadowgraph), 1/2 = deflection_x/y in mrad (schlieren).

## Feeding an OpenFOAM field (`foam_shadowgraph.py`)
Helper to go from a reconstructed OpenFOAM case → `grid.bin` → `gpuShadow`. Two fidelity knobs:

- **`method`** — how the mesh field lands on the uniform grid:
  - `scatter` (default): bin cells by centre. Order-independent, and it **keeps detached
    droplets** — each liquid cell's value survives.
  - `resample` (`vtkResampleToImage`): smooth interpolation, but a 1–2-cell droplet averages
    with surrounding gas and drops **below α=0.5 → it vanishes**. Use for a clean coherent-core
    figure; use `scatter` for atomization/breakup fidelity.
  - Both are fast at *native* grid resolution. Upsampling the grid is slow and adds no real
    detail — render the **image** at high `res` instead (the kernel trilinear-upsamples the
    volume for free), so image resolution is decoupled from grid cost.

- **`presmooth`** (Gaussian σ, cells) — the kernel locates the α=0.5 surface *implicitly*
  (trilinear threshold-crossing) with the normal from the local gradient — **not** a
  marching-cubes surface with averaged normals like the CPU tracer. A light σ≈1 restores that
  **smoothed iso-surface + smoothed normals**; larger σ smooths small droplets away (→ the
  `resample` look).

- **schlieren scaling** — `outmode` 1/2 give deflection (mrad). Map the magnitude through
  `log1p(mag/A)` rather than a linear clip so a strong local bloom doesn't saturate and its
  internal structure still reads.

## Validation
Eikonal deflection matches the analytic linear-gradient case to 0.19% and a closed-form
Gaussian phase object to 7e-4; hybrid reduces exactly to sharp (gas off) and to the gas
bend (liquid off). See `gpu_validate.py`.
