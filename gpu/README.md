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

## Validation
Eikonal deflection matches the analytic linear-gradient case to 0.19% and a closed-form
Gaussian phase object to 7e-4; hybrid reduces exactly to sharp (gas off) and to the gas
bend (liquid off). See `gpu_validate.py`.
