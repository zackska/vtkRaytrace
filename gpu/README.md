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
- Args: `grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1 srcAng knife cutoff sgain Lbg]`

| `outmode` | output | units |
|---|---|---|
| 0 | transmittance (**shadowgraph**) — rays deviated beyond `ACCdeg` read dark | 0–1 |
| 1 / 2 | mean deflection ε_x / ε_y | mrad |
| 3 | **classical schlieren** — knife edge at the focal plane | 0–1 |
| 4 / 5 / 6 | **BOS** apparent background displacement: \|d\|, dx, dy | mm |

Extra optional args: `srcAng` diffuse-source half-angle (deg); `knife` 1 = cuts ε_x
(vertical edge), 2 = cuts ε_y; `cutoff` fraction of the source image passing at zero
deflection (0.5 = the usual half-cut); `sgain` = f2/a in 1/rad (schlieren sensitivity);
`Lbg` BOS background distance in grid length units.

## Classical schlieren (`outmode 3`)
A source image of height `a` is cut by a knife so a fraction `cutoff` passes undeflected;
a ray deflected by ε is displaced `f2·ε` at the knife plane, so

```
T = clamp(cutoff + (f2/a)·ε, 0, 1)
```

which reproduces the standard result that schlieren contrast is `ΔI/I0 = f2·ε/a`. The image
therefore responds to the **first** derivative of refractive index along the knife-cut
direction, where shadowgraphy responds to the second. Saturation at 0 or 1 is physical — the
knife fully blocking or fully passing the source image.

## Background-oriented schlieren (`outmode 4/5/6`)
A patterned background a distance `Lbg` behind the test section appears displaced by
`Lbg·tan(ε)`. Rays that cannot reach the background — turned back by TIR at an interface, or
extinguished in opaque liquid — are flagged **invalid in the kernel** and returned as `NaN`,
which is what a real rig sees (the cross-correlation simply fails there). Without that flag,
grazing-incidence rays drive `tan(ε)` toward infinity and produce nonsense.

`outmode 4/5/6` give the **ground-truth** displacement, which no experiment can measure
directly. `bos_correlate.py` closes that gap: it generates a speckle background, warps it by
the displacement field, and recovers the result by windowed FFT cross-correlation with
subpixel Gaussian peak fitting — i.e. the same processing a real BOS rig applies. It
reproduces window-scale smoothing, peak-locking, and correlation dropout, and includes the
**loss-of-pairs correction** (normalising by the window autocorrelation); without that the
recovered displacement is systematically ~8% low at 4–6 px shifts.

### Measurable envelope — read this before trusting a BOS render
A correlation window tracks displacements only between roughly the noise floor (~0.1 px) and
about a quarter of the window (8 px for a 32 px window). On interface-resolved sprays the
deflection distribution spans **~4 orders of magnitude** (median ~0.02 mrad in the plume,
p99 ~90 mrad at the interface), so **no single standoff `Lbg` puts the whole field in range**:

| standoff | outcome on an evaporating jet |
|---|---|
| 300 mm (typical gas-phase BOS) | interface displacements reach 10²–10³ px — unrecoverable |
| 20 mm | interfaces measurable, but the plume displaces ~0.02 px — below the noise floor |

BOS is well suited to the smooth gas-phase gradients it was designed for, and structurally
ill-suited to sprays with resolved interfaces. Choose `Lbg` for the part of the flow you
actually want to measure, and state it.

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

`validate_schlieren_bos.py` covers the newer outmodes against closed-form answers on a linear
index ramp, where ε = (dn/dx)·L exactly: deflection and BOS displacement both to ~1.6%
(ray-marching discretisation, and BOS inherits it 1:1 as it should), and schlieren
transmission matching `clamp(cutoff + (f2/a)·ε, 0, 1)` across sensitivities from 0 to
20000 /rad — linear then correctly saturating. The correlation emulator is checked by
recovering known uniform shifts: 0–8 px to better than 0.15 px.

### Target and camera model

`bos_correlate.py` generates an idealised Gaussian speckle background. Two further scripts
model what a real rig records, so you can ask how the optics and the target — rather than the
flow — limit the measurement.

`bos_realistic_target.py` builds a **printed** target (hard-edged dots, not Gaussian) and
puts both exposures through a camera: optical PSF, shot noise, read noise and quantisation,
with **independent** noise on each frame. A real rig cannot reuse one exposure, and that
matters: measured on a helium-jet field, Gaussian speckle falls to ~81 % valid vectors at
0.32 px RMS under a realistic camera, while a printed target holds ~99 % at 0.27 px — hard
ink edges carry more high-frequency content, so the correlation peak stands clear of the
noise floor. Defocus dominates everything else; light level and bit depth barely register.

`bos_window_sweep.py` sweeps interrogation window size against dot density. The usual
"32 px window" advice is only half the story: at 7.6 dots/window the same 32 px window gave
3.7 px RMS with a +0.9 px bias, and at 38 dots/window it gave 0.115 px. Across 8–64 px
windows the **density mattered more than the window size**, so if you need finer windows the
fix is more dots, not smaller boxes.

Both take a deflection field as an `.npz` containing `dx`, `dy` (apparent background
displacement, mm) and `mm_per_px` — i.e. what `gpuShadow` outmodes 5/6 produce:

```bash
python3 gpu/bos_realistic_target.py field.npz
python3 gpu/bos_window_sweep.py     field.npz
```

Numbers quoted above came from a single jet at one instant; treat them as illustrative of the
trends rather than as calibration constants.
