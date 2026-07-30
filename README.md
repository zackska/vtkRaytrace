# vtkRaytrace

A C++ ray tracer for simulating diffuse back-lit illumination through transmissive objects, built on the [VTK](https://vtk.org/) library. Originally developed to model the **diffuse light imaging utility** used to characterise high-pressure fuel sprays during PhD research at Chalmers.

The tracer reads a triangulated surface mesh (STL), converts it to a set of density iso-contours, builds an OBB tree for fast intersection, and shoots rays through the geometry to produce a synthetic shadowgram / transmission image.

**See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for how the pieces fit together** — the
two tracers, the physics chain from field to detector, the data contracts, and where fidelity
is won and lost.

## Files

### CPU tracer (`ReadSTL`)

| File | Purpose |
|---|---|
| `main.cpp` | Entry point — loads the input, configures the tracer, runs it |
| `vtkRaytrace.cpp/.h` | Core class: mesh setup, ray casting, intersection logic |
| `meshTools.cpp/.h` | VTK helpers for mesh manipulation and visualisation |
| `geometry.h` | Vector / geometric primitive utilities |
| `dirent.h` | POSIX `dirent` polyfill for Windows builds |
| `CMakeLists.txt` | CMake build script (produces `ReadSTL` executable) |

### GPU renderer (`gpu/`)

| File | Purpose |
|---|---|
| `gpu/gpuShadow.cu` | CUDA kernel — one thread per pixel, `sharp`/`eikonal`/`hybrid` ray march, 7 output modes |
| `gpu/foam_shadowgraph.py` | OpenFOAM case → `grid.bin` (`scatter`/`resample`, `presmooth`) |
| `gpu/bos_correlate.py` | Speckle warp + windowed FFT cross-correlation with loss-of-pairs correction |
| `gpu/bos_realistic_target.py` | Printed target + camera model (PSF, shot/read noise, quantisation) |
| `gpu/bos_window_sweep.py` | Interrogation window size × dot density sweep |
| `gpu/bos_design_envelope.py` | **Rig sizing before building** — deflection distribution → measurable-fraction vs standoff and window |
| `gpu/optics_sweep.py` | Finite-aperture / defocus sweep (f/32 → f/2.8), blur traced not post-applied |
| `gpu/bos_velocimetry.py` | Whether BOS yields velocity — correlates deflection fields between instants |
| `gpu/correction_transfer.py` | Experiment-only vs CFD-informed bias correction, scored on held-out frames |
| `gpu/gpu_validate.py` | Eikonal/hybrid validation against closed-form cases |
| `gpu/validate_schlieren_bos.py` | Schlieren + BOS validation on an analytic index ramp |
| `gpu/README.md` | Modes, arguments, output modes, and the measurable-envelope caveats |

### Study harnesses (`analysis/`)

Not library code — each carries assumptions from the study it was written for. See
[`analysis/README.md`](analysis/README.md).

| File | Purpose |
|---|---|
| `analysis/ladder_strip.py` | 5-panel common-crop shadowgraph strip across a Weber-number ladder |
| `analysis/sigma_gap.py` | Instantaneous singular-value diagnostics (δ, r99, spectral entropy) on shadowgraph sets |

## Build

Requires CMake ≥ 3.12 and **VTK 9**. `CMakeLists.txt` requests only the modules actually used
rather than calling a bare `find_package(VTK)` — a bare call pulls in every module present,
which breaks when building against ParaView's bundled VTK. Classes are factory-registered via
`vtk_module_autoinit`, which VTK ≥ 9 requires for the IO and Rendering factories to
initialise. OpenMP is used if found.

```bash
mkdir build && cd build
cmake ..
cmake --build .
./ReadSTL [input.stl] [nPixels] [absorptionCoeff]
```

Arguments are optional and positional: input path, square image resolution (default 400), and
absorption coefficient in 1/length (default 80). With no arguments a hard-coded path in
`main.cpp` is used — edit `inputFilename` or pass `argv[1]`.

Behaviour is further set by environment variable:

| Variable | Effect |
|---|---|
| `VTKRT_VOLUME=1` | read a single ASCII `.vts` structured grid (cell array `rho_s`) instead of an STL surface |
| `VTKRT_INTERFACE` | `sharp` = one refracting surface at the α=0.5 iso-level; `diffuse` (default) = nested iso-index shells |
| `VTKRT_NCONTOURS` | shell count in the `diffuse` case (default 10) |
| `VTKRT_ISOVALUE` | iso-level in `rho_s` units for the `sharp` case (default: field midpoint) |
| `VTKRT_OUT` | output filename, so concurrent processes don't collide |
| `VTKRT_FRAME_BOUNDS` | `"xmin xmax ymin ymax zmin zmax"` — fixes the bounding box so an image *sequence* shares one framing |

## GPU renderer

`gpu/` contains **gpuShadow**, a CUDA reimplementation that supersedes the VTK tracer for
volumetric fields: one thread per pixel, rays marched through a 3D texture with hardware
trilinear sampling, ~300 Mrays/s (a frame in 3-5 ms against ~15-30 s here). It adds
**classical schlieren** and **background-oriented schlieren (BOS)** alongside shadowgraphy,
and ships analytic validation for each. See `gpu/README.md`.

The two are complementary rather than redundant: this VTK tracer works from a triangulated
surface (STL), gpuShadow from a sampled volume.

## References

Provenance for every implemented method — Gladstone–Dale, the Sharma RK4 graded-index ray
trace, Settles for the knife-edge relation, Raffel for BOS and the defocus geometry, Adrian for
the correlation processing, and the datasets — is in
[`docs/REFERENCES.md`](docs/REFERENCES.md), with DOIs resolved against Crossref.

## Citing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21701176.svg)](https://doi.org/10.5281/zenodo.21701176)

Archived on Zenodo. There are **two DOIs** and the distinction matters:

| DOI | Meaning |
|---|---|
| [`10.5281/zenodo.21701176`](https://doi.org/10.5281/zenodo.21701176) | **Concept** — always resolves to the latest release. **Cite this** unless you need to pin a version. |
| [`10.5281/zenodo.21701177`](https://doi.org/10.5281/zenodo.21701177) | **Version** — v0.1.0 exactly. Use when a paper's results must be reproducible against specific code. |

```bibtex
@software{falgout_cfd_optical_diagnostics,
  author    = {Falgout, Zachary},
  title     = {{cfd-optical-diagnostics}: synthetic shadowgraph, schlieren
               and BOS rendering from CFD fields},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21701176},
  url       = {https://doi.org/10.5281/zenodo.21701176}
}
```

`CITATION.cff` carries the same metadata machine-readably — GitHub renders a *Cite this
repository* button from it.

If you use a specific method rather than the software as a whole, cite the original work from
[`docs/REFERENCES.md`](docs/REFERENCES.md) as well.

## License

MIT - see `LICENSE`.
