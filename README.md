# vtkRaytrace

A C++ ray tracer for simulating diffuse back-lit illumination through transmissive objects, built on the [VTK](https://vtk.org/) library. Originally developed to model the **diffuse light imaging utility** used to characterise high-pressure fuel sprays during PhD research at Chalmers.

The tracer reads a triangulated surface mesh (STL), converts it to a set of density iso-contours, builds an OBB tree for fast intersection, and shoots rays through the geometry to produce a synthetic shadowgram / transmission image.

## Files

| File | Purpose |
|---|---|
| `main.cpp` | Entry point — loads an STL, configures the tracer, runs it |
| `vtkRaytrace.cpp/.h` | Core class: mesh setup, ray casting, intersection logic |
| `meshTools.cpp/.h` | VTK helpers for mesh manipulation and visualisation |
| `geometry.h` | Vector / geometric primitive utilities |
| `dirent.h` | POSIX `dirent` polyfill for Windows builds |
| `CMakeLists.txt` | CMake build script (produces `ReadSTL` executable) |

## Build

Requires CMake ≥ 3.0 and a VTK installation (tested with VTK 7–8).

```bash
mkdir build && cd build
cmake ..
cmake --build .
./ReadSTL
```

The STL path is currently hard-coded in `main.cpp` (look for `inputFilename`) — edit it before running.

## GPU renderer

`gpu/` contains **gpuShadow**, a CUDA reimplementation that supersedes the VTK tracer for
volumetric fields: one thread per pixel, rays marched through a 3D texture with hardware
trilinear sampling, ~300 Mrays/s (a frame in 3-5 ms against ~15-30 s here). It adds
**classical schlieren** and **background-oriented schlieren (BOS)** alongside shadowgraphy,
and ships analytic validation for each. See `gpu/README.md`.

The two are complementary rather than redundant: this VTK tracer works from a triangulated
surface (STL), gpuShadow from a sampled volume.

## License

MIT - see `LICENSE`.
