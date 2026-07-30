# References

Provenance for the methods implemented here. Every DOI below was resolved against Crossref;
where a DOI could not be verified it is omitted rather than guessed.

## Refractive index

**Gladstone, J. H. & Dale, T. P.** (1863) *Researches on the refraction, dispersion, and
sensitiveness of liquids.* Philosophical Transactions of the Royal Society of London **153**,
317–343. [10.1098/rstl.1863.0014](https://doi.org/10.1098/rstl.1863.0014)

The linear relation `n − 1 = K·ρ` used throughout, with `K` composition-weighted for mixtures.
Note that `docs/ARCHITECTURE.md` specifies **Lorentz–Lorenz** with per-species refractivity;
Gladstone–Dale is what is implemented, and the gap is recorded in the
*Implemented vs specified* table.

## Ray propagation in a graded-index medium

**Sharma, A., Kumar, D. V. & Ghatak, A. K.** (1982) *Tracing rays through graded-index media:
a new method.* Applied Optics **21**(6), 984–987.
[10.1364/AO.21.000984](https://doi.org/10.1364/AO.21.000984)

The RK4 formulation integrated by `marchRay()` in `mode eikonal`. Snell refraction with total
internal reflection (`mode sharp`) and Beer–Lambert attenuation are textbook and not
separately cited.

## Schlieren and shadowgraphy

**Settles, G. S.** (2001) *Schlieren and Shadowgraph Techniques: Visualizing Phenomena in
Transparent Media.* Springer.
[10.1007/978-3-642-56640-0](https://doi.org/10.1007/978-3-642-56640-0)

Source for the knife-edge relation implemented as `outmode 3`,
`T = clamp(cutoff + (f₂/a)·ε, 0, 1)`, and for the distinction that schlieren responds to the
first derivative of refractive index while shadowgraphy responds to the second.

## Background-oriented schlieren

**Raffel, M.** (2015) *Background-oriented schlieren (BOS) techniques.* Experiments in Fluids
**56**, 60. [10.1007/s00348-015-1927-5](https://doi.org/10.1007/s00348-015-1927-5)

Source for `d = L_bg · tan ε`, and for the focusing argument behind `gpu/optics_sweep.py`: a
real rig focuses on the *background*, so the density object is deliberately defocused and the
combined blur adds in quadrature, `d_Σ = √(d_d² + d_i²)`.

## Correlation processing

**Adrian, R. J.** (2005) *Twenty years of particle image velocimetry.* Experiments in Fluids
**39**, 159–169. [10.1007/s00348-005-0991-7](https://doi.org/10.1007/s00348-005-0991-7)

Background for the windowed FFT cross-correlation, subpixel Gaussian peak fitting,
peak-locking and loss-of-pairs behaviour emulated in `gpu/bos_correlate.py`. The loss-of-pairs
correction (normalising by the window autocorrelation) matters: without it the recovered
displacement is systematically ~8% low at 4–6 px shifts.

## Singular-value coherence diagnostic

**Zamani Ashtiani, S. & Fukami, K.** (2026) *Data-Driven Time-Dependent Bases for Turbulent
Airfoil Wake–Extreme Gust Interactions.* AIAA Journal.
[10.2514/1.J067033](https://doi.org/10.2514/1.J067033)

Source of the relative singular-value gap `δ(t) = (σ₁ − σ₂)/σ₁` used in
`analysis/sigma_gap.py`. What is implemented is the paper's *diagnostic*, not its algorithm —
their time-dependent bases solve evolution equations to avoid storing the history of large 3-D
fields, which is unnecessary at image sizes where an instantaneous SVD costs microseconds. See
`analysis/README.md` for the two ways δ did **not** transfer to shadowgraphy.

## Experimental dataset

**Liang, C., Johansen, C. T. & Linne, M. A.** (2022) *Characteristics of sprays produced by
coaxial non-swirling and swirling air-blast atomizers.* Physics of Fluids **34**, 093606.
[10.1063/5.0107480](https://doi.org/10.1063/5.0107480)

The high-speed shadowgraph dataset `analysis/sigma_gap.py` was developed against, including
the FWI / SWI / Bag / Fiber regime classification. The data are open at Edinburgh DataShare,
[10.7488/ds/3459](https://doi.org/10.7488/ds/3459).

## Toolkits

- **VTK** — the CPU tracer is built on it. <https://vtk.org/>
- **CUDA** — `gpu/gpuShadow.cu`; build with `nvcc -O3 -arch=sm_86 -cudart static`.
