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

**Keane, R. D. & Adrian, R. J.** (1990) *Optimization of particle image velocimeters. I.
Double pulsed systems.* Measurement Science and Technology **1**, 1202.
[10.1088/0957-0233/1/11/013](https://doi.org/10.1088/0957-0233/1/11/013)

**Scharnowski, S., Sciacchitano, A. & Kähler, C. J.** (2019) *On the universality of Keane &
Adrian's valid detection probability in PIV.* Measurement Science and Technology **30**.
[10.1088/1361-6501/aafe9d](https://doi.org/10.1088/1361-6501/aafe9d)

The valid-detection criterion — roughly six effective particle images per interrogation
window, conditional on five further non-dimensional parameters including in-plane displacement
below about a quarter of the window. `gpu/bos_window_sweep.py` **reproduces** this rather than
extending it; the README records which of those conditions the helium-jet field violates and
why a "density beats window size" reading of the sweep is not supportable.

## Setup design: sensitivity and geometric blur

**Schmidt, B. E., Bathel, B. F., Grauer, S. J., Hargather, M. J., Heineck, J. T. & Raffel, M.**
(2025) *Twenty-Five Years of Background-Oriented Schlieren: Advances and Novel Applications.*
AIAA Journal **63**(12), 5028–5058.
[10.2514/1.J065669](https://doi.org/10.2514/1.J065669) —
author copy free at [NTRS 20240014988](https://ntrs.nasa.gov/citations/20240014988)

Section IV is the basis of `gpu/bos_setup_optimise.py`: the sensitivity factor
`S = f z_D/(z_D + z_A − f)` (Eqs 22–23), the field-of-view constraint on focal length
(Eq 24), the circle of confusion `CoC = f² z_D/(f# z_A (z_A + z_D − f))` (Eq 25) and its
*linear* proportionality to `S`, the result that sensitivity is maximised at
`z_A/z_B = 0.5`, and the ~50 % of-CoC minimum resolvable feature size (after Schwarz &
Braukmann). Note the review's literature cut-off is late 2024.

**Rajendran, L. K., Bane, S. P. M. & Vlachos, P. P.** (2019) *Dot tracking methodology for
background-oriented schlieren (BOS).* Experiments in Fluids **60**(11).
[10.1007/s00348-019-2793-3](https://doi.org/10.1007/s00348-019-2793-3)

The prior art for synthetic BOS rendering as a setup-design aid — the review cites it as the
ray-tracing tool for predicting geometric blur from setup parameters, and calls such methods
"computationally intensive". Its successor from the same group is MIRAGE
([10.1088/1361-6501/ae4f0a](https://doi.org/10.1088/1361-6501/ae4f0a), 2026), which postdates
the review.

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

**Liang, Y., Johansen, L. C. & Linne, M.** (2022) *Breakup of a laminar liquid jet by coaxial
non-swirling and swirling air streams.* Physics of Fluids **34**(9), 093606.
[10.1063/5.0100456](https://doi.org/10.1063/5.0100456)

The high-speed shadowgraph dataset `analysis/sigma_gap.py` was developed against, including
the FWI / SWI / Bag / Fiber regime classification and the non-dimensional breakup lengths.
The data are open at Edinburgh DataShare,
[10.7488/ds/3459](https://doi.org/10.7488/ds/3459).

Atomizer geometry, quoted because downstream comparisons depend on it: liquid tube
**4 mm** inner / 5 mm outer diameter, gas tube **10 mm** inner / 14 mm outer, liquid Reynolds
number held at **480** for every flow case. The aerodynamic Weber number is
`We_A = ρ_g U² D_l / σ_l` with `U` the gas–liquid relative velocity and `D_l` the *liquid tube
inner diameter*; breakup length and first-droplet location are both normalised by `D_l`.

Note the companion paper, same authors and year, is a **different** reference and is easily
confused with this one: *Characteristics of sprays produced by coaxial non-swirling and
swirling air–water jets with high aerodynamic Weber numbers*, Physics of Fluids **34**(10),
103604, [10.1063/5.0107480](https://doi.org/10.1063/5.0107480). Earlier revisions of this file
cited a conflation of the two (this DOI, the other's article number and title).

## Toolkits

- **VTK** — the CPU tracer is built on it. <https://vtk.org/>
- **CUDA** — `gpu/gpuShadow.cu`; build with `nvcc -O3 -arch=sm_86 -cudart static`.
