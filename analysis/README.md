# analysis/

Study harnesses, not library code. These consume the renderers in `gpu/` (or plain image
data) to answer one specific question, and each carries assumptions from the study it was
written for. Read the header before reusing.

| Script | What it does | Study-specific assumptions |
|---|---|---|
| `ladder_strip.py` | Renders a 5-panel shadowgraph comparison strip across a Weber-number ladder, cropped to a *common* box so panels are like-for-like | `CASES` is hard-coded to `We{16,18,20,22,25}_S0_L2`; expects each case staged as a directory containing `system/blockMeshDict`, `case.foam` and reconstructed time dirs holding `alpha.water` |
| `sigma_gap.py` | Instantaneous singular-value diagnostics on a shadowgraph image set — the relative gap δ = (σ₁−σ₂)/σ₁, the rank for 99% spectral energy, and normalised spectral entropy | Expects the `frames_real.npz` layout (`X`, `meta`, `H`, `W`) built from the Liang (2022) coaxial-swirl dataset; the `AIR2WE` map is calibrated for **S=0 only** |

```bash
python3 analysis/ladder_strip.py  /path/to/staged/cases     # writes ladder_strip.png there
python3 analysis/sigma_gap.py     /path/to/frames_real.npz  # writes *_sigma_gap.npy alongside
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
