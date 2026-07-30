#!/usr/bin/env python3
"""Instantaneous singular-value-gap diagnostic on Liang's experimental shadowgraphs.

Taken from Zamani Ashtiani & Fukami, AIAA J 2026 (10.2514/1.J067033), Fig. 7:
    delta(t) = (sigma_1(t) - sigma_2(t)) / sigma_1(t)
computed from the SVD of an INSTANTANEOUS snapshot. delta -> 1 means energy is
concentrated in one coherent structure; delta -> 0 means it is spread across many modes,
i.e. multiscale and disorganised. They use it to rank extreme-gust cases by coherence.

WHAT THIS IS AND IS NOT. This is the paper's DIAGNOSTIC, not its ALGORITHM. Full
time-dependent bases solves evolution equations for (Phi, T, Psi) to avoid storing full
history of huge 3-D fields. These are 100x80 images: the instantaneous SVD is ~microseconds,
so the machinery solves a problem we do not have. What transfers is the metric.

TWO VARIANTS, because the choice matters and the paper's field (vorticity) has no DC:
  delta_raw  : SVD of the frame as-is. sigma_1 is then dominated by the bright static
               background and the nozzle, so delta saturates near 1 and says little.
  delta_fluc : SVD after subtracting the per-condition temporal mean frame. This removes
               the static background and leaves the fluctuating spray, which is the
               closer analogue of their vorticity field. PRIMARY.

Also reported per frame:
  r99   : number of modes for 99% of squared-singular-value energy (complexity/rank)
  Hn    : normalised spectral entropy of sigma^2 (0 = one mode, 1 = flat spectrum)

Data: frames_real.npz, 96 conditions (8 swirl x 4 air x 3 samplings) x 40 frames,
meta = (S, air, sampling). At S=0 the air->We_A map is exact and validated:
air 3/4/5/8/10 -> We_A 9/16/25/64/100.
"""
import os, sys
import json
import numpy as np

FRAMES = sys.argv[1] if len(sys.argv) > 1 else 'frames_real.npz'
Z = np.load(FRAMES, allow_pickle=True)
X, META = Z['X'], Z['meta']
H, W = int(Z['H']), int(Z['W'])
AIR2WE = {3: 9, 4: 16, 5: 25, 8: 64, 10: 100}


def spectrum_stats(img):
    s = np.linalg.svd(img, compute_uv=False)
    s = s[s > 0]
    if s.size < 2:
        return np.nan, np.nan, np.nan
    e = s ** 2
    p = e / e.sum()
    delta = (s[0] - s[1]) / s[0]
    r99 = int(np.searchsorted(np.cumsum(p), 0.99) + 1)
    hn = float(-(p * np.log(p + 1e-300)).sum() / np.log(len(p)))
    return delta, r99, hn


conds = np.unique(META, axis=0)
rows = []
for c in conds:
    sel = np.all(META == c, axis=1)
    frames = X[sel].reshape(-1, H, W)
    mean_frame = frames.mean(axis=0)
    dr, df, rr, hh = [], [], [], []
    for f in frames:
        a, _, _ = spectrum_stats(f)
        b, r, h = spectrum_stats(f - mean_frame)
        dr.append(a); df.append(b); rr.append(r); hh.append(h)
    S, air, samp = float(c[0]), float(c[1]), float(c[2])
    rows.append(dict(S=S, air=air, samp=samp, n=len(frames),
                     d_raw=float(np.nanmean(dr)),
                     d_fluc=float(np.nanmean(df)), d_fluc_sd=float(np.nanstd(df)),
                     r99=float(np.nanmean(rr)), Hn=float(np.nanmean(hh))))

np.save(os.path.splitext(FRAMES)[0] + '_sigma_gap.npy', rows, allow_pickle=True)

print(f"frames {X.shape[0]}  image {H}x{W}  conditions {len(conds)}\n")
print("=== S = 0 (air->We_A map exact and validated here) ===")
print(f"{'air':>4} {'We_A':>5} {'samp':>5} {'d_raw':>7} {'d_fluc':>8} {'d_sd':>7} {'r99':>6} {'Hn':>6}")
s0 = [r for r in rows if r['S'] == 0]
for r in sorted(s0, key=lambda r: (r['air'], r['samp'])):
    we = AIR2WE.get(int(r['air']), None)
    print(f"{int(r['air']):>4} {str(we) if we else '?':>5} {int(r['samp']):>5} "
          f"{r['d_raw']:>7.4f} {r['d_fluc']:>8.4f} {r['d_fluc_sd']:>7.4f} "
          f"{r['r99']:>6.1f} {r['Hn']:>6.3f}")

print("\n--- S=0 collapsed over the 3 samplings ---")
print(f"{'We_A':>5} {'d_fluc':>8} {'+/-':>7} {'r99':>6} {'Hn':>6}")
for air in sorted({r['air'] for r in s0}):
    g = [r for r in s0 if r['air'] == air]
    we = AIR2WE.get(int(air), '?')
    print(f"{str(we):>5} {np.mean([r['d_fluc'] for r in g]):>8.4f} "
          f"{np.std([r['d_fluc'] for r in g]):>7.4f} "
          f"{np.mean([r['r99'] for r in g]):>6.1f} {np.mean([r['Hn'] for r in g]):>6.3f}")

print("\n=== swirl dependence, collapsed over air and sampling ===")
print(f"{'S':>5} {'d_fluc':>8} {'r99':>6} {'Hn':>6}")
for S in sorted({r['S'] for r in rows}):
    g = [r for r in rows if r['S'] == S]
    print(f"{S:>5.1f} {np.mean([r['d_fluc'] for r in g]):>8.4f} "
          f"{np.mean([r['r99'] for r in g]):>6.1f} {np.mean([r['Hn'] for r in g]):>6.3f}")

# monotonicity of the S=0 trend, which is the clean anchor
airs = sorted({r['air'] for r in s0})
mu = [np.mean([r['d_fluc'] for r in s0 if r['air'] == a]) for a in airs]
try:
    from scipy.stats import spearmanr
    rho, pv = spearmanr(airs, mu)
    print(f"\nS=0 monotonicity of d_fluc vs air (=We_A): Spearman rho={rho:+.3f} p={pv:.3f}")
except Exception as e:
    print(f"\n(scipy unavailable: {e})")
