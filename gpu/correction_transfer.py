#!/usr/bin/env python3
"""Does a CFD-informed correction beat an experiment-only one? (Zack's framing)

Two classes of correction are available to someone doing BOS velocimetry:

  A. EXPERIMENT-ONLY — calibrate the recovered field against something measurable. The
     minimal realistic version is a single-point calibration: one probe/LDV reading on the
     axis, used to set one global scale factor. No simulation needed.

  B. CFD-INFORMED — run a representative simulation, extract the spatial bias profile
     beta(x) = U_true(x) / U_bos(x) which is NOT measurable, and transfer it to the
     experiment.

B is only worth anything if it TRANSFERS. Zack's earlier point is the risk: the relation
between the optical weighting and the velocity varies with flow and time, so a beta(x)
derived from one realisation may not apply to another.

This tests transfer in time, which is the part the present data can settle: the 19
frame-pairs are split into an early half (the "representative CFD") and a late half (the
"experiment"). beta(x) is fitted on the early half only and applied to the late half, and
both corrections are scored against the late half's true centreline velocity.

Scoring on the SAME half that fitted the correction would be circular, which is the whole
reason for the split.
"""
import os, sys
import numpy as np

CASE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
z = np.load(os.path.join(CASE, 'velo_result.npz'))
US = z['US']; xg = z['xg']; yg = z['yg']; b = z['bounds']
mmpx = float(z['mm_per_px']); NXv = z['U_centreline'].size
Ucl = z['U_centreline']

npairs = US.shape[2]
half = npairs // 2
print(f"  {npairs} frame-pairs; representative = pairs 1-{half}, experiment = pairs {half+1}-{npairs}")

# axis row of the correlation grid
y_mm = b[2] * 1e3 + yg * mmpx
axis_row = int(np.argmin(np.abs(y_mm - 0.5 * (b[2] + b[3]) * 1e3)))
x_corr = b[0] * 1e3 + xg * mmpx                      # mm
x_vol = np.linspace(b[0], b[1], NXv) * 1e3
Ucl_at_corr = np.interp(x_corr, x_vol, Ucl)          # accessible truth on the correlation grid

rep = np.nanmedian(US[axis_row, :, :half], axis=1)
exp = np.nanmedian(US[axis_row, :, half:], axis=1)

ok = np.isfinite(rep) & np.isfinite(exp) & (Ucl_at_corr > 2.0)   # inside the jet
print(f"  usable stations along the axis: {ok.sum()} of {ok.size}")

# --- A: experiment-only, one probe station, one global scale ---
probe = np.where(ok)[0][len(np.where(ok)[0]) // 2]                # a single mid-jet station
sA = Ucl_at_corr[probe] / exp[probe]
corrA = exp * sA

# --- B: CFD-informed, full beta(x) from the representative half ---
beta = np.full_like(rep, np.nan)
beta[ok] = Ucl_at_corr[ok] / rep[ok]
corrB = exp * beta

def rms(a, t, m):
    d = a[m] - t[m]
    return float(np.sqrt(np.nanmean(d ** 2)))

e_raw = rms(exp, Ucl_at_corr, ok)
e_A = rms(corrA, Ucl_at_corr, ok)
e_B = rms(corrB, Ucl_at_corr, ok)
print(f"\n  RMS error against true centreline over the EXPERIMENT half [m/s]:")
print(f"    uncorrected                       {e_raw:7.3f}")
print(f"    A  experiment-only (1 probe)      {e_A:7.3f}   ({100*(1-e_A/e_raw):+.0f}% vs raw)")
print(f"    B  CFD-informed beta(x)           {e_B:7.3f}   ({100*(1-e_B/e_raw):+.0f}% vs raw)")
print(f"\n  B beats A by {100*(1-e_B/e_A):+.0f}% " +
      ("-> CFD correction transfers" if e_B < e_A else "-> CFD correction does NOT transfer"))

# --- does beta drift between the halves? (Zack's objection, measured) ---
beta_exp = np.full_like(exp, np.nan)
beta_exp[ok] = Ucl_at_corr[ok] / exp[ok]
drift = np.abs(beta_exp[ok] - beta[ok]) / np.abs(beta[ok])
print(f"\n  beta(x) drift between halves: median {100*np.median(drift):.1f}%, "
      f"p90 {100*np.percentile(drift,90):.1f}%")
print(f"  mean beta: representative {np.nanmean(beta[ok]):.3f}, experiment {np.nanmean(beta_exp[ok]):.3f}")

np.savez(os.path.join(CASE, 'transfer_result.npz'),
         x_corr=x_corr, ok=ok, rep=rep, exp=exp, Ucl=Ucl_at_corr,
         corrA=corrA, corrB=corrB, beta=beta, beta_exp=beta_exp,
         e_raw=e_raw, e_A=e_A, e_B=e_B, probe=probe,
         U_rho=z['U_rho'], U_gradrho=z['U_gradrho'], U_centreplane=z['U_centreplane'],
         x_vol=x_vol, bounds=b)
print("  wrote transfer_result.npz")
