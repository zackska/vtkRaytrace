#!/usr/bin/env python3
"""Optimise a full BOS optical setup, not just the standoff -- sensitivity vs geometric blur.

WHY THIS SUPERSEDES bos_design_envelope.py. That script sweeps the background distance alone
and counts a pixel "measurable" if its displacement lands in [noise floor, win/4]. Two things
are wrong with treating that as a design answer, both from Schmidt et al., "Twenty-Five Years
of Background-Oriented Schlieren", AIAA J 63(12) 2025 (10.2514/1.J065669), Sec. IV, which
follows Schwarz & Braukmann:

1. DISPLACEMENT IS NOT SET BY STANDOFF ALONE. The image-plane displacement is

       delta_v = S . eps        with    S = f . z_D / (z_D + z_A - f)          (Eqs 22, 23)

   for focal length f, lens-to-object z_A and object-to-background z_D. The naive
   d = L_bg . eps ignores the (z_D + z_A - f) denominator and OVERESTIMATES the recorded
   shift -- by ~35% at the geometry used for the helium-jet renders here.

2. SENSITIVITY AND BLUR ARE THE SAME KNOB. The circle of confusion of the schlieren object

       CoC = f^2 z_D / (f# z_A (z_A + z_D - f))                                (Eq 25)

   is LINEARLY proportional to S. You cannot raise sensitivity without proportionally blurring
   the object, at fixed f#, field of view and sensor size. Only f# breaks the coupling. The
   minimum resolvable flow feature is about 50% of the CoC mapped into object space
   (Schwarz & Braukmann). So chasing the coverage peak by pushing the background away blurs
   out the very structure that made the tail measurable -- a trade the standoff-only sweep
   cannot see.

Two further results from the same section, worth having explicitly:

  * For a fixed field of view and total setup length z_B = z_A + z_D, sensitivity S is
    MAXIMISED at z_A / z_B = 0.5 -- object halfway between camera and background. The review
    notes this is rarely stated explicitly.
  * Increasing z_B at fixed z_A/z_B needs longer focal length and gives higher S. Longer
    benches win; moving only the background does not.

The focal length follows from the field-of-view constraint (Eq 24):

       f = L_sens . z_A / (L_FoV + L_sens)

USAGE
    python3 gpu/bos_setup_optimise.py field.npz [L_bg_of_field_mm]

`field.npz` holds `dx`, `dy` (mm) and `mm_per_px`, i.e. gpuShadow outmodes 5/6. The standoff
argument is only used to convert the rendered displacement back to deflection angle.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- camera / optics, stated so the numbers are reproducible -------------------
L_SENS = 28.0          # mm, sensor width (1400 px at 20 um)
PIX = 0.020            # mm, pixel pitch
L_FOV = 160.0          # mm, required field of view (the jet domain)
WIN = 32               # px, interrogation window
NOISE_PX = 0.1         # px, subpixel correlation noise floor
FNUMS = [2.8, 5.6, 11.0, 22.0]
Z_A = 1000.0           # mm, lens-to-object (fixed for the S sweep)
Z_B_SET = [500.0, 1000.0, 2000.0, 4000.0]   # mm, total setup lengths for the ratio panel

field = sys.argv[1] if len(sys.argv) > 1 else 'field.npz'
lbg_ref = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0

z = np.load(field)
mm_per_px = float(z['mm_per_px'])
d_mm = np.sqrt(np.nan_to_num(z['dx']) ** 2 + np.nan_to_num(z['dy']) ** 2)
eps = (d_mm[d_mm > 0] / lbg_ref)                 # deflection magnitude, rad
eps.sort()
pct = {q: float(np.percentile(eps, q)) * 1e3 for q in (50, 90, 99, 99.9)}

f_len = L_SENS * Z_A / (L_FOV + L_SENS)          # Eq 24
M = f_len / (Z_A - f_len)                        # magnification
print(f"field {field}: {eps.size:,} px | FoV {L_FOV:.0f} mm, sensor {L_SENS:.0f} mm "
      f"-> f = {f_len:.1f} mm, M = {M:.4f}")
print(f"  |eps| median {pct[50]:.3f}  p90 {pct[90]:.3f}  p99 {pct[99]:.3f}  "
      f"p99.9 {pct[99.9]:.3f} mrad")


def S_of(zD, zA=Z_A, f=None):
    f = f_len if f is None else f
    return f * zD / (zD + zA - f)                # Eq 23


def coc_obj(S, fnum, zA=Z_A, f=None):
    """Circle of confusion mapped into object space [mm]; CoC = S.f/(f# zA)."""
    f = f_len if f is None else f
    return (S * f / (fnum * zA)) / M


def coverage_S(S):
    dpx = S * eps / PIX
    return float(np.mean((dpx >= NOISE_PX) & (dpx <= WIN / 4.0)))


# --- naive vs correct displacement, at the rendered standoff -------------------
S_ref = S_of(lbg_ref)
naive_px = np.percentile(eps, 99) * lbg_ref / mm_per_px
corr_px = S_ref * np.percentile(eps, 99) / PIX
print(f"\n  at z_D = {lbg_ref:.0f} mm: S = {S_ref:.2f} mm")
print(f"    p99 shift, naive  d = L_bg.eps : {naive_px:.2f} px")
print(f"    p99 shift, Eq 22  d = S.eps    : {corr_px:.2f} px "
      f"({100*(naive_px/corr_px-1):+.0f}% naive error)")

ZD = np.logspace(np.log10(20), np.log10(20000), 300)
S_sweep = S_of(ZD)
cov = np.array([coverage_S(s) for s in S_sweep])
ibest = int(np.argmax(cov))
print(f"\n  coverage peaks at {100*cov[ibest]:.1f}% : S = {S_sweep[ibest]:.1f} mm, "
      f"z_D = {ZD[ibest]:.0f} mm")
print(f"  min resolvable feature there:")
for fn in FNUMS:
    print(f"    f/{fn:<5} {coc_obj(S_sweep[ibest], fn)*0.5:7.3f} mm")

fig, (axR, axS) = plt.subplots(1, 2, figsize=(13.4, 5.6), dpi=150)

# ---- Panel A: the z_A/z_B = 0.5 optimum --------------------------------------
r = np.linspace(0.02, 0.98, 400)
for zB, c in zip(Z_B_SET, ['#1b7f5f', '#b9741a', '#7a4fa3', '#a9432c']):
    zA = r * zB
    fz = L_SENS * zA / (L_FOV + L_SENS)
    axR.plot(r, S_of(zB - zA, zA, fz), color=c, lw=1.9, label=f'$z_B$ = {zB:.0f} mm')
axR.axvline(0.5, color='#41474c', ls='--', lw=1.0)
axR.text(0.51, 0.04, ' $z_A/z_B$ = 0.5\n optimum (Schmidt et al. §IV)', fontsize=8.2,
         color='#41474c', transform=axR.get_xaxis_transform(), va='bottom')
axR.set_xlabel('$z_A / z_B$   (lens-to-object / total setup length)')
axR.set_ylabel('sensitivity  $S$  [mm]')
axR.set_title('A · sensitivity peaks with the object halfway,\nand longer benches win',
              fontsize=10.5)
axR.legend(fontsize=8, frameon=False, loc='upper left')
axR.grid(alpha=0.18, lw=0.6)

# ---- Panel B: coverage vs blur, the real trade --------------------------------
axS.plot(S_sweep, 100 * cov, color='#2c5f8a', lw=2.2, label='measurable fraction')
axS.plot([S_sweep[ibest]], [100 * cov[ibest]], 'o', color='#2c5f8a', ms=6)
axS.set_xscale('log')
axS.set_xlabel('sensitivity  $S$  [mm]')
axS.set_ylabel('measurable fraction of field  [%]', color='#2c5f8a')
axS.tick_params(axis='y', labelcolor='#2c5f8a')
axS.set_ylim(0, 100)
axS.grid(alpha=0.18, lw=0.6)

ax2 = axS.twinx()
for fn, c in zip(FNUMS, ['#a9432c', '#b9741a', '#7a4fa3', '#1b7f5f']):
    ax2.plot(S_sweep, 0.5 * coc_obj(S_sweep, fn), color=c, lw=1.5, ls='--',
             label=f'f/{fn:g}')
ax2.set_yscale('log')
ax2.set_ylabel('min resolvable feature  [mm]  (dashed)', color='#41474c')
ax2.tick_params(axis='y', labelcolor='#41474c')
ax2.legend(fontsize=7.6, frameon=False, loc='upper left', title='blur, by $f_\\#$',
           title_fontsize=7.6)
axS.set_title('B · raising sensitivity blurs the object in proportion —\n'
              'only $f_\\#$ breaks the coupling', fontsize=10.5)

fig.suptitle('BOS setup optimisation: sensitivity and geometric blur are one knob '
             '(Schmidt et al., AIAA J 2025, §IV)', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
base = os.path.splitext(field)[0]
fig.savefig(base + '_setup_optimise.png', facecolor='white', bbox_inches='tight')
np.savez(base + '_setup_optimise.npz', zD=ZD, S=S_sweep, coverage=cov,
         fnums=np.array(FNUMS),
         feature_mm=np.array([0.5 * coc_obj(S_sweep, fn) for fn in FNUMS]),
         f_len=f_len, M=M, z_A=Z_A, pix=PIX, L_sens=L_SENS, L_fov=L_FOV, win=WIN,
         naive_p99_px=naive_px, correct_p99_px=corr_px)
print(f"\nwrote {base}_setup_optimise.png and .npz")
