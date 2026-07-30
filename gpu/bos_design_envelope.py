#!/usr/bin/env python3
"""Size a BOS rig BEFORE building it, from a predicted deflection field.

THE PROBLEM THIS SOLVES. A correlation window tracks displacements only between a subpixel
noise floor (~0.1 px) and roughly a quarter of the window. Which part of a flow lands inside
that band depends on the background standoff L_bg, because d = L_bg.tan(eps). So L_bg is a
design decision -- and in a real experiment you cannot make it well, because you do not know
the deflection distribution until you have measured, by which point L_bg is already baked
into the data. The distribution is also the hard part to guess: the median is predictable from
a similarity estimate, but the TAIL is set by the smallest resolved structure, which is
precisely what you are trying to measure. Practice therefore bootstraps -- multi-pass
interrogation with a coarse predictor (Keane & Adrian's conditions; see docs/REFERENCES.md),
or a bracketing shot at short standoff to read the histogram before moving the background out.

A CFD-driven synthetic render sidesteps the loop: it gives the distribution before the rig
exists. This script turns that distribution into the two plots you actually need to choose
hardware.

  Panel A  the ground-truth displacement field |d| -- what the ray trace says the background
           actually moves by, which is the quantity a real rig can never see directly.
  Panel B  WHERE the field is unmeasurable at a chosen standoff: below the subpixel noise
           floor, inside the band, or past the win/4 correlation ceiling. This is the panel
           that makes the trade physical -- on a jet the out-of-range pixels are the shear
           layer and the near-nozzle fronts, i.e. the part worth measuring.
  Panel C  the deflection distribution against the usable band for several standoffs,
           showing directly that no single L_bg covers a wide-dynamic-range flow.
  Panel D  measurable fraction vs L_bg for a range of window sizes -- the design curve,
           with the best achievable coverage marked. On flows like this it peaks well
           below 100%, and that ceiling is the honest thing to know in advance.

Spatial panels use aspect='equal' so 1 mm in x is 1 mm in y -- the field is 160x64 mm, so
they render 2.5:1. Squaring the axes would misrepresent the geometry.

USAGE
    python3 gpu/bos_design_envelope.py field.npz [L_bg_of_field_mm]

`field.npz` holds `dx`, `dy` (apparent background displacement, mm) and `mm_per_px` -- i.e.
what gpuShadow outmodes 5/6 produce. The second argument is the standoff the field was
rendered at, needed to convert displacement back to deflection angle (default 300 mm).
Writes <field>_design_envelope.png and _design_envelope.npz beside the input.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NOISE_PX = 0.1          # subpixel correlation noise floor
WINS = [16, 32, 48, 64, 96]
LBG_MM = np.logspace(np.log10(10), np.log10(3000), 220)
SHOW_LBG = [20, 100, 300, 1000]

field = sys.argv[1] if len(sys.argv) > 1 else 'field.npz'
lbg_ref = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0

z = np.load(field)
mm_per_px = float(z['mm_per_px'])
d_mm = np.sqrt(np.nan_to_num(z['dx']) ** 2 + np.nan_to_num(z['dy']) ** 2)
eps = (d_mm[d_mm > 0] / lbg_ref) * 1e3          # deflection magnitude, mrad
eps.sort()

pct = {q: float(np.percentile(eps, q)) for q in (50, 90, 99, 99.9)}
print(f"field {field}: {eps.size:,} deflected pixels, {mm_per_px:.4f} mm/px, "
      f"rendered at L_bg = {lbg_ref:.0f} mm")
print(f"  |eps|  median {pct[50]:.3f}  p90 {pct[90]:.3f}  p99 {pct[99]:.3f}  "
      f"p99.9 {pct[99.9]:.3f}  max {eps[-1]:.3f}  mrad")
print(f"  dynamic range median -> p99.9 : {pct[99.9]/pct[50]:.0f}x")


def coverage(lbg, win):
    """Fraction of the field inside [noise floor, win/4] at this standoff."""
    d_px = eps * 1e-3 * lbg / mm_per_px
    return float(np.mean((d_px >= NOISE_PX) & (d_px <= win / 4.0)))


cov = {w: np.array([coverage(L, w) for L in LBG_MM]) for w in WINS}
best = {w: (LBG_MM[int(np.argmax(cov[w]))], float(cov[w].max())) for w in WINS}
print("\n  best achievable coverage per window:")
for w in WINS:
    L, c = best[w]
    print(f"    win {w:>3} px : {100*c:5.1f} % of the field, at L_bg = {L:6.0f} mm")

MASK_LBG, MASK_WIN = 300.0, 32
d_px_ref = d_mm / mm_per_px * (MASK_LBG / lbg_ref)
cat = np.full(d_px_ref.shape, 1, dtype=np.int8)          # 1 = measurable
cat[d_px_ref < NOISE_PX] = 0                             # 0 = below noise floor
cat[d_px_ref > MASK_WIN / 4.0] = 2                       # 2 = past correlation ceiling
fr = [float(np.mean(cat == k)) for k in (0, 1, 2)]
print(f"\n  spatial breakdown at L_bg = {MASK_LBG:.0f} mm, {MASK_WIN} px window:")
print(f"    below noise floor  {100*fr[0]:5.1f} %")
print(f"    measurable         {100*fr[1]:5.1f} %")
print(f"    past win/4 ceiling {100*fr[2]:5.1f} %")

extent = [0, d_mm.shape[1] * mm_per_px, 0, d_mm.shape[0] * mm_per_px]

fig = plt.figure(figsize=(12.4, 13.4), dpi=150)
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.45], hspace=0.26, wspace=0.22)
axT = fig.add_subplot(gs[0, :])
axM = fig.add_subplot(gs[1, :])
axA = fig.add_subplot(gs[2, 0])
axB = fig.add_subplot(gs[2, 1])

# ---------------- Panel A: ground-truth displacement field ----------------
imT = axT.imshow(d_mm, origin='lower', extent=extent, aspect='equal',
                 cmap='magma', vmin=0, vmax=float(np.percentile(d_mm, 99.5)))
axT.set_title(f'A · ground-truth displacement |d| at $L_{{bg}}$ = {lbg_ref:.0f} mm — '
              f'what no experiment can read directly', fontsize=10.5)
axT.set_xlabel('x (mm)'); axT.set_ylabel('y (mm)')
# reserve identical colourbar width on BOTH spatial axes so A and B render the same size
# and can be compared pixel-for-pixel; B's slot is created then hidden.
from mpl_toolkits.axes_grid1 import make_axes_locatable
capT = make_axes_locatable(axT).append_axes('right', size='1.8%', pad=0.10)
cbT = fig.colorbar(imT, cax=capT)
cbT.set_label('|d|  [mm]', fontsize=9); cbT.ax.tick_params(labelsize=8)

# ---------------- Panel B: where it is unmeasurable ----------------
from matplotlib.colors import ListedColormap
cmap = ListedColormap(['#cfd8e0', '#1b7f5f', '#a9432c'])
axM.imshow(cat, origin='lower', extent=extent, aspect='equal', cmap=cmap, vmin=-0.5, vmax=2.5,
           interpolation='nearest')
axM.set_title(f'B · where the measurement fails at $L_{{bg}}$ = {MASK_LBG:.0f} mm, '
              f'{MASK_WIN} px window', fontsize=10.5)
axM.set_xlabel('x (mm)'); axM.set_ylabel('y (mm)')
capM = make_axes_locatable(axM).append_axes('right', size='1.8%', pad=0.10)
capM.axis('off')
from matplotlib.patches import Patch
axM.legend(handles=[Patch(fc='#cfd8e0', ec='none', label=f'below noise floor  {100*fr[0]:.0f}%'),
                    Patch(fc='#1b7f5f', ec='none', label=f'measurable  {100*fr[1]:.0f}%'),
                    Patch(fc='#a9432c', ec='none', label=f'past win/4 ceiling  {100*fr[2]:.0f}%')],
           loc='upper right', fontsize=8.0, frameon=True, framealpha=0.93, ncol=1,
           borderpad=0.55, handlelength=1.4, labelspacing=0.45,
           edgecolor='#8a8f94')

# ---------------- Panel C: distribution vs usable bands ----------------
axA.hist(eps, bins=np.logspace(np.log10(max(eps[0], 1e-4)), np.log10(eps[-1]), 90),
         color='#9fb3c2', edgecolor='none')
axA.set_xscale('log')
axA.set_xlabel('deflection |ε|  [mrad]')
axA.set_ylabel('pixels')
axA.set_title('C · the flow spans more range than one standoff covers', fontsize=10.5)
ymax = axA.get_ylim()[1]
cols = ['#1b7f5f', '#b9741a', '#7a4fa3', '#a9432c']
WIN_A = 32
for c, L in zip(cols, SHOW_LBG):
    lo = NOISE_PX * mm_per_px / L * 1e3
    hi = (WIN_A / 4.0) * mm_per_px / L * 1e3
    frac = coverage(L, WIN_A)
    yb = ymax * (0.95 - 0.085 * SHOW_LBG.index(L))
    axA.plot([lo, hi], [yb] * 2, color=c, lw=3.2, solid_capstyle='butt', zorder=4)
    axA.plot([lo, hi], [yb] * 2, '|', color=c, ms=7, zorder=4)
    axA.text(lo / 1.35, yb, f'{100*frac:.0f}%', color=c, fontsize=8.4, va='center',
             ha='right', fontweight='bold')
    axA.text(hi * 1.25, yb, f'$L_{{bg}}$={L:.0f}', color=c, fontsize=8.0, va='center')
for q, lbl in [(50, 'median'), (90, 'p90'), (99, 'p99'), (99.9, 'p99.9')]:
    axA.axvline(pct[q], color='#41474c', lw=0.9, ls=':')
    axA.text(pct[q], ymax * 0.035, lbl, rotation=90, fontsize=7.4, color='#41474c',
             ha='right', va='bottom')
axA.text(0.02, 0.02, f'{WIN_A} px window · band = [{NOISE_PX} px, win/4]',
         transform=axA.transAxes, fontsize=7.8, color='#5a6066', style='italic')

# ---------------- Panel D: the design curve ----------------
for w, c in zip(WINS, ['#1b7f5f', '#b9741a', '#7a4fa3', '#a9432c', '#2c5f8a']):
    axB.plot(LBG_MM, 100 * cov[w], color=c, lw=1.9, label=f'{w} px window')
    L, cc = best[w]
    axB.plot([L], [100 * cc], 'o', color=c, ms=5, zorder=5)
axB.set_xscale('log')
axB.set_xlabel('background standoff  $L_{bg}$  [mm]')
axB.set_ylabel('measurable fraction of field  [%]')
axB.set_ylim(0, 100)
axB.set_title('D · the design curve — choose $L_{bg}$ before building', fontsize=10.5)
axB.grid(alpha=0.18, lw=0.6)
axB.legend(fontsize=8, frameon=False, loc='lower center', ncol=2)
wbest = max(WINS, key=lambda w: best[w][1])
axB.annotate(f'ceiling {100*best[wbest][1]:.0f}% — no standoff\nmeasures the whole field',
             xy=(best[wbest][0], 100 * best[wbest][1]),
             xytext=(0.045, 0.845), textcoords='axes fraction', fontsize=8.4,
             color='#41474c', ha='left', va='top',
             arrowprops=dict(arrowstyle='->', color='#8a8f94', lw=0.9,
                             shrinkA=6, shrinkB=4,
                             connectionstyle='arc3,rad=0.16'))
axB.axhline(100 * best[wbest][1], color='#a9432c', lw=0.8, ls='--', alpha=0.55)

fig.suptitle('BOS design envelope from a predicted deflection field — '
             'the distribution a real rig cannot know until after it has measured',
             fontsize=11)

base = os.path.splitext(field)[0]
fig.savefig(base + '_design_envelope.png', facecolor='white', bbox_inches='tight')
np.savez(base + '_design_envelope.npz', lbg_mm=LBG_MM,
         coverage=np.array([cov[w] for w in WINS]), wins=np.array(WINS),
         percentiles=np.array([[q, pct[q]] for q in (50, 90, 99, 99.9)]),
         max_mrad=float(eps[-1]), mm_per_px=mm_per_px, lbg_ref_mm=lbg_ref,
         noise_px=NOISE_PX, mask_lbg_mm=MASK_LBG, mask_win_px=MASK_WIN,
         mask_fractions=np.array(fr))
print(f"\nwrote {base}_design_envelope.png and .npz")
