#!/usr/bin/env python3
"""How fine can the BOS interrogation window go, and what limits it?

Two things shrink together when you shrink the window, and they are easy to conflate:
  1. the number of speckle dots inside it   -> correlation quality (the usual >=8-10 rule)
  2. the maximum trackable shift, ~win/4 px -> dynamic range ceiling

So a window sweep at FIXED speckle density tests both at once and cannot separate them.
This sweeps window size at three dot densities, which does separate them: if small windows
recover once the density is raised, the limit was pattern; if they do not, it is the
correlation itself.

Truth is the window-averaged displacement, so each window size is scored against its own
correct answer rather than against a fixed reference.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bos_correlate import synthesise

def main():
    if len(sys.argv) < 2:
        sys.exit('usage: bos_window_sweep.py <deflection_field.npz>\n'
                 '  npz must contain dx, dy (apparent displacement in mm) and mm_per_px.')
    z = np.load(sys.argv[1])
    dx, dy, mm_per_px = z['dx'], z['dy'], float(z['mm_per_px'])

    true_px = np.sqrt(np.nan_to_num(dx) ** 2 + np.nan_to_num(dy) ** 2) / mm_per_px
    fin = true_px[true_px > 0]
    print(f"field: {dx.shape[1]}x{dx.shape[0]} px, {mm_per_px:.4f} mm/px")
    print(f"true |d|: median {np.median(fin):.3f} px, p99 {np.percentile(fin,99):.2f} px, "
          f"max {fin.max():.2f} px\n")

    WINS = [8, 12, 16, 24, 32, 48, 64]
    DENS = [0.06, 0.15, 0.30]

    rows = []      # persisted; printing alone loses the numbers the moment the shell scrolls
    print(f"{'win':>4} {'step':>5} {'density':>8} {'dots/win':>9} {'limit px':>9} "
          f"{'valid %':>8} {'RMS px':>8} {'bias px':>8} {'vectors':>9}")
    for dens in DENS:
        for win in WINS:
            step = max(win // 2, 2)
            dots = dens * win * win / (np.pi * 1.6 ** 2)
            res = synthesise(dx, dy, mm_per_px, win=win, step=step, seed=5, density=dens)
            u, v = res['u_mm'] / mm_per_px, res['v_mm'] / mm_per_px
            tu, tv = res['tu_mm'] / mm_per_px, res['tv_mm'] / mm_per_px
            rec = np.sqrt(u ** 2 + v ** 2); tru = np.sqrt(tu ** 2 + tv ** 2)
            lim = win / 4.0
            # COMMON evaluation band. Scoring each window on its own tru<win/4 mask is a
            # selection effect: a small window is then graded only on small, easy shifts and
            # looks artificially accurate. Every window is instead scored on the same band,
            # set by the strictest ceiling in the sweep (smallest win/4).
            BAND = min(WINS) / 4.0
            ok = np.isfinite(rec) & (tru < BAND)
            err = rec[ok] - tru[ok]
            rms = float(np.sqrt(np.mean(err ** 2))) if err.size else float('nan')
            bias = float(np.mean(err)) if err.size else float('nan')
            rows.append(dict(win=win, step=step, density=dens, dots_per_win=float(dots),
                             limit_px=lim, valid_frac=float(res['valid']), rms_px=rms,
                             bias_px=bias, n_vectors=int(rec.size), band_px=BAND))
            print(f"{win:>4} {step:>5} {dens:>8.2f} {dots:>9.1f} {lim:>9.1f} "
                  f"{100*res['valid']:>8.1f} {rms:>8.3f} {bias:>+8.3f} {rec.size:>9,}",
                  flush=True)
        print()


if __name__ == '__main__':
    main()

    out = os.path.splitext(CACHE)[0] + '_window_sweep.npz'
    np.savez(out, rows=np.array(rows, dtype=object),
             mm_per_px=mm_per_px, max_true_px=float(fin.max()),
             median_true_px=float(np.median(fin)), p99_true_px=float(np.percentile(fin, 99)))
    print(f"\nwrote {out}")
    print("NOTE: interpret against Keane & Adrian (10.1088/0957-0233/1/11/013) -- see README.")
