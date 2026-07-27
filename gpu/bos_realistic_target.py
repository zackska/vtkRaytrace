#!/usr/bin/env python3
"""Does synthetic Gaussian speckle flatter the BOS result? Test against a realistic target.

The idealised pipeline in bos_correlate.py generates Gaussian dots, warps them, and
correlates. Two things about that are kinder than reality:

  1. GAUSSIAN dots have a smooth, compact spatial spectrum. A PRINTED target is closer to
     binary -- hard ink edges -- so it carries much more high-frequency content, which the
     camera then band-limits through its PSF and pixel aperture. The correlation peak shape
     is set by that spectrum, so this is not cosmetic.
  2. Reference and distorted frames are generated from the SAME noiseless array, so the only
     difference between them is the displacement. A real rig takes two exposures with
     INDEPENDENT shot and read noise, which decorrelates them slightly even at zero flow.

This builds a printed-target model with the physical chain in the right order --
print -> optical warp -> camera PSF -> pixel integration -> shot noise -> read noise ->
quantisation -- with independent noise on each exposure, and scores it against the same
ground truth.
"""
import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bos_correlate import correlate, speckle, warp


def printed_target(shape, density=0.30, radius_px=1.6, seed=0, edge=0.15):
    """Near-binary printed dots: hard-edged discs, not Gaussian.

    `edge` is the ink-edge softness in pixels at the print stage (a real printer/plotter
    does not produce a perfect step). The optical PSF is applied later, separately.
    """
    rng = np.random.default_rng(seed)
    H, W = shape
    img = np.zeros((H, W), np.float32)
    n = int(density * H * W / (np.pi * radius_px ** 2))
    ys = rng.uniform(0, H, n); xs = rng.uniform(0, W, n)
    r = int(np.ceil(radius_px + 3 * edge + 2))
    yy, xx = np.mgrid[0:H, 0:W]
    for y0, x0 in zip(ys, xs):
        iy, ix = int(y0), int(x0)
        y1, y2 = max(iy - r, 0), min(iy + r + 1, H)
        x1, x2 = max(ix - r, 0), min(ix + r + 1, W)
        if y2 <= y1 or x2 <= x1:
            continue
        d = np.hypot(yy[y1:y2, x1:x2] - y0, xx[y1:y2, x1:x2] - x0)
        # smoothstep edge: ~1 inside the disc, ~0 outside, soft over `edge` px
        img[y1:y2, x1:x2] = np.maximum(img[y1:y2, x1:x2],
                                       np.clip((radius_px - d) / max(edge, 1e-6), 0, 1))
    return np.clip(img, 0, 1)


def camera(img, psf_px=0.8, full_well=8000.0, read_e=12.0, bits=8, rng=None):
    """Optical PSF -> shot noise -> read noise -> quantisation. Independent per exposure."""
    if rng is None:
        rng = np.random.default_rng()
    x = gaussian_filter(img.astype(np.float32), psf_px)          # optics + pixel aperture
    e = np.clip(x, 0, 1) * full_well                              # electrons
    e = rng.poisson(e).astype(np.float32)                         # shot noise
    e += rng.normal(0.0, read_e, e.shape)                         # read noise
    y = np.clip(e / full_well, 0, 1)
    levels = 2 ** bits - 1
    return (np.round(y * levels) / levels).astype(np.float32)     # quantise


def run_case(dx, dy, mm_per_px, kind, win=32, step=16, density=0.30, seed=11,
             psf_px=0.8, full_well=8000.0, read_e=12.0, bits=8, noisy=True):
    dx_px = np.nan_to_num(dx) / mm_per_px
    dy_px = np.nan_to_num(dy) / mm_per_px
    if kind == 'gaussian':
        base = speckle(dx.shape, density=density, seed=seed)
    else:
        base = printed_target(dx.shape, density=density, seed=seed)
    warped = warp(base, dx_px, dy_px)
    if noisy:
        r1 = np.random.default_rng(seed * 7 + 1)
        r2 = np.random.default_rng(seed * 7 + 2)          # INDEPENDENT noise per exposure
        ref = camera(base, psf_px, full_well, read_e, bits, r1)
        dis = camera(warped, psf_px, full_well, read_e, bits, r2)
    else:
        ref, dis = base, warped
    u, v, xg, yg = correlate(ref, dis, win=win, step=step)
    tu = np.full_like(u, np.nan); tv = np.full_like(u, np.nan)
    for jj, y0 in enumerate((yg - win / 2).astype(int)):
        for ii, x0 in enumerate((xg - win / 2).astype(int)):
            tu[jj, ii] = np.nanmean(dx_px[y0:y0 + win, x0:x0 + win])
            tv[jj, ii] = np.nanmean(dy_px[y0:y0 + win, x0:x0 + win])
    rec = np.sqrt(u ** 2 + v ** 2); tru = np.sqrt(tu ** 2 + tv ** 2)
    BAND = 2.0
    ok = np.isfinite(rec) & (tru < BAND)
    err = rec[ok] - tru[ok]
    return dict(valid=float(np.isfinite(u).mean()),
                rms=float(np.sqrt(np.mean(err ** 2))) if err.size else np.nan,
                bias=float(np.mean(err)) if err.size else np.nan,
                contrast=float(ref.std()))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: bos_realistic_target.py <deflection_field.npz>\n'
                 '  npz must contain dx, dy (apparent displacement in mm) and mm_per_px,\n'
                 '  as written by gpuShadow outmodes 5/6.')
    z = np.load(sys.argv[1])
    dx, dy, mmpx = z['dx'], z['dy'], float(z['mm_per_px'])
    print(f"{'target':<12} {'noise':<7} {'psf px':>7} {'bits':>5} {'contrast':>9} "
          f"{'valid %':>8} {'RMS px':>8} {'bias px':>8}")
    cases = [
        ('gaussian', dict(noisy=False),                       'none'),
        ('printed',  dict(noisy=False),                       'none'),
        ('gaussian', dict(noisy=True),                        'full'),
        ('printed',  dict(noisy=True),                        'full'),
        ('printed',  dict(noisy=True, psf_px=1.5),            'full'),
        ('printed',  dict(noisy=True, full_well=1500.0),      'full'),
        ('printed',  dict(noisy=True, bits=12),               'full'),
    ]
    for kind, kw, lab in cases:
        r = run_case(dx, dy, mmpx, kind, **kw)
        print(f"{kind:<12} {lab:<7} {kw.get('psf_px',0.8):>7.1f} {kw.get('bits',8):>5} "
              f"{r['contrast']:>9.4f} {100*r['valid']:>8.1f} {r['rms']:>8.3f} "
              f"{r['bias']:>+8.3f}", flush=True)
