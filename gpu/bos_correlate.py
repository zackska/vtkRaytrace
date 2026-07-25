"""Emulate what a real BOS rig actually measures, rather than the ground-truth deflection.

gpuShadow outputs the exact apparent displacement of the background, d = Lbg*tan(eps).
No experiment can measure that directly. A real BOS rig:
  1. photographs a speckled background through the undisturbed medium   -> reference image
  2. photographs it again through the flow, where each ray lands displaced -> distorted image
  3. recovers displacement by CROSS-CORRELATING interrogation windows between the two

Step 3 is the lossy part: the recovered field is smoothed at the window scale, suffers
peak-locking toward integer pixels, and drops out entirely where the pattern is sheared
beyond recognition. This module reproduces all three by doing the actual correlation, so
a synthetic BOS image can be compared with an experimental one on equal terms.
"""
import numpy as np

def speckle(shape, density=0.06, radius=1.6, seed=0):
    """Random dot background, the usual BOS target."""
    rng = np.random.default_rng(seed)
    H, W = shape
    img = np.zeros((H, W), np.float32)
    n = int(density*H*W/(np.pi*radius**2))
    ys = rng.uniform(0, H, n); xs = rng.uniform(0, W, n)
    yy, xx = np.mgrid[0:H, 0:W]
    # accumulate gaussian dots (vectorised over a small stamp for speed)
    r = int(np.ceil(3*radius))
    for y0, x0 in zip(ys, xs):
        iy, ix = int(y0), int(x0)
        y1, y2 = max(iy-r, 0), min(iy+r+1, H); x1, x2 = max(ix-r, 0), min(ix+r+1, W)
        if y2 <= y1 or x2 <= x1: continue
        sy = yy[y1:y2, x1:x2]-y0; sx = xx[y1:y2, x1:x2]-x0
        img[y1:y2, x1:x2] += np.exp(-(sy*sy+sx*sx)/(2*radius*radius))
    return np.clip(img, 0, 1)

def warp(ref, dx_px, dy_px):
    """Distorted image: background feature at (x,y) appears shifted by (dx,dy)."""
    from scipy.ndimage import map_coordinates
    H, W = ref.shape
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    d_x = np.nan_to_num(dx_px); d_y = np.nan_to_num(dy_px)
    return map_coordinates(ref, [yy-d_y, xx-d_x], order=1, mode='nearest').astype(np.float32)

def _subpix(c):
    """3-point Gaussian peak fit -> subpixel offset from the correlation peak."""
    j, i = np.unravel_index(np.argmax(c), c.shape)
    if not (0 < i < c.shape[1]-1 and 0 < j < c.shape[0]-1): return None
    def g(a, b, d):
        a, b, d = max(a, 1e-9), max(b, 1e-9), max(d, 1e-9)
        den = 2*(np.log(a)-2*np.log(b)+np.log(d))
        return 0.0 if abs(den) < 1e-12 else (np.log(a)-np.log(d))/den
    return (i + g(c[j, i-1], c[j, i], c[j, i+1]) - c.shape[1]//2,
            j + g(c[j-1, i], c[j, i], c[j+1, i]) - c.shape[0]//2)

def correlate(ref, dist, win=32, step=16, peak_ratio=1.15):
    """Windowed FFT cross-correlation, as a PIV/BOS processor does.

    Returns (u, v, x, y) on the vector grid. Windows whose correlation peak is not
    clearly above the second peak are rejected (NaN) -- the real dropout mechanism.
    """
    H, W = ref.shape
    ys = np.arange(0, H-win+1, step); xs = np.arange(0, W-win+1, step)
    u = np.full((len(ys), len(xs)), np.nan); v = np.full_like(u, np.nan)
    wnd = np.hanning(win)[:, None]*np.hanning(win)[None, :]
    # autocorrelation of the window: the per-shift effective overlap weight
    wcorr = np.fft.fftshift(np.real(np.fft.ifft2(np.abs(np.fft.fft2(wnd))**2)))
    wcorr = wcorr/wcorr.max()
    for jj, y0 in enumerate(ys):
        for ii, x0 in enumerate(xs):
            a = ref[y0:y0+win, x0:x0+win]*wnd
            b = dist[y0:y0+win, x0:x0+win]*wnd
            a = a-a.mean(); b = b-b.mean()
            if a.std() < 1e-6 or b.std() < 1e-6: continue
            c = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(b)*np.conj(np.fft.fft2(a)))))
            # Loss-of-pairs correction: the windowed correlation is biased toward zero
            # because the overlapping area shrinks with shift. Normalise by the
            # autocorrelation of the window function (standard PIV practice) -- without
            # this the recovered displacement is systematically LOW (~8% at 4-6 px).
            c = c/np.maximum(wcorr, 1e-6)
            pk = c.max()
            if pk <= 0: continue
            # second-peak test for validity (standard BOS/PIV dropout criterion)
            cc = c.copy(); j0, i0 = np.unravel_index(np.argmax(cc), cc.shape)
            cc[max(j0-2, 0):j0+3, max(i0-2, 0):i0+3] = 0
            if cc.max() > 0 and pk/cc.max() < peak_ratio: continue
            s = _subpix(c)
            if s is None: continue
            u[jj, ii], v[jj, ii] = s
    return u, v, xs+win/2, ys+win/2

def synthesise(dx_mm, dy_mm, mm_per_px, win=32, step=16, seed=0, density=0.06):
    """Full pipeline: true displacement (mm) -> synthetic BOS measurement.

    Returns dict with the reference/distorted images, the recovered displacement in mm
    on the vector grid, and the true displacement averaged over the same windows.
    """
    dx_px = np.nan_to_num(dx_mm)/mm_per_px
    dy_px = np.nan_to_num(dy_mm)/mm_per_px
    ref = speckle(dx_mm.shape, density=density, seed=seed)
    dist = warp(ref, dx_px, dy_px)
    u, v, xg, yg = correlate(ref, dist, win=win, step=step)
    # window-averaged truth, so recovered and true are compared like for like
    tu = np.full_like(u, np.nan); tv = np.full_like(u, np.nan)
    for jj, y0 in enumerate((yg-win/2).astype(int)):
        for ii, x0 in enumerate((xg-win/2).astype(int)):
            tu[jj, ii] = np.nanmean(dx_px[y0:y0+win, x0:x0+win])
            tv[jj, ii] = np.nanmean(dy_px[y0:y0+win, x0:x0+win])
    return dict(ref=ref, dist=dist,
                u_mm=u*mm_per_px, v_mm=v*mm_per_px,
                tu_mm=tu*mm_per_px, tv_mm=tv*mm_per_px,
                xg=xg, yg=yg,
                valid=float(np.isfinite(u).mean()))
