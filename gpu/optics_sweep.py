#!/usr/bin/env python3
"""Effect of finite aperture / defocused flow on BOS accuracy.

Real BOS focuses on the BACKGROUND, so the density object is deliberately out of focus
(Raffel 2015). Every ray reaching a given pixel converges on the same background point but
crosses the flow at a DIFFERENT transverse position, so the recorded deflection is an
average over a disc at the object plane. Everything we rendered until now used
naper=1, apR=0 — a pinhole, all-in-focus — which makes our recovered fields sharper than
any real rig could achieve.

gpuShadow already supports this properly: `naper` rays per pixel converging at focal plane
`zf` with aperture radius `apR`. Here zf is placed at the background plane and apR is set
from a stated lens geometry, so the defocus is traced rather than approximated by blurring.

Geometry (stated so the numbers are reproducible):
    lens-to-object      Z_A = 1.00 m
    object-to-background Z_D = 0.30 m   (= L_bg used throughout)
    focal length         f   = 50 mm
    aperture diameter    d_A = f / N
    blur disc at object  D_obj = d_A * Z_D / (Z_A + Z_D)      [Raffel Eq. 6 geometry]

In the tracer's coordinates rays start at z=0 over radius apR and converge at zf, so the
bundle radius at the jet mid-plane z_j is apR*(1 - z_j/zf). apR is chosen to reproduce D_obj.
"""
import os, struct, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bos_correlate import speckle, warp, correlate
from bos_realistic_target import printed_target, camera

W = os.environ.get('GPUSHADOW_DIR', HERE) + '/'
B = W + 'gpuShadow'
CASE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()   # OpenFOAM case dir
K_AIR, K_HE = 2.39e-4, 1.96e-4
NX, NY, NZ = 600, 170, 170
RES = 1000
LBG = 0.30
WIN, STEP = 32, 16
Z_A, F_LENS = 1.00, 0.050
NAPER = 64

import vtk
from vtk.util.numpy_support import vtk_to_numpy


def find_ug(o):
    if o.IsA('vtkUnstructuredGrid') and o.GetNumberOfCells() > 0:
        return o
    if hasattr(o, 'GetNumberOfBlocks'):
        for i in range(o.GetNumberOfBlocks()):
            b = o.GetBlock(i)
            if b is not None:
                r = find_ug(b)
                if r is not None:
                    return r
    return None


# one representative frame is enough for an optics sweep
t = sorted([d for d in os.listdir(CASE) if d.startswith('0.030')
            and os.path.isfile(os.path.join(CASE, d, 'rho'))], key=float)[10]
foam = os.path.join(CASE, 'case.foam'); open(foam, 'a').close()
rd = vtk.vtkPOpenFOAMReader(); rd.SetFileName(foam); rd.SetCaseType(1)
rd.EnableAllCellArrays(); rd.Update(); rd.UpdateTimeStep(float(t)); rd.Update()
ug = find_ug(rd.GetOutput()); b = ug.GetBounds()
c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(ug); c2p.Update()
rs = vtk.vtkResampleToImage(); rs.SetInputDataObject(c2p.GetOutput()); rs.UseInputBoundsOff()
rs.SetSamplingBounds(*b); rs.SetSamplingDimensions(NX, NY, NZ); rs.Update()
pd = rs.GetOutput().GetPointData()
msk = vtk_to_numpy(pd.GetArray('vtkValidPointMask')) > 0
rho = np.nan_to_num(vtk_to_numpy(pd.GetArray('rho')))
He = np.clip(np.nan_to_num(vtk_to_numpy(pd.GetArray('He'))), 0, 1)
rho[~msk] = float(rho[msk].max()); He[~msk] = 0.0
rho = rho.reshape(NZ, NY, NX); He = He.reshape(NZ, NY, NX)
Kmix = He * K_HE + (1 - He) * K_AIR
vol = (rho * Kmix / K_AIR).astype(np.float32)
mm_per_px = (b[1] - b[0]) * 1e3 / RES
z_jet = 0.5 * (b[4] + b[5])
z_bg = b[5] + LBG
print(f'  frame t={t}   mm/px {mm_per_px:.4f}   jet mid-plane z={z_jet:.4f} m   '
      f'background z={z_bg:.4f} m', flush=True)

with open(W + 'op.bin', 'wb') as f:
    f.write(struct.pack('3i', NX, NY, NZ))
    f.write(struct.pack('6f', *[float(x) for x in b]))
    vol.tofile(f)


def deflect(naper, apR):
    out = []
    for mode in (5, 6):
        p = subprocess.run([B, W + 'op.bin', W + 'op_o.bin', 'eikonal', str(RES), '0.9', '0.33',
                            f'{K_AIR:.6e}', '0', '1.0', str(mode), str(naper), f'{apR:.6e}',
                            f'{z_bg:.6f}', '', '0', '0', '2', '0.5', '1000', str(LBG)],
                           capture_output=True, text=True, timeout=1800)
        if p.returncode != 0:
            raise RuntimeError(p.stderr[-300:])
        with open(W + 'op_o.bin', 'rb') as fh:
            rx, ry = struct.unpack('2i', fh.read(8))
            out.append(np.nan_to_num(np.fromfile(fh, np.float32).reshape(ry, rx)))
    return out


MAPS = {}
CASES = [('pinhole (what we used)', None), ('f/32', 32.0), ('f/16', 16.0),
         ('f/8', 8.0), ('f/2.8', 2.8)]
ref_t = printed_target((int(RES * (b[3] - b[2]) / (b[1] - b[0])), RES), density=0.30, seed=7)
results = []
truth = None
for lab, N in CASES:
    if N is None:
        apR, D_obj = 0.0, 0.0
        dx, dy = deflect(1, 0.0)
    else:
        dA = F_LENS / N
        D_obj = dA * LBG / (Z_A + LBG)                     # blur disc at the object [m]
        apR = 0.5 * D_obj / (1.0 - z_jet / z_bg)           # tracer start-disc radius
        dx, dy = deflect(NAPER, apR)
    mag = np.sqrt(dx ** 2 + dy ** 2)
    if truth is None:
        truth = mag.copy()                                  # pinhole = the underlying field
    # image it and recover, with the full camera model
    warped = warp(ref_t, dx / mm_per_px, dy / mm_per_px)
    im0 = camera(ref_t, 0.8, 8000.0, 12.0, 8, np.random.default_rng(11))
    im1 = camera(warped, 0.8, 8000.0, 12.0, 8, np.random.default_rng(12))
    u, v, xg, yg = correlate(im0, im1, win=WIN, step=STEP)
    rec = np.sqrt(u ** 2 + v ** 2) * mm_per_px
    # window-average the true field for a like-for-like score
    tw = np.full_like(rec, np.nan)
    for jj, y0 in enumerate((yg - WIN / 2).astype(int)):
        for ii, x0 in enumerate((xg - WIN / 2).astype(int)):
            tw[jj, ii] = np.nanmean(truth[y0:y0 + WIN, x0:x0 + WIN])
    ok = np.isfinite(rec) & (tw < WIN / 4 * mm_per_px)
    rms = float(np.sqrt(np.nanmean((rec[ok] - tw[ok]) ** 2)))
    # loss of fine structure: high-frequency energy retained vs pinhole
    hf = float(np.nanstd(np.diff(mag, axis=1)) / np.nanstd(np.diff(truth, axis=1)))
    results.append((lab, N, D_obj * 1e3, D_obj * 1e3 / mm_per_px, 100 * np.isfinite(u).mean(),
                    rms, float(mag.max()), hf))
    MAPS[lab] = mag
    print(f'  {lab:22} blur {D_obj*1e3:6.3f} mm ({D_obj*1e3/mm_per_px:5.1f} px)  '
          f'valid {100*np.isfinite(u).mean():5.1f}%  RMS {rms:.4f} mm  '
          f'peak|d| {mag.max():.3f} mm  HF retained {hf:.3f}', flush=True)

np.savez_compressed(os.path.join(CASE, 'optics_sweep.npz'),
         results=np.array(results, dtype=object),
         bounds=np.array(b), mm_per_px=mm_per_px, win=WIN,
         **{f'map{i}': MAPS[l] for i, (l, _) in enumerate(CASES)})
print('\n  saved optics_sweep.npz')
