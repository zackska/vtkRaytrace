#!/usr/bin/env python3
"""BOS velocimetry on the helium jet: what does it actually measure?

Design note, after Zack's two corrections.

  1. My first proposal was to score the recovered field against a |grad rho|-weighted
     velocity. That weighting was an intuition ("gradients make the signal"), not a
     derivation. The transport equation disagrees: with n-1 = K.rho, the column mass
     m = int(rho dz) obeys  dm/dt + div(m U) = 0  with U = int(rho u dz)/int(rho dz),
     i.e. the MASS-weighted projected velocity. So |grad rho| was arbitrary.

  2. More decisively: BOTH weighted quantities need the full 3-D rho and u fields, which
     an experimentalist does not have. A reference you cannot compute is not a reference.
     So they are demoted here from "truth" to "explanation", and clearly labelled
     simulation-only.

The comparison is therefore against what an experiment CAN access:
     - centreplane velocity  (what planar PIV on the symmetry plane would give)
     - centreline velocity   (a probe / LDV on the axis)
     - bulk inlet velocity   (known from the flow rate)
and the deliverable is the BIAS between BOS-recovered velocity and those, since that is
the only correction an experimentalist could actually apply.

Pipeline: ray-trace each frame -> BOS deflection -> warp ONE fixed speckle target by each
frame's deflection -> cross-correlate consecutive frames -> apparent structure velocity.
That is exactly what a real rig does; nothing here uses the 3-D field except the labelled
diagnostic curves.
"""
import os, re, struct, subprocess, sys, glob
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bos_correlate import speckle, warp, correlate

W = os.environ.get('GPUSHADOW_DIR', HERE) + '/'
B = W + 'gpuShadow'
CASE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()   # OpenFOAM case dir
K_AIR, K_HE = 2.39e-4, 1.96e-4
NX, NY, NZ = 600, 170, 170
RES = 1000
LBG = 0.30
WIN, STEP = 32, 16
OUT = os.path.join(CASE, 'bos_velocimetry.png')


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


times = sorted([d for d in os.listdir(CASE)
                if re.match(r'^0\.030', d) and os.path.isfile(os.path.join(CASE, d, 'rho'))],
               key=float)
print(f'  {len(times)} frames, dt = {float(times[1])-float(times[0]):.3e} s', flush=True)

foam = os.path.join(CASE, 'case.foam'); open(foam, 'a').close()
rd = vtk.vtkPOpenFOAMReader(); rd.SetFileName(foam); rd.SetCaseType(1)
rd.EnableAllCellArrays(); rd.Update()

frames = []
for t in times:
    rd.UpdateTimeStep(float(t)); rd.Update()
    ug = find_ug(rd.GetOutput())
    b = ug.GetBounds()
    c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(ug); c2p.Update()
    rs = vtk.vtkResampleToImage()
    rs.SetInputDataObject(c2p.GetOutput()); rs.UseInputBoundsOff()
    rs.SetSamplingBounds(*b); rs.SetSamplingDimensions(NX, NY, NZ); rs.Update()
    pd = rs.GetOutput().GetPointData()
    msk = vtk_to_numpy(pd.GetArray('vtkValidPointMask')) > 0
    rho = np.nan_to_num(vtk_to_numpy(pd.GetArray('rho')))
    He = np.clip(np.nan_to_num(vtk_to_numpy(pd.GetArray('He'))), 0, 1)
    U = np.nan_to_num(vtk_to_numpy(pd.GetArray('U')))
    rho_amb = float(rho[msk].max()); rho[~msk] = rho_amb; He[~msk] = 0.0
    U[~msk] = 0.0
    rho = rho.reshape(NZ, NY, NX); He = He.reshape(NZ, NY, NX)
    Ux = U[:, 0].reshape(NZ, NY, NX)
    Kmix = He * K_HE + (1 - He) * K_AIR
    frames.append(dict(t=float(t), rho=rho, Ux=Ux, rho_eff=(rho * Kmix / K_AIR).astype(np.float32),
                       bounds=b))
    print(f'    t={t}  rho {rho.min():.3f}..{rho.max():.3f}  Ux max {Ux.max():.1f} m/s', flush=True)

b = frames[0]['bounds']
mm_per_px = (b[1] - b[0]) * 1e3 / RES


def deflect(vol):
    with open(W + 'vv.bin', 'wb') as f:
        f.write(struct.pack('3i', NX, NY, NZ))
        f.write(struct.pack('6f', *[float(x) for x in b]))
        vol.tofile(f)
    out = []
    for mode in (5, 6):
        p = subprocess.run([B, W + 'vv.bin', W + 'vv_o.bin', 'eikonal', str(RES), '0.9', '0.33',
                            f'{K_AIR:.6e}', '0', '1.0', str(mode), '1', '0', '0.003', '', '0',
                            '0', '2', '0.5', '1000', str(LBG)],
                           capture_output=True, text=True, timeout=900)
        if p.returncode != 0:
            raise RuntimeError(p.stderr[-300:])
        with open(W + 'vv_o.bin', 'rb') as fh:
            rx, ry = struct.unpack('2i', fh.read(8))
            out.append(np.fromfile(fh, np.float32).reshape(ry, rx))
    return out


# CORRECTION. The first version warped a FIXED target by each frame's deflection and
# correlated consecutive warped images. That is wrong: the target never moves, so the
# undistorted speckle dominates the correlation and it returns ~0 shift. What moves between
# frames is the DISTORTION, not the pattern carrying it — and the garbage it produced
# (mean beta ~ 8800) is what exposed it.
#
# Correct pipeline: the deflection field d(x,y,t) IS the image in which structures move.
# Recover it per instant (that is ordinary BOS), then correlate d between instants.
DX, DY = [], []
for i, fr in enumerate(frames):
    dx, dy = deflect(fr['rho_eff'])
    DX.append(np.nan_to_num(dx)); DY.append(np.nan_to_num(dy))
    print(f'    deflected frame {i+1}/{len(frames)}', flush=True)
DX = np.dstack(DX); DY = np.dstack(DY)
np.savez_compressed(os.path.join(CASE, 'deflection_series.npz'),
                    DX=DX, DY=DY, bounds=np.array(b), mm_per_px=mm_per_px,
                    times=np.array([f['t'] for f in frames]))
print(f'  saved deflection series {DX.shape}', flush=True)

dt = frames[1]['t'] - frames[0]['t']
us = []
for i in range(DX.shape[2] - 1):
    a = np.sqrt(DX[:, :, i] ** 2 + DY[:, :, i] ** 2)
    c = np.sqrt(DX[:, :, i+1] ** 2 + DY[:, :, i+1] ** 2)
    rng = np.nanpercentile(a[a > 0], 99) if np.any(a > 0) else 1.0
    u, v, xg, yg = correlate(np.clip(a / rng, 0, 1), np.clip(c / rng, 0, 1), win=WIN, step=STEP)
    us.append(u * mm_per_px * 1e-3 / dt)          # px -> m/s
US = np.dstack(us)
Ubos = np.nanmedian(US, axis=2)
print(f'  correlated {len(us)} deflection-field pairs, dt = {dt:.3e} s', flush=True)

# ---- accessible references ----
mid = NZ // 2
rho0 = frames[0]['rho']; Ux0 = frames[0]['Ux']
x_mm = np.linspace(b[0], b[1], NX) * 1e3
U_centreplane = Ux0[mid, :, :]                     # planar PIV on the symmetry plane
U_centreline = Ux0[mid, NY // 2, :]                # probe on the axis
# ---- simulation-only diagnostics (cannot be measured) ----
U_rho = (rho0 * Ux0).sum(axis=0) / rho0.sum(axis=0)
g = np.gradient(rho0, axis=2)
wgt = np.abs(g) + 1e-12
U_gradrho = (wgt * Ux0).sum(axis=0) / wgt.sum(axis=0)

np.savez_compressed(os.path.join(CASE, 'velo_result.npz'),
                    Ubos=Ubos, US=US, xg=xg, yg=yg, x_mm=x_mm, dt=dt, mm_per_px=mm_per_px,
                    U_centreplane=U_centreplane, U_centreline=U_centreline,
                    U_rho=U_rho, U_gradrho=U_gradrho, bounds=np.array(b))
print('  wrote velo_result.npz', flush=True)
