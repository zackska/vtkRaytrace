"""Ladder strip, cropped to the actual liquid extent.

The parent script (ladder_strip_latest.py) hard-codes set_xlim(0,70) mm, which was right
for the developed MPLIC runs but wastes ~95% of the canvas now that the PLIC rebuilds are
only at t~0.02 and have penetrated ~3 mm. Same physics, same gpuShadow surface-Snell
render -- only the framing changes.

A COMMON crop is used across all five panels so the comparison stays like-for-like: the
crop box is the union of every case's liquid bounding box plus a margin. Panels are laid
out as columns because the cropped aspect is tall, and aspect='equal' is kept so the
morphology is not distorted.
"""
import os, sys, re, struct, subprocess, time
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get('GPUSHADOW_DIR', os.path.join(HERE, '..', 'gpu')) + '/'
B = W + 'gpuShadow'
BASE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()   # dir of staged cases
CASES = ['We16_S0_L2', 'We18_S0_L2', 'We20_S0_L2', 'We22_S0_L2', 'We25_S0_L2']
LABELS = ['We=16\n(Liang: SWI)', 'We=18', 'We=20', 'We=22', 'We=25\n(Liang: bag)']
NX, NY, NZ = 600, 150, 150
RES = 1000; DN = '0.33'; KGD = '2.3e-4'
MARGIN_MM = 1.5


def fug(o):
    if o.IsA('vtkUnstructuredGrid') and o.GetNumberOfCells() > 0:
        return o
    if hasattr(o, 'GetNumberOfBlocks'):
        for i in range(o.GetNumberOfBlocks()):
            b = o.GetBlock(i)
            if b is not None and fug(b) is not None:
                return fug(b)
    return None


imgs, exts, tlabs, fracs = [], [], [], []
for C in CASES:
    d = f'{BASE}/{C}'; t0 = time.time()
    try:
        s = open(f'{d}/system/blockMeshDict').read()
        v = re.findall(r'\(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+)\)', s)
        LX = max(float(a) for a, b, c in v if 0 < float(a) < 1)
        DT = max(float(b) for a, b, c in v if 0 < float(b) < 1)
        ts = sorted(float(fn) for fn in os.listdir(d)
                    if re.match(r'^0\.[0-9]+$', fn) and os.path.exists(f'{d}/{fn}/alpha.water'))
        t = ts[-1]
        r = vtk.vtkPOpenFOAMReader(); r.SetFileName(f'{d}/case.foam'); r.SetCaseType(1)
        r.EnableAllCellArrays(); r.Update(); r.UpdateTimeStep(t); r.Update()
        ug = fug(r.GetOutput())
        c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(ug); c2p.Update()
        rs = vtk.vtkResampleToImage(); rs.SetInputDataObject(c2p.GetOutput()); rs.UseInputBoundsOff()
        rs.SetSamplingBounds(0, LX, 0, DT, 0, DT); rs.SetSamplingDimensions(NX, NY, NZ); rs.Update()
        al = np.clip(np.nan_to_num(vtk_to_numpy(
            rs.GetOutput().GetPointData().GetArray('alpha.water')).reshape(NZ, NY, NX)), 0, 1).astype(np.float32)
        fr = float((al > 0.5).mean())
        f = open(W + 'lz.bin', 'wb'); f.write(struct.pack('3i', NX, NY, NZ))
        f.write(struct.pack('6f', 0, LX, 0, DT, 0, DT)); al.tofile(f); f.close()
        subprocess.run([B, W + 'lz.bin', W + 'lzo.bin', 'sharp', str(RES), '8', DN, KGD, '300', '1.33', '0'],
                       check=True, stdout=subprocess.DEVNULL, timeout=200)
        fp = open(W + 'lzo.bin', 'rb'); rx, ry = struct.unpack('2i', fp.read(8))
        im = np.fromfile(fp, np.float32).reshape(ry, rx)
        imgs.append(im); exts.append((LX, DT)); tlabs.append(t); fracs.append(fr)
        print(f'  {C} t={t} frac={fr:.4f} {time.time()-t0:.0f}s', flush=True)
    except Exception as e:
        print(f'  {C} FAIL {e}', flush=True)
        imgs.append(None); exts.append((0.07, 0.024)); tlabs.append(None); fracs.append(None)

# ---- common crop box: union of liquid bounding boxes, in mm ----
x0s, x1s, y0s, y1s = [], [], [], []
for im, (LX, DT) in zip(imgs, exts):
    if im is None:
        continue
    dark = im < 0.5                       # shadowgraph: liquid blocks light
    if not dark.any():
        continue
    ys, xs = np.where(dark)
    x0s.append(xs.min() / im.shape[1] * LX * 1e3); x1s.append(xs.max() / im.shape[1] * LX * 1e3)
    y0s.append(ys.min() / im.shape[0] * DT * 1e3); y1s.append(ys.max() / im.shape[0] * DT * 1e3)
CX0 = max(0.0, min(x0s) - MARGIN_MM); CX1 = max(x1s) + MARGIN_MM
CY0 = max(0.0, min(y0s) - MARGIN_MM); CY1 = max(y1s) + MARGIN_MM
print(f'  common crop: x {CX0:.2f}-{CX1:.2f} mm, y {CY0:.2f}-{CY1:.2f} mm', flush=True)

fig, axs = plt.subplots(1, 5, figsize=(13, 6.2), dpi=150, constrained_layout=True)
for ax, im, (LX, DT), lab, tl, fr in zip(axs, imgs, exts, LABELS, tlabs, fracs):
    if im is not None:
        ax.imshow(im, origin='lower', extent=[0, LX * 1e3, 0, DT * 1e3],
                  aspect='equal', cmap='gray', vmin=0, vmax=1)
    ax.set_xlim(CX0, CX1); ax.set_ylim(CY0, CY1)
    sub = f'{lab}\nt = {tl:.4f} s' if tl is not None else f'{lab}\n(no frame)'
    if fr is not None:
        sub += f'\nliquid {fr*1e3:.2f}‰'
    ax.set_title(sub, fontsize=9)
    ax.set_xlabel('x (mm)', fontsize=8)
    ax.tick_params(labelsize=7)
axs[0].set_ylabel('y (mm)', fontsize=8)
for ax in axs[1:]:
    ax.set_yticklabels([])
fig.suptitle('SWI → bag regime ladder, S=0, PLIC rebuilds — gpuShadow surface-Snell (8°), '
             'common crop', fontsize=11)
fig.savefig(os.path.join(BASE, 'ladder_strip.png'), facecolor='white')
print('ZOOM STRIP DONE', flush=True)
