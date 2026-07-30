#!/usr/bin/env python3
"""Centre-plane density and helium mass fraction from the completed v2 LES at t = 0.150 s.

Reads the DECOMPOSED case directly (vtkPOpenFOAMReader CaseType 0), so no reconstructPar and
no serial polyMesh are needed -- which is why the S3 pull could skip the 500 MB serial mesh.

The plane is extracted by resampling with SamplingDimensions (NX, NY, 1) and zmin = zmax =
z_mid, which probes a single plane rather than building a full 3-D image and slicing it.

Spatial panels use aspect='equal'.
"""
import os
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
T = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
NX, NY = 1200, 520
OUT = os.path.join(CASE, 'density_plane.png')

foam = os.path.join(CASE, 'case.foam')
open(foam, 'a').close()

rd = vtk.vtkPOpenFOAMReader()
rd.SetFileName(foam)
rd.SetCaseType(0)                    # 0 = DECOMPOSED_CASE
rd.EnableAllCellArrays()
rd.Update()
rd.UpdateTimeStep(T)
rd.Update()


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


ug = find_ug(rd.GetOutput())
if ug is None:
    raise SystemExit('no unstructured grid found — check CaseType / time value')
b = ug.GetBounds()
ncells = ug.GetNumberOfCells()
zmid = 0.5 * (b[4] + b[5])
print(f'cells {ncells:,}   bounds x {b[0]:.4f}..{b[1]:.4f}  y {b[2]:.4f}..{b[3]:.4f}  '
      f'z {b[4]:.4f}..{b[5]:.4f} m')
print(f'centre plane at z = {zmid:.5f} m')

c2p = vtk.vtkCellDataToPointData()
c2p.SetInputData(ug)
c2p.Update()

rs = vtk.vtkResampleToImage()
rs.SetInputDataObject(c2p.GetOutput())
rs.UseInputBoundsOff()
rs.SetSamplingBounds(b[0], b[1], b[2], b[3], zmid, zmid)
rs.SetSamplingDimensions(NX, NY, 1)
rs.Update()
pd = rs.GetOutput().GetPointData()

msk = vtk_to_numpy(pd.GetArray('vtkValidPointMask')).reshape(NY, NX) > 0
rho = np.nan_to_num(vtk_to_numpy(pd.GetArray('rho'))).reshape(NY, NX)
He = np.clip(np.nan_to_num(vtk_to_numpy(pd.GetArray('He'))).reshape(NY, NX), 0, 1)
rho[~msk] = np.nan
He[~msk] = np.nan

ext = [b[0] * 1e3, b[1] * 1e3, b[2] * 1e3, b[3] * 1e3]   # mm
print(f'rho  {np.nanmin(rho):.4f} .. {np.nanmax(rho):.4f} kg/m3')
print(f'He   {np.nanmin(He):.4f} .. {np.nanmax(He):.4f}')

fig, axs = plt.subplots(2, 1, figsize=(12.6, 8.4), dpi=150)
for ax, F, cmap, lab, ttl, vmin, vmax in [
        (axs[0], rho, 'viridis', r'$\rho$  [kg/m$^3$]',
         'centre-plane density — the field that bends the light',
         float(np.nanpercentile(rho, 0.2)), float(np.nanpercentile(rho, 99.8))),
        (axs[1], He, 'magma', r'$Y_{He}$  [-]',
         'centre-plane helium mass fraction — the composition driving it', 0.0, 1.0)]:
    im = ax.imshow(F, origin='lower', extent=ext, aspect='equal', cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_title(ttl, fontsize=10.5)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='1.8%', pad=0.10)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(lab, fontsize=9)
    cb.ax.tick_params(labelsize=8)

fig.suptitle(f'Helium jet into air — v2 LES (WALE, nonUnityLewisEddyDiffusivity Sct 0.7), '
             f'{ncells/1e6:.2f} M cells, t = {T} s', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(OUT, facecolor='white', bbox_inches='tight')
print(f'wrote {OUT}')

# persist the plane so the design-envelope figure can mount it as its first panel
npz = os.path.splitext(OUT)[0] + '.npz'
np.savez_compressed(npz, rho=rho.astype(np.float32), He=He.astype(np.float32),
                    extent_mm=np.array(ext), t=T, ncells=ncells)
print(f'wrote {npz}')
