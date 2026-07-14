#!/usr/bin/env python3
"""
foam_shadowgraph.py — feed an OpenFOAM field to the gpuShadow renderer.

Pipeline:  OpenFOAM case  ->  structured grid.bin  ->  gpuShadow  ->  image
The grid step has two fidelity knobs (see FIDELITY NOTES at bottom):
  * method   : 'scatter' (cell-center binning; preserves detached droplets) or
               'resample' (vtkResampleToImage; smooth, but drops sub-threshold droplets)
  * presmooth: Gaussian sigma (cells) on the field before upload = the smoothed 0.5
               iso-surface + smoothed refraction normals on the GPU (0 = off).

Requires: vtk, numpy, scipy, and the compiled ./gpuShadow binary.
"""
import numpy as np, struct, subprocess, os
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.ndimage import distance_transform_edt, gaussian_filter

def _first_ug(o):
    if o.IsA('vtkUnstructuredGrid') and o.GetNumberOfCells() > 0:
        return o
    if hasattr(o, 'GetNumberOfBlocks'):
        for i in range(o.GetNumberOfBlocks()):
            b = o.GetBlock(i)
            if b is not None:
                r = _first_ug(b)
                if r is not None:
                    return r
    return None

def read_case(case_foam, time):
    """Return the unstructured grid at `time` (reconstructed case)."""
    r = vtk.vtkPOpenFOAMReader(); r.SetFileName(case_foam); r.SetCaseType(0)
    r.EnableAllCellArrays(); r.Update(); r.UpdateTimeStep(float(time)); r.Update()
    return _first_ug(r.GetOutput())

def field_to_grid(ug, field, bounds, dims, method='scatter', presmooth=0.0, cache=None):
    """
    Sample a cell field onto a uniform structured grid.
      bounds = (X0,X1,Y0,Y1,Z0,Z1);  dims = (NX,NY,NZ)
      method 'scatter' bins cells by centre (order-independent, keeps droplets);
             'resample' uses vtkResampleToImage (smooth interpolation).
      presmooth: Gaussian sigma in cells applied to the grid (smoothed iso-surface).
      cache: optional dict; on 'scatter' the bin indices (static mesh) are stored/reused.
    Returns float32 array shaped (NZ,NY,NX), C-order (k,j,i).
    """
    X0,X1,Y0,Y1,Z0,Z1 = bounds; NX,NY,NZ = dims
    if method == 'resample':
        c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(ug); c2p.Update()
        rs = vtk.vtkResampleToImage(); rs.SetInputDataObject(c2p.GetOutput()); rs.UseInputBoundsOff()
        rs.SetSamplingBounds(X0,X1,Y0,Y1,Z0,Z1); rs.SetSamplingDimensions(NX,NY,NZ); rs.Update()
        g = np.clip(vtk_to_numpy(rs.GetOutput().GetPointData().GetArray(field)).reshape(NZ,NY,NX),0,1).astype(np.float32)
    else:  # scatter
        if cache is not None and 'IJK' in cache:
            I,J,K,FI,sel = cache['IJK']
        else:
            cc = vtk.vtkCellCenters(); cc.SetInputData(ug); cc.Update()
            P = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
            sel = (P[:,0]>=X0)&(P[:,0]<X1)&(P[:,1]>=Y0)&(P[:,1]<Y1)&(P[:,2]>=Z0)&(P[:,2]<Z1)
            I = np.clip(((P[sel,0]-X0)/(X1-X0)*NX).astype(int),0,NX-1)
            J = np.clip(((P[sel,1]-Y0)/(Y1-Y0)*NY).astype(int),0,NY-1)
            K = np.clip(((P[sel,2]-Z0)/(Z1-Z0)*NZ).astype(int),0,NZ-1)
            cnt = np.zeros((NZ,NY,NX)); np.add.at(cnt,(K,J,I),1)
            FI = distance_transform_edt(~(cnt>0), return_distances=False, return_indices=True)  # nearest-fill graded gaps
            if cache is not None: cache['IJK'] = (I,J,K,FI,sel)
        vals = vtk_to_numpy(ug.GetCellData().GetArray(field))
        g = np.zeros((NZ,NY,NX),np.float32); g[K,J,I] = vals[sel]; g = g[tuple(FI)].astype(np.float32)
    if presmooth > 0:
        g = gaussian_filter(g, presmooth).astype(np.float32)
    return g

def write_grid(path, grid, bounds):
    NZ,NY,NX = grid.shape; X0,X1,Y0,Y1,Z0,Z1 = bounds
    with open(path,'wb') as f:
        f.write(struct.pack('3i',NX,NY,NZ)); f.write(struct.pack('6f',X0,X1,Y0,Y1,Z0,Z1)); grid.tofile(f)

def render(gpushadow, grid_bin, out_bin, mode='sharp', res=1600, acc=8.0, dn=0.33,
           K=0.0, absorb=80.0, nliq=1.33, outmode=0, extra=None):
    """Call the gpuShadow binary; return the image as (RESY,RESX) float32."""
    a = [gpushadow, grid_bin, out_bin, mode, str(res), str(acc), str(dn), str(K), str(absorb), str(nliq), str(outmode)]
    if extra: a += [str(x) for x in extra]
    subprocess.run(a, check=True)
    with open(out_bin,'rb') as f:
        rx,ry = struct.unpack('2i', f.read(8)); return np.fromfile(f, np.float32).reshape(ry,rx)

# ------------------------------------------------------------------ FIDELITY NOTES
# scatter vs resample : resample smoothly interpolates the mesh, so a 1-2 cell droplet
#   averages with surrounding gas and falls below alpha=0.5 -> it vanishes. scatter keeps
#   each cell's value -> droplets survive. Use scatter for atomization/breakup fidelity;
#   resample for a clean coherent-core figure. Both are fast at NATIVE grid resolution
#   (upsampling the grid is what was slow, and adds no real detail -- render the IMAGE at
#   high `res` instead; gpuShadow trilinear-upsamples the volume for free).
# presmooth : gpuShadow finds the alpha=0.5 surface implicitly (trilinear threshold crossing)
#   with the normal from the local gradient -- NOT a marching-cubes surface with averaged
#   normals like the CPU tracer. A light Gaussian (sigma~1 cell) on the grid restores that
#   smoothed iso-surface + smoothed normals; larger sigma smooths droplets away.
# schlieren : outmode 1/2 give ray deflection (mrad). For a plume with a strong local bloom,
#   map the deflection magnitude through log1p(mag/A) (not a linear clip) so the bloom's
#   internal structure reads instead of saturating.
