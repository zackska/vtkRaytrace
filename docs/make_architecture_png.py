#!/usr/bin/env python3
"""Render docs/architecture.png — a slide-usable copy of the Mermaid diagram in
docs/ARCHITECTURE.md.

The Mermaid block in ARCHITECTURE.md is the SOURCE OF TRUTH: it renders natively on GitHub
and diffs cleanly. This script exists only because a raster copy is handy for talks, and
because no mermaid CLI is installed on the build box. If you change the pipeline, change both.

Dashed edges/boxes mark paths that are SPECIFIED BUT NOT YET IMPLEMENTED (see the
"Implemented vs specified" table in ARCHITECTURE.md) — currently the experimental-image
branch into the correlator, and the density inversion.

    python3 docs/make_architecture_png.py

Layout note: grouping is carried by box FILL COLOUR plus the legend, not by swimlane
rectangles -- band labels kept colliding with the vertical spine. Arrows never cross a box.
Look at the render after editing.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK = '#1b1d1e'
GREY = '#8a8f94'
SUB = '#41474c'
LBL = '#5a6066'
C_IO = '#e8eef2'
C_PHYS = '#dfe8e0'
C_PROC = '#e6e2ee'
C_VAL = '#f0e4e4'
C_OFF = '#f1f1f1'

fig, ax = plt.subplots(figsize=(12.6, 12.8), dpi=150)
ax.set_xlim(0, 100)
ax.set_ylim(0, 126)
ax.axis('off')

SPINE_X, SPINE_W = 33.0, 32.0
CX = SPINE_X + SPINE_W / 2


def box(x, y, w, h, title, body='', fc=C_IO, ts=9.6, bs=7.3, mono=False, dim=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.0',
                                linewidth=1.0, edgecolor=GREY, facecolor=fc, zorder=2,
                                linestyle=(0, (4, 3)) if dim else '-'))
    tc = '#7b8085' if dim else INK
    bc = '#8d9296' if dim else SUB
    if body:
        ax.text(x + w / 2, y + h * 0.70, title, ha='center', va='center', zorder=4,
                fontsize=ts, color=tc, fontweight='bold')
        ax.text(x + w / 2, y + h * 0.29, body, ha='center', va='center', zorder=4,
                fontsize=bs, color=bc, linespacing=1.55,
                family='monospace' if mono else 'sans-serif')
    else:
        ax.text(x + w / 2, y + h / 2, title, ha='center', va='center', zorder=4,
                fontsize=ts, color=tc, fontweight='bold')


def arrow(x1, y1, x2, y2, dashed=False, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=11,
                                 linewidth=1.0, color='#5f6469', zorder=3,
                                 linestyle=(0, (4, 3)) if dashed else '-',
                                 connectionstyle=f'arc3,rad={rad}'))


# ------------------------------------------------------------------- spine
box(SPINE_X, 116.0, SPINE_W, 7.0, 'LES field', 'rho, T, Y_i', C_IO, mono=True)
box(SPINE_X, 104.5, SPINE_W, 8.6, 'Refractive index',
    'Lorentz-Lorenz\nper-species refractivity', C_PHYS)
box(SPINE_X, 91.0, SPINE_W, 10.0, 'Interpolant',
    'tricubic, C1-continuous\n(trilinear breaks grad n)', C_PHYS)
box(SPINE_X, 77.5, SPINE_W, 10.0, 'Ray integrator',
    'd/dz(n dx/dz) = dn/dx\nRK4, backward from pixel', C_PHYS)
box(SPINE_X, 66.0, SPINE_W, 8.0, 'Aperture sampling', 'Monte Carlo, N rays/pixel', C_PHYS)
box(SPINE_X, 54.0, SPINE_W, 8.6, 'Background-plane', 'intersection + sampling', C_PHYS)
box(SPINE_X, 42.5, SPINE_W, 8.0, 'Synthetic image pair', 'flow-on / flow-off', C_IO)
box(SPINE_X, 31.0, SPINE_W, 8.0, 'Cross-correlation', 'PIV-style windows', C_PROC)
box(SPINE_X, 19.5, SPINE_W, 8.0, 'Displacement field', 'dx, dy', C_PROC, mono=True)
box(SPINE_X, 8.0, SPINE_W, 8.0, 'Comparison layer',
    'vs. measured BOS displacement', C_VAL, bs=7.0)

for y1, y2 in [(116.0, 113.4), (104.5, 101.2), (91.0, 87.8), (77.5, 74.2),
               (66.0, 62.9), (54.0, 50.7), (42.5, 39.2), (31.0, 27.8), (19.5, 16.2)]:
    arrow(CX, y1, CX, y2)

# ------------------------------------------------------------------- side inputs
box(74, 77.5, 24, 10.0, 'Optical config', 'M, Z_D, aperture, f', C_IO, mono=True)
arrow(73.8, 82.5, 65.3, 82.5)

box(74, 54.0, 24, 8.6, 'Background pattern', 'speckle / dot array', C_IO)
arrow(73.8, 58.3, 65.3, 58.3)

box(74, 31.0, 24, 8.0, 'Experimental BOS', 'images', C_OFF, dim=True)
arrow(73.8, 35.0, 65.3, 35.0, dashed=True)

# ------------------------------------------------------------------- side output
box(2, 19.5, 25, 9.6, 'Density inversion', 'Poisson, ill-posed\ndiagnostic only',
    C_OFF, dim=True, bs=6.9)
arrow(32.8, 23.5, 27.3, 23.9, dashed=True)

# ------------------------------------------------------------------- legend
LEG = [('Input / configuration', C_IO), ('Physics / optics', C_PHYS),
       ('Processing', C_PROC), ('Validation', C_VAL), ('Not yet implemented', C_OFF)]
lx = 2.0
for name, col in LEG:
    ax.add_patch(FancyBboxPatch((lx, 3.4), 2.6, 2.2,
                                boxstyle='round,pad=0.15,rounding_size=0.5',
                                linewidth=0.9, edgecolor=GREY, facecolor=col,
                                linestyle=(0, (3, 2)) if col == C_OFF else '-', zorder=2))
    ax.text(lx + 3.4, 4.5, name, ha='left', va='center', fontsize=7.4, color=SUB)
    lx += 3.4 + len(name) * 1.02 + 3.2

ax.text(2, 0.9, 'Dashed = specified but not yet implemented: the experimental-image branch into '
                'the correlator (so nothing is yet compared against a measurement) and the '
                'density inversion.',
        ha='left', va='center', fontsize=7.3, color=LBL, style='italic')
ax.text(98, 122.5, 'cfd-optical-diagnostics', ha='right', va='center',
        fontsize=8.0, color=GREY, fontweight='bold')
ax.text(98, 120.6, 'docs/ARCHITECTURE.md is the source of truth', ha='right', va='center',
        fontsize=6.8, color=GREY)
ax.text(98, 118.2, 'shipped interpolant is trilinear, not tricubic —', ha='right',
        va='center', fontsize=6.8, color='#a08a8a')
ax.text(98, 116.6, 'see "Implemented vs specified"', ha='right', va='center',
        fontsize=6.8, color='#a08a8a')

fig.savefig('/home/zack/projects/vtkRaytrace/docs/architecture.png',
            facecolor='white', bbox_inches='tight')
print('wrote docs/architecture.png')
