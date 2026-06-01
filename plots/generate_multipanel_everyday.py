#!/usr/bin/env python3
"""
Everyday Materials Multipanel Fermi Surface Generator
=====================================================
Generates 4x2 and 2x4 grid images showing the evolution of the 
Fermi surface across 8 different energy levels for Copper, Iron, and Graphite.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import Normalize
import skimage.measure

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'Images', 'Everyday_Materials_Multipanel')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Material Definitions ───────────────────────────────────────────────────
MATERIALS = {
    'copper_gold_silver': {
        'name': 'Copper / Gold / Silver (FCC)',
        'dispersion': lambda kx, ky, kz: -4.0 * (
            np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky) + 
            np.cos(2*np.pi*ky)*np.cos(2*np.pi*kz) + 
            np.cos(2*np.pi*kz)*np.cos(2*np.pi*kx)
        ),
        'cmap': 'managua_r',
        # Fixed range requested by user
        'sweep_energies': np.linspace(-3.0, 3.0, 8)
    },
    'iron_tungsten': {
        'name': 'Iron / Tungsten (BCC)',
        'dispersion': lambda kx, ky, kz: -8.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi*ky) * np.cos(2*np.pi*kz),
        'cmap': 'cividis',
        # Fixed range requested by user
        'sweep_energies': np.linspace(-3.0, 3.0, 8)
    },
    'graphite': {
        'name': 'Graphite (Layered Hexagonal)',
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(4*np.pi*kx) + 2.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi * np.sqrt(3) * ky))
            - 2.0 * 0.5 * np.cos(2*np.pi*kz)
        ),
        'cmap': 'bone',
        # Fixed range requested by user
        'sweep_energies': np.linspace(-3.0, 3.0, 8)
    }
}

# ── Color Mapping (Flat / Unshaded) ────────────────────────────────────────
def compute_numerical_gradient_magnitude(dispersion_fn, kx, ky, kz, dk=1e-5):
    dE_dx = (dispersion_fn(kx + dk, ky, kz) - dispersion_fn(kx - dk, ky, kz)) / (2 * dk)
    dE_dy = (dispersion_fn(kx, ky + dk, kz) - dispersion_fn(kx, ky - dk, kz)) / (2 * dk)
    dE_dz = (dispersion_fn(kx, ky, kz + dk) - dispersion_fn(kx, ky, kz - dk)) / (2 * dk)
    mag = np.sqrt(dE_dx**2 + dE_dy**2 + dE_dz**2)
    mag[mag < 1e-8] = 1e-8
    return mag

def compute_flat_colors(velocities, cmap_name='inferno'):
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=velocities.min(), vmax=velocities.max())
    return cmap(norm(velocities)).copy()

def render_grid(mat_id, mat_info, meshes, rows, cols):
    BG = '#07070f'
    # Figsize proportional to grid dimensions
    fig = plt.figure(figsize=(cols * 4, rows * 4), facecolor=BG)
    
    for i, (verts_tri, EF) in enumerate(meshes):
        ax = fig.add_subplot(rows, cols, i + 1, projection='3d', facecolor=BG)
        ax.set_proj_type('persp')
        
        if verts_tri is not None:
            centroids = verts_tri.mean(axis=1)
            velocities = compute_numerical_gradient_magnitude(
                mat_info['dispersion'], centroids[:, 0], centroids[:, 1], centroids[:, 2]
            )
            face_colors = compute_flat_colors(velocities, mat_info['cmap'])
            poly = Poly3DCollection(verts_tri.tolist(), zsort='average', 
                                    facecolors=face_colors, edgecolors='none')
            ax.add_collection3d(poly)
            
        ax.set_box_aspect([1, 1, 1])
        span = 0.5
        ax.set_xlim(-span, span)
        ax.set_ylim(-span, span)
        ax.set_zlim(-span, span)
        ax.set_axis_off()
        ax.view_init(elev=25, azim=45)
        
        ax.set_title(f"E = {EF:+.2f} eV", color='#ffaa66', fontsize=14, pad=0)
        
    fig.suptitle(f"{mat_info['name']} - Energy Evolution", color='#ffffff', fontsize=22, y=0.98)
    plt.tight_layout()
    # Adjust top to make room for the suptitle
    plt.subplots_adjust(top=0.92) 
    
    out_path = os.path.join(OUTPUT_DIR, f"multipanel_{mat_id}_{rows}x{cols}.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print(f"  Saved {rows}x{cols} grid to {out_path}")

# ── Generation ─────────────────────────────────────────────────────────────
def main():
    print(f"Saving multipanel images to {OUTPUT_DIR}")
    npts_mesh = 50  # Optimized resolution for subplots (fast and clean)
    kx = np.linspace(0, 1, npts_mesh + 1)
    ky = np.linspace(0, 1, npts_mesh + 1)
    kz = np.linspace(0, 1, npts_mesh + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    for mat_id, mat_info in MATERIALS.items():
        print(f"\nProcessing {mat_info['name']}...")
        E_mesh = mat_info['dispersion'](KX - 0.5, KY - 0.5, KZ - 0.5)
        
        # 1. Precompute all 8 meshes
        meshes = []
        for EF in mat_info['sweep_energies']:
            grid = E_mesh - EF
            try:
                verts, faces, _, _ = skimage.measure.marching_cubes(
                    grid, level=0.0, spacing=(1.0 / npts_mesh, 1.0 / npts_mesh, 1.0 / npts_mesh)
                )
                verts -= 0.5  # Center at origin
                verts_tri = verts[faces]
                meshes.append((verts_tri, EF))
            except (ValueError, RuntimeError):
                # Empty surface (e.g. energy outside band)
                meshes.append((None, EF))
                
        # 2. Render grids
        render_grid(mat_id, mat_info, meshes, rows=4, cols=2)
        render_grid(mat_id, mat_info, meshes, rows=2, cols=4)

if __name__ == "__main__":
    main()
