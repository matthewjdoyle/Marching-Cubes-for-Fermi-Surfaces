#!/usr/bin/env python3
"""
Real-Life Materials Fermi Surface Generator
===========================================
Plots static 3D isosurfaces for specific, physical energy levels 
of real-world materials and topological models.
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
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'Images', 'Materials')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Material Definitions ───────────────────────────────────────────────────
# t = 1.0 is used as the base energy scale (eV) for the cuprate and ruthenate.
MATERIALS = {
    'ybco_cuprate': {
        'name': 'YBCO-like Cuprate (High-Tc SC)',
        # t = 1.0, t' = -0.3, t'' = 0.1, tz = 0.05
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky))
            - 4.0 * (-0.3) * np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky)
            - 2.0 * 0.1 * (np.cos(4*np.pi*kx) + np.cos(4*np.pi*ky))
            - 2.0 * 0.05 * np.cos(2*np.pi*kz)
        ),
        'cmap': 'plasma',
        'energies': [-1.70, -1.60, -1.40] # -1.6 is the Van Hove singularity (Lifshitz transition)
    },
    'sr2ruo4_gamma': {
        'name': 'Sr2RuO4 Gamma Band (p-wave SC candidate)',
        # t = 1.0, t' = 0.4, tz = 0.05
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky))
            - 4.0 * 0.4 * np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky)
            - 2.0 * 0.05 * np.cos(2*np.pi*kz)
        ),
        'cmap': 'viridis',
        'energies': [1.40, 1.60, 1.80] # 1.6 is the Van Hove singularity
    },
    'weyl_semimetal': {
        'name': 'Minimal Weyl Semimetal (Topological)',
        # Conduction band of a minimal model with Weyl nodes. M = 2.0
        'dispersion': lambda kx, ky, kz: np.sqrt(
            np.sin(2*np.pi*kx)**2 + np.sin(2*np.pi*ky)**2 + 
            (2.0 - np.cos(2*np.pi*kx) - np.cos(2*np.pi*ky) - np.cos(2*np.pi*kz))**2
        ),
        'cmap': 'magma',
        'energies': [0.15, 0.70, 1.50] # 0.0 are the Weyl nodes, pockets grow and merge as E increases
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

# ── Generation ─────────────────────────────────────────────────────────────
def main():
    print(f"Saving material static images to {OUTPUT_DIR}")
    npts_mesh = 120  # Very high resolution for clean Fermi arcs and pockets
    kx = np.linspace(0, 1, npts_mesh + 1)
    ky = np.linspace(0, 1, npts_mesh + 1)
    kz = np.linspace(0, 1, npts_mesh + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    for mat_id, mat_info in MATERIALS.items():
        print(f"\nProcessing {mat_info['name']}...")
        E_mesh = mat_info['dispersion'](KX - 0.5, KY - 0.5, KZ - 0.5)
        
        for EF in mat_info['energies']:
            print(f"  Generating plot for E = {EF} eV...")
            grid = E_mesh - EF
            try:
                verts, faces, _, _ = skimage.measure.marching_cubes(
                    grid, level=0.0, spacing=(1.0 / npts_mesh, 1.0 / npts_mesh, 1.0 / npts_mesh)
                )
                verts -= 0.5  # Center at origin
                verts_tri = verts[faces]
                
                centroids = verts_tri.mean(axis=1)
                velocities = compute_numerical_gradient_magnitude(
                    mat_info['dispersion'], centroids[:, 0], centroids[:, 1], centroids[:, 2]
                )
                face_colors = compute_flat_colors(velocities, mat_info['cmap'])

                BG = '#07070f'
                fig = plt.figure(figsize=(9, 9), facecolor=BG)
                ax = fig.add_subplot(111, projection='3d', facecolor=BG)
                ax.set_proj_type('persp')

                poly = Poly3DCollection(verts_tri.tolist(), zsort='average', 
                                        facecolors=face_colors, edgecolors='none')
                ax.add_collection3d(poly)
                
                ax.set_box_aspect([1, 1, 1])
                span = 0.5
                ax.set_xlim(-span, span)
                ax.set_ylim(-span, span)
                ax.set_zlim(-span, span)
                ax.set_axis_off()
                
                # Use a nice isometric view
                ax.view_init(elev=25, azim=45)
                
                fig.text(0.5, 0.92, f"Fermi Energy: {EF:+.2f} eV", 
                         color='#ffaa66', fontsize=16, ha='center', fontweight='bold')
                
                # Title text for the material
                fig.text(0.5, 0.88, f"{mat_info['name']}", 
                         color='#aaaaaa', fontsize=12, ha='center')
                
                out_path = os.path.join(OUTPUT_DIR, f"material_{mat_id}_{EF:+.2f}eV.png")
                plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor=BG)
                plt.close(fig)
                
            except (ValueError, RuntimeError) as e:
                print(f"    Failed to generate surface at E = {EF} eV: {e}")

if __name__ == "__main__":
    main()
