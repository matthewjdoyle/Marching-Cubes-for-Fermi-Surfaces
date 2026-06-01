#!/usr/bin/env python3
"""
Static Fermi Surface Generator
==============================
Plots static 3D isosurfaces for specific, interesting energy levels 
(Van Hove singularities, topological transitions) for each lattice.
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
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'Images', 'Static')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Lattice Definitions ────────────────────────────────────────────────────
LATTICES = {
    'sc': {
        'name': 'Simple Cubic (SC)',
        'dispersion': lambda kx, ky, kz: -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz)),
        'cmap': 'inferno',
        'energies': [-2.0, 0.0, 2.0]
    },
    'sc_nnn': {
        'name': 'Simple Cubic with NNN',
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz))
            - 4.0 * (-0.2) * (np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky) + 
                              np.cos(2*np.pi*ky)*np.cos(2*np.pi*kz) + 
                              np.cos(2*np.pi*kz)*np.cos(2*np.pi*kx))
        ),
        'cmap': 'plasma',
        'energies': [-1.5, 0.4, 4.0]
    },
    'sc_anisotropic': {
        'name': 'Anisotropic Simple Cubic',
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky)) - 2.0 * 0.15 * np.cos(2*np.pi*kz)
        ),
        'cmap': 'cividis',
        'energies': [-2.0, 0.0, 2.0]
    },
    'bcc': {
        'name': 'Body-Centered Cubic (BCC)',
        'dispersion': lambda kx, ky, kz: -8.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi*ky) * np.cos(2*np.pi*kz),
        'cmap': 'viridis',
        'energies': [-4.0, 0.0, 4.0]
    },
    'fcc': {
        'name': 'Face-Centered Cubic (FCC)',
        'dispersion': lambda kx, ky, kz: -4.0 * (
            np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky) + 
            np.cos(2*np.pi*ky)*np.cos(2*np.pi*kz) + 
            np.cos(2*np.pi*kz)*np.cos(2*np.pi*kx)
        ),
        'cmap': 'magma',
        'energies': [-8.0, -4.0, 0.0]
    },
    'hexagonal': {
        'name': 'Simple Hexagonal',
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(4*np.pi*kx) + 2.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi * np.sqrt(3) * ky))
            - 2.0 * 0.5 * np.cos(2*np.pi*kz)
        ),
        'cmap': 'coolwarm',
        'energies': [-3.0, -1.0, 1.0]
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
    print(f"Saving static images to {OUTPUT_DIR}")
    npts_mesh = 100  # High resolution for static images
    kx = np.linspace(0, 1, npts_mesh + 1)
    ky = np.linspace(0, 1, npts_mesh + 1)
    kz = np.linspace(0, 1, npts_mesh + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    for lat_id, lat_info in LATTICES.items():
        print(f"\nProcessing {lat_info['name']}...")
        E_mesh = lat_info['dispersion'](KX - 0.5, KY - 0.5, KZ - 0.5)
        
        for EF in lat_info['energies']:
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
                    lat_info['dispersion'], centroids[:, 0], centroids[:, 1], centroids[:, 2]
                )
                face_colors = compute_flat_colors(velocities, lat_info['cmap'])

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
                ax.view_init(elev=30, azim=45)
                
                fig.text(0.5, 0.92, f"Fermi Energy: {EF:+.2f} eV", 
                         color='#ffaa66', fontsize=16, ha='center', fontweight='bold')
                
                out_path = os.path.join(OUTPUT_DIR, f"static_{lat_id}_{EF:+.2f}eV.png")
                plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor=BG)
                plt.close(fig)
                
            except (ValueError, RuntimeError) as e:
                print(f"    Failed to generate surface at E = {EF} eV: {e}")

if __name__ == "__main__":
    main()
