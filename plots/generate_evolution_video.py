#!/usr/bin/env python3
"""
Fermi Surface Evolution Video Generator
=======================================
Generates a 3D video of the Fermi surface evolution from E = -4.0 eV to E = +4.0 eV.
- Video settings: 5 seconds, 30 FPS, 150 frames.

Output:
    Images/Rendered/fermi_surface_evolution.mp4
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')                          # Headless backend
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import skimage.measure

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
OUTPUT_DIR   = os.path.join(PROJECT_ROOT, "Images", "Rendered")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH  = os.path.join(OUTPUT_DIR, "fermi_surface_evolution.mp4")

# ── Grid & Energy Dispersion Setup ─────────────────────────────────────────
NPTS = 100 
kx = np.linspace(0, 1, NPTS + 1)
ky = np.linspace(0, 1, NPTS + 1)
kz = np.linspace(0, 1, NPTS + 1)
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

# Analytic tight-binding dispersion for simple cubic lattice:
# E(k) = -2*cos(kx) - 2*cos(ky) - 2*cos(kz)
E = -2.0 * (np.cos(2.0 * np.pi * (KX - 0.5)) + 
            np.cos(2.0 * np.pi * (KY - 0.5)) + 
            np.cos(2.0 * np.pi * (KZ - 0.5)))

# ── Shading Helper ─────────────────────────────────────────────────────────
def compute_shading_colors(verts_tri, cmap_name='inferno'):
    """
    Computes Fermi-velocity-colored faces with outward-aligned diffuse shading.
    """
    centroids = verts_tri.mean(axis=1)           # (N, 3)

    # Face normals for shading
    v1 = verts_tri[:, 1, :] - verts_tri[:, 0, :]
    v2 = verts_tri[:, 2, :] - verts_tri[:, 0, :]
    normals = np.cross(v1, v2)
    
    # Align normals to point outwards from BZ center (Gamma / origin)
    dots = np.sum(normals * centroids, axis=1, keepdims=True)
    normals = np.where(dots < 0, -normals, normals)
    
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms < 1e-14] = 1.0
    normals /= norms

    # Simple cubic analytical dispersion gradient (Fermi velocity speed)
    # centroids is centered at origin: [-0.5, 0.5]
    grad_x = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 0])
    grad_y = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 1])
    grad_z = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 2])
    grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    grad_mag[grad_mag < 1e-8] = 1e-8

    # Map color to Fermi velocity values
    color_vals = grad_mag
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=color_vals.min(), vmax=color_vals.max())
    face_colors = cmap(norm(color_vals)).copy()   # (N, 4) RGBA

    # Premium diffuse shading: light from top-front-left
    light = np.array([0.6, 0.5, 1.0])
    light /= np.linalg.norm(light)
    shading = np.clip(normals @ light, 0, 1)
    
    # Mix ambient (0.35) + diffuse (0.65)
    face_colors[:, :3] *= (0.35 + 0.65 * shading[:, np.newaxis])

    return face_colors

# ── Animation Setup ────────────────────────────────────────────────────────
FPS = 30
DURATION = 8  # seconds (increased from 5 to 8)
TOTAL_FRAMES = FPS * DURATION  # 240 frames (increased from 150)
ENERGIES = np.linspace(-6.0, 6.0, TOTAL_FRAMES)  # expanded from [-4.0, 4.0] to [-6.0, 6.0]

BG = '#07070f'
fig = plt.figure(figsize=(9, 9), facecolor=BG)
ax  = fig.add_subplot(111, projection='3d', facecolor=BG)
ax.set_proj_type('persp')

# Render static overlay text once
# The only text shown should be the "Fermi Energy: value eV" line.
subtitle_text = fig.text(0.5, 0.92, "Fermi Energy: -- eV", 
                         color='#ffaa66', fontsize=14, ha='center', fontweight='bold')

# Initialize empty list to track current polygon collection
poly_container = []

def update(frame_idx):
    EF = ENERGIES[frame_idx]
    
    # Clear previous surface polygons
    for p in poly_container:
        p.remove()
    poly_container.clear()

    # Isosurface values grid
    grid = E - EF

    # Extract isosurface using marching cubes (level=0.0 means E = EF)
    try:
        verts, faces, _, _ = skimage.measure.marching_cubes(
            grid, level=0.0, spacing=(1.0 / NPTS, 1.0 / NPTS, 1.0 / NPTS)
        )
        verts -= 0.5  # Center Gamma at origin
        verts_tri = verts[faces]
        
        # Compute face colors and shading
        face_colors = compute_shading_colors(verts_tri, 'inferno')

        # Add new polygons
        poly = Poly3DCollection(verts_tri.tolist(), zsort='average', 
                                facecolors=face_colors, edgecolors='none')
        ax.add_collection3d(poly)
        poly_container.append(poly)
    except (RuntimeError, ValueError):
        # Marching cubes can fail if the energy is completely outside the band (no surface)
        pass

    # Ensure equal aspect ratios and correct boundaries
    ax.set_box_aspect([1, 1, 1])
    span = 0.5
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_zlim(-span, span)

    # Completely hide grid, labels, axes lines, panes
    ax.set_axis_off()

    # Camera viewing angle
    ax.view_init(elev=22, azim=35)

    # Update text label
    subtitle_text.set_text(f"Fermi Energy: {EF:+.2f} eV")
    
    # Simple console progress log
    if (frame_idx + 1) % 10 == 0 or frame_idx == 0:
        print(f"  Frame {frame_idx + 1}/{TOTAL_FRAMES} rendered ({EF:+.2f} eV)")

# ── Write Video ────────────────────────────────────────────────────────────
def main():
    print(f"\nGenerating Fermi Surface Evolution Video ({TOTAL_FRAMES} frames)...")
    print(f"Target: {OUTPUT_PATH}")

    # Set up FFMpegWriter for high-quality H.264 MP4 export
    metadata = dict(title='Fermi Surface Evolution', artist='M Doyle')
    writer = animation.FFMpegWriter(fps=FPS, metadata=metadata, codec='libx264',
                                    extra_args=['-pix_fmt', 'yuv420p', '-crf', '18'])

    with writer.saving(fig, OUTPUT_PATH, dpi=160):
        for f in range(TOTAL_FRAMES):
            update(f)
            writer.grab_frame()
            
    plt.close(fig)
    print(f"\nSuccess! Video saved successfully to:\n{OUTPUT_PATH}\n")

if __name__ == '__main__':
    main()
