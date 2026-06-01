#!/usr/bin/env python3
"""

Generates high-resolution Density of States (DOS) plots and premium 3D 
marching-cubes evolution videos for multiple crystallographic lattices:
1. Simple Cubic (SC) — Nearest-Neighbor only
2. Simple Cubic with Next-Nearest-Neighbor (SC-NNN) hopping (t' = -0.2 eV)
3. Body-Centered Cubic (BCC) tight-binding lattice
4. Face-Centered Cubic (FCC) tight-binding lattice


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

# Plotting parameters
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'axes.linewidth': 1.0,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.5,
})

# ── Lattice Specifications ─────────────────────────────────────────────────
LATTICES = {
    'sc': {
        'name': 'Simple Cubic (SC)',
        'dispersion': lambda kx, ky, kz: -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz)),
        'grad': lambda kx, ky, kz: (
            4.0 * np.pi * np.sin(2*np.pi*kx),
            4.0 * np.pi * np.sin(2*np.pi*ky),
            4.0 * np.pi * np.sin(2*np.pi*kz)
        ),
        'cmap': 'inferno',
        'xlim': (-7.5, 7.5),
        'ylim': (0.0, 0.16),
        'xticks': [-7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5],
        'yticks': [0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150]
    },
    'sc_nnn': {
        'name': 'Simple Cubic with NNN ($t\' = -0.2$ eV)',
        # E(k) = -2*t*(cos_x + cos_y + cos_z) - 4*t'*(cos_x*cos_y + cos_y*cos_z + cos_z*cos_x)
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz))
            - 4.0 * (-0.2) * (np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky) + 
                              np.cos(2*np.pi*ky)*np.cos(2*np.pi*kz) + 
                              np.cos(2*np.pi*kz)*np.cos(2*np.pi*kx))
        ),
        'grad': lambda kx, ky, kz: (
            4.0*np.pi*np.sin(2*np.pi*kx) + 8.0*np.pi*(-0.2)*np.sin(2*np.pi*kx)*(np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz)),
            4.0*np.pi*np.sin(2*np.pi*ky) + 8.0*np.pi*(-0.2)*np.sin(2*np.pi*ky)*(np.cos(2*np.pi*kx) + np.cos(2*np.pi*kz)),
            4.0*np.pi*np.sin(2*np.pi*kz) + 8.0*np.pi*(-0.2)*np.sin(2*np.pi*kz)*(np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky))
        ),
        'cmap': 'plasma',
        'xlim': (-7.5, 7.5),
        'ylim': (0.0, 0.18),
        'xticks': [-7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5],
        'yticks': [0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18]
    },
    'sc_anisotropic': {
        'name': 'Anisotropic Simple Cubic (Layered)',
        # E(k) = -2*t_xy*(cos_x + cos_y) - 2*t_z*cos_z (t_xy = 1.0, t_z = 0.15)
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky)) - 2.0 * 0.15 * np.cos(2*np.pi*kz)
        ),
        'grad': lambda kx, ky, kz: (
            4.0*np.pi*np.sin(2*np.pi*kx),
            4.0*np.pi*np.sin(2*np.pi*ky),
            4.0*np.pi*0.15*np.sin(2*np.pi*kz)
        ),
        'cmap': 'cividis',
        'xlim': (-4.5, 4.5),
        'ylim': (0.0, 0.35),
        'xticks': [-4.5, -3.0, -1.5, 0.0, 1.5, 3.0, 4.5],
        'yticks': [0.00, 0.07, 0.14, 0.21, 0.28, 0.35]
    },
    'bcc': {
        'name': 'Body-Centered Cubic (BCC)',
        # E(k) = -8*t*cos(2pi*kx)*cos(2pi*ky)*cos(2pi*kz)
        'dispersion': lambda kx, ky, kz: -8.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi*ky) * np.cos(2*np.pi*kz),
        'grad': lambda kx, ky, kz: (
            16.0 * np.pi * np.sin(2*np.pi*kx) * np.cos(2*np.pi*ky) * np.cos(2*np.pi*kz),
            16.0 * np.pi * np.cos(2*np.pi*kx) * np.sin(2*np.pi*ky) * np.cos(2*np.pi*kz),
            16.0 * np.pi * np.cos(2*np.pi*kx) * np.cos(2*np.pi*ky) * np.sin(2*np.pi*kz)
        ),
        'cmap': 'viridis',
        'xlim': (-9.5, 9.5),
        'ylim': (0.0, 0.20),
        'xticks': [-8.0, -4.0, 0.0, 4.0, 8.0],
        'yticks': [0.00, 0.04, 0.08, 0.12, 0.16, 0.20]
    },
    'fcc': {
        'name': 'Face-Centered Cubic (FCC)',
        # E(k) = -4*t*(cos(2pi*kx)*cos(2pi*ky) + cos(2pi*ky)*cos(2pi*kz) + cos(2pi*kz)*cos(2pi*kx))
        'dispersion': lambda kx, ky, kz: -4.0 * (
            np.cos(2*np.pi*kx)*np.cos(2*np.pi*ky) + 
            np.cos(2*np.pi*ky)*np.cos(2*np.pi*kz) + 
            np.cos(2*np.pi*kz)*np.cos(2*np.pi*kx)
        ),
        'grad': lambda kx, ky, kz: (
            8.0 * np.pi * np.sin(2*np.pi*kx) * (np.cos(2*np.pi*ky) + np.cos(2*np.pi*kz)),
            8.0 * np.pi * np.sin(2*np.pi*ky) * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*kz)),
            8.0 * np.pi * np.sin(2*np.pi*kz) * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky))
        ),
        'cmap': 'magma',
        'xlim': (-13.5, 5.5),
        'ylim': (0.0, 0.24),
        'xticks': [-12.0, -8.0, -4.0, 0.0, 4.0],
        'yticks': [0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24]
    },
    'hexagonal': {
        'name': 'Simple Hexagonal',
        # E(k) = -2*t_xy*(cos(4pi*kx) + 2*cos(2pi*kx)*cos(2pi*sqrt(3)*ky)) - 2*t_z*cos(2pi*kz)
        # With t_xy = 1.0, t_z = 0.5
        'dispersion': lambda kx, ky, kz: (
            -2.0 * (np.cos(4*np.pi*kx) + 2.0 * np.cos(2*np.pi*kx) * np.cos(2*np.pi * np.sqrt(3) * ky))
            - 2.0 * 0.5 * np.cos(2*np.pi*kz)
        ),
        'grad': lambda kx, ky, kz: (
            8.0*np.pi*np.sin(4*np.pi*kx) + 8.0*np.pi*np.sin(2*np.pi*kx)*np.cos(2*np.pi * np.sqrt(3) * ky),
            8.0*np.sqrt(3)*np.pi*np.cos(2*np.pi*kx)*np.sin(2*np.pi * np.sqrt(3) * ky),
            2.0*np.pi*np.sin(2*np.pi*kz)
        ),
        'cmap': 'coolwarm',
        'xlim': (-8.5, 5.5),
        'ylim': (0.0, 0.18),
        'xticks': [-8.0, -4.0, 0.0, 4.0],
        'yticks': [0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18]
    }
}

# ── Outward Shading Helper ─────────────────────────────────────────────────
def compute_shading_colors(verts_tri, velocities, cmap_name='inferno'):
    centroids = verts_tri.mean(axis=1)           # (N, 3)

    v1 = verts_tri[:, 1, :] - verts_tri[:, 0, :]
    v2 = verts_tri[:, 2, :] - verts_tri[:, 0, :]
    normals = np.cross(v1, v2)
    
    # Align normals to point outwards from center
    dots = np.sum(normals * centroids, axis=1, keepdims=True)
    normals = np.where(dots < 0, -normals, normals)
    
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms < 1e-14] = 1.0
    normals /= norms

    # Map color to velocity values
    color_vals = velocities
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=color_vals.min(), vmax=color_vals.max())
    face_colors = cmap(norm(color_vals)).copy()   # (N, 4) RGBA

    # Shadows removed as requested
    # light = np.array([0.6, 0.5, 1.0])
    # light /= np.linalg.norm(light)
    # shading = np.clip(normals @ light, 0, 1)
    # face_colors[:, :3] *= (0.35 + 0.65 * shading[:, np.newaxis])

    return face_colors

# ── Primary Engine Process ─────────────────────────────────────────────────
def process_lattice(lat_id, lat_info):
    print(f"\n==================================================================")
    print(f"PROCESSING LATTICE: {lat_info['name']}")
    print(f"==================================================================")

    # ── 1. high-resolution DOS Grid Evaluation (100^3 grid) ──
    print("Evaluating dispersion on a high-resolution 100^3 k-mesh...")
    npts = 100
    kx = np.linspace(0, 1, npts + 1)
    ky = np.linspace(0, 1, npts + 1)
    kz = np.linspace(0, 1, npts + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    # Centered BZ fractional coordinates [-0.5, 0.5]
    E_grid = lat_info['dispersion'](KX - 0.5, KY - 0.5, KZ - 0.5)
    E_flat = E_grid.flatten()

    E_min, E_max = E_flat.min(), E_flat.max()
    print(f"  Energy range: {E_min:+.4f} eV to {E_max:+.4f} eV")

    # Histogram-based DOS
    num_bins = 200
    counts, bin_edges = np.histogram(E_flat, bins=num_bins, range=(E_min, E_max), density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Normalize DOS
    histogram_dos = counts * (E_max - E_min) / num_bins * (npts**3 / len(E_flat))

    # ── 2. Mesh-Based DOS Surface Integrals (100^3 grid, 40 energy points) ──
    print("Performing marching cubes sweep to compute surface-integral DOS...")
    npts_mesh = 100
    kx_m = np.linspace(0, 1, npts_mesh + 1)
    ky_m = np.linspace(0, 1, npts_mesh + 1)
    kz_m = np.linspace(0, 1, npts_mesh + 1)
    KXM, KYM, KZM = np.meshgrid(kx_m, ky_m, kz_m, indexing='ij')

    E_mesh = lat_info['dispersion'](KXM - 0.5, KYM - 0.5, KZM - 0.5)

    # Sweep energies slightly inside band edges to prevent marching cubes from rendering empty spaces
    sweep_energies = np.linspace(E_min * 0.98, E_max * 0.98, 40)
    mesh_energies = []
    mesh_dos_vals = []

    for EF in sweep_energies:
        grid = E_mesh - EF
        try:
            verts, faces, _, _ = skimage.measure.marching_cubes(
                grid, level=0.0, spacing=(1.0 / npts_mesh, 1.0 / npts_mesh, 1.0 / npts_mesh)
            )
            verts -= 0.5  # Shift to center Gamma at origin
            verts_tri = verts[faces]
            centroids = verts_tri.mean(axis=1)

            # Triangle Area
            v1 = verts_tri[:, 1, :] - verts_tri[:, 0, :]
            v2 = verts_tri[:, 2, :] - verts_tri[:, 0, :]
            cross_prod = np.cross(v1, v2)
            areas = 0.5 * np.linalg.norm(cross_prod, axis=1)

            # Gradient from analytic definition
            grad_x, grad_y, grad_z = lat_info['grad'](centroids[:, 0], centroids[:, 1], centroids[:, 2])
            grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
            
            # Robust Denominator Safeguard
            grad_mag[grad_mag < 1e-8] = 1e-8

            # Surface integral DOS: sum( area / |grad E| )
            dos_val = np.sum(areas / grad_mag)

            mesh_energies.append(EF)
            mesh_dos_vals.append(dos_val)
        except (ValueError, RuntimeError):
            pass

    mesh_energies = np.array(mesh_energies)
    mesh_dos_vals = np.array(mesh_dos_vals)

    # Reconcile scaling of histogram to match surface integrals
    if len(mesh_dos_vals) > 0:
        scale_factor = mesh_dos_vals[len(mesh_dos_vals)//2] / histogram_dos[len(histogram_dos)//2]
        histogram_dos *= scale_factor

    # ── 3. Plot Minimalist DOS Figure (no title, no annotations) ──
    print("Plotting minimalist DOS figure...")
    fig1, ax1 = plt.subplots(figsize=(8.5, 6.5), dpi=300, facecolor='white')
    ax1.set_facecolor('white')
    ax1.grid(True, linestyle='--', color='#d3d3d3', linewidth=0.5, alpha=0.7)

    # Analytical curve
    ax1.plot(bin_centers, histogram_dos, color='black', linewidth=0.8, 
             label=f'Analytical DOS ({lat_info["name"]})')

    # Mesh-based surface integral points
    if len(mesh_dos_vals) > 0:
        ax1.scatter(mesh_energies, mesh_dos_vals, color='red', s=12, 
                    label='Marching Cubes Surface Integral', zorder=3)

    # Styling limits, tick coordinates, and inward ticks
    ax1.set_xlim(lat_info['xlim'])
    ax1.set_ylim(lat_info['ylim'])
    ax1.set_xticks(lat_info['xticks'])
    ax1.set_yticks(lat_info['yticks'])

    ax1.tick_params(direction='in', top=True, right=True, length=6, width=1.0)
    ax1.tick_params(axis='x', pad=8)
    ax1.tick_params(axis='y', pad=8)

    ax1.set_xlabel('E (eV)', labelpad=12, fontsize=12)
    ax1.set_ylabel('D(E)', labelpad=12, fontsize=12)
    ax1.legend(loc='upper right', framealpha=0.9, edgecolor='#cccccc')

    dos_out_path = os.path.join(OUTPUT_DIR, f"fermi_dos_plot_{lat_id}.png")
    plt.savefig(dos_out_path, bbox_inches='tight', dpi=300)
    plt.close(fig1)
    print(f"  Saved DOS plot to: {dos_out_path}")

    # ── 4. Generate Premium 3D Fermi Surface Sweep Video ──
    FPS = 30
    DURATION = 8  # seconds
    TOTAL_FRAMES = FPS * DURATION  # 240 frames
    
    # Use optimized 64^3 mesh for fluid real-time 3D marching cubes video generation
    npts_mesh_vid = 64
    kx_v = np.linspace(0, 1, npts_mesh_vid + 1)
    ky_v = np.linspace(0, 1, npts_mesh_vid + 1)
    kz_v = np.linspace(0, 1, npts_mesh_vid + 1)
    KXV, KYV, KZV = np.meshgrid(kx_v, ky_v, kz_v, indexing='ij')
    E_mesh_vid = lat_info['dispersion'](KXV - 0.5, KYV - 0.5, KZV - 0.5)

    # Expand energy sweep range slightly outside band edges for full birth-death sweep
    sweep_energies_vid = np.linspace(E_min * 1.02, E_max * 1.02, TOTAL_FRAMES)

    viewing_angles = [(22, 35), (45, 125), (0, 0)]
    for angle_idx, (cam_elev, cam_azim) in enumerate(viewing_angles):
        BG = '#07070f'
        fig = plt.figure(figsize=(9, 9), facecolor=BG)
        ax  = fig.add_subplot(111, projection='3d', facecolor=BG)
        ax.set_proj_type('persp')

        # Render static overlay text once (minimalist style: ONLY the Fermi Energy line is shown)
        subtitle_text = fig.text(0.5, 0.92, "Fermi Energy: -- eV", 
                                 color='#ffaa66', fontsize=14, ha='center', fontweight='bold')

        poly_container = []

        def update_frame(frame_idx):
            EF = sweep_energies_vid[frame_idx]
            
            # Clear previous polygons
            for p in poly_container:
                p.remove()
            poly_container.clear()

            # Isosurface grid values
            grid = E_mesh_vid - EF

            try:
                verts, faces, _, _ = skimage.measure.marching_cubes(
                    grid, level=0.0, spacing=(1.0 / npts_mesh_vid, 1.0 / npts_mesh_vid, 1.0 / npts_mesh_vid)
                )
                verts -= 0.5  # Center Gamma at origin
                verts_tri = verts[faces]
                
                # Compute shading colors using aligned normals
                # Gradient from analytic definition
                centroids = verts_tri.mean(axis=1)
                grad_x, grad_y, grad_z = lat_info['grad'](centroids[:, 0], centroids[:, 1], centroids[:, 2])
                grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
                grad_mag[grad_mag < 1e-8] = 1e-8
                
                face_colors = compute_shading_colors(verts_tri, grad_mag, lat_info['cmap'])

                poly = Poly3DCollection(verts_tri.tolist(), zsort='average', 
                                        facecolors=face_colors, edgecolors='none')
                ax.add_collection3d(poly)
                poly_container.append(poly)
            except (ValueError, RuntimeError):
                pass

            # Equal aspect ratio and bounds
            ax.set_box_aspect([1, 1, 1])
            span = 0.5
            ax.set_xlim(-span, span)
            ax.set_ylim(-span, span)
            ax.set_zlim(-span, span)

            # Hide axes, grid, labels, panes
            ax.set_axis_off()

            # Camera viewpoint
            ax.view_init(elev=cam_elev, azim=cam_azim)

            # Update energy text dynamically
            subtitle_text.set_text(f"Fermi Energy: {EF:+.2f} eV")
            
            if (frame_idx + 1) % 20 == 0 or frame_idx == 0:
                print(f"    Frame {frame_idx + 1}/{TOTAL_FRAMES} rendered ({EF:+.2f} eV)")

    # Set up FFMpegWriter for high-quality MP4 H.264
        video_out_path = os.path.join(OUTPUT_DIR, f"fermi_surface_evolution_{lat_id}_angle{angle_idx+1}.mp4")
        print(f"Generating Fermi Surface Evolution Video ({TOTAL_FRAMES} frames) for angle {angle_idx+1}...")
        
        metadata = dict(title=f'Fermi Surface Evolution — {lat_info["name"]}', artist='M Doyle')
        writer = animation.FFMpegWriter(fps=FPS, metadata=metadata, codec='libx264',
                                        extra_args=['-pix_fmt', 'yuv420p', '-crf', '18'])

        with writer.saving(fig, video_out_path, dpi=160):
            for f in range(TOTAL_FRAMES):
                update_frame(f)
                writer.grab_frame()
                
        plt.close(fig)
        print(f"  Saved sweep video to: {video_out_path}")

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("\n==================================================================")
    print("STARTING SOLID STATE PHYSICS TB-ENGINE SWEEP PROCESS")
    print("==================================================================")

    for lat_id, lat_info in LATTICES.items():
        process_lattice(lat_id, lat_info)
        
    print("\n==================================================================")
    print("ALL RUNS COMPLETED SUCCESSFULLY!")
    print("All final plots and evolution videos are stored in: " + OUTPUT_DIR)
    print("==================================================================\n")

if __name__ == '__main__':
    main()
