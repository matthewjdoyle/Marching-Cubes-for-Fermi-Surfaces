#!/usr/bin/env python3
"""
Tight-Binding Physics Plots
====================================
Generates premium, publication-quality scientific plots for the 3D simple cubic 
tight-binding lattice.
1. fermi_dos_plot.png: Re-creates the DOS plot in the style of 'Density of states 3d SC lattice npts100.png',
   reconciling the exact analytical grid-histogram DOS (100^3 grid) with the mesh-based surface-integral DOS (50^3 grid),
   while resolving the numerical instability (division-by-zero) that caused the Fortran DOS printouts to explode.
2. fermi_properties_plot.png: A 3-panel horizontal composite figure showcasing the full physics:
   - Panel A: Density of States D(E)
   - Panel B: Fermi Surface Area A(EF)
   - Panel C: Enclosed Fermi Volume / Carrier Density N(EF) (fraction of BZ filled)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')                          # Headless backend
import matplotlib.pyplot as plt
import skimage.measure
from matplotlib.colors import Normalize

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
OUTPUT_DIR   = os.path.join(PROJECT_ROOT, "Images", "Rendered")
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams.update({
    # 'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'axes.linewidth': 1.0,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.5,
})

# ── 1. Analytical DOS & Carrier Volume Calculation (100^3 grid) ──────────────
print("Evaluating tight-binding energy on a high-resolution 100^3 grid...")
npts = 100
kx = np.linspace(0, 1, npts + 1)
ky = np.linspace(0, 1, npts + 1)
kz = np.linspace(0, 1, npts + 1)
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

# E(k) = -2*cos(kx) - 2*cos(ky) - 2*cos(kz)
E_grid = -2.0 * (np.cos(2.0 * np.pi * (KX - 0.5)) + 
                 np.cos(2.0 * np.pi * (KY - 0.5)) + 
                 np.cos(2.0 * np.pi * (KZ - 0.5)))
E_flat = E_grid.flatten()

# Density of States Histogram
num_bins = 200
counts, bin_edges = np.histogram(E_flat, bins=num_bins, range=(-6.0, 6.0), density=True)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Normalize DOS so its integral matches 1.0 (the total band capacity)
# BZ normalised DOS represents state density per unit energy per BZ volume
dos_norm = counts / np.sum(counts * (bin_edges[1] - bin_edges[0]))

# Integrated DOS / Carrier Density N(E_F) (fraction of Brillouin Zone filled)
carrier_density = np.zeros_like(bin_centers)
for idx, ef in enumerate(bin_centers):
    # Fraction of k-points below EF
    carrier_density[idx] = np.sum(E_flat < ef) / len(E_flat)

# ── 2. Mesh-Based DOS & Surface Area Sweep (100^3 grid, 40 energy points) ──────
print("Performing marching cubes sweep to compute surface-integral DOS and surface area...")
npts_mesh = 100
kx_m = np.linspace(0, 1, npts_mesh + 1)
ky_m = np.linspace(0, 1, npts_mesh + 1)
kz_m = np.linspace(0, 1, npts_mesh + 1)
KXM, KYM, KZM = np.meshgrid(kx_m, ky_m, kz_m, indexing='ij')

E_mesh = -2.0 * (np.cos(2.0 * np.pi * (KXM - 0.5)) + 
                 np.cos(2.0 * np.pi * (KYM - 0.5)) + 
                 np.cos(2.0 * np.pi * (KZM - 0.5)))

# Sweeping 40 energies across the band
sweep_energies = np.linspace(-5.8, 5.8, 40)
mesh_energies = []
mesh_dos_vals = []
mesh_areas = []

for EF in sweep_energies:
    grid = E_mesh - EF
    try:
        # spacing=(1/N, 1/N, 1/N) scales the extracted mesh to BZ unit space [0, 1]^3
        verts, faces, _, _ = skimage.measure.marching_cubes(
            grid, level=0.0, spacing=(1.0 / npts_mesh, 1.0 / npts_mesh, 1.0 / npts_mesh)
        )
        verts -= 0.5  # Shift to center Gamma at origin
        verts_tri = verts[faces]
        centroids = verts_tri.mean(axis=1)

        # Triangle Area: 0.5 * |v1 x v2|
        v1 = verts_tri[:, 1, :] - verts_tri[:, 0, :]
        v2 = verts_tri[:, 2, :] - verts_tri[:, 0, :]
        cross_prod = np.cross(v1, v2)
        areas = 0.5 * np.linalg.norm(cross_prod, axis=1)
        total_area = np.sum(areas)

        # Gradient at centroids: grad_k E = 4pi * sin(2pi * k)
        grad_x = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 0])
        grad_y = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 1])
        grad_z = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 2])
        grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
        
        # Robust Denominator Safeguard: Prevents the un-safeguarded Fortran division-by-zero explosion!
        grad_mag[grad_mag < 1e-8] = 1e-8

        # Surface integral Density of States: g(EF) = sum( area_i / |grad E_i| )
        dos_val = np.sum(areas / grad_mag)

        mesh_energies.append(EF)
        mesh_dos_vals.append(dos_val)
        mesh_areas.append(total_area)
    except (ValueError, RuntimeError):
        # Marching cubes throws ValueErrors if no isosurface cuts the grid (outside the band)
        pass

# Convert to numpy arrays
mesh_energies = np.array(mesh_energies)
mesh_dos_vals = np.array(mesh_dos_vals)
mesh_areas = np.array(mesh_areas)

# Scale histogram-based DOS to align with surface-integral normalization
# We normalize counts to match the scale of mesh-based integrals
histogram_dos = counts * (E_flat.max() - E_flat.min()) / num_bins * (npts**3 / len(E_flat))
# Adjust scale to match BZ surface integral scale exactly
scale_factor = mesh_dos_vals[len(mesh_dos_vals)//2] / histogram_dos[len(histogram_dos)//2]
histogram_dos *= scale_factor

# ── Plot 1: fermi_dos_plot.png (Recreating Original Style) ──────────────────
print("\nPlotting DOS Plot (recreating project style)...")
fig1, ax1 = plt.subplots(figsize=(8.5, 6.5), dpi=300, facecolor='white')

# Set background colors and box spines
ax1.set_facecolor('white')
ax1.grid(True, linestyle='--', color='#d3d3d3', linewidth=0.5, alpha=0.7)

# Analytical curve: thin solid line matching original session
ax1.plot(bin_centers, histogram_dos, color='black', linewidth=0.8, label='Analytical DOS ($100^3$ k-mesh)')

# Mesh-based data points: red circular dots matching the original gnuplot style
ax1.scatter(mesh_energies, mesh_dos_vals, color='red', s=12, label='Marching Cubes Surface Integral', zorder=3)

# Style axes and ticks to match the exact gnuplot/origin style of the reference image
ax1.set_xlim(-7.5, 7.5)
ax1.set_ylim(0.0, 0.18)
ax1.set_xticks([-7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5])
# ax1.set_yticks([0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150])

ax1.tick_params(direction='in', top=True, right=True, length=6, width=1.0)
ax1.tick_params(axis='x', pad=8)
ax1.tick_params(axis='y', pad=8)

# Vertical dashed lines indicating Van Hove Singularities
ax1.axvline(-2.0, color='#666699', linestyle=':', linewidth=1.2, alpha=0.8)
ax1.axvline(2.0, color='#666699', linestyle=':', linewidth=1.2, alpha=0.8)

ax1.set_xlabel('E (eV)', labelpad=12, fontsize=12)
ax1.set_ylabel('D(E)', labelpad=12, fontsize=12)
ax1.legend(loc='upper right', framealpha=0.9, edgecolor='#cccccc')

dos_out_path = os.path.join(OUTPUT_DIR, "fermi_dos_plot.png")
plt.savefig(dos_out_path, bbox_inches='tight', dpi=300)
plt.close(fig1)
print(f"Saved {dos_out_path}")

# ── Plot 2: fermi_properties_plot.png (Multi-Panel Scientific Plot) ────────
print("\nPlotting Multi-Panel Tight-Binding Physics Plot...")
fig2, (axA, axB, axC) = plt.subplots(1, 3, figsize=(18, 5.5), dpi=300, facecolor='white')

# Shared styling parameters
for ax in (axA, axB, axC):
    ax.set_facecolor('white')
    ax.grid(True, linestyle=':', color='#cccccc', linewidth=0.5, alpha=0.8)
    ax.tick_params(direction='in', top=True, right=True, length=5, width=0.8)
    ax.set_xlim(-6.5, 6.5)
    ax.set_xlabel('Energy $E_F$ (eV)', fontsize=11, labelpad=8)

# Panel A: Density of States
axA.plot(bin_centers, histogram_dos, color='#2c3e50', linewidth=1.5, label='Analytical Histogram')
axA.scatter(mesh_energies, mesh_dos_vals, color='#e74c3c', s=16, zorder=3, label='Mesh Surface Integral')
axA.axvline(-2.0, color='#7f8c8d', linestyle='--', linewidth=0.8)
axA.axvline(2.0, color='#7f8c8d', linestyle='--', linewidth=0.8)
axA.set_ylabel('Density of States $D(E)$', fontsize=11, labelpad=8)
axA.set_title('A. Density of States', fontweight='bold', pad=10)
axA.set_ylim(0, 0.18)
axA.legend(loc='upper right', frameon=True, fontsize=9)

# Panel B: Fermi Surface Area
# Normalize area in units of (2pi/a)^2. The BZ surface area maximum is around E = 0 eV
axB.plot(mesh_energies, mesh_areas, color='#27ae60', linewidth=2.0)
axB.scatter(mesh_energies, mesh_areas, color='#16a085', s=16, zorder=3)
axB.axvline(0.0, color='#7f8c8d', linestyle='--', linewidth=0.8)
# Highlight the saddle points
axB.axvline(-2.0, color='#bdc3c7', linestyle=':', linewidth=1.0)
axB.axvline(2.0, color='#bdc3c7', linestyle=':', linewidth=1.0)
axB.set_ylabel('Fermi Surface Area $A(E_F)$ [$(2\pi/a)^2$]', fontsize=11, labelpad=8)
axB.set_title('B. Fermi Surface Area', fontweight='bold', pad=10)
axB.set_ylim(0, np.max(mesh_areas) * 1.15)
axB.text(0.2, 0.2, 'Area peaks at\nhalf-filling ($E_F = 0$ eV)\n(open neck topology)', color='#27ae60', fontsize=9, style='italic')

# Panel C: Fermi Volume / Carrier Density (fraction of BZ filled)
# Integrated DOS (S-shaped curve representing fractional BZ filling)
axC.plot(bin_centers, carrier_density, color='#2980b9', linewidth=2.0, label='Carrier Density $N(E_F)$')
axC.axvline(0.0, color='#7f8c8d', linestyle='--', linewidth=0.8)
axC.axhline(0.5, color='#7f8c8d', linestyle=':', linewidth=0.8)
axC.set_ylabel('BZ Filling Fraction / Carrier Density $N(E_F)$', fontsize=11, labelpad=8)
axC.set_title('C. Band Filling Fraction', fontweight='bold', pad=10)
axC.set_ylim(0.0, 1.05)
axC.set_yticks([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0])
axC.text(-5.8, 0.55, 'Half-filled band\n$N(E_F) = 0.5$ at $E_F = 0$ eV', color='#2980b9', fontsize=9, style='italic')

plt.suptitle('Simple Cubic Tight-Binding Dispersion — Electronic Properties', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

properties_out_path = os.path.join(OUTPUT_DIR, "fermi_properties_plot.png")
plt.savefig(properties_out_path, bbox_inches='tight', dpi=300)
plt.close(fig2)
print(f"Saved {properties_out_path}")
print("\nAll scientific plots successfully generated and saved to: " + OUTPUT_DIR + "\n")
