import numpy as np
import matplotlib.pyplot as plt
import skimage.measure
from matplotlib.colors import Normalize

# 1. Histogram-based DOS on 100^3 grid (extremely fast and smooth!)
print("Calculating histogram-based DOS...")
npts_hist = 100
kx = np.linspace(0, 1, npts_hist + 1)
ky = np.linspace(0, 1, npts_hist + 1)
kz = np.linspace(0, 1, npts_hist + 1)
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

E_hist = -2.0 * (np.cos(2.0 * np.pi * (KX - 0.5)) + 
                 np.cos(2.0 * np.pi * (KY - 0.5)) + 
                 np.cos(2.0 * np.pi * (KZ - 0.5)))

# Flat array of energies
E_flat = E_hist.flatten()

# Histogram
num_bins = 200
counts, bin_edges = np.histogram(E_flat, bins=num_bins, range=(-6.0, 6.0), density=True)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

print("Calculating mesh-based DOS at discrete energies...")
# Let's test mesh-based DOS at a few energies on a 50^3 grid
npts_mesh = 50
kx_m = np.linspace(0, 1, npts_mesh + 1)
ky_m = np.linspace(0, 1, npts_mesh + 1)
kz_m = np.linspace(0, 1, npts_mesh + 1)
KXM, KYM, KZM = np.meshgrid(kx_m, ky_m, kz_m, indexing='ij')

E_mesh = -2.0 * (np.cos(2.0 * np.pi * (KXM - 0.5)) + 
                 np.cos(2.0 * np.pi * (KYM - 0.5)) + 
                 np.cos(2.0 * np.pi * (KZM - 0.5)))

test_energies = [-4.0, -2.0, 0.0, 2.0, 4.0]
for EF in test_energies:
    grid = E_mesh - EF
    try:
        verts, faces, _, _ = skimage.measure.marching_cubes(
            grid, level=0.0, spacing=(1.0 / npts_mesh, 1.0 / npts_mesh, 1.0 / npts_mesh)
        )
        verts -= 0.5
        verts_tri = verts[faces]
        centroids = verts_tri.mean(axis=1)

        # Areas of triangles
        v1 = verts_tri[:, 1, :] - verts_tri[:, 0, :]
        v2 = verts_tri[:, 2, :] - verts_tri[:, 0, :]
        cross_prod = np.cross(v1, v2)
        areas = 0.5 * np.linalg.norm(cross_prod, axis=1)

        # Gradient at centroids
        # grad E = 4pi * sin(2pi * k)
        grad_x = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 0])
        grad_y = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 1])
        grad_z = 4.0 * np.pi * np.sin(2.0 * np.pi * centroids[:, 2])
        grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
        grad_mag[grad_mag < 1e-10] = 1e-10

        # DOS sum
        # In BZ units, the surface integral matches the BZ-normalized DOS
        # Volume of BZ in kx,ky,kz is 1.0, so the BZ-normalized DOS is:
        # g(E) = (1 / V_BZ) * sum(area / |grad E|)
        dos_val = np.sum(areas / grad_mag)
        
        # In BZ units, V_BZ = 1.0
        # Let's see: the DOS value from the histogram at this energy is:
        hist_idx = np.abs(bin_centers - EF).argmin()
        hist_dos = counts[hist_idx]
        
        # Print comparison
        print(f"E = {EF:+.2f} eV: Mesh DOS = {dos_val:.4f}, Histogram DOS = {hist_dos:.4f}")
    except Exception as e:
        print(f"E = {EF:+.2f} eV: Failed ({e})")
