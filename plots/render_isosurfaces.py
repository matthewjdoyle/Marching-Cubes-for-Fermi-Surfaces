#!/usr/bin/env python3
"""
Fermi Surface Renderer
======================
Reads the *isosurf.dn gnuplot pm3d triangle files produced by the MC33
Fortran code and renders them as beautiful 3D isosurface images.

Output format of isosurf.dn (per triangle):
    p1.x  p1.y  p1.z
    p1.x  p1.y  p1.z    <- duplicate (gnuplot pm3d scan-line format)
                         <- blank line
    p2.x  p2.y  p2.z
    p3.x  p3.y  p3.z
                         <- blank line
                         <- blank line (end of triangle block)

Coordinates are in fractional BZ units [0, 1], BZ centre at (0.5, 0.5, 0.5).
This script shifts them to [-0.5, 0.5] to centre Gamma at the origin.

Usage:
    python render_isosurfaces.py
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')                          # headless backend
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR   = os.path.join(PROJECT_ROOT, "Images", "Rendered")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Parser ─────────────────────────────────────────────────────────────────
def parse_isosurf_dn(filepath):
    """
    Parse a gnuplot pm3d isosurf.dn file into a list of triangles and velocities.
    Each triangle is a (3, 3) numpy array: rows = vertices, cols = x,y,z.
    Each velocity is the average velocity of the three vertices of that triangle.

    Handles:
      - Windows CRLF line endings (Python text mode strips CR automatically)
      - Fortran D-exponent notation  (e.g. 0.260000D+00 → 0.260000E+00)
    """
    with open(filepath, 'r') as f:
        raw = f.read()

    # Fortran writes D exponents; replace with E so Python can parse them
    raw = raw.replace('D+', 'E+').replace('D-', 'E-')
    # Normalise any remaining stray \r
    raw = raw.replace('\r', '')

    triangles = []
    velocities = []
    # Each triangle block ends with a double blank line → split on \n\n\n
    blocks = raw.split('\n\n\n')

    for block in blocks:
        block = block.strip()
        if not block:
            continue
        # Within a block:  "p1\np1\n\np2\np3"
        halves = [h.strip() for h in block.split('\n\n') if h.strip()]
        if len(halves) < 2:
            continue
        try:
            p1_lines  = [l for l in halves[0].split('\n') if l.strip()]
            p23_lines = [l for l in halves[1].split('\n') if l.strip()]
            if len(p1_lines) < 1 or len(p23_lines) < 2:
                continue
            
            p1_raw = [float(x) for x in p1_lines[0].split()]
            p2_raw = [float(x) for x in p23_lines[0].split()]
            p3_raw = [float(x) for x in p23_lines[1].split()]
            
            if len(p1_raw) >= 3 and len(p2_raw) >= 3 and len(p3_raw) >= 3:
                p1 = np.array(p1_raw[:3])
                p2 = np.array(p2_raw[:3])
                p3 = np.array(p3_raw[:3])
                
                # Filter out degenerate triangles (which have zero area and skew color ranges)
                d12 = np.linalg.norm(p1 - p2)
                d23 = np.linalg.norm(p2 - p3)
                d31 = np.linalg.norm(p3 - p1)
                if d12 < 1e-5 or d23 < 1e-5 or d31 < 1e-5:
                    continue
                
                triangles.append(np.array([p1, p2, p3]))
                
                # Check for 4th column (velocity); default to average z-coordinate (p[:, 2]) if not present
                v1 = p1_raw[3] if len(p1_raw) > 3 else p1[2]
                v2 = p2_raw[3] if len(p2_raw) > 3 else p2[2]
                v3 = p3_raw[3] if len(p3_raw) > 3 else p3[2]
                velocities.append((v1 + v2 + v3) / 3.0)
                
        except (ValueError, IndexError):
            continue

    return triangles, velocities



def energy_from_filename(filepath):
    """Extract Fermi energy from a filename like '-2.000isosurf.dn' or 'fermi_surface_neg2.000_isosurf.dn'."""
    base = os.path.basename(filepath)
    prefix = base.replace('_isosurf.dn', '').replace('isosurf.dn', '').strip()
    if 'fermi_surface_' in prefix:
        prefix = prefix.replace('fermi_surface_', '')
    if 'neg' in prefix:
        prefix = prefix.replace('neg', '-')
    try:
        return float(prefix)
    except ValueError:
        return None


# ── Renderer ───────────────────────────────────────────────────────────────
CMAPS = {
    -4.0: 'plasma',
    -2.0: 'inferno',
     0.0: 'magma',
     2.0: 'viridis',
     4.0: 'cividis',
}
DEFAULT_CMAP = 'plasma'

# Light source for shading
_LS = LightSource(azdeg=225, altdeg=45)


def compute_face_colors(verts, velocities, cmap_name):
    """
    Colour faces by their Fermi velocity (or parsed scalar value), then apply diffuse shading via face normals.
    verts : (N, 3, 3) array
    velocities : (N,) array
    Returns RGBA array (N, 4).
    """
    centroids = verts.mean(axis=1)           # (N, 3)

    # Face normals for shading
    v1 = verts[:, 1, :] - verts[:, 0, :]
    v2 = verts[:, 2, :] - verts[:, 0, :]
    normals = np.cross(v1, v2)
    
    # Align normals to point outwards from BZ center (Gamma / origin)
    dots = np.sum(normals * centroids, axis=1, keepdims=True)
    normals = np.where(dots < 0, -normals, normals)
    
    norms   = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms < 1e-14] = 1.0
    normals /= norms


    # Map colour to parsed velocity values
    color_vals = velocities
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=color_vals.min(), vmax=color_vals.max())
    face_colors = cmap(norm(color_vals)).copy()   # (N, 4) RGBA

    # Diffuse shading: light from upper-left
    light = np.array([0.6, 0.5, 1.0])
    light /= np.linalg.norm(light)
    shading = np.clip(normals @ light, 0, 1)
    # Mix ambient (0.35) + diffuse (0.65)
    face_colors[:, :3] *= (0.35 + 0.65 * shading[:, np.newaxis])

    return face_colors


def render(triangles, velocities, energy, output_path, cmap_name=DEFAULT_CMAP):
    """Render a single Fermi surface and save to output_path."""
    if not triangles:
        print(f"  [skip] No triangles for E={energy:+.2f}")
        return

    # Stack to (N, 3, 3) and centre BZ at origin
    verts = np.array(triangles)     # (N, 3, 3)
    verts -= 0.5                    # shift from [0,1] to [-0.5, 0.5]

    face_colors = compute_face_colors(verts, np.array(velocities), cmap_name)

    # ── Figure ──────────────────────────────────────────────────────────
    BG = '#07070f'
    fig = plt.figure(figsize=(9, 9), facecolor=BG)
    ax  = fig.add_subplot(111, projection='3d', facecolor=BG)
    ax.set_proj_type('persp')
    ax.set_box_aspect([1, 1, 1]) # Ensure equal aspect ratio

    poly = Poly3DCollection(verts.tolist(), zsort='average')
    poly.set_facecolor(face_colors)
    poly.set_edgecolor('none')
    ax.add_collection3d(poly)

    # ── Axis limits & labels ────────────────────────────────────────────
    all_pts = verts.reshape(-1, 3)
    span    = max(np.abs(all_pts).max() * 1.15, 0.35)
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_zlim(-span, span)

    tick_col = '#555577'
    ax.set_xlabel('$k_x \ (2\pi/a)$', color='#9999bb', labelpad=10, fontsize=11)
    ax.set_ylabel('$k_y \ (2\pi/a)$', color='#9999bb', labelpad=10, fontsize=11)
    ax.set_zlabel('$k_z \ (2\pi/a)$', color='#9999bb', labelpad=10, fontsize=11)
    ax.tick_params(colors=tick_col, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(tick_col)

    # Pane styling
    pane_col = '#0e0e1c'
    edge_col = '#1a1a2e'
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = True
        pane.set_facecolor(pane_col)
        pane.set_edgecolor(edge_col)
    ax.grid(True, color='#1c1c30', linewidth=0.5, linestyle='--')

    # ── Title ────────────────────────────────────────────────────────────
    sign   = '+' if energy >= 0 else ''
    n_tris = len(triangles)
    title  = (f'Fermi Surface  —  $E_F = {sign}{energy:.1f}$ eV\n'
              f'Colored by Fermi Velocity  ·  MC33  ·  {n_tris:,} triangles')
    ax.set_title(title, color='#ccccee', fontsize=12, pad=16)

    # Nice viewing angle
    ax.view_init(elev=22, azim=35)

    plt.tight_layout(pad=0.5)
    plt.savefig(output_path, dpi=160, bbox_inches='tight',
                facecolor=BG, edgecolor='none')
    plt.close(fig)
    print(f"  Saved  {os.path.basename(output_path)}   ({n_tris:,} triangles)")


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print(f"\nSearching for isosurf.dn files in:  {DATA_DIR}\n")
    files = sorted(glob.glob(os.path.join(DATA_DIR, "*isosurf.dn")))

    if not files:
        print("No *isosurf.dn files found — run the Fortran code first.")
        return

    for filepath in files:
        energy = energy_from_filename(filepath)
        label  = f"{energy:+.2f}" if energy is not None else "unknown"
        print(f"Parsing  {os.path.basename(filepath)}   (E = {label} eV)...")

        triangles, velocities = parse_isosurf_dn(filepath)
        print(f"  Found {len(triangles):,} triangles")

        cmap_name = CMAPS.get(energy, DEFAULT_CMAP)

        # Build output filename
        safe_label = label.replace('+', 'p').replace('-', 'm').replace('.', '_')
        out_name   = f"fermi_surface_E{safe_label}.png"
        out_path   = os.path.join(OUTPUT_DIR, out_name)

        render(triangles, velocities, energy if energy is not None else 0.0,
               out_path, cmap_name)

    print(f"\nAll images saved to: {OUTPUT_DIR}\n")


if __name__ == '__main__':
    main()
