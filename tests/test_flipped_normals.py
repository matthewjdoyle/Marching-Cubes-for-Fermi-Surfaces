import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

filepath = "data/fermi_surface_neg2.000_isosurf.dn"
print(f"Reading file: {filepath}")

# Parse triangles
with open(filepath, 'r') as f:
    raw = f.read().replace('D+', 'E+').replace('D-', 'E-').replace('\r', '')

triangles = []
for block in raw.split('\n\n\n'):
    block = block.strip()
    if not block:
        continue
    halves = [h.strip() for h in block.split('\n\n') if h.strip()]
    if len(halves) < 2:
        continue
    try:
        p1_lines  = [l for l in halves[0].split('\n') if l.strip()]
        p23_lines = [l for l in halves[1].split('\n') if l.strip()]
        p1 = np.array([float(x) for x in p1_lines[0].split()])
        p2 = np.array([float(x) for x in p23_lines[0].split()])
        p3 = np.array([float(x) for x in p23_lines[1].split()])
        triangles.append(np.array([p1, p2, p3]))
    except:
        continue

verts = np.array(triangles) - 0.5
centroids = verts.mean(axis=1)           # (N, 3)

# Shading logic with fixed normals!
def compute_face_colors_fixed(verts, cmap_name):
    centroids = verts.mean(axis=1)           # (N, 3)

    v1 = verts[:, 1, :] - verts[:, 0, :]
    v2 = verts[:, 2, :] - verts[:, 0, :]
    normals = np.cross(v1, v2)
    
    # Flip normals that point inwards!
    # A normal points inwards if its dot product with the centroid is negative.
    dots = np.sum(normals * centroids, axis=1, keepdims=True)
    normals = np.where(dots < 0, -normals, normals)
    
    norms   = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms < 1e-14] = 1.0
    normals /= norms

    color_vals = centroids[:, 2]
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

face_colors = compute_face_colors_fixed(verts, 'inferno')

BG = '#07070f'
fig = plt.figure(figsize=(9, 9), facecolor=BG)
ax  = fig.add_subplot(111, projection='3d', facecolor=BG)
ax.set_proj_type('persp')

poly = Poly3DCollection(verts.tolist(), zsort='average')
poly.set_facecolor(face_colors)
poly.set_edgecolor('none')
ax.add_collection3d(poly)

span = 0.5
ax.set_xlim(-span, span)
ax.set_ylim(-span, span)
ax.set_zlim(-span, span)
ax.view_init(elev=22, azim=35)

os.makedirs("Images/Tests", exist_ok=True)
plt.savefig("Images/Tests/test_flipped_normals.png", dpi=160, bbox_inches='tight',
            facecolor=BG, edgecolor='none')
plt.close(fig)
print("Saved Images/Tests/test_flipped_normals.png")
