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

# Shading logic from render_isosurfaces.py
def compute_face_colors(verts, cmap_name):
    centroids = verts.mean(axis=1)           # (N, 3)

    v1 = verts[:, 1, :] - verts[:, 0, :]
    v2 = verts[:, 2, :] - verts[:, 0, :]
    normals = np.cross(v1, v2)
    norms   = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms < 1e-14] = 1.0
    normals /= norms

    color_vals = centroids[:, 2]
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=color_vals.min(), vmax=color_vals.max())
    face_colors = cmap(norm(color_vals)).copy()   # (N, 4) RGBA

    light = np.array([0.6, 0.5, 1.0])
    light /= np.linalg.norm(light)
    shading = np.clip(normals @ light, 0, 1)
    face_colors[:, :3] *= (0.35 + 0.65 * shading[:, np.newaxis])

    return face_colors, normals

face_colors, normals = compute_face_colors(verts, 'inferno')

# Camera vector for elev=22, azim=35
# In matplotlib, azim is degrees counter-clockwise from the -y axis, elev is degrees above xy plane.
# So the camera position vector is:
elev, azim = np.radians(22), np.radians(35)
# Camera unit vector
v_cam = np.array([
    np.cos(elev) * np.cos(azim),
    np.cos(elev) * np.sin(azim),
    np.sin(elev)
])

# Let's check which way the normals point
dot_cam = np.sum(normals * v_cam, axis=1)

# Try culling faces pointing away from camera.
# If normals point outwards, then faces pointing towards the camera have dot_cam > 0.
# If normals point inwards, then faces pointing towards the camera have dot_cam < 0.
# We will do both and see which one is correct!

for name, mask in [
    ("outward", dot_cam > -0.05), # Keep faces pointing towards camera (allow slight overlap at edges)
    ("inward", dot_cam < 0.05)
]:
    verts_culled = verts[mask]
    colors_culled = face_colors[mask]

    BG = '#07070f'
    fig = plt.figure(figsize=(9, 9), facecolor=BG)
    ax  = fig.add_subplot(111, projection='3d', facecolor=BG)
    ax.set_proj_type('persp')

    poly = Poly3DCollection(verts_culled.tolist(), zsort='average', facecolors=colors_culled, edgecolors='none')
    ax.add_collection3d(poly)

    span = 0.5
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_zlim(-span, span)

    ax.view_init(elev=22, azim=35)

    os.makedirs("Images/Tests", exist_ok=True)
    plt.savefig(f"Images/Tests/test_cull_{name}.png", dpi=160, bbox_inches='tight',
                facecolor=BG, edgecolor='none')
    plt.close(fig)
    print(f"Saved Images/Tests/test_cull_{name}.png with {len(verts_culled)} triangles")
