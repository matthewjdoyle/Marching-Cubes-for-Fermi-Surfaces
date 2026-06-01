# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:10:34 2019

@author: matth
"""
import os
import glob

# ── Path configuration ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

# Try to find a kmesh.isosurf file in the data/ directory first, otherwise fallback to local
kmesh_files = sorted(glob.glob(os.path.join(DATA_DIR, "*kmesh.isosurf")))
if kmesh_files:
    filepath = kmesh_files[-1]  # Use the latest generated file
    print(f"Loading mesh data from: {filepath}")
else:
    filepath = "kmesh.isosurf"
    print(f"No mesh data found in {DATA_DIR}. Falling back to current working directory: {filepath}")

with open(filepath, "r") as f:
    lines = f.readlines()

kx = []
ky = []
kz = []
kweight = []

for x in lines:
    parts = x.split()
    if len(parts) >= 4:
        try:
            kx.append(float(parts[0]))
            ky.append(float(parts[1]))
            kz.append(float(parts[2]))
            kweight.append(float(parts[3]))
        except ValueError:
            continue

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.tri as mtri

tri = mtri.Triangulation(kx, kweight)
tri2 = tri.triangles

fig = plt.figure(figsize = (12,9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(kx, ky, kz, s = 15, depthshade = True, c = kweight)
#ax.plot_trisurf(kx, ky, kz, 
 #               cmap='Spectral', lw=1)
#ax.set_title("E = %f" % e)
#ax.set_axis_off()
plt.show()
