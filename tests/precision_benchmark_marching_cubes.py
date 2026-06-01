import os
import glob
import subprocess
import numpy as np
import skimage.measure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def compute_flat_colors(velocities, cmap_name='jet', vmin=None, vmax=None):
    cmap = plt.get_cmap(cmap_name)
    if vmin is None:
        vmin = velocities.min()
    if vmax is None:
        vmax = velocities.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    return cmap(norm(velocities)).copy()

def parse_isosurf_dn(filepath):
    with open(filepath, 'r') as f:
        raw = f.read()
    raw = raw.replace('D+', 'E+').replace('D-', 'E-')
    raw = raw.replace('\r', '')
    triangles = []
    velocities = []
    blocks = raw.split('\n\n\n')
    for block in blocks:
        block = block.strip()
        if not block:
            continue
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
                triangles.append(np.array([p1, p2, p3]))
                
                # Check for 4th column (velocity); default to average z-coordinate if not present
                v1 = p1_raw[3] if len(p1_raw) > 3 else p1[2]
                v2 = p2_raw[3] if len(p2_raw) > 3 else p2[2]
                v3 = p3_raw[3] if len(p3_raw) > 3 else p3[2]
                velocities.append((v1 + v2 + v3) / 3.0)
        except (ValueError, IndexError):
            continue
    return triangles, velocities



def exact_energy(verts):
    kx, ky, kz = verts[:, 0], verts[:, 1], verts[:, 2]
    # Evaluate the exact tight binding analytical function. Both Python and Fortran assume range is [0, 1] and shifted by 0.5
    return -2.0 * (np.cos(2*np.pi*(kx-0.5)) + np.cos(2*np.pi*(ky-0.5)) + np.cos(2*np.pi*(kz-0.5)))

def get_python_mae(npts):
    kx = np.linspace(0, 1, npts + 1)
    ky = np.linspace(0, 1, npts + 1)
    kz = np.linspace(0, 1, npts + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    E_grid = -2.0 * (np.cos(2*np.pi*(KX-0.5)) + np.cos(2*np.pi*(KY-0.5)) + np.cos(2*np.pi*(KZ-0.5)))
    
    try:
        verts, faces, _, _ = skimage.measure.marching_cubes(
            E_grid, level=0.0, spacing=(1.0 / npts, 1.0 / npts, 1.0 / npts)
        )
        
        energies = exact_energy(verts)
        mae = np.mean(np.abs(energies - 0.0))
        return mae
    except Exception as e:
        print(f"Python error: {e}")
        return np.nan

def modify_and_compile_fortran_precision(npts, src_dir):
    f90_path = os.path.join(src_dir, "driver_modular_mc33.f90")
    bench_f90_path = os.path.join(src_dir, "_driver_modular_mc33_precision_gen.f90")
    exe_path = os.path.join(src_dir, "mc33_modular_precision.exe")
    
    with open(f90_path, "r") as f:
        content = f.read()

    # Modify npts and NE
    content = content.replace("npts=100", f"npts={npts}")
    content = content.replace("NE = 1", "NE = 1")
    
    # Change EF to 0.0 directly to guarantee we test the 0.0 isosurface
    content = content.replace("ef = cmplx(-4.0d0 + dble(i-1)*2.0d0, 0.0d0, 8)", "ef = cmplx(0.0d0, 0.0d0, 8)")
    
    with open(bench_f90_path, "w") as f:
        f.write(content)
        
    compile_cmd = f'gfortran -J . -O3 dispersion_interface_mod.f90 mc33_core_mod.f90 brillouin_zone_mod.f90 fermi_bisection_mod.f90 "{bench_f90_path}" -o "{exe_path}"'
    subprocess.run(compile_cmd, shell=True, check=True, cwd=src_dir)
    return exe_path

def get_fortran_mae(npts, src_dir):
    exe_path = modify_and_compile_fortran_precision(npts, src_dir)
    project_root = os.path.dirname(src_dir)
    data_dir = os.path.join(project_root, "data")
    
    # Clear any old isosurf files in data/
    old_files = glob.glob(os.path.join(data_dir, "*kmesh.isosurf"))
    for f in old_files:
        os.remove(f)
        
    subprocess.run(f'"{exe_path}"', shell=True, capture_output=True, text=True, cwd=project_root)
    
    # Find the generated isosurf file in data/
    new_files = glob.glob(os.path.join(data_dir, "*kmesh.isosurf"))
    if not new_files:
        print("Fortran error: isosurf file not found")
        return np.nan
        
    iso_file = new_files[0]
    
    try:
        data = np.loadtxt(iso_file)
        if len(data) == 0:
            return np.nan
            
        # The first 3 columns are kx, ky, kz. 
        verts = data[:, :3]
        energies = exact_energy(verts)
        mae = np.mean(np.abs(energies - 0.0))
        return mae
    except Exception as e:
        print(f"Fortran error reading file: {e}")
        return np.nan

def main():
    mesh_sizes = [10, 20, 30, 40, 50, 60, 80, 100]
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    src_dir = os.path.join(project_root, "src")
    plots_dir = os.path.join(project_root, "Images", "Plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    python_maes = []
    fortran_maes = []
    
    for npts in mesh_sizes:
        print(f"Benchmarking precision for mesh size: {npts}^3")
        
        py_mae = get_python_mae(npts)
        python_maes.append(py_mae)
        print(f"  Python MAE:  {py_mae:.4e} eV")
        
        f_mae = get_fortran_mae(npts, src_dir)
        fortran_maes.append(f_mae)
        print(f"  Fortran MAE: {f_mae:.4e} eV")
        
    # Plotting
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.5,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.5,
    })
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    ax.plot(mesh_sizes, python_maes, 'o-', color='#1f77b4', linewidth=2, markersize=8, label='skimage')
    ax.plot(mesh_sizes, fortran_maes, 's-', color='#d62728', linewidth=2, markersize=8, label='Fortran')
    
    ax.set_xlabel('Mesh Size ($N$ for $N^3$ grid)', labelpad=10, fontweight='bold')
    ax.set_ylabel('Mean Absolute Error (eV)', labelpad=10, fontweight='bold')
    ax.set_title('Precision Comparison: Marching Cubes Algorithm\nAnalytical Error vs Target Fermi Surface', pad=15, fontweight='bold')
    
    ax.grid(True, linestyle='--', color='#cccccc')
    ax.legend(framealpha=0.9, edgecolor='#999999', fancybox=True, loc='center right')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_ylim(bottom=1e-13, top=1e-1)

    # Custom ticks for X
    ax.set_xticks(mesh_sizes)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    
    ax.tick_params(direction='in', length=6, width=1.5, top=True, right=True, which='major')
    ax.tick_params(direction='in', length=3, width=1.0, top=True, right=True, which='minor')
    
    plot_path = os.path.join(plots_dir, "benchmark_precision_marching_cubes.png")
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"\\nPrecision plot saved to: {plot_path}")
    
    for n in [4, 8, 16, 64]:
        print(f"\nGenerating comparison plot for N={n}...")
        plot_skimage_vs_fortran(n, src_dir, plots_dir)

def plot_skimage_vs_fortran(npts, src_dir, plots_dir):
    kx = np.linspace(0, 1, npts + 1)
    ky = np.linspace(0, 1, npts + 1)
    kz = np.linspace(0, 1, npts + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    E_grid = -2.0 * (np.cos(2*np.pi*(KX-0.5)) + np.cos(2*np.pi*(KY-0.5)) + np.cos(2*np.pi*(KZ-0.5)))
    
    py_verts, py_faces, _, _ = skimage.measure.marching_cubes(
        E_grid, level=0.0, spacing=(1.0 / npts, 1.0 / npts, 1.0 / npts)
    )
    py_verts -= 0.5  # Center at origin
    py_verts_tri = py_verts[py_faces]
    
    # Filter degenerate triangles (which have zero area and skew color ranges)
    py_non_deg = []
    for t in py_verts_tri:
        d12 = np.linalg.norm(t[0] - t[1])
        d23 = np.linalg.norm(t[1] - t[2])
        d31 = np.linalg.norm(t[2] - t[0])
        if not (d12 < 1e-5 or d23 < 1e-5 or d31 < 1e-5):
            py_non_deg.append(t)
    py_verts_tri_clean = np.array(py_non_deg)
    
    # Python (skimage) analytical velocities evaluated at centroids of clean triangles
    py_centroids = py_verts_tri_clean.mean(axis=1)
    grad_x = 4.0 * np.pi * np.sin(2.0 * np.pi * py_centroids[:, 0])
    grad_y = 4.0 * np.pi * np.sin(2.0 * np.pi * py_centroids[:, 1])
    grad_z = 4.0 * np.pi * np.sin(2.0 * np.pi * py_centroids[:, 2])
    py_grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    py_grad_mag[py_grad_mag < 1e-8] = 1e-8
    
    project_root = os.path.dirname(src_dir)
    data_dir = os.path.join(project_root, "data")
    exe_path = modify_and_compile_fortran_precision(npts, src_dir)
    
    # Remove old isosurf.dn and kmesh.isosurf files in data_dir
    for f in glob.glob(os.path.join(data_dir, "*isosurf.dn")) + glob.glob(os.path.join(data_dir, "*kmesh.isosurf")):
        try:
            os.remove(f)
        except Exception:
            pass
            
    subprocess.run(f'"{exe_path}"', shell=True, capture_output=True, text=True, cwd=project_root)
    
    new_files = glob.glob(os.path.join(data_dir, "*isosurf.dn"))
    if not new_files:
        print("Error: Fortran isosurf.dn output not found for comparison plot.")
        return
        
    triangles, velocities = parse_isosurf_dn(new_files[0])
    if not triangles:
        print("Error: Could not parse triangles from Fortran output.")
        return
        
    f_verts_tri = np.array(triangles)
    f_verts_tri -= 0.5  # Center at origin
    
    # Filter degenerate triangles for Fortran
    f_non_deg = []
    for t in f_verts_tri:
        d12 = np.linalg.norm(t[0] - t[1])
        d23 = np.linalg.norm(t[1] - t[2])
        d31 = np.linalg.norm(t[2] - t[0])
        if not (d12 < 1e-5 or d23 < 1e-5 or d31 < 1e-5):
            f_non_deg.append(t)
    f_verts_tri_clean = np.array(f_non_deg)
    
    # Fortran (MC33) analytical velocities evaluated at centroids of clean triangles
    f_centroids = f_verts_tri_clean.mean(axis=1)
    grad_x = 4.0 * np.pi * np.sin(2.0 * np.pi * f_centroids[:, 0])
    grad_y = 4.0 * np.pi * np.sin(2.0 * np.pi * f_centroids[:, 1])
    grad_z = 4.0 * np.pi * np.sin(2.0 * np.pi * f_centroids[:, 2])
    f_grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    f_grad_mag[f_grad_mag < 1e-8] = 1e-8
    
    # Define shared limits based on clean non-degenerate triangles
    vmin = min(py_grad_mag.min(), f_grad_mag.min())
    vmax = max(py_grad_mag.max(), f_grad_mag.max())
    if vmax - vmin < 1e-4:
        vmin -= 0.5
        vmax += 0.5
    
    py_colors = compute_flat_colors(py_grad_mag, 'jet', vmin=vmin, vmax=vmax)
    f_colors = compute_flat_colors(f_grad_mag, 'jet', vmin=vmin, vmax=vmax)
    
    BG = '#07070f'
    fig = plt.figure(figsize=(18, 9), facecolor=BG)
    
    # Python (skimage) Panel
    ax1 = fig.add_subplot(121, projection='3d', facecolor=BG)
    ax1.set_proj_type('persp')
    poly1 = Poly3DCollection(py_verts_tri_clean.tolist(), zsort='average', facecolors=py_colors, edgecolors='none', shade=False)
    ax1.add_collection3d(poly1)
    
    ax1.set_box_aspect([1, 1, 1])
    span = 0.5
    ax1.set_xlim(-span, span)
    ax1.set_ylim(-span, span)
    ax1.set_zlim(-span, span)
    ax1.set_axis_off()
    ax1.view_init(elev=30, azim=45)
    ax1.set_title('scikit-image (Python)', color='#ccccee', fontsize=14, pad=10, fontweight='bold')
    
    # Add colorbar for Python panel
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'), norm=norm)
    sm.set_array([])
    cb1 = fig.colorbar(sm, ax=ax1, shrink=0.55, pad=0.01)
    cb1.ax.yaxis.set_tick_params(color='#9999bb', labelcolor='#9999bb')
    cb1.outline.set_edgecolor('#555577')
    cb1.set_label('Fermi Velocity', color='#9999bb', fontsize=10, labelpad=8)
    
    # Fortran (MC33) Panel
    ax2 = fig.add_subplot(122, projection='3d', facecolor=BG)
    ax2.set_proj_type('persp')
    poly2 = Poly3DCollection(f_verts_tri_clean.tolist(), zsort='average', facecolors=f_colors, edgecolors='none', shade=False)
    ax2.add_collection3d(poly2)
    
    ax2.set_box_aspect([1, 1, 1])
    ax2.set_xlim(-span, span)
    ax2.set_ylim(-span, span)
    ax2.set_zlim(-span, span)
    ax2.set_axis_off()
    ax2.view_init(elev=30, azim=45)
    ax2.set_title('MC33 Algorithm (Fortran)', color='#ccccee', fontsize=14, pad=10, fontweight='bold')
    
    # Add colorbar for Fortran panel
    cb2 = fig.colorbar(sm, ax=ax2, shrink=0.55, pad=0.01)
    cb2.ax.yaxis.set_tick_params(color='#9999bb', labelcolor='#9999bb')
    cb2.outline.set_edgecolor('#555577')
    cb2.set_label('Fermi Velocity', color='#9999bb', fontsize=10, labelpad=8)
    
    fig.text(0.5, 0.95, f"Fermi Surface Reconstruction Comparison (Mesh: {npts}^3)", 
             color='#ffaa66', fontsize=18, ha='center', fontweight='bold')
    
    plot_path = os.path.join(plots_dir, f"skimage_vs_fortran_N{npts}.png")
    plt.savefig(plot_path, bbox_inches='tight', dpi=300, facecolor=BG)
    plt.close()
    print(f"\nComparison plot saved to: {plot_path}")

    # Clean up generated precision benchmark files
    print("\nCleaning up temporary benchmark files...")
    temp_patterns = [
        os.path.join(src_dir, "_driver_modular_mc33_precision_gen.f90"),
        os.path.join(src_dir, "mc33_modular_precision.exe"),
        os.path.join(data_dir, "*isodown*"),
        os.path.join(data_dir, "*isosurf*"),
        os.path.join(data_dir, "DOS.txt"),
    ]
    for pattern in temp_patterns:
        for f in glob.glob(pattern):
            try:
                os.remove(f)
            except Exception:
                pass

if __name__ == '__main__':
    main()
