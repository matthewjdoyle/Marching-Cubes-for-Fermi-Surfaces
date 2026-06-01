import os
import time
import subprocess
import numpy as np
import skimage.measure
import matplotlib.pyplot as plt

def modify_and_compile_fortran(npts, src_dir):
    f90_path = os.path.join(src_dir, "driver_modular_mc33.f90")
    bench_f90_path = os.path.join(src_dir, "_driver_modular_mc33_benchmark_gen.f90")
    exe_path = os.path.join(src_dir, "mc33_modular_benchmark.exe")
    
    with open(f90_path, "r") as f:
        content = f.read()

    # Modify npts and NE
    content = content.replace("npts=100", f"npts={npts}")
    content = content.replace("NE = 1", "NE = 1")
    
    # Timing variables are already present in driver_modular_mc33.f90
    
    loop_end = """
END DO
call cpu_time(t_end)
print *, "FORTRAN_TIME:", t_end - t_start
"""
    content = content.replace("END DO\ncall cpu_time(t_end)", loop_end)

    with open(bench_f90_path, "w") as f:
        f.write(content)
        
    # Compile
    compile_cmd = f'gfortran -J . -O3 dispersion_interface_mod.f90 mc33_core_mod.f90 brillouin_zone_mod.f90 fermi_bisection_mod.f90 "{bench_f90_path}" -o "{exe_path}"'
    subprocess.run(compile_cmd, shell=True, check=True, cwd=src_dir)
    return exe_path

def run_fortran_benchmark(exe_path, src_dir):
    project_root = os.path.dirname(src_dir)
    result = subprocess.run(f'"{exe_path}"', shell=True, capture_output=True, text=True, cwd=project_root)
    for line in result.stdout.splitlines():
        if "FORTRAN_TIME:" in line:
            return float(line.split("FORTRAN_TIME:")[1].strip())
    raise RuntimeError(f"Fortran time not found in output: {result.stdout}")

def time_python_marching_cubes(npts):
    kx = np.linspace(0, 1, npts + 1)
    ky = np.linspace(0, 1, npts + 1)
    kz = np.linspace(0, 1, npts + 1)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    E_grid = -2.0 * (np.cos(2*np.pi*(KX-0.5)) + np.cos(2*np.pi*(KY-0.5)) + np.cos(2*np.pi*(KZ-0.5)))
    
    t0 = time.perf_counter()
    grid = E_grid - 0.0
    try:
        verts, faces, _, _ = skimage.measure.marching_cubes(
            grid, level=0.0, spacing=(1.0 / npts, 1.0 / npts, 1.0 / npts)
        )
    except Exception:
        pass
    t1 = time.perf_counter()
    
    return t1 - t0

def main():
    mesh_sizes = [10, 20, 30, 40, 50, 60, 80, 100]
    repeats = 10
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    src_dir = os.path.join(project_root, "src")
    plots_dir = os.path.join(project_root, "Images", "Plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    fortran_means = []
    fortran_stds = []
    python_means = []
    python_stds = []
    
    for npts in mesh_sizes:
        print(f"Benchmarking mesh size: {npts}^3")
        
        # Python
        py_times = []
        for _ in range(repeats):
            py_times.append(time_python_marching_cubes(npts))
        py_mean = np.mean(py_times)
        py_std = np.std(py_times)
        python_means.append(py_mean)
        python_stds.append(py_std)
        print(f"  Python:  {py_mean:.4f} ± {py_std:.4f} s")
        
        # Fortran
        exe_path = modify_and_compile_fortran(npts, src_dir)
        fortran_times = []
        for _ in range(repeats):
            fortran_times.append(run_fortran_benchmark(exe_path, src_dir))
        f_mean = np.mean(fortran_times)
        f_std = np.std(fortran_times)
        fortran_means.append(f_mean)
        fortran_stds.append(f_std)
        print(f"  Fortran: {f_mean:.4f} ± {f_std:.4f} s")
        
    # Plotting
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.5,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.5,
    })
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    ax.errorbar(mesh_sizes, python_means, yerr=python_stds, fmt='o-', 
                color='#1f77b4', linewidth=2, capsize=5, capthick=2, markersize=8, label='skimage')
                
    ax.errorbar(mesh_sizes, fortran_means, yerr=fortran_stds, fmt='s-', 
                color='#d62728', linewidth=2, capsize=5, capthick=2, markersize=8, label='Fortran')
    
    ax.set_xlabel('Mesh Size ($N$ for $N^3$ grid)', labelpad=10, fontweight='bold')
    ax.set_ylabel('Execution Time (seconds)', labelpad=10, fontweight='bold')
    ax.set_title('Marching Cubes Algorithm Performance\nPython vs Fortran', pad=15, fontweight='bold')
    
    ax.grid(True, linestyle='--', color='#cccccc')
    ax.legend(framealpha=0.9, edgecolor='#999999', fancybox=True, loc='upper left')
    
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylim(bottom=1e-4)  # Set a reasonable lower limit for log scale
    
    # Custom ticks
    ax.set_xticks(mesh_sizes)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    
    ax.tick_params(direction='in', length=6, width=1.5, top=True, right=True, which='major')
    ax.tick_params(direction='in', length=3, width=1.0, top=True, right=True, which='minor')
    
    plot_path = os.path.join(plots_dir, "benchmark_marching_cubes.png")
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"\\nBenchmark plot saved to: {plot_path}")

    # Clean up generated benchmark files
    import glob
    project_root = os.path.dirname(src_dir)
    data_dir = os.path.join(project_root, "data")
    print("\nCleaning up temporary benchmark files...")
    temp_patterns = [
        os.path.join(src_dir, "_driver_modular_mc33_benchmark_gen.f90"),
        os.path.join(src_dir, "mc33_modular_benchmark.exe"),
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
