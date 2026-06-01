#!/usr/bin/env python3
"""
Multi-Lattice Density of States (DOS) Plot
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless rendering backend
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "Images", "Plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

fig_width_inch = 7.2
fig_height_inch = 6.0

plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.4,
    'grid.alpha': 0.4,
    'figure.dpi': 300,
})

def load_dos_data(lattice_name):
    """
    Parses the energy and DOS columns from the Fortran-generated output files.
    """
    filepath = os.path.join(DATA_DIR, f"DOS_{lattice_name}.txt")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Missing dataset for {lattice_name} at: {filepath}")
    
    data = np.loadtxt(filepath)
    energies = data[:, 0]
    dos = data[:, 1]
    return energies, dos

def create_4panel_plot():
    fig, axs = plt.subplots(2, 2, figsize=(fig_width_inch, fig_height_inch), sharex=False, sharey=False)
    plt.subplots_adjust(wspace=0.25, hspace=0.30)
    
    # Curated physical colours (non-generic, tailored for premium solid-state publication)
    colours = {
        'sc': '#2c3e6b',   # Deep slate blue
        'bcc': '#b83b30',  # Crimson red
        'fcc': '#306b34',  # Forest green
        'hex': '#b87730'   # Amber orange
    }

    # ── Panel 1: Simple Cubic (SC) ──────────────────────────────────────────
    ax = axs[0, 0]
    e_sc, d_sc = load_dos_data('sc')
    
    # Plot Fortran MC33 solver data
    ax.plot(e_sc, d_sc, color=colours['sc'], linewidth=1.2, marker='o', markersize=2.5,
            markevery=1, zorder=2)
    
    # Setup inward ticks and round-number boundaries
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.grid(True, linestyle=':', color='#cccccc')
    ax.set_xlim(-6.0, 6.0)
    ax.set_ylim(0.0, 0.20)
    ax.set_xticks([-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0])
    ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20])
    ax.set_ylabel(r'$D(E)$')
    ax.set_xlabel(r'$E$ (eV)')
    ax.text(-5.5, 0.18, '(a) Simple Cubic (SC)', fontweight='bold', fontsize=9)

    # ── Panel 2: Body-Centred Cubic (BCC) ───────────────────────────────────
    ax = axs[0, 1]
    e_bcc, d_bcc = load_dos_data('bcc')
    
    # Plot Fortran MC33 solver data
    ax.plot(e_bcc, d_bcc, color=colours['bcc'], linewidth=1.2, marker='o', markersize=2.5,
            markevery=1, zorder=2)
    
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.grid(True, linestyle=':', color='#cccccc')
    ax.set_xlim(-8.0, 8.0)
    ax.set_ylim(0.0, 0.50)  # Increased max y limit from 0.40 to 0.50
    ax.set_xticks([-8.0, -4.0, 0.0, 4.0, 8.0])
    ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])  # Starts and ends on round boundaries
    ax.set_ylabel(r'$D(E)$')
    ax.set_xlabel(r'$E$ (eV)')
    ax.text(-7.3, 0.45, '(b) Body-Centred Cubic (BCC)', fontweight='bold', fontsize=9)

    # ── Panel 3: Face-Centred Cubic (FCC) ───────────────────────────────────
    ax = axs[1, 0]
    e_fcc, d_fcc = load_dos_data('fcc')
    
    # Plot Fortran MC33 solver data
    ax.plot(e_fcc, d_fcc, color=colours['fcc'], linewidth=1.2, marker='o', markersize=2.5,
            markevery=1, zorder=2)
    
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.grid(True, linestyle=':', color='#cccccc')
    ax.set_xlim(-12.0, 4.0)
    ax.set_ylim(0.0, 0.30)
    ax.set_xticks([-12.0, -8.0, -4.0, 0.0, 4.0])
    ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax.set_ylabel(r'$D(E)$')
    ax.set_xlabel(r'$E$ (eV)')
    ax.text(-11.0, 0.27, '(c) Face-Centred Cubic (FCC)', fontweight='bold', fontsize=9)

    # ── Panel 4: Simple Hexagonal (HEX) ─────────────────────────────────────
    ax = axs[1, 1]
    e_hex, d_hex = load_dos_data('hex')
    
    # Plot Fortran MC33 solver data
    ax.plot(e_hex, d_hex, color=colours['hex'], linewidth=1.2, marker='o', markersize=2.5,
            markevery=1, zorder=2)
    
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.grid(True, linestyle=':', color='#cccccc')
    ax.set_xlim(-8.0, 4.0)
    ax.set_ylim(0.0, 0.30)  # Increased max y limit from 0.20 to 0.30
    ax.set_xticks([-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0])
    ax.set_yticks([0.0, 0.1, 0.2, 0.3])  # Starts and ends on round boundaries
    ax.set_ylabel(r'$D(E)$')
    ax.set_xlabel(r'$E$ (eV)')
    ax.text(-7.3, 0.27, '(d) Simple Hexagonal (HEX)', fontweight='bold', fontsize=9)

    # Save output to both standard high-DPI PNG and vector PDF for publication
    png_path = os.path.join(OUTPUT_DIR, "multi_lattice_dos_4panel.png")
    pdf_path = os.path.join(OUTPUT_DIR, "multi_lattice_dos_4panel.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    
    print(f"Nature-quality 4-panel DOS plot successfully saved to:")
    print(f"  PNG: {png_path}")
    print(f"  PDF: {pdf_path}")

if __name__ == "__main__":
    create_4panel_plot()
