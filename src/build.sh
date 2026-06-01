#!/bin/bash
echo "Building MC33 Modular Benchmark..."

# Ensure we are compiling in the correct dependency order
gfortran -O3 \
    dispersion_interface_mod.f90 \
    mc33_core_mod.f90 \
    brillouin_zone_mod.f90 \
    fermi_bisection_mod.f90 \
    driver_modular_mc33.f90 \
    -o mc33_modular_benchmark.exe

if [ $? -eq 0 ]; then
    echo "Build successful: mc33_modular_benchmark.exe"
else
    echo "Build failed!"
    exit 1
fi
