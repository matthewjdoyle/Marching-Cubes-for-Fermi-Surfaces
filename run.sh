#!/bin/bash
# Linux Terminal Script to compile and run the Fortran MC33 engine from the project root.
# Usage: chmod +x run_fortran.sh && ./run_fortran.sh

# Colors for nice output
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${CYAN}=========================================${NC}"
echo -e "${CYAN}  Fortran MC33 Fermi Surface Engine      ${NC}"
echo -e "${CYAN}=========================================${NC}"

# 1. Check if gfortran is installed
if ! command -v gfortran &> /dev/null; then
    echo -e "${RED}Error: gfortran is not installed or not in your PATH. Please install gcc-gfortran.${NC}"
    exit 1
fi

# 2. Compile the Fortran codebase
echo -e "${YELLOW}Compiling Fortran files...${NC}"
SRC_DIR="src"
EXE_PATH="src/mc33_modular"

# gfortran command compiling in correct dependency order
# -J src/ outputs compiler module files (.mod) directly to src/ to keep root clean
gfortran -O3 \
    "$SRC_DIR/dispersion_interface_mod.f90" \
    "$SRC_DIR/mc33_core_mod.f90" \
    "$SRC_DIR/brillouin_zone_mod.f90" \
    "$SRC_DIR/fermi_bisection_mod.f90" \
    "$SRC_DIR/driver_modular_mc33.f90" \
    -J "$SRC_DIR" \
    -o "$EXE_PATH"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Compilation successful: $EXE_PATH${NC}"
else
    echo -e "${RED}Error: Compilation failed!${NC}"
    exit 1
fi

# 3. Run the compiled executable inside the src folder
echo -e "${YELLOW}Running MC33 Modular Solver...${NC}"
mkdir -p "$SRC_DIR/data"
cd "$SRC_DIR" || exit 1
./mc33_modular
cd ..

# 4. Clean up and move generated data files to data/ directory
echo -e "${YELLOW}Organizing output files to data/...${NC}"
mkdir -p data

moved_count=0
for f in "$SRC_DIR"/data/*isosurf* "$SRC_DIR"/data/*isodown* "$SRC_DIR"/data/DOS.txt; do
    if [ -f "$f" ]; then
        mv -f "$f" data/
        ((moved_count++))
    fi
done

echo -e "${GREEN}Success! Moved $moved_count output datasets to data/ folder.${NC}"
echo -e "${GREEN}Fortran execution complete.${NC}"
