# Windows PowerShell Script to compile and run the Fortran MC33 engine from the project root.
# Usage: .\run_fortran.ps1

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "  Fortran MC33 Fermi Surface Engine  " -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

# 1. Check if gfortran is installed
if (!(Get-Command gfortran -ErrorAction SilentlyContinue)) {
    Write-Error "gfortran is not installed or not in your PATH. Please install a Fortran compiler (e.g., MSYS2/MinGW-w64)."
    Exit 1
}

# 2. Compile the Fortran codebase
Write-Host "Compiling Fortran files..." -ForegroundColor Yellow
$srcDir = "src"
$exePath = "src\mc33_modular.exe"

# gfortran command compiling in correct dependency order
# -J src/ outputs compiler module files (.mod) directly to src/ to keep root clean
& gfortran -O3 `
    "$srcDir\dispersion_interface_mod.f90" `
    "$srcDir\mc33_core_mod.f90" `
    "$srcDir\brillouin_zone_mod.f90" `
    "$srcDir\fermi_bisection_mod.f90" `
    "$srcDir\driver_modular_mc33.f90" `
    -J "$srcDir" `
    -o "$exePath"

if ($LASTEXITCODE -ne 0) {
    Write-Error "Compilation failed!"
    Exit 1
}
Write-Host "Compilation successful: $exePath" -ForegroundColor Green

# 3. Run the compiled executable inside the src folder
Write-Host "Running MC33 Modular Solver..." -ForegroundColor Yellow
New-Item -ItemType Directory -Force -Path "$srcDir\data" | Out-Null
Push-Location $srcDir
& ".\mc33_modular.exe"
Pop-Location

# 4. Clean up and move generated data files to data/ directory
Write-Host "Organizing output files to data/..." -ForegroundColor Yellow
New-Item -ItemType Directory -Force -Path "data" | Out-Null

$movedFilesCount = 0
Get-ChildItem -Path "$srcDir\data" -Filter "*isosurf*" | ForEach-Object {
    Move-Item -Force $_.FullName -Destination "data\"
    $movedFilesCount++
}
Get-ChildItem -Path "$srcDir\data" -Filter "*isodown*" | ForEach-Object {
    Move-Item -Force $_.FullName -Destination "data\"
    $movedFilesCount++
}
if (Test-Path "$srcDir\data\DOS.txt") {
    Move-Item -Force "$srcDir\data\DOS.txt" -Destination "data\"
    $movedFilesCount++
}

Write-Host "Success! Moved $movedFilesCount output datasets to data/ folder." -ForegroundColor Green
Write-Host "Fortran execution complete." -ForegroundColor Green
