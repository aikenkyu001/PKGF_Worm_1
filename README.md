# PKGF-Worm: Deterministic Connectome Dynamics Model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19344285.svg)](https://doi.org/10.5281/zenodo.19344285)

PKGF-Worm is a research project that maps the full connectome of *C. elegans* (302 neurons) onto the manifold structure of Parallel Key Geometric Flow (PKGF), a differential geometric theory of intelligence, to simulate deterministic neural and kinetic dynamics.

This model demonstrates the emergence of **sustained non-equilibrium attractors with goal-directed bias** driven solely by geometric necessity on a manifold, without reliance on probabilistic frameworks or statistical optimization.

## Core Concepts
- **Deterministic Emergence**: Spontaneous movement (undulation) emerges from infinitesimal initial asymmetries and geometric connectivity without the use of random number generators.
- **Cross-Platform Synchronization**: The project features a Python implementation for research and prototyping and a Fortran 95 implementation for high-precision, rigorous computation, designed to maintain bit-level synchronization.
- **Geometric Intelligence**: Intelligence is redefined as a process of seeking stable attractors within warped manifolds, quantified by the numerical **Intelligence Metric** ($\mathcal{I}$).

## Directory Structure
```text
/PKGF_Worm/
├── Docs/               # Academic papers (EN/JP), research plans, and experimental reports.
├── Scripts/            # Simulation source code.
│   ├── pkgf_worm_unified.py    # Unified Python simulator (NumPy-based).
│   ├── main.f90                # Fortran main entry point.
│   ├── pkgf_worm_engine_mod.f90 # PKGF flow equation engine.
│   ├── pkgf_worm_node_mod.f90   # Node (neuron) definitions.
│   └── ...                     # Mathematical utilities and constant modules.
├── References/         # Connectome data (JSON) and theoretical references.
├── Logs/               # Execution logs for cross-implementation behavior comparison.
└── README.md           # This file.
```

## Build and Execution

### 1. Python Implementation
A pure NumPy-based implementation with minimal dependencies.
```bash
cd Scripts
python pkgf_worm_unified.py
```
Detailed internal states and coordinates are output to `Logs/worm_log_python.txt`.

### 2. Fortran Implementation
Uses `gfortran` for rigorous numerical stability and reproducibility.
```bash
cd Scripts
# Compile modules and the main program
gfortran -c pkgf_constants.f90 pkgf_math_utils.f90 pkgf_connectome_data.f90
gfortran -c pkgf_worm_node_mod.f90 pkgf_environment_mod.f90 pkgf_worm_engine_mod.f90
gfortran main.f90 *.o -o pkgf_worm_fortran

# Run
./pkgf_worm_fortran
```
Generates `Logs/worm_log_fortran.txt` upon completion.

### 3. Verification of Determinism
By comparing logs from Python and Fortran, you can verify that the system follows identical "thought and movement processes" across different environments.
```bash
# Use plot_results.py to visualize trajectory bifurcation and Lyapunov-like sensitivity
python plot_results.py
```

## Key Scientific Documentation
For detailed mathematical formulations and experimental results, please refer to:
- [Academic Paper (English)](Docs/PKGF_Worm_Academic_Paper_en.md)
- [学術論文 (Japanese)](Docs/PKGF_Worm_Academic_Paper_jp.md)

## Author
Fumio Miyata (2026)
