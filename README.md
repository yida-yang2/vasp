# vasp
Convergence.f90

Description

convergence.f90 is a Fortran program designed to extract energy, force, stress, and computational time from VASP's OUTCAR file. It provides a detailed convergence analysis of structural relaxation calculations.

Features

Extracts total energy from OUTCAR.

Computes the maximum atomic force for each step.

Extracts and converts stress tensor to GPa.

Reads and records computational time for each iteration.

Checks convergence based on the EDIFFG criterion.

Outputs results in a structured format.

Compilation

To compile the program, use a Fortran compiler such as gfortran:

 gfortran -o convergence convergence.f90

Usage

Run the program in a directory containing OUTCAR:

 ./convergence

Output Format

The program prints the following table to the terminal:

 Step |   E (eV)  |  dE (eV)   | Fmax (eV/\AA) | Smax (GPa) | Time (s) | Average drift | Convergence |
-------------------------------------------------------------------------------------------------------------
   1   -1234.567   --------          0.0234       1.5678       120.45        0.00001         YES
   2   -1234.890   -0.323            0.0123       0.9876       110.32        0.00002          NO
EDIFFG = -0.01
Where:

E (eV): Total energy.

dE (eV): Energy difference between steps.

Fmax (eV/\AA): Maximum force.

Smax (GPa): Maximum stress.

Time (s): Computation time per step.

Average drift: Average atomic drift.

Convergence: "YES" if force is below EDIFFG, otherwise "NO".

EDIFFG : Convergence standard 

Dependencies

Fortran compiler (e.g., gfortran).

OUTCAR file from VASP.

Notes

The program assumes OUTCAR follows a standard VASP format.

Author

Yida Yang
