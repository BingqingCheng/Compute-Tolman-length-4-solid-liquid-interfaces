# Compute-Tolman-length-4-solid-liquid-interfaces
The LAMMPS input files and data analysis scripts used in the paper 

"Communication: Computing the Tolman length for solid-liquid interfaces"

Bingqing Cheng and Michele Ceriotti

The Journal of Chemical Physics 148 (23), 231102

## 100-template/110-template
for running interface pinning simulations of LJ system with a solid-liquid interface along 100/110 lattice direction.

## \*-solid-template/\*-liquid-template
for running simulations of bulk solid or bulk liquid LJ systems, in order to get the the amplitudes of fluctuations in the bulk phases.

## quick-direct-ft.py
A python script to extract the amplitudes of Fourier modes of the whole simulation box.

## average-ft.ipynb
Obtain interface stiffness from the amplitudes of Fourier modes.

## q6/m6
Using the postprocessing functionality of PLUMED, compute Q6 and locally-averaged Q6 order parameters of each atom in the system. 

