# TimeESR v1.0.2
An STM-ESR code solving a QME in the time domain. TimeESR is given as is and subject to gnu 3.0 CopyLeft. 
Enjoy! TimeESR team, February 17, 2024.

## Intended Use
The code is intended for computing ESR in an STM junction, where the driving is produced by the modulation of the electron hopping integrals. It models one electronic orbital coupled to a complex system of spins. It is intended for (i) obtaining continuous-wave ESR spectra (our Floquet code is better suited for this), (ii) doing quantum operations driving spins with short pulses, and (iii) using multiple driving frequencies as is done in state-of-the-art STM-ESR.

## Compilation Instructions
1. Prerequisities: LAPACK, BLAS, Fortran compiler to run the program. Gnuplot, Python, xdg-open, and qutip are recommended if using the provided plot and movie scripts for any example in the `./example` folder. 
2. Review `./src/make.sys` to make sure your compiler and compilation options are set properly. 
3. Go to `./src` and execute `make`.

## Examples
The directory `./example` contains a few examples to get started. For each, execute `./run.sh` in the example directory, which cleans the `./src` folder, makes the TimeESR exectuible, runs the program, and then makes plots. Output plots can include population of states, spin projections, and current as a function of time. Post-processed plots of representative Bloch spheres, and movies of appropriate Bloch sphere projections as a function of time, are created for some examples. A full description of examples can be found in `Examples.md`.

## Input Files
For both input files, see the provided examples for descriptions of each option.
1. `H_QD.input`: Descriptors for the Hamiltonian, with the first spin is a driven $S=1/2$ electron. The first site has four states: single-occupied up spin $z$ direction, single-occupied down spin $z$ direction, zero occupation, and double occupation. Other interactions, such as additional spins possessing magnetic anisotropies, local magnetic fields, and pair-exchange interactions, can be added.
2. `TimeESR.input`: Options for the QME-related subroutines which generate the time-depenedent outputs.

## Output Files
After running the code, if one enables default output options, the following output files will be generated. 
Note: for reported eigenvectors $\bf{V}$ and Hamiltonians, components are given in pairs of $\Re(\bf{V})$ $\Im(\bf{V})$. Components are ordered by the origional occupation basis unless truncated: $\ket{\uparrow}$, $\ket{\downarrow}$, $\ket{\emptyset}$, $\ket{2}$. For more than one spin, the basis is ordered by the tensor product. For example, for two spin 1/2 objects (where the first is the driven electron site) the order of the basis is $\ket{\uparrow\uparrow}$, $\ket{\uparrow\downarrow}$, $\ket{\downarrow\uparrow}$, $\ket{\downarrow\downarrow}$, $\ket{\emptyset\uparrow}$, $\ket{\emptyset\downarrow}$, $\ket{2 \uparrow}$, $\ket{2 \downarrow}$.
1. `C.dat`: Current (pA) as a function of time (ns)
2. `EIGENVECT.dat`: Eigenvectors, listed in order of energy (GHz) reported by diagonalization and relative energy to ground state
3. `ESR.dat`: DC "measured" current (pA), which uses an averaging scheme to determine the value
4. `Hamiltonian.output`: Diagonalized Hamiltonian for future runs
5. `Hamiltonian_Eigen.dat`: Eigenvectors, listed in order of energy (meV) reported by diagonalization and relative energy to ground state
6. `PD_HAMIL.dat`: Pre-diagonalized Hamiltonian components (GHz)
7. `POP_AVE.dat`: DC "measured" populations, which uses an averaging scheme to determine the value
8. `POPULATIONS.dat`: Populations (diagonal elements) of the density matrix as a function of time (ns). The data format is: time, rho_11, rho_22., etc.
9. `RHO.dat`: Full density matrix
10. `Spin_distribution.dat`: Average spin values $\langle S_i \rangle$, for each direction $i$, for each eigenstate. A list of excitation energies (in GHZ, and meV), the value of the $S^2$ operator in units of $\hbar$, and the corresponding s quantum number for each state is also given.
11. `SpinDynamics.dat`: Several spin components as a function of time (ns). The data format is: time, ($S^x_i$, $S^y_i$, $S^z_i$) for site 1, ($S^x_i$, $S^y_i$, $S^z_i$) for site 2, etc. until site N, $S^{2}$, $\sqrt(\sqrt(1+4S^{2})-1)$, trace of the real part of the density matrix

## Troubleshooting Tips
- The total duration of the running interval must correspond to the sum of all intervals of all pulses
- When you set up the driving pulses, make sure that frequency and amplitude of the pulse correspond (this is important for ENDOR experiments or other with several frequencies)
- If you have an overall magnetic field, this goes in `H_QD.input` and is added to each local magnetic field, there is no key or entry for an overall magnetic field
- Make sure that the number of pairs with exchange interactions match the number of spins divided by two (the code will crash if not)
- Check if the key is `.true.` or `.false.` for reading the eigenstates of the previous run. Depending on what you are looping through, you may not want to read again the same eigenstates (basically you will be doing the "same" run if nothing changed in TimeESR.input)
