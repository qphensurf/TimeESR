# TimeESR v1.0.1
An STM-ESR code solving a QME in the time domain. TimeESR is given as is and subject to gnu 3.0 CopyLeft. 
Enjoy! TimeESR team, October 8, 2023.

## Intended Use
The code is intended for computing ESR in an STM junction, where the driving is produced by the modulation of the electron hopping integrals. It models one electronic orbital coupled to a complex system of spins. It is intended for (i) obtaining continuous-wave ESR spectra (our Floquet code is better suited for this), (ii) doing quantum operations driving spins with short pulses, and (iii) using multiple driving frequencies as is done in state-of-the-art STM-ESR.

## Compilation Instructions
1. Prerequisities: LAPACK, BLAS, Fortran compiler to run the program. Gnuplot, Python, xdg-open, and qutip are recommended if using the provided plot and movie scripts for any example in the `./example` folder. 
2. Review `./src/make.sys` to make sure your compiler and compilation options are set properly. 
3. Go to `./src` and execute `make`.

## Examples
The directory `./example` contains a few examples to get started. For each, execute `./run.sh` in the example directory, which cleans the `./src` folder, makes the TimeESR exectuible, runs the program, and then makes plots. Output plots include population of states, spin projections, and current as a function of time. A post-processed plot of a representative Bloch sphere, and a movie corresponding with an appropriate Bloch sphere projection as a function of time, are also given.
1. CW_ESR: A model simulation of a continous-wave ESR experiment, using a single spin 1/2 transport site and several DC measurements with several different driving frequencies which range between 16.9 to 17.1 GHz. The current vs driving frequency plot is then computed.
2. ONE_SPIN: Driving a spin 1/2 transport site from the prepared |0> state, taken to be the ground state with an applied B field in the x direction with strength 0.5 T (i.e., a |down> state in the x basis), to a (|0> + |1>)/Sqrt[2] state, i.e., a Hadamard operation on a single spin 1/2 qubit. The calculated driving frequency is 13.9953 Ghz for this operation. A rotating frame of reference is used in the post-processing step to remove the inherent larmor frequency of the system, and the system slowly decoheres after reaching the final state.
3. TWO_SPINS: Driving a spin 1/2 transport site exchange coupled to another spin 1/2 with strength J = -5.0 GHz (isotropically applied) and an applied B field in the x direction with strength 0.5 T. The first four states are, in the x basis of both spins with the transport electron first, | down down >, (| up down > + | down up >)/Sqrt[2], (| up down > - | down up >)/Sqrt[2], and | up up >. The intent is to transition all populations in the | down down > state to the | up up > state. This is achieved by driving the system from the ground state to the first excited state at 19.0 GHz, then transitioning from the first excited state to the third excited state with a drive at 9.0 GHz. Again, the system slowly decoheres after reaching the final state.
4. FIVE_SPINS: A model simulation of a DC sweeped-voltage experiment for a spin 1/2 transport site coupled to a spin system of the following geometry: transport site - spin 1 - spin 1 - spin 1/2, with another spin 1/2 site exchange coupled to the first spin 1 site. This represents a radical porphyrin (Open Acess https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202105906) except that there are two metallic centers (two Fe atoms, each with S=1). A driving frequency of 50.0 GHz is used, with bias values swept from -15.0 to 15.0 mV, and then the current vs voltage plot is computed.

## Input Files
For both input files, see the provided examples for descriptions of each option.
1. `H_QD.input`: Descriptors for the Hamiltonian, with the first spin is a driven s=1/2 electron. The first site has four states: single-occupied up spin z direction, single-occupied down spin z direction, zero occupation, and double occupation. Other interactions, such as additional spins possessing magnetic anisotropies, local magnetic fields, and pair-exchange interactions, can be added.
2. `TimeESR.input`: Options for the QME-related subroutines which generate the time-depenedent outputs.

## Output Files
After running the code, if one enables default output options, the following output files will be generated. 
Note: for reported eigenvectors and Hamiltonians, components are given in pairs of re(V) im(V). Components are ordered by the origional occupation basis unless truncated: [up], [down], [0], [2]. For more than one spin, the basis is ordered by the tensor product. For example, for two spin 1/2 objects (where the first is the driven electron site) the order of the basis is [up up], [up down], [down up], [down down], [0 up], [0 down], [2 up], [2 down].
1. `C.dat`: Current (pA) as a function of time (ns)
2. `EIGENVECT.dat`: Eigenvectors, listed in order of energy (GHz) reported by diagonalization and relative energy to ground state
3. `ESR.dat`: DC "measured" current (pA), which uses an averaging scheme to determine the value
4. `Hamiltonian.output`: Diagonalized Hamiltonian for future runs
5. `Hamiltonian_Eigen.dat`: Eigenvectors, listed in order of energy (meV) reported by diagonalization and relative energy to ground state
6. `PD_HAMIL.dat`: Pre-diagonalized Hamiltonian components (GHz)
7. `POP_AVE.dat`: DC "measured" populations, which uses an averaging scheme to determine the value
8. `POPULATIONS.dat`: Populations (diagonal elements) of the density matrix as a function of time (ns). The data format is: time, rho_11, rho_22., etc.
9. `RHO.dat`: Full density matrix
10. `Spin_distribution.dat`: Average spin values <S_i> for each direction i for each eigenstate. A list of excitation energies (in GHZ, and meV), the value of the S² operator in units of hbar, and the corresponding s quantum number for each state is also given.
11. `SpinDynamics.dat`: Several spin components as a function of time (ns). The data format is: time, (Sx, Sy, Sz) for site 1, (Sx, Sy, Sz) for site 2, etc. until site N, spin2, sqrt(sqrt(1+4*spin2)-1), trace of the real part of the density matrix

## Troubleshooting Tips
- The total duration of the running interval must correspond to the sum of all intervals of all pulses
- When you set up the driving pulses, make sure that frequency and amplitude of the pulse correspond (this is important for ENDOR experiments or other with several frequencies)
- If you have an overall magnetic field, this goes in `H_QD.input` and is added to each local magnetic field, there is no key or entry for an overall magnetic field
- Make sure that the number of pairs with exchange interactions match the number of spins divided by two (the code will crash if not)
- Check if the key is `.true.` or `.false.` for reading the eigenstates of the previous run. Depending on what you are looping through, you may not want to read again the same eigenstates (basically you will be doing the "same" run if nothing changed in TimeESR.input)
