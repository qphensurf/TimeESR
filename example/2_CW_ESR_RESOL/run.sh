#!/bin/bash

rundir=$(pwd)
progdir="../../src"
./clean.sh

cd $progdir
progdir=$(pwd)
make clean
make

cd $rundir
mkdir rough
mkdir fine

# Perform the simulations on a rough time grid

cd rough
# Write the Hamiltonian input file
cat >H_QD.input <<@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
1                  ! Number of spins (molecules or sites)
------------------- Spin properties: Spin 1 ---------------------
0.5                ! First spin must be 0.5 for transport
0 0 0 0            ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
.60827625 0.0 0.0  ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Exchange interactions -----------------------
0                  ! Number of exchange connected pairs
------------------- Electronic interactions ---------------------
-10.0              ! eps_QD electronic level (meV)
100.0000           ! Hubbard U (meV)
------------------- Output options ------------------------------
Hamiltonian.output ! Name of the output file with H
4                  ! Number of states to print in Spin_distribution.dat
.true.             ! .true. Write pre-diagonalized Hamiltonian to PD_HAMIL.dat
.true.             ! .true. Write all eigenvectors to EIGENVECT.dat
*****************************************************************
***************** END of INPUT FILE H_QD.input ******************
*****************************************************************
@

# Perform a DC simulation for a series of frequencies in order to produce a CW-ESR-like plot
for i in `LC_ALL=C seq 16.8 0.04 17.2`
do
cat >TimeESR.input <<@@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
100000       ! Number of points used in the time propagation
0.0   150.0  ! Initial and final time (ns)
---------------Pulse definition block----------------------------
1            ! Number of pulses (must contain zero pulses too)
1            ! Maximum number of frequencies
0.0   150.0  ! 1 - times for pulse (ns)
1.0          ! 1 - F1 amplitude
$i           ! 1 - F1 pulse frequency (GHz)
0.0          ! 1 - phase shift (radians)
---------------electrode set-up-----------------------------------
0.00500      ! gamma_R_0= 2*pi*W_R_0*W_R_0*rho (meV)
0.00010      ! gamma_L_0 2*pi*W_L_0*W_L_0*rho (meV)
0.00000      ! gamma_R_1= 2*pi*W_R_0*W_R_1*rho (meV)
0.00015      ! gamma_L_1= 2*pi*W_L_0*W_L_1*rho (meV)
100.000      ! Cutoff for integral in Lambshift (meV) 
0.00500      ! Broadening of Green's function (meV)
500000       ! Number of points in I11 and I21 (see Manual)
---------------bias, temperature, spin polarization---------------
1            ! Number of bias pulses
15.0         ! 1 - Right electrode bias (mV)
-15.0        ! 1 - Left electrode bias (mV)
150.0        ! 1 - Duration of bias pulse (ns)
0.5          ! Temperature (K)
0.0          ! Spin polarization of electrode R, between -1 and 1
1.0          ! Spin polarization of electrode L, between -1 and 1
0            ! Current measurement: 0 is left and 1 is right electrode
---------------output file names----------------------------------
.true.       ! .true. write POPULATIONS.dat
.false.      ! .true. write RHO.dat
C.dat        ! file name for the time-dependent current
S.dat        ! file name for the current Fourier transform 
ESR.dat      ! file DC current
---------------read previous run (DANGEROUS)---------------------
.true.       ! Keep to false to avoid reading previous runs
---------------other options-------------------------------------
.true.       ! Plot time evolution of spins
.false.      ! if .true. reduce states to open plus a few closed channels
4            ! Redimension number - must include empty states for non-zero current
*****************************************************************
*************** END of INPUT FILE TimeESR.input *****************
*****************************************************************
@@
echo "Frequency (GHz)::" $i
$progdir/TimeESR.x
E=`tail -1 ESR.dat` ; echo $i $E >> SpectraESR.dat
cp POPULATIONS.dat POP$i.dat
cp SpinDynamics.dat SP$i.dat
done

# Perform the simulations on a fine time grid

cd ../fine
# Write the Hamiltonian input file
cat >H_QD.input <<@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
1                  ! Number of spins (molecules or sites)
------------------- Spin properties: Spin 1 ---------------------
0.5                ! First spin must be 0.5 for transport
0 0 0 0            ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
.60827625 0.0 0.0  ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Exchange interactions -----------------------
0                  ! Number of exchange connected pairs
------------------- Electronic interactions ---------------------
-10.0              ! eps_QD electronic level (meV)
100.0000           ! Hubbard U (meV)
------------------- Output options ------------------------------
Hamiltonian.output ! Name of the output file with H
4                  ! Number of states to print in Spin_distribution.dat
.true.             ! .true. Write pre-diagonalized Hamiltonian to PD_HAMIL.dat
.true.             ! .true. Write all eigenvectors to EIGENVECT.dat
*****************************************************************
***************** END of INPUT FILE H_QD.input ******************
*****************************************************************
@

# Perform a DC simulation for a series of frequencies in order to produce a CW-ESR-like plot
for i in `LC_ALL=C seq 16.8 0.04 17.2` 
do
cat >TimeESR.input <<@@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
500000       ! Number of points used in the time propagation
0.0   150.0  ! Initial and final time (ns)
---------------Pulse definition block----------------------------
1            ! Number of pulses (must contain zero pulses too)
1            ! Maximum number of frequencies
0.0   150.0  ! 1 - times for pulse (ns)
1.0          ! 1 - F1 amplitude
$i           ! 1 - F1 pulse frequency (GHz)
0.0          ! 1 - phase shift (radians)
---------------electrode set-up-----------------------------------
0.00500      ! gamma_R_0= 2*pi*W_R_0*W_R_0*rho (meV)
0.00010      ! gamma_L_0 2*pi*W_L_0*W_L_0*rho (meV)
0.00000      ! gamma_R_1= 2*pi*W_R_0*W_R_1*rho (meV)
0.00015      ! gamma_L_1= 2*pi*W_L_0*W_L_1*rho (meV)
100.000      ! Cutoff for integral in Lambshift (meV) 
0.00500      ! Broadening of Green's function (meV)
500000       ! Number of points in I11 and I21 (see Manual)
---------------bias, temperature, spin polarization---------------
1            ! Number of bias pulses
15.0         ! 1 - Right electrode bias (mV)
-15.0        ! 1 - Left electrode bias (mV)
150.0        ! 1 - Duration of bias pulse (ns)
0.5          ! Temperature (K)
0.0          ! Spin polarization of electrode R, between -1 and 1
1.0          ! Spin polarization of electrode L, between -1 and 1
0            ! Current measurement: 0 is left and 1 is right electrode
---------------output file names----------------------------------
.true.       ! .true. write POPULATIONS.dat
.false.      ! .true. write RHO.dat
C.dat        ! file name for the time-dependent current
S.dat        ! file name for the current Fourier transform 
ESR.dat      ! file DC current
---------------read previous run (DANGEROUS)---------------------
.true.       ! Keep to false to avoid reading previous runs
---------------other options-------------------------------------
.true.       ! Plot time evolution of spins
.false.      ! if .true. reduce states to open plus a few closed channels
4            ! Redimension number - must include empty states for non-zero current
*****************************************************************
*************** END of INPUT FILE TimeESR.input *****************
*****************************************************************
@@
echo "Frequency (GHz)::" $i
$progdir/TimeESR.x
E=`tail -1 ESR.dat` ; echo $i $E >> SpectraESR.dat
cp POPULATIONS.dat POP$i.dat
cp SpinDynamics.dat SP$i.dat
done

cd ..
# Plot the results
python plot_SpectraESR_compare.py rough/SpectraESR.dat fine/SpectraESR.dat rough fine Compare.png
xdg-open Compare.png &

