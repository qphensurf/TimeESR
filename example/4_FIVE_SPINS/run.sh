#!/bin/bash

rundir=$(pwd)
progdir="../../src"
./clean.sh

cd $progdir
make clean
make

cd $rundir

# Write the Hamiltonian input file
cat >H_QD.input <<@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
5                  ! Number of spins (molecules or sites)
------------------- Spin properties: Spin 1 ---------------------
0.5                ! First spin must be 0.5 for transport
0 0 0 0            ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
0.0 0.0 0.0        ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Spin properties: Spin 2 ---------------------
1.0                ! Second Spin
1.334 0 0 0        ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
0.0 0.0 0.0        ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Spin properties: Spin 3 ---------------------
1.0                ! Third Spin
1.334 0 0 0        ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
0.0 0.0 0.0        ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Spin properties: Spin 4 ---------------------
0.5                ! Fourth Spin
0 0 0 0            ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
0.0 0.0 0.0        ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Spin properties: Spin 5 ---------------------
0.5                ! Fifth Spin
0 0 0 0            ! Stephen Coefficients: B20 B22 B40 B44
0 0 1              ! Axis for Stephen Operators
0.0 0.0 0.0        ! Local magnetic field (T)
2.0 2.0 2.0        ! Local gyromagnetic factors
------------------- Exchange interactions -----------------------
4                  ! Number of exchange connected pairs
1  2               ! Connected pair
-1000 -1000 -1000  ! Jx, Jy, Jz (GHz)
2  3               ! Connected pair
3750 3750 3750     ! Jx, Jy, Jz (GHz)
3  4               ! Connected pair
-6800 -6800 -6800  ! Jx, Jy, Jz (GHz)
2  5               ! Connected pair
-6800 -6800 -6800  ! Jx, Jy, Jz (GHz)
------------------- Electronic interactions ---------------------
30.0               ! eps_QD electronic level (meV)
0.0000             ! Hubbard U (meV)
------------------- Output options ------------------------------
Hamiltonian.output ! Name of the output file with H
20                 ! Number of states to print in Spin_distribution.dat
.false.            ! .true. Write pre-diagonalized Hamiltonian to PD_HAMIL.dat
.false.            ! .true. Write all eigenvectors to EIGENVECT.dat
*****************************************************************
***************** END of INPUT FILE H_QD.input ******************
*****************************************************************
@

# Perform a DC simulation for a series of negative frequencies in order to produce a CW-ESR-like plot
for i in `seq 15. -0.5 0.5`
do
cat >TimeESR.input <<@@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
10000        ! Number of points used in the time propagation
0.0   0.2    ! Initial and final time (ns)
---------------Pulse definition block----------------------------
1            ! Number of pulses (must contain zero pulses too)
1            ! Maximum number of frequencies
0.0   0.2    ! 1 - times for pulse (ns)
1.0          ! 1 - F1 amplitude
50.0         ! 1 - F1 pulse frequency (GHz)
0.0          ! 1 - phase shift (radians)
---------------electrode set-up-----------------------------------
100.00       ! gamma_R_0= 2*pi*W_R_0*W_R_0*rho (meV)
0.1000       ! gamma_L_0 2*pi*W_L_0*W_L_0*rho (meV)
0.0000       ! gamma_R_1= 2*pi*W_R_0*W_R_1*rho (meV)
0.0000       ! gamma_L_1= 2*pi*W_L_0*W_L_1*rho (meV)
100.00       ! Cutoff for integral in Lambshift (meV) 
1.0000       ! Broadening of Green's function (meV)
500000       ! Number of points in I11 and I21 (see Manual)
---------------bias, temperature, spin polarization---------------
1            ! Number of bias pulses
0            ! 1 - Right electrode bias (mV)
-$i          ! 1 - Left electrode bias (mV)
0.2          ! 1 - Duration of bias pulse (ns)
5.0          ! Temperature (K)
0.0          ! Spin polarization of electrode R, between -1 and 1
0.0          ! Spin polarization of electrode L, between -1 and 1
0            ! Current measurement: 0 is left and 1 is right electrode
---------------output file names----------------------------------
.true.       ! .true. write POPULATIONS.dat
.false.      ! .true. write RHO.dat
C.dat        ! file name for the time-dependent current
S.dat        ! file name for the current Fourier transform 
DC.d at      ! file DC current
---------------read previous run (DANGEROUS)---------------------
.true.       ! Keep to false to avoid reading previous runs
---------------other options-------------------------------------
.true. ! Plot time evolution of spins
.true. ! if .true. reduce states to open plus a few closed channels
12     ! Redimension number - must include empty states for non-zero current
*****************************************************************
*************** END of INPUT FILE TimeESR.input *****************
*****************************************************************
@@
echo "Bias (meV)::" -$i
$progdir/TimeESR.x
b=`echo "$i" |bc`
E=`tail -1 DC.dat` ; echo -$b $E >> Current_DC.dat
cp POPULATIONS.dat POPn$i.dat
cp SpinDynamics.dat SPn$i.dat
done

# Do the same for the positive frequencies
for i in `seq 0.5 0.5 15`
do 
cat >TimeESR.input <<@@@
*****************************************************************
******** Please keep this FORMAT, including seperators **********
*****************************************************************
10000        ! Number of points used in the time propagation
0.0   0.2    ! Initial and final time (ns)
---------------Pulse definition block----------------------------
1            ! Number of pulses (must contain zero pulses too)
1            ! Maximum number of frequencies
0.0   0.2    ! 1 - times for pulse (ns)
1.0          ! 1 - F1 amplitude
50.0         ! 1 - F1 pulse frequency (GHz)
0.0          ! 1 - phase shift (radians)
---------------electrode set-up-----------------------------------
100.00       ! gamma_R_0= 2*pi*W_R_0*W_R_0*rho (meV)
0.1000       ! gamma_L_0 2*pi*W_L_0*W_L_0*rho (meV)
0.0000       ! gamma_R_1= 2*pi*W_R_0*W_R_1*rho (meV)
0.0000       ! gamma_L_1= 2*pi*W_L_0*W_L_1*rho (meV)
100.00       ! Cutoff for integral in Lambshift (meV) 
1.0000       ! Broadening of Green's function (meV)
500000       ! Number of points in I11 and I21 (see Manual)
---------------bias, temperature, spin polarization---------------
1            ! Number of bias pulses
0            ! 1 - Right electrode bias (mV)
$i           ! 1 - Left electrode bias (mV)
0.2          ! 1 - Duration of bias pulse (ns)
5.0          ! Temperature (K)
0.0          ! Spin polarization of electrode R, between -1 and 1
0.0          ! Spin polarization of electrode L, between -1 and 1
0            ! Current measurement: 0 is left and 1 is right electrode
---------------output file names----------------------------------
.true.       ! .true. write POPULATIONS.dat
.false.      ! .true. write RHO.dat
C.dat        ! file name for the time-dependent current
S.dat        ! file name for the current Fourier transform 
DC.d at      ! file DC current
---------------read previous run (DANGEROUS)---------------------
.true.       ! Keep to false to avoid reading previous runs
---------------other options-------------------------------------
.true. ! Plot time evolution of spins
.true. ! if .true. reduce states to open plus a few closed channels
12     ! Redimension number - must include empty states for non-zero current
*****************************************************************
*************** END of INPUT FILE TimeESR.input *****************
*****************************************************************
@@@
echo "Bias (meV)::" $i
$progdir/TimeESR.x
b=`echo "$i" |bc`
E=`tail -1 DC.dat` ; echo $b $E >> Current_DC.dat
cp POPULATIONS.dat POP$i.dat
cp SpinDynamics.dat SP$i.dat
done

# Plot the results
gnuplot CurrentESR.plot
xdg-open CurrentESR.png &


