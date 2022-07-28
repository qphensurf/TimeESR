# TimESR
STM-ESR code solving a QME in the time domain
TimeESR.f90 code is given as is and subjected to gnu 3.0 CopyLeft.

The code is intended for computing ESR in an STM junction, where the driving is produced by the modulation of the electron hopping integrals.

V0.5 only contains one electronic orbital that is coupled to a complex system of spins.

V0.5 is intended for (i) obtaining continuous-wave ESR spectra (our Floquet code is better suited for this) (ii) doing quantum operations driving spins with short pulses (iii) using multiple frequencies for extraordinary operations as is currently done in STM setups already

In the directory ./EXAMPLE you have some scripts and input files. There is a particular example for a "very" long pulse in order to obtain the continuous-wave ESR spectra (CW).

The script "script_bash_for_TimESR" removes old files, performs a loop on driving frequencies, creates two new input files as it performs the loop, runs the TimeESR.out executable that has to be in the path and then copies some of the output files with names corresponding to each frequency.

It is easy to modify the loop, add more loops, etc, adding/modifying the input values.

The present results of running "script_bash_for_TimESR" have being plotted in the document TimeESR.pdf of the DOCS directory in this distribution.

The two input files are

(1) TimeESR.input ! contains data for running the ESR code

(2) H_QD.input ! sets up the Hamiltonian corresponding to the "impurity", where the impurity is a single electronic level (with 4 states: donw, up, zero and double, corresponding to the single, empty or double electron occupations) Hunds-coupled to a "Spin Hamiltonian" (meaning any sequence of spins with anisotropies, local magnetic fields and pair-exchange interactions).

The files are commented (in the bash script joined here) so you can now what entry is what.

Troubleshooting:: Possible mistakes when running the code: 1- The total duration of the running interval must correspond to the sum of all intervals of all pulses 2- When you set up the driving pulses, make sure that frequency and amplitude of the pulse correspond (this is important for ENDOR experiments or other with several frequencies) 3- If you have an overall magnetic field, this goes in the H_QD.input and is added to each local magnetic field, there is no key or entry for an overall magnetic field 4- Make sure that the number of pairs with exchange interactions match the number of spins divided by two (the code will crash if not) 5- Check if the key is .true. or .false. for reading the impurity's eigenstates of the previous run. Depending on what you are looping on you may not want to read again the same eigenstates (basically you will be doing the "same" run if nothing changed in TimeESR.input)

Enjoy! Nicolas Lorente_ July 27, 2022
