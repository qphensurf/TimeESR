program TimeESR

! Time-dependent master equation
! based on the theory by Galvez et al
! Phys. Rev. B Physical Review B 104 (24), 245435 (2021)

!
! gnu licence 3.0 (c) J. Reina Galvez, E. D. Switzer, C. Wolf & N. Lorente
!
!  Release version 1.0.1
!  2023 October 7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the main file and calls in:
Use declarations !all variables, constants etc are defined here
Use OpenFiles !all units of files in open statements are here
USE info ! Pretty information output for the user
Use io !(input/output) the reading and writing routines are here
Use timer
Use H_QD !(Hamiltonian Quantum Dot) contains the Hamiltonian and solution
Use QME !(Quantum Master Equation) contains rates, Runge-Kutta, and aux routines
Use SpinTools ! calcualtion of the spin vector for each spin site (or "molecule")
Use Transport ! computation of electron current
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

! 1. Start environment, print info
   CALL environment_start("TimeESR")

! 2. Read input files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          input run's values                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call clock ('STEP 1:: Reading INPUT for the QME ', 1)

     call reading_input (Nfreq, Ntime, t0, t1, Amplitude_seq, t_seq, stept,  &
       Freq_seq, Phase_seq, time, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, &
       N_int, GammaC, Cutoff, redimension, Nd, Nbias, &
       bias_R, bias_L, bias_time, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, output_file, output_fourier, output_ESR, runs, population, density_matrix, spindyn)
     call clock ('Finished reading INPUT for the QME ', 2)

! 3. Solve the Hamiltonian and redimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the spin+orbital system                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call clock ('STEP 2:: Entering Hamiltonian ', 1)

     call Hamiltonian (Ndim, lambda, Delta, Eigenvalues, H, Ss, runs, Nm, hx, hy, hz)


     if (redimension) then

    write (*,*) ' '
    write (*,*) 'Redimensioning from',Ndim,'to'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Since everything is on the Dimensions of H_QD, Ndim, and the eignestates 
! are ordered by increasing eigenvalues, we can apply an energy cutoff (do
! not confuse with Cutoff that refers to the integrals in the selfenergies)
! the energy cutoff is internal to redimensioning and it is the
! applied bias plus ten times the GammaC smearing of the Fermi function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    call redimensioning (Ndim, Delta, bias_R, bias_L, GammaC)
     Ndim = Nd
    write (*,*) 'New dimension',Ndim
    write (*,*) ' '

! Ndim value is DIFFERENT now... and substantially smaller...

    endif

     call clock ('Finished working on the Hamiltonian ', 2)

! 3. Solve the QME and 4. Electronic current
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the time-dependent                            !
!                         master equation                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call initial_population (Ndim, Ntime, Eigenvalues, rho, Temperature)

     call clock ('STEP 3:: Runge Kutta solver ', 1)
      call RungeKutta (Ndim, Ntime, stept, Delta, rho, NFreq,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Nbias, &
         bias_R, bias_L, bias_time, Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude_seq, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         N_int, GammaC, Cutoff, time, Temperature)
     call clock ('Finished Runge Kutta solver ', 2)
        
     call clock ('STEP 4:: Electronic current ', 1)
      call Current (Ndim, Ntime, rho, NFreq, Nbias, bias_time, &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude_seq, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         Temperature, Electrode, curr)
     call clock ('Finished electronic current ', 2)

      if (spindyn) then
        call SpinDynamics (Nm, Ndim, Ntime, time, Ss, H, rho, hx, hy, hz)
      endif

! 5. Print outputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    output the results                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call clock ('STEP 5:: Writing output ', 1)
     call writing_output (Ntime, time, curr, Ndim, rho, population, density_matrix, & 
          Freq_seq, output_file, output_fourier, output_ESR)
     call clock ('Finished writing output ', 2)

! 6. End environment
     call environment_end( )
      

     stop
end program TimeESR
