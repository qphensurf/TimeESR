program TimeESR

! Time-dependent master equation
! based on the theory by Galvez et al
! Phys. Rev. B Physical Review B 104 (24), 245435 (2021)

!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!  Version 1.0 Including Lamb shift with cutoff
!  July 2022

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the main file and calls in:
Use declarations !all variables, constants etc are defined here
Use OpenFiles !all units of files in open statements are here
Use io !(input/output) the reading and writing routines are here
Use H_QD !(Hamiltonian Quantum Dot) contains the Hamiltonian and solution
Use QME !(Quantum Master Equation) contains rates, Runge-Kutta, and aux routines
Use SpinTools ! calcualtion of the spin vector for each spin site (or "molecule")
Use Transport ! computation of electron current
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          input run's values                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call reading_input (Nfreq, Ntime, t0, t1, Amplitude_seq, t_seq, stept,  &
       Freq_seq, Phase_seq, time, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, &
       N_int, GammaC, Cutoff, &
       bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, output_file, output_fourier, output_ESR, runs, population, spindyn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the spin+orbital system                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call Hamiltonian (Ndim, lambda, Delta, Eigenvalues, H, Ss, runs, Nm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    calculate the rates for the                         !
!                    master equation for a time array                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the time-dependent                            !
!                    master equation for a time array                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call initial_population (Ndim, Ntime, Eigenvalues, rho, Temperature)

      call RungeKutta (Ndim, Ntime, stept, Delta, rho, NFreq,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         bias_R, bias_L, Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude_seq, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         N_int, GammaC, Cutoff, time, Temperature)
        
      call Current (Ndim, Ntime, rho, stept, Delta, NFreq,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         bias_R, bias_L, Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude_seq, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         N_int, GammaC, Cutoff, time, Temperature, Electrode, curr)

      if (spindyn) then
        call SpinDynamics (Nm, Ndim, Ntime, time, Ss, H, rho)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    output the results                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call writing_output (Ntime, time, curr, Ndim, rho, population, & 
          Freq_seq, output_file, output_fourier, output_ESR)

      

     stop
end program TimeESR
