module io
Use declarations
Use OpenFiles
CONTAINS
   subroutine reading_input ( Nfreq, Ntime, t0, t1, Amplitude_seq, t_seq, stept,  &
       Freq_seq, Phase_seq, time, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, &
       bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, output_file, output_fourier, output_ESR, runs, population, spindyn)
   implicit none
   integer :: Ninterval, Nfreq, Ntime, n, i, Electrode
   real (q), allocatable :: t0 (:), t1(:), Amplitude_seq (:,:)
   real (q), allocatable :: Freq_seq (:,:), Phase_seq (:)
   character ( len = 100 ) :: output_file, output_fourier, output_ESR
   real (q) :: t_initial, t_final, stept
   real (q), allocatable :: time (:)
   integer, allocatable ::  t_seq (:)
   real (q) :: Spin_polarization_R, Spin_polarization_L, Temperature 
   real (q) :: gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1
   real (q) :: bias_R, bias_L
   real (q) :: tolerance 
   logical :: runs, population, spindyn

   tolerance = 1.e-6 !nanoseconds

   open (unit_input, file='TimeESR.input', status='old')

! READ the pulse data

       read (unit_input, *)  ! separator
       read (unit_input, *)  ! separator
       read (unit_input, *)  ! separator
       read (unit_input, *) Ntime ! number of times for Runge-Kutta
       read (unit_input, *) t_initial, t_final ! Initial and final time in nanoseconds
         if (t_final < t_initial) then
              print *, ' '
              print *, ' ERROR: t_final must be larger than t_initial'
              print *, ' Please change them in TimeESR.input'
              print *, 'STOP.'
              print *, ' '
              stop
         endif

       read (unit_input, *)  ! separator
       read (unit_input,*) Ninterval ! Number of pulses (must contain zero pulses too)
       read (unit_input,*) Nfreq ! Maximum number of frequencies if it is a flat pulse
                                 ! you can still write Nfreq=1 in so far as you enter
                                 ! a 0 at the frequency line
                                 ! the rational for this is to be able to write
                                 ! several frequencies in the same pulse
                                 ! and be as flexible as possible with the pulse sequence

! allocates
       allocate (t0 (Ninterval), t1(Ninterval)) 
       allocate (Amplitude_seq (Ninterval, Nfreq))
       allocate (Freq_seq (Ninterval, Nfreq))
       allocate (Phase_seq (Ninterval))
       allocate (t_seq (Ntime))
       allocate (time (Ntime))

       
! definition of intervals in the input file (nanoseconds!!)
   do n = 1, Ninterval
       read (unit_input,*) t0(n), t1(n) ! in nanoseconds
         if (t1(n)<t0(n)) then
              print *, ' '
              print *, ' ERROR: t1 must be larger than t0'
              print *, ' Please change them in TimeESR.input'
              print *, 'STOP.'
              print *, ' '
              stop
         endif
     do i = 1, Nfreq
       read (unit_input,*) Amplitude_seq (n,i) ! between 0 and 1 (or more but it is relative)
                                             ! it can also be negative
                                             ! one Amplityde per frequency, even if it zero
       read (unit_input,*) Freq_seq (n, i) !Frequencies for interval n, value i. In GHz
                                           ! These are frequencies, we multiply times 2 pi
                                           ! and use them as Omega (radial freq) in the code.
     enddo
       read (unit_input,*) Phase_seq (n) ! phase shift in radians for interval n
   enddo
       read (unit_input, *)  ! separator
       read (unit_input, *) gamma_R_0 ! in meV : gamma_R_0= 2*pi*W_R_0*W_R_0*rho
       read (unit_input, *) gamma_L_0
       read (unit_input, *) gamma_R_1 ! in meV: gamma_R_1= 2*pi*W_R_0*W_R_1*rho
       read (unit_input, *) gamma_L_1
       read (unit_input, *)  ! separator
       read (unit_input, *) bias_R ! mV right electrode
       read (unit_input, *) bias_L ! mV left electrode
       read (unit_input, *) Temperature
       read (unit_input, *) Spin_polarization_R !  between -1 and 1
       read (unit_input, *) Spin_polarization_L !  between -1 and 1
       read (unit_input, *) Electrode !  0 is left and 1 is right
       read (unit_input, *)  ! separator
       read (unit_input, *) population ! .true. write POPULATIONS.dat
       read (unit_input, *) output_file !name of the file for the time-dependent current
       read (unit_input, *) output_fourier !name of the file for frequency-dependent current
       read (unit_input, *) output_ESR !name of the file for ESR signal
       read (unit_input, *)  ! separator
       read (unit_input, *) runs !if .true. read previous runs if they exist
       read (unit_input, *) spindyn !if .true. compute spin time evolution

! Then gamma_R_1/gamma_R_0 is exactly "the driving"
       print *, ' '
       print *, ' The driving is:'
       print *, '         Left electrode  =', 100*gamma_L_1/gamma_L_0, 'in %'
       print *, '         Right electrode =', 100*gamma_R_1/gamma_R_0, 'in %'
       print *, ' '
       print *, ' The applied bias is:'
       print *, '         Left electrode  =', bias_L, ' in mV'
       print *, '         Right electrode =', bias_R, ' in mV'
       print *, ' '


       call atomic_units (t0, t1, t_initial, t_final, &
         gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1,Temperature, bias_R, &
         bias_L, Freq_seq, tolerance) ! Freq_seq transformed from frequency to radial freq
          

! Check everything is correct
!   you need to have the same span in pulses as in time
         if (abs(t0(1)-t_initial) > tolerance) then   ! initial times must agree
              print *, ' ' 
              print *, ' ERROR: the initial time must be the same as the first interval time.'
              print *, ' Please change them in TimeESR.input'
              print *, 'STOP.'
              print *, ' '
              stop
         endif
        if (abs(t1(Ninterval)-t_final) > tolerance) then   ! final times must agree
              print *, ' '
              print *, ' ERROR: the final time must be the same as the last interval time.'
              print *, ' Please change them in TimeESR.input'
              print *, 'STOP.'
              print *, ' '
              stop
         endif


! creation of the time arrays for Runge Kutta and for pulse generation
         stept = (t_final - t_initial) / (Ntime-1)
    do i=1, Ntime
         time (i) = t_initial + (i-1) * stept ! in atomic units
    enddo

! t_seq is an array of indices that maps the time into the corresponding pulse
    
    do i=1, Ntime
     do n = 1, Ninterval
       if ((time (i)- t0(n)) > 0 .and. (t1 (n) - time (i)) > 0 ) then
                             t_seq(i) = n
       endif
     enddo
    enddo

   return

   end subroutine reading_input 
!
! Write output file
!
   subroutine writing_output (Ntime, time, curr, Ndim, rho, population, &
          Freq_seq, output_file, output_fourier, output_ESR)
   Use declfft
   Use algebra
   implicit none
   integer :: i, j, Ntime, Ndim
   real (q) :: suma, sumb, weight, stept, periode
   real (q) :: time (:), curr (:), Freq_seq (:,:)
   real (q) :: t1, t2, t3, DC_current
   complex (qc) :: a (3,3)
   real (q), parameter :: pA = .662371573E10_q !a.u. to pA
   real (q), allocatable :: frequencies (:)
   real (q), allocatable :: matrix (:)
   complex (qc), allocatable :: matrix_out (:), rho (:,:,:)
   character ( len = 100 ) :: output_file, output_fourier, output_ESR
   logical :: population

   if (population) then ! only if tag population is set to .true.
! Diagonal of density matrix or populations as a function of time
! name is hardwired to POPULATIONS.dat

     print *, ' '
     print *, 'Diagonal POPULATIONS written in POPULATIONS.dat '
     print *, ' '

      open (unit_output, file='POPULATIONS.dat')

      do i=1,Ntime
      write(unit_output,*) time(i)*time_unit, (dble(rho(j,j,i)), j = 1, Ndim)
      enddo

      close (unit_output)
   endif


! time-dependent current
   open (unit_output, file=output_file)

    print *, ' '
    print *, 'Output written in ',output_file,'in nanoseconds and picoAmps.'
    print *, 'The calculation is finished.'
    print *, ' '
    do i = 1, Ntime
       write (unit_output,*) time (i)*time_unit, curr (i)*pA !nanoseconds, picoAmp
    enddo

    close (unit_output)

!! Fourier transform for spectra
!  open (unit_fourier, file=output_fourier)
!
!
!  allocate (matrix (Ntime))
!  allocate (matrix_out (Ntime/2 + 1))
!  allocate (frequencies (Ntime/2 + 1))
!
!  matrix = curr*pA !in picoAmps
!
!  call dfftw_plan_dft_r2c_1d (plan,Ntime,matrix,matrix_out,FFTW_ESTIMATE)
!  call dfftw_execute_dft_r2c (plan,matrix,matrix_out)
!  call dfftw_destroy_plan(plan)
!
!  do i = 1, Ntime/2
!     frequencies (i) = (i-1)/(time(Ntime)-time(1))
!  enddo
!
!  frequencies = frequencies/time_unit !from a.u. to GHz
!  
!  do i = 1, Ntime/2 
!     write (unit_fourier,*) frequencies (i), dble(matrix_out (i) / Ntime), Aimag(matrix_out (i) / Ntime) ! GHz, picoAmps
!  enddo
!
!  close (unit_fourier)
!
! ESR-signal (DC component of the current) as a function of driving frequency component in the given pulse)
!
  open (unit_output, file=output_ESR)

! FIRST option/ average almost brute force
! this is a costly and not precise way
!   periode=2000*2._q*pi_d/Freq_seq (1,1)
!   stept = (time(2)-time(1))/periode
!   if (periode + time(1000) > time (Ntime)) then
!     print *, ''
!     print *, 'ERROR: periode+time(1000) bigger than end time'
!     print *, ' cannot compute average for CW calculation '
!     stop
!   endif

!   suma = curr(1000)*pA*stept*0.5
!   do i=1001, Ntime
!   if (time(i)-time(1000) < periode) then
!   suma = suma+curr(i)*pA*stept
!   sumb = curr(i)*pA*stept*0.5
!   endif
!   enddo
!   suma = suma - sumb


! write in output_ESR
!     write (unit_output, *) suma !DC

! SECOND option fit a cosine at the end of the evolution
! Take three "random" times
   periode=2._q*pi_d/Freq_seq (1,1)
   t1=time(Ntime)
   t2=time(Ntime)-0.5*periode
   t3=time(Ntime)-periode
   a(1,1)=cos (Freq_seq(1,1)*t1)
   a(1,2)=-sin (Freq_seq(1,1)*t1)
   a(1,3)=1._q
   a(2,1)=cos (Freq_seq(1,1)*t2)
   a(2,2)=cos (Freq_seq(1,1)*t2)
   a(2,3)=1._q
   a(3,1)=cos (Freq_seq(1,1)*t3)
   a(3,2)=-sin (Freq_seq(1,1)*t3)
   a(3,3)=1._q

   call inversion (a)

   DC_current= dble(a(3,1)+a(3,2)+a(3,3))*curr (Ntime)*pA
      write (unit_output, *) DC_current

  close (unit_output)

   return
   end subroutine writing_output 

!
! change into atomic units using the parameters of declarations.f90
!

   subroutine atomic_units (t0, t1, t_initial, t_final, &
         gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1,Temperature, bias_R, &
         bias_L, Freq_seq, tolerance)
       implicit none
       integer, parameter :: q = SELECTED_REAL_KIND(10)
       real (q), parameter :: pi_d = 3.14159265358979323846_q
       real (q), parameter :: time_unit = 2.4188843266E-8_q ! nanoseconds
       real (q), parameter :: Hartree = 27211.6_q ! meV
       real (q), parameter :: BohrMagneton=5.7883818066E-5_q !eV/T
   real (q) :: t0 (:), t1(:)
   real (q) :: Freq_seq (:,:)
   real (q) :: t_initial, t_final
   real (q) :: Temperature 
   real (q) :: gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1
   real (q) :: bias_R, bias_L
   real (q) :: tolerance 
   logical :: runs

       t0(:)=t0(:)/time_unit; t1(:)=t1(:)/time_unit
       t_initial=t_initial/time_unit; t_final=t_final/time_unit
       tolerance = tolerance/time_unit
       gamma_R_0 = gamma_R_0/Hartree
       gamma_L_0 = gamma_L_0/Hartree
       gamma_R_1 = gamma_R_1/Hartree
       gamma_L_1 = gamma_L_1/Hartree
       bias_R = bias_R/Hartree
       bias_L = bias_L/Hartree
       Temperature = Temperature * 25.852_q / (Hartree*300._q)
       Freq_seq (:,:) = Freq_seq (:,:)*2.*pi_d*time_unit

   end subroutine atomic_units

end module io
