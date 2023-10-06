module io
Use declarations
Use OpenFiles
CONTAINS
   subroutine reading_input ( Nfreq, Ntime, t0, t1, Amplitude_seq, t_seq, stept,  &
       Freq_seq, Phase_seq, time, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, &
       N_int, GammaC, Cutoff, redimension, Nd, Nbias, &
       bias_R, bias_L, bias_time, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, output_file, output_fourier, output_ESR, runs, population, density_matrix, spindyn)
   implicit none
   integer :: Ninterval, Nfreq, Ntime, n, i, Electrode, N_int
   real (q), allocatable :: t0 (:), t1(:), Amplitude_seq (:,:)
   real (q), allocatable :: Freq_seq (:,:), Phase_seq (:)
   character ( len = 100 ) :: output_file, output_fourier, output_ESR
   real (q) :: t_initial, t_final, stept
   real (q), allocatable :: time (:), b_time (:), j_time (:)
   integer, allocatable ::  t_seq (:), bias_time (:), exch_time (:)
   real (q) :: Spin_polarization_R, Spin_polarization_L, Temperature 
   real (q) :: gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, Cutoff, GammaC
   real (q), allocatable :: bias_R (:), bias_L (:)
   real (q) :: tolerance
   real (q) :: btime_sum 
   logical :: runs, population, spindyn, redimension, density_matrix
   integer :: Nd, Nbias

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
       allocate (bias_time (Ntime))
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
       read (unit_input, *) Cutoff !  For convergence in Lambshift (meV)
       read (unit_input, *) gammaC ! Broadening of Green's function (meV)
        read (unit_input, *) N_int ! Number of points in I11 and I21 (see Manual)
       read (unit_input, *)  ! separator
       read (unit_input, *) Nbias ! Number of bias pulses
       allocate (bias_R(Nbias), bias_L(Nbias), b_time(Nbias))
        do n = 1, Nbias
            read (unit_input, *) bias_R (n) ! mV right electrode
            read (unit_input, *) bias_L (n) ! mV left electrode
            read (unit_input, *) b_time (n) ! ns duration of bias pulse
        enddo
       read (unit_input, *) Temperature
       read (unit_input, *) Spin_polarization_R !  between -1 and 1
       read (unit_input, *) Spin_polarization_L !  between -1 and 1
       read (unit_input, *) Electrode !  0 is left and 1 is right
       read (unit_input, *)  ! separator
       read (unit_input, *) population ! .true. write POPULATIONS.dat
       read (unit_input, *) density_matrix ! .true. write RHO.dat
       read (unit_input, *) output_file !name of the file for the time-dependent current
       read (unit_input, *) output_fourier !name of the file for frequency-dependent current
       read (unit_input, *) output_ESR !name of the file for ESR signal
       read (unit_input, *)  ! separator
       read (unit_input, *) runs !if .true. read previous runs if they exist
       read (unit_input, *)  ! separator
       read (unit_input, *) spindyn !if .true. compute spin time evolution
       read (unit_input, *) redimension !if .true. reduce states to open plus a few closed channels
       read (unit_input, *) Nd  !new dimension

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


       call atomic_units (t0, t1, t_initial, t_final, GammaC, Cutoff, &
         gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1,Temperature, bias_R, &
         bias_L, b_time, Freq_seq, tolerance) ! Freq_seq transformed from frequency to radial freq
          

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
         
         if (Freq_seq (1,1) .eq. 0._q) then   ! must need at least a non-zero frequency for DC Current and DC Populations
              print *, ' '
              print *, ' ERROR: non-zero frequency needed for first frequency in order for DC current and DC Pop to be calculated.'
              print *, ' Please change this in TimeESR.input'
              print *, 'STOP.'
              print *, ' '
              stop
         endif


! creation of the time arrays for Runge Kutta and for pulse generation
         stept = (t_final - t_initial) / (Ntime-1)
    do i=1, Ntime
         time (i) = t_initial + (i-1) * stept ! in atomic units
    enddo

!initialization of t_seq for keeping ifort from complaining
         t_seq = 1
! t_seq is an array of indices that maps the time into the corresponding pulse
    
    do i=1, Ntime
     do n = 1, Ninterval
       if ((time (i)- t0(n)) > 0 .and. (t1 (n) - time (i)) > 0 ) then
                             t_seq(i) = n
       endif
     enddo
    enddo

!initialization of bias_time for keeping ifort from complaining
         bias_time = 1
! bias_time is an array of indices that maps the time into the corresponding bias pulse
    
    do i=1, Ntime
     btime_sum = 0._q
     do n = 2, Nbias
       
       btime_sum = btime_sum + b_time(n-1)
     
       if ((time (i)- btime_sum) > 0 .and. ((btime_sum + b_time (n)) - time (i)) > 0 ) then
                              bias_time(i) = n
       endif
     enddo
    enddo
   return

   end subroutine reading_input 
!
! Write output file
!
   subroutine writing_output (Ntime, time, curr, Ndim, rho, population, density_matrix, &
          Freq_seq, output_file, output_fourier, output_ESR)
   Use declfft
   Use algebra
   implicit none
   integer :: p, i, j, k, Ntime, Ndim, indexperiod
   integer :: Nfreq, Ninterval
   real (q) :: suma, sumb, weight, stept, period 
   real (q) :: time (:), curr (:), Freq_seq (:,:)
   real (q) :: t1, t2, t3, DC_current
   complex (qc) :: a (3,3)
   real (q), parameter :: pA = .662371573E10_q !a.u. to pA
   real (q), allocatable :: frequencies (:)
   real (q), allocatable :: matrix (:), array (:)
   complex (qc), allocatable :: matrix_out (:), rho (:,:,:), DC_pop (:), DC_rho (:,:)
   character ( len = 100 ) :: output_file, output_fourier, output_ESR
   logical :: population, density_matrix

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
   end if
   
   
  !    open (unit_output, file='COHERENCES.dat')
  !    do i=1,Ntime
  !    write(unit_output,'(4g14.4)') time(i)*time_unit, dble(rho(1,2,i)), &
  !     &  dble(rho(1,3,i)), dble(rho(2,3,i))
  !    enddo
  !    close (unit_output)
   
   if (density_matrix) then   
      open (unit_output, file='RHO.dat')

      do i=1,Ntime
      write(unit_output,*) time(i)*time_unit, real(rho(:,:,i),q), imag(rho(:,:,i))
      enddo

      close (unit_output)
      
   end if


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

! Determine whether we have one or two frequencies 

      allocate (DC_pop (Ndim)) !DC_pop contains the "average" population
      allocate (DC_rho (Ndim, Ntime)) !DC_rho contains the "average" diagonal density matrix elements
      do i = 1, Ndim
        DC_rho (i,:) = 0.5* rho(i,i,:)
      enddo
! to compare with qutip calculations


      Nfreq = size (Freq_seq(1,:))
      if (Nfreq == 2) then
      Ninterval = size (Freq_seq(:,1))

! combination of two sines gives a low freq difference of two freq

      period =2._q*pi_d/abs(Freq_seq (Ninterval,1)-Freq_seq (Ninterval,2))
      stept = (time(2)-time(1))
      indexperiod = int(period/stept)

! Compute DC current as the mean of the integral over a number of periods
! take 20 periods from the last point in time (Ntime)
! we have to be careful to have reached the long-time regime
! and we have to be careful to have enough points so as to minimize
! mismatch in frequencies and integrate as many oscillations as possible
      p=20

        if (p*indexperiod > int(Ntime*0.3)) then
             print *, ' '
             print *, 'WARNING: the integration time to get the DC current '
             print *, 'is bigger than one third of the total time'
             print *, 'the populations are probably not stable yet!'
             print *, 'EITHER: go into io.f90 and reduce p'
             print *, 'OR: increase Ntime'
             print *, ' '
        endif

        DC_current = 0.5* curr(Ntime-p*indexperiod)
        DC_pop (:) = 0.5* DC_rho(:,Ntime-p*indexperiod)
        
      do i = Ntime-p*indexperiod+1, Ntime-1
        DC_current = DC_current + curr(i)
        DC_pop (:) = DC_pop (:) + DC_rho (:,i)
      enddo
        DC_current = DC_current + 0.5* curr(Ntime)
        DC_pop (:) = DC_pop (:) + 0.5* DC_rho (:,Ntime)

! mean value
      DC_current = DC_current / (p*indexperiod)
      DC_pop (:) = DC_pop (:) / (p*indexperiod)  
      

      else ! I assume one only freq. Probably no sense for more than two Freq
      
      ! Starting frequency of zero doesn't work here, so make sure that we catch this in the input
      period =2._q*pi_d/Freq_seq (1,1)
      stept = (time(2)-time(1))
      indexperiod = int(period/stept)

! take five periods, this is hardwired and can
! have consequence if the driving is not sinusoidal
! or the period is to long (when doing DC calculations...)
! Comment out just the max-min for sinusoidal currents
!     p=5
!     l=0
!     allocate (array (p*indexperiod+1))
!     array = 0
!     do i = Ntime-p*indexperiod, Ntime
!        l = l+1
!        array (l) = curr (i)
!     enddo
!     DC_current= 0.5*(MAXVAL(array)+MINVAL(array))
!     deallocate (array)
!
! INSTEAD
! We apply general averaging over p periods
      p=20

        if (p*indexperiod > int(Ntime*0.3)) then
             print *, ' '
             print *, 'WARNING: the integration time to get the DC current '
             print *, 'is bigger than one third of the total time'
             print *, 'the populations are probably not stable yet!'
             print *, 'EITHER: go into io.f90 and reduce p'
             print *, 'OR: increase Ntime.'
             print *, ' '
        endif
! Integration follows:
        DC_current = 0.5* curr(Ntime-p*indexperiod)
        DC_pop (:) = 0.5* DC_rho(:,Ntime-p*indexperiod)
        
      do i = Ntime-p*indexperiod+1, Ntime-1
        DC_current = DC_current + curr(i)
        DC_pop (:) = DC_pop (:) + DC_rho (:,i)
      enddo
        DC_current = DC_current + 0.5* curr(Ntime)
        DC_pop (:) = DC_pop (:) + 0.5* DC_rho (:,Ntime)
! mean value
      DC_current = DC_current / (p*indexperiod)
      DC_pop (:) = DC_pop (:) / (p*indexperiod) 

      endif !close if on number of frequencies

      if (5*stept > period) then
         print *, 'WARNING!'
         print *, 'stept:',stept
         print *, 'has to be at least 5 times SMALLER than the period:',period
         print *, 'if you are running a time-depedent calculation this'
         print *, 'looks bad (except for issues on how to define the period) ...'
         print *, 'And this affects the calculation of the DC current.'
         print *, 'we continue but there might be a memory error if'
         print *, 'int(period/stept) is not an integer!'
         print *, ' '
      endif

  open (unit_output, file=output_ESR)

      write (unit_output, *) DC_current*pA !picoAmp units

  close (unit_output)
  open (unit_output, file='POP_AVE.dat')

      write (unit_output, *) (real(DC_pop (i)), i=1,Ndim) ! occupation units

  close (unit_output)

   return
   end subroutine writing_output 

!
! change into atomic units using the parameters of declarations.f90
!

   subroutine atomic_units (t0, t1, t_initial, t_final, GammaC, Cutoff,&
         gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1,Temperature, bias_R, &
         bias_L, b_time, Freq_seq, tolerance)
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
   real (q) :: bias_R (:), bias_L (:), b_time (:)
   real (q) :: tolerance, Cutoff, gammaC 
   logical :: runs

       t0(:)=t0(:)/time_unit; t1(:)=t1(:)/time_unit
       t_initial=t_initial/time_unit; t_final=t_final/time_unit
       tolerance = tolerance/time_unit
       gamma_R_0 = gamma_R_0/Hartree
       gamma_L_0 = gamma_L_0/Hartree
       gamma_R_1 = gamma_R_1/Hartree
       gamma_L_1 = gamma_L_1/Hartree
       gammaC = gammaC/Hartree
       Cutoff=Cutoff/Hartree
       bias_R = bias_R/Hartree
       bias_L = bias_L/Hartree
       b_time = b_time/time_unit
       Temperature = Temperature * 25.852_q / (Hartree*300._q)
       Freq_seq (:,:) = Freq_seq (:,:)*2.*pi_d*time_unit

   end subroutine atomic_units
   

  subroutine print_2darray_complex16( array, filename, factor )
  
    implicit none
  	
    complex(q), intent(in) :: array(:,:)
    character(len=*), intent(in) :: filename
    real(q), intent(in) :: factor
    integer :: i, j, size1, size2
    size1 = size(array,1)
    size2 = size(array,2)
  
    open (501, file=filename)
    do i = 1, size1
      write (501,*)  (Real(array(i,j))*factor, aimag(array(i,j))*factor, j=1, size2)
    enddo
    close (501)
  
  end subroutine

     
  subroutine print_2darray_abscomplex16( array, filename, factor )
  
    implicit none
  	
    complex(q), intent(in) :: array(:,:)
    character(len=*), intent(in) :: filename
    real(q), intent(in) :: factor
    integer :: i, j, size1, size2
    size1 = size(array,1)
    size2 = size(array,2)
  
    open (501, file=filename)
    do i = 1, size1
      write (501,*)  (sqrt((Real(array(i,j))**2.0)+(aimag(array(i,j))**2.0))*factor, j=1, size2)
    enddo
    close (501)
  
  end subroutine

end module io
