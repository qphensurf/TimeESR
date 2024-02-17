module QME
Use declarations
Use timer
CONTAINS
     subroutine rates (Ndim, NFreq, Ntime, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a,       &
         ufermiR_a, ufermiL_a, Temperature, G)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
     integer :: Ndim, Ntime, NFreq, Nbias
! Output: the Rates called G (:,:,:,:,:,:) here
     complex (qc) :: G (:,:,:,:,:,:) ! for QME
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc) :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! Only used in this subroutine
     integer :: v, l, j, u, n
     complex (qc) :: g0p_up, g0m_up, g1p_up, g1m_up
     complex (qc) :: g0p_dn, g0m_dn, g1p_dn, g1m_dn, pulse
     complex (qc) :: Lvluj, Ljulv 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   
   
     open(666, file='rates', status='unknown')
     open(667, file='rates_R', status='unknown')
     open(668, file='rates_L', status='unknown')    
     open(669, file='rates_R_major', status='unknown')
     open(670, file='rates_L_major', status='unknown')
     open(671, file='lambda', status='unknown')
     open(672, file='gamma_pos_R', status='unknown')
     open(673, file='gamma_neg_R', status='unknown')
     open(674, file='gamma_pos_L', status='unknown')
     open(675, file='gamma_neg_L', status='unknown')
 
     do n=1,Nbias
     do j=1,Ndim
     do u=1,Ndim

      fermiR = fermiR_a(j,u,n)
      fermiL = fermiL_a(j,u,n)
      ufermiR = ufermiR_a(j,u,n)
      ufermiL = ufermiL_a(j,u,n)
      
! Right electrode
      
      g0p_up = gamma_R_0*fermiR*(1._q+Spin_polarization_R)*0.5_q
      g0m_up = gamma_R_0*ufermiR*(1._q+Spin_polarization_R)*0.5_q
      g0p_dn = gamma_R_0*fermiR*(1._q-Spin_polarization_R)*0.5_q
      g0m_dn = gamma_R_0*ufermiR*(1._q-Spin_polarization_R)*0.5_q

      write(672,*) j, u, n, g0p_up*Hartree/GHz, g0p_dn*Hartree/GHz ! Units in GHz
      write(673,*) j, u, n, g0m_up*Hartree/GHz, g0m_dn*Hartree/GHz ! Units in GHz

! Left electrode

      g1p_up = gamma_L_0*fermiL*(1._q+Spin_polarization_L)*0.5_q
      g1m_up = gamma_L_0*ufermiL*(1._q+Spin_polarization_L)*0.5_q
      g1p_dn = gamma_L_0*fermiL*(1._q-Spin_polarization_L)*0.5_q
      g1m_dn = gamma_L_0*ufermiL*(1._q-Spin_polarization_L)*0.5_q
      
      write(674,*) j, u, n, g1p_up*Hartree/GHz, g1p_dn*Hartree/GHz ! Units in GHz
      write(675,*) j, u, n, g1m_up*Hartree/GHz, g1m_dn*Hartree/GHz ! Units in GHz

        do v=1,Ndim
        do l=1,Ndim
            
            if (j .eq. 1 .and. u .eq. 1) write(671,*) v, l, lambda (v,l,1), lambda (v,l,2)
        
            Lvluj = lambda (v,l,1)*conjg(lambda(u,j,1))*g0p_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g0p_dn
            Ljulv = lambda (j,u,1)*conjg(lambda(l,v,1))*g0m_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g0m_dn

! Right electrode
            G (v,l,j,u,1,n) = 0.5_q*(Lvluj + Ljulv)

            Lvluj = lambda (v,l,1)*conjg(lambda(u,j,1))*g1p_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g1p_dn
            Ljulv = lambda (j,u,1)*conjg(lambda(l,v,1))*g1m_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g1m_dn

! Left electrode
            G (v,l,j,u,2,n) = 0.5_q*(Lvluj + Ljulv)
            
            if ( (abs(REAL(G(v,l,j,u,1,n)*Hartree/GHz,q)) >= eps4 ) &
                 .or. (abs(IMAG(G(v,l,j,u,1,n)*Hartree/GHz)) >= eps4) &
                 .or. (abs(REAL(G(v,l,j,u,2,n)*Hartree/GHz,q)) >= eps4) &
                 .or. (abs(IMAG(G(v,l,j,u,2,n)*Hartree/GHz)) >= eps4) ) then
               write(666,*) v, l, j, u, n, G(v,l,j,u,1,n)*Hartree/GHz, G(v,l,j,u,2,n)*Hartree/GHz ! Units in GHz
            end if

            if ( (abs(REAL(G(v,l,j,u,1,n)*Hartree/GHz,q)) >= eps2 ) &
                 .or. (abs(IMAG(G(v,l,j,u,1,n)*Hartree/GHz)) >= eps2) ) then 
               write(667,*) v, l, j, u, n, G(v,l,j,u,1,n)*Hartree/GHz ! Units in GHz
            end if


            if ( (abs(REAL(G(v,l,j,u,2,n)*Hartree/GHz,q)) >= eps2 ) &
                 .or. (abs(IMAG(G(v,l,j,u,2,n)*Hartree/GHz)) >= eps2) ) then 
               write(668,*) v, l, j, u, n, G(v,l,j,u,2,n)*Hartree/GHz ! Units in GHz
            end if


            if ( (abs(REAL(G(v,l,j,u,1,n)*Hartree/GHz,q)) >= eps1 ) &
                 .or. (abs(IMAG(G(v,l,j,u,1,n)*Hartree/GHz)) >= eps1) ) then
               write(669,*) v, l, j, u, n, G(v,l,j,u,1,n)*Hartree/GHz ! Units in GHz
            end if


            if ( (abs(REAL(G(v,l,j,u,2,n)*Hartree/GHz,q)) >= eps1 ) &
                 .or. (abs(IMAG(G(v,l,j,u,2,n)*Hartree/GHz)) >= eps1) ) then
               write(670,*) v, l, j, u, n, G(v,l,j,u,2,n)*Hartree/GHz ! Units in GHz
            end if


        enddo
        enddo
      enddo
     enddo
     enddo
     
     close(666)
     close(667)
     close(668)
     close(669)
     close(670)
     close(671)
     close(672)
     close(673)
     close(674)
     close(675)    

     return

     end subroutine rates 
!
! For the Current NOW
!
     subroutine ratesC (Ndim, NFreq, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a, ufermiR_a, ufermiL_a, &
         Temperature, Electrode,  GC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
     integer :: Ndim, NFreq, Nbias
     integer :: Electrode
! Output: the Rates called GC (:,:,:,:,:) here
     complex (qc) :: GC (:,:,:,:,:) ! for current
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc) :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! Only used in this subroutine
     integer :: v, l, j, u, n, i, m
     complex (qc) :: g0pa_up, g0ma_up, g1pa_up, g1ma_up
     complex (qc) :: g0pa_dn, g0ma_dn, g1pa_dn, g1ma_dn
     complex (qc) :: Lvluja, Ljulva
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if( Electrode == 1) then
      open(666, file='ratesC_R', status='unknown')
     else
      open(666, file='ratesC_L', status='unknown')
     end if
               
     do n =1,Nbias
     do j=1,Ndim
     do u=1,Ndim

! Static contribution
! Speed up if we create a Fermi array, we leave it as a function
! but it slows down the code (only two evaluation per pair u,j )


      fermiR = fermiR_a(j,u,n)
      fermiL = fermiL_a(j,u,n)
      ufermiR = ufermiR_a(j,u,n)
      ufermiL = ufermiL_a(j,u,n)
      

!Right electrode
      g0pa_up = Electrode*gamma_R_0*fermiR*(1._q+Spin_polarization_R)*0.5_q
      g0ma_up = Electrode*gamma_R_0*ufermiR*(1._q+Spin_polarization_R)*0.5_q
      g0pa_dn = Electrode*gamma_R_0*fermiR*(1._q-Spin_polarization_R)*0.5_q
      g0ma_dn = Electrode*gamma_R_0*ufermiR*(1._q-Spin_polarization_R)*0.5_q

! Left electrode
     g1pa_up = (1-Electrode)*gamma_L_0*fermiL*(1._q+Spin_polarization_L)*0.5_q
     g1ma_up = (1-Electrode)*gamma_L_0*ufermiL*(1._q+Spin_polarization_L)*0.5_q
     g1pa_dn = (1-Electrode)*gamma_L_0*fermiL*(1._q-Spin_polarization_L)*0.5_q
     g1ma_dn = (1-Electrode)*gamma_L_0*ufermiL*(1._q-Spin_polarization_L)*0.5_q

        do v=1,Ndim
        do l=1,Ndim

!Right electrode
            Lvluja = lambda (v,l,1)*conjg(lambda(u,j,1))*g0pa_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g0pa_dn
            Ljulva = lambda (j,u,1)*conjg(lambda(l,v,1))*g0ma_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g0ma_dn

            GC (v,l,j,u,n) = 0.5_q*(Lvluja - Ljulva)

! Left electrode
            Lvluja = lambda (v,l,1)*conjg(lambda(u,j,1))*g1pa_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g1pa_dn
            Ljulva = lambda (j,u,1)*conjg(lambda(l,v,1))*g1ma_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g1ma_dn
              
              
            GC (v,l,j,u,n) = GC (v,l,j,u,n) + 0.5_q*(Lvluja - Ljulva) ! This has the right admixture of electrodes

            if( (abs(REAL(GC(v,l,j,u,n)*Hartree/GHz,q)) >= eps2) &
                 .or. (abs(IMAG(GC(v,l,j,u,n)*Hartree/GHz)) >= eps2) ) then
               write(666,*) v, l, j, u, n, GC(v,l,j,u,n)*Hartree/GHz ! Units in GHz
            end if
        enddo
        enddo
      enddo
     enddo
     enddo
     
     close(666)
         
     return

     end subroutine ratesC 
!
! initial population for the propagation
!
     subroutine initial_population (Ndim, Ntime, Eigenvalues, rho, Temperature)
     implicit none
     integer, intent (in) :: Ndim, Ntime
     real (q), intent (in) :: Temperature
     real (q), intent (in) :: Eigenvalues (:)
     complex (qc), allocatable, intent (out) :: rho (:,:,:)
     real (q) :: Z, expval, maxarg, ln_val, sum_exp
     real (q), allocatable :: arg(:)
     integer :: i
     logical, allocatable :: truncate_flag(:)
     logical :: report_truncate

     ! Overflow and underflow errors are common in calculation of partition functions
     ! so we utilize a log sum exp (LSE) trick which is commonly used in data science

     allocate (rho (Ndim, Ndim, Ntime))
     allocate(arg(size(Eigenvalues)))
     allocate(truncate_flag(size(Eigenvalues)))
     
     Z = 0._q; rho = zero;
     truncate_flag = .FALSE.
     report_truncate = .FALSE.
     arg = 0._q
   
     ! LSE: Get exponential arguments
     do i=1, Ndim 
      arg(i) = -(Eigenvalues(i)-Eigenvalues(1))/Temperature 

        ! Checking for underflow errors
        if ( arg(i) < log(tiny(1.0_q)) ) then
          truncate_flag(i) = .TRUE.
          report_truncate = .FALSE. ! Changed from true to get rid of extra messages
        end if
     end do
     
     if(report_truncate) then
      
       write(*,*) "Large negative exponential argument(s) found for initial population calculation"
       write(*,*) "The exponent of this argument will be truncated to zero"
       write(*,*) "(state number, eigenvalue, temperature, exp argument, &
          & and largest negative exponent argument allowed):"
       do i=1, Ndim
         if (truncate_flag(i)) then
            write(*,*) i, Eigenvalues(i)-Eigenvalues(1), Temperature, arg(i), log(tiny(1.0_q))
         end if
       end do
       
     end if

     ! LSE: Find maximum value
     maxarg = maxval(arg)
     
     ! LSE: Sum up all exponentials
     sum_exp = 0._q 
     do i =1, Ndim
      if ( arg(i) .ge. log(tiny(1.0_q))) then
         sum_exp = sum_exp + exp(arg(i) - maxarg)
      end if
     end do
     
     ! LSE: Log that sum, and add back the maximum amount that was removed
     ln_val = maxarg + log(sum_exp)
     
     ! Get Z from LSE
     Z = exp(ln_val)

     ! Now do the populations
     do i = 1, Ndim
      if ( arg(i) .ge. log(tiny(1.0_q)) ) then
         rho (i,i,1) = rho (i,i,1) + exp ( arg(i) )/Z
      end if
     enddo

!      write(*,*) Temperature, Z,Eigenvalues
!       
!      do i=1,ndim
!         write(*,*)
!         write(*,*) (rho(i,j,1),j=1,ndim)
!         write(*,*)  
!      enddo

     deallocate(arg)
     
     return

     end subroutine initial_population 
!
! Master equation solver
! RK4 method
!
     subroutine RungeKutta (Ndim, Ntime, stept, Delta, rho, NFreq,  &
      lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Nbias, &
      bias_R, bias_L, bias_time, Spin_polarization_R, Spin_polarization_L,  &
      t_seq, Amplitude, Freq_seq, Phase_seq,  &
      fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
      N_int, GammaC, Cutoff, time, Temperature)

     implicit none
     real (q), intent (in) :: stept
     real (q), intent (in) :: Delta (:,:), time (:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature, Cutoff, GammaC
     real (q), intent (in):: bias_R (:) , bias_L (:)
     real (q) :: Spin_polarization_R, Spin_polarization_L
     complex (qc), intent (in):: lambda (:,:,:)
     integer :: Ndim, Ntime, NFreq, N_int, Nbias
     real (q) :: Pulse  !for time i and i+1
     complex (qc), allocatable:: Gstatic (:,:,:,:,:,:)
     complex (qc), allocatable:: G (:,:,:,:,:)
     complex (qc) :: rho (:,:,:)
     complex (qc) :: Temp1, Temp2 !buffers for changing names...
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
     integer, intent (in):: t_seq (:), bias_time (:) !sequence of pulses gives where time is
     real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
     real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
     real (q), intent (in):: Phase_seq (:) ! sequence of pulses
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc), allocatable:: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc), allocatable:: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! internal
     integer :: n, l, j, i, u, v, m, nb
     real (q) :: half
     complex (qc), allocatable :: k1(:), k2(:), k3(:), k4(:), D(:), P(:)
     logical :: print_flag

      
     allocate (k1(Ndim*Ndim), k2(Ndim*Ndim), k3(Ndim*Ndim), k4(Ndim*Ndim))
     allocate (D(Ndim*Ndim), P (Ndim*Ndim))
     allocate (Gstatic(Ndim,Ndim,Ndim,Ndim,2,Nbias))
     allocate (G(Ndim,Ndim,Ndim,Ndim,2))
     allocate (fermiR_a(Ndim,Ndim,Nbias))
     allocate (fermiL_a(Ndim,Ndim,Nbias))
     allocate (ufermiR_a(Ndim,Ndim,Nbias))
     allocate (ufermiL_a(Ndim,Ndim,Nbias))

! Bias intergal
     print_flag = .FALSE.
     do n = 1, Nbias
! I11 and I21 from the Manual are honored here:
     do j=1,Ndim
     do u=1,Ndim
      call ExtendedFermiIntegral ( Delta (j,u), bias_R (n), Temperature, Cutoff, GammaC, N_int, fermiR, n, j, u, 'R', print_flag )
      call ExtendedFermiIntegral ( Delta (j,u), bias_L (n), Temperature, Cutoff, GammaC, N_int,  fermiL, n, j, u, 'L', print_flag )
      call ExtendeduFermiIntegral ( Delta (j,u), bias_R (n), Temperature, Cutoff, GammaC, N_int, ufermiR)
      call ExtendeduFermiIntegral ( Delta (j,u), bias_L (n), Temperature, Cutoff, GammaC, N_int,  ufermiL)
!The idea is to create an array and pass it to the rates routine
      fermiR_a(j,u,n) = fermiR  / pi_d !important pi factor, see Manual
      fermiL_a(j,u,n) = fermiL  / pi_d
      ufermiR_a(j,u,n) = ufermiR  / pi_d
      ufermiL_a(j,u,n) = ufermiL  / pi_d

     enddo
     enddo
     enddo

       call rates (Ndim, NFreq, Ntime, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a,   &
         ufermiR_a, ufermiL_a, Temperature, Gstatic)

! loop on time
      call clock ('Entering time loop in RK after computing static rates', 2)
       do i = 1, Ntime-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate rates for this time i, and  next one i+1
        n= t_seq(i) ! index of the pulse interval that contains time (i)
        nb= bias_time (i) ! same idea but for the bias pulse
! Pulse sequence
          Pulse = 0._q
        do m= 1, Nfreq
          Pulse = Pulse + Amplitude (n,m)*cos(Freq_seq(n,m)*time(i)+Phase_seq(n))
        enddo

! Add pulse on G, I first create a temporal value to add left and right electrode
! and then I crash the left and right G by the first and second times (ugly, I KNOW)
          
        G (:,:,:,:,1) = Gstatic (:,:,:,:,1,nb)*(1._q+Pulse*gamma_R_1/gamma_R_0) + &
               Gstatic (:,:,:,:,2,nb)*(1._q+Pulse*gamma_L_1/gamma_L_0) ! The driving is the ratio gamma_R_1/gamma_R_0

! Repeat for time i+1
        n= t_seq(i+1) ! index of the pulse interval that contains time (i)
! Pulse sequence
          Pulse  = 0._q
        do m= 1, Nfreq
          Pulse = Pulse  + Amplitude (n,m)*cos(Freq_seq(n,m)*time(i+1)+Phase_seq(n))
        enddo

        G (:,:,:,:,2) = Gstatic (:,:,:,:,1,nb)*(1._q+Pulse*gamma_R_1/gamma_R_0) + &
               Gstatic (:,:,:,:,2,nb)*(1._q+Pulse*gamma_L_1/gamma_L_0) ! The driving is the ratio gamma_R_1/gamma_R_0
! G(:,:,:,:,1) is for time i and G(:,:,:,:,2) is for time i+1
! end generate rates for this time i, and i+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initial step
! Respects detailed balance because it is thermal
          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo


       half = 0
       call matrix_rate (D, G, Delta, P, half)
       k1 = stept*P

       D = D + 0.5_q*k1
       half = 0.5_q
       call matrix_rate (D, G, Delta, P, half)
       k2 = stept*P
       
        n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo

       D = D + 0.5_q*k2
       half = 0.5_q
       call matrix_rate (D, G, Delta, P, half)
       k3 = stept*P
       
         n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo
       D = D + k3
       half = 1._q
       call matrix_rate (D, G, Delta, P, half)
       k4 = stept*P

       D = k1+2._q*k2+2._q*k3+k4

       
          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       rho (l,j,i+1) = rho (l,j,i) + D (n) /6._q
       enddo
       enddo

       enddo

     deallocate (k1, k2, k3, k4)
     deallocate (D,P)
     deallocate (Gstatic)
     deallocate (G)

       return
     end subroutine RungeKutta

!
! The rate equation is written as D'=M.D
!
! it is:
!  
! 
     subroutine matrix_rate (D, G, Delta, P, half)
     implicit none
     complex (qc), intent (in):: G (:,:,:,:,:)
     complex (qc), intent (in) :: D (:)
     complex (qc), intent (out) :: P(:)
     real (q), intent(in) :: Delta (:,:)
     real (q) :: half
     complex (qc) :: rate
     integer ::  m
     integer :: n, l, j, r, s

          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1

              P (n) = ui*Delta (l,j)*D(n)

! contribution on rho (v,u)
             m = 0
          do v = 1, Ndim
          do u = 1, Ndim
             m = m+1

              rate = half*(G (v,l,j,u,2)-G (v,l,j,u,1))+G (v,l,j,u,1)
              rate = rate + conjg(half*(G (u,j,l,v,2)-G (u,j,l,v,1))+G (u,j,l,v,1))

              P (n) = P(n) + rate * D(m)

! contribution on rho (l,u)
         
             r = (l-1)*Ndim + u

             rate = half*(G (j,v,v,u,2) - G (j,v,v,u,1)) + G (j,v,v,u,1)

             P(n) = P (n) - rate * D (r)
              
! contribution on rho (u, j)

            s = (u-1)*Ndim + j
              
            rate = half*(G (l,v,v,u,2) - G (l,v,v,u,1)) + G (l,v,v,u,1)

            P(n) = P (n) - conjg(rate) * D (s)

          enddo
          enddo


       enddo
       enddo

     return

     end subroutine matrix_rate
!
! Fermi occupation function
!
!    function Fermi (e, T)
!      implicit none
!      real (q) :: Fermi, e, T
!         if ( e > 0._q ) then
!            Fermi = exp(-e/T)/(exp(-e/T)+1._q)
!         else
!           Fermi = 1._q/(exp(e/T)+1._q)
!         endif
!      return
!    end function Fermi

!
! Fermi occupation function with limitors to avoid underflow/overflows
!
    function Fermi (e, T)
    
      implicit none
      real (q) :: Fermi
      real (q), intent(in) :: e, T
      real (q) :: beta
      
      Fermi = 0._q
      beta = 0._q
      
      beta = 1._q / T
      Fermi = beta * e
      
      ! Take care of underflow and overflow errors here, which can occur within some distance from tiny and huge
      if(Fermi-4._q .lt. log(tiny(1._q))) then
#ifdef __DEBUG
         write(*,*) "Large negative beta*energy value found"
         write(*,*) "(beta, energy, beta*energy, largest negative value allowed):"
         write(*,*) beta, e, Fermi, log(tiny(1._q))
         write(*,*) "Fermi function will default to the value of 1"
#endif
         Fermi = 1._q
      elseif (Fermi+4._q .gt. log(huge(1._q))) then
#ifdef __DEBUG
         write(*,*) "Large positive beta*energy value found"
         write(*,*) "(beta, energy, beta*energy, largest positive value allowed):"
         write(*,*) beta, e, Fermi, log(huge(1._q))
         write(*,*) "Fermi function will default to the value of 0"
#endif
         Fermi = 0._q
      else
         Fermi = 1._q/(1._q + exp(Fermi))
      end if
      
    end function


! Calculation of energy integration of rates involving the Fermi function
      subroutine ExtendedFermiIntegral ( D, V, T, Cutoff, GammaC, N,  fermiA, i_n, i_j, i_u, bias_dir, print_flag )
      implicit none
      real (q) :: D, V, T, Cutoff, GammaC
      real (q) :: e, step_e
      integer :: i, N
      complex (qc):: fermiA
      logical :: truncate_flag
      integer :: i_n, i_j, i_u
      character(len=*) :: bias_dir
      logical :: print_flag
!fermiA is Integral I11 of the Manual
! Trapeze-integration (the best among the better)

      step_e = 2._q*Cutoff/(N-1._q)
      e = -Cutoff
      fermiA = 0.5_q*Fermi(e-V,T) / (e - D + ui*GammaC)
      truncate_flag = .FALSE.
      
      do i = 2, N-1
         e= -Cutoff + (i-1._q)*step_e       
         
         ! Found an issue with dividing a small e-D quantity into a very small Fermi
         ! By including GammaC, Fortran performs a multiplication of the numerator by
         ! (e-D,-GammaC) and the denominator becomes |(e-D,GammaC)|^2. The first step
         ! may underflow the real register, and so we handle it the best we can.
         if(Fermi(e-V,T) .ne. 0._q) then
            if((log(Fermi(e-V,T)) + log(abs(e-D)))-6._q .ge. log(tiny(1._q))) then
               fermiA = fermiA + Fermi(e-V,T) / (e - D + ui*GammaC)
            else
               truncate_flag = .TRUE.
            end if
         end if
        
      enddo
      
      if(truncate_flag) then
         if (print_flag) then
            write(*,*) "Truncated ExtendedFermiIntegral, &
                        for bias direction, number, and state combination: "
            print_flag = .FALSE.
         end if
         !write(*,*) trim(bias_dir), i_n, i_j, i_u ! Changed to get rid of extra messages
      end if
      
      e = Cutoff
      fermiA = fermiA + 0.5_q*Fermi(e-V,T) / (e - D + ui*GammaC)
      fermiA = step_e*ui*fermiA

      return
      end subroutine ExtendedFermiIntegral
      
! Calculation of energy integration of rates involving 1-Fermi function
      subroutine ExtendeduFermiIntegral ( D, V,  T, Cutoff, GammaC, N, ufermiA)
      implicit none
      real (q) :: D, V, T, Cutoff, GammaC
      real (q) :: e, step_e
      integer :: i, N
      complex (qc):: ufermiA
!ufermiA is Integral I21 of the Manual
! Trapeze-integration (the best among the better)

      step_e = 2._q*Cutoff/(N-1._q)
      e= -Cutoff 
      ufermiA=0.5_q*(1._q-Fermi (e-V, T)) / (e+D-ui*GammaC)

      do i = 2, N-1
      e= -Cutoff + (i-1._q)*step_e
      ufermiA=ufermiA+(1._q-Fermi (e-V, T)) / (e+D-ui*GammaC)
      enddo
      e = Cutoff
      ufermiA=ufermiA+0.5_q*(1._q-Fermi (e-V, T)) / (e+D-ui*GammaC)

      ufermiA = -step_e*ui*ufermiA

      return
      end subroutine ExtendeduFermiIntegral

     
end module QME
