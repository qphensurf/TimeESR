module Transport
Use declarations
Use QME
CONTAINS
      subroutine Current (Ndim, Ntime, rho, NFreq,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         Temperature, Electrode, curr)

     implicit none
     integer, intent (in) :: Ndim, Ntime, NFreq, Electrode
     complex (qc), intent (in) :: rho (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
     complex (qc), intent (in):: lambda (:,:,:)
     complex (qc), allocatable ::  GC (:,:,:,:)
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
     integer, intent (in):: t_seq (:) !sequence of pulses gives where time is
     real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
     real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
     real (q), intent (in):: Phase_seq (:) ! sequence of pulses
     real (q) :: Pulse  !for time i
! Computed in ExtendedFermiIntegral called in RungeKutta
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:), fermiL_a(:,:)
     complex (qc) :: ufermiR_a(:,:), ufermiL_a(:,:)
     real (q), allocatable :: curr (:)
! internal
      integer :: l,j,u,n,m

      allocate (curr (Ntime))
      allocate (GC(Ndim,Ndim,Ndim,Ndim))


         call  ratesC (Ndim, NFreq, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a, ufermiR_a, ufermiL_a, &
         Temperature, Electrode,  GC)



      curr = 0._q

       do i = 1, Ntime

        n= t_seq(i) ! index of the pulse interval that contains time (i)
! Pulse sequence
          Pulse = 0._q
        do m= 1, Nfreq
          Pulse = Pulse + Amplitude (n,m)*cos(Freq_seq(n,m)*time(i)+Phase_seq(n))
        enddo



         do l = 1, Ndim
         do u = 1, Ndim
         do j = 1, Ndim
      curr (i) = curr (i) +    &
            (rho (l,u,i)*GC(l,j,j,u)+conjg(rho (l,u,i))*conjg(GC(l,j,j,u))*  &
            (1+Pulse*((1-Electrode)*gamma_L_1/gamma_L_0+Electrode*gamma_R_1/gamma_R_0)))
         enddo
         enddo
         enddo
       enddo

      deallocate (GC)
      
      return 
      end subroutine Current 

end module Transport
