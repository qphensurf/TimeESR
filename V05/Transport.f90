module Transport
Use declarations
Use QME
CONTAINS
      subroutine Current (Ndim, Ntime, rho, stept, Delta, NFreq,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         bias_R, bias_L, Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
         N_int, GammaC, Cutoff, time, Temperature, Electrode, curr)

     implicit none
     integer, intent (in) :: Ndim, Ntime, NFreq, N_int, Electrode
     complex (qc), intent (in) :: rho (:,:,:)
     real (q), intent (in) :: stept
     real (q), intent (in) :: Delta (:,:), time (:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature, Cutoff, GammaC
     real (q), intent (in):: bias_R, bias_L, Spin_polarization_R, Spin_polarization_L
     complex (qc), intent (in):: lambda (:,:,:)
     complex (qc), allocatable ::  GC (:,:,:,:)
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
     integer, intent (in):: t_seq (:) !sequence of pulses gives where time is
     real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
     real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
     real (q), intent (in):: Phase_seq (:) ! sequence of pulses
! Computed in ExtendedFermiIntegral called in RungeKutta
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:), fermiL_a(:,:)
     complex (qc) :: ufermiR_a(:,:), ufermiL_a(:,:)
     real (q), allocatable :: curr (:)
! internal
      integer :: l,j,u

      allocate (curr (Ntime))
      allocate (GC(Ndim,Ndim,Ndim,Ndim))





      curr = 0._q

       do i = 1, Ntime

         call  ratesC (i, Ndim, NFreq, Ntime, lambda, Delta, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, t_seq, Amplitude, Freq_seq, Phase_seq,  &
         fermiR_a, fermiL_a, ufermiR_a, ufermiL_a, &
         time, Temperature, Electrode,  GC)

         do l = 1, Ndim
         do u = 1, Ndim
         do j = 1, Ndim
      curr (i) = curr (i) +    &
            real(rho (l,u,i)*(GC(l,j,j,u)+conjg(GC(u,j,j,l))))
         enddo
         enddo
         enddo
       enddo

      
      return 
      end subroutine Current 

end module Transport
