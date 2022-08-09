module SpinTools
Use OpenFiles
Use declarations, only: q, qc, zero, time_unit
CONTAINS
! Plot how the spins evolve in time
    subroutine SpinDynamics(Nm, Ndim, Ntime, Time, Ss, H, rho)
     implicit none
     integer :: Nm, Ndim, Ntime, i, j, m, j1, j2
     real (q), intent (in) :: Time (:)
     real (q) :: sumrule
     complex (qc) :: Sx(Nm,Ntime), Sy(Nm,Ntime), Sz(Nm,Ntime)
     complex (qc) :: Spinx(Ndim,Ndim), Spiny(Ndim, Ndim), Spinz(Ndim, Ndim)
     complex (qc), intent (in) :: Ss (:,:,:,:), H (:,:), rho(:,:,:)

     Sx = zero; Sy=zero; Sz=zero

     open (unit_spinoutput,  file='SpinDynamics.dat')


     do m = 1, Nm ! molecules or spin sites
     spinX = zero; spinY=zero; spinZ=zero
     do i = 1, Ndim
     do j = 1, Ndim
! computation of the Spin
      do j1=1, Ndim
      do j2=1, Ndim
      spinX (i,j)=spinX (i,j)+conjg(H(j1,i))*Ss (m, 1, j1, j2)*H(j2, j)
      spinY (i,j)=spinY (i,j)+conjg(H(j1,i))*Ss (m, 2, j1, j2)*H(j2, j)
      spinZ (i,j)=spinZ (i,j)+conjg(H(j1,i))*Ss (m, 3, j1, j2)*H(j2, j)
      enddo
      enddo

! Trace over density matrix rho

      Sx (m,:) = Sx (m,:) + spinX (i,j)*rho(j,i,:)
      Sy (m,:) = Sy (m,:) + spinY (i,j)*rho(j,i,:)
      Sz (m,:) = Sz (m,:) + spinZ (i,j)*rho(j,i,:)


     enddo
     enddo
     enddo

      do j = 1, Ntime
      sumrule  = 0
       do i = 1, Ndim
      sumrule  = sumrule  + real (rho(i,i,j))
       enddo
      write (unit_spinoutput,*) Time(j)*time_unit , &
    &  (real(Sx(i,j)), real(Sy(i,j)), real(Sz(i,j)), &
    &   sqrt(abs(Sx(i,j))**2+abs(Sy(i,j))**2+abs(Sz(i,j))**2), i=1, Nm), sumrule ! nanoseconds
      enddo

      close (unit_spinoutput)
      return

    end subroutine SpinDynamics
end module SpinTools
