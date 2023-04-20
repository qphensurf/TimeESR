module SpinTools
Use OpenFiles
Use declarations, only: q, qc, zero, time_unit
CONTAINS
! Plot how the spins evolve in time
    subroutine SpinDynamics(Nm, Ndim, Ntime, Time, Ss, H, rho, hx, hy, hz)
     implicit none
     integer :: Nm, Ndim, Ntime, i,ii, j, m, j1, j2,j3,j4, Ndim_old
     real (q), intent (in) :: Time (:), hx(:), hy(:), hz(:)
     real (q) :: sumrule, mod_h
     complex (qc) :: Sx(Nm,Ntime),Sy(Nm,Ntime),Sz(Nm,Ntime),Sh(Nm,Ntime),spin2(Ndim,Ndim)&
     &,spin2_time(Ntime)
     complex (qc) :: Spinx(Ndim,Ndim,Nm),Spiny(Ndim,Ndim,Nm),Spinz(Ndim,Ndim,Nm)
     complex (qc), intent (in) :: Ss (:,:,:,:), H (:,:), rho(:,:,:)

     Sx = zero; Sy=zero; Sz=zero; spin2=zero

     open (unit_spinoutput,  file='SpinDynamics.dat')

    Ndim_old=size(H(1,:))
     
     do i = 1, Ndim
     do ii = 1, Ndim
      do j = 1, Nm
       do j4 = 1, Nm
        do j1=1, Ndim_old
         do j2=1, Ndim_old
          do j3=1, Ndim_old
               spin2(i,ii)=spin2(i,ii)+conjg(H(j1,i))*(Ss (j, 1, j1, j3)*Ss (j4, 1, j3, j2)+  &
      &        Ss (j, 2, j1, j3)*Ss (j4, 2, j3, j2)+ Ss (j, 3, j1, j3)*Ss (j4, 3, j3, j2) &
      &        )*H(J2, ii)
          enddo
         enddo
        enddo
       enddo
      enddo
     spin2_time(:)=spin2_time(:)+spin2(i,ii)*rho(ii,i,:)
     enddo
     enddo
     
    do m=1,Nm ! loop over the molecules
     spinX = zero; spinY=zero; spinZ=zero
     do i = 1, Ndim
     do j = 1, Ndim
! computation of the Spin
      do j1=1, Ndim_old
      do j2=1, Ndim_old
      spinX (i,j,m)=spinX (i,j,m)+conjg(H(j1,i))*Ss (m, 1, j1, j2)*H(j2, j)
      spinY (i,j,m)=spinY (i,j,m)+conjg(H(j1,i))*Ss (m, 2, j1, j2)*H(j2, j)
      spinZ (i,j,m)=spinZ (i,j,m)+conjg(H(j1,i))*Ss (m, 3, j1, j2)*H(j2, j)
      
      enddo
      enddo
! print*,spinZ (i,j,m),i,j,m ,rho(j,i,1),H(j1,i)
      
      ! Trace over density matrix rho

      Sx (m,:) = Sx (m,:) + spinX (i,j,m)*rho(j,i,:)
      Sy (m,:) = Sy (m,:) + spinY (i,j,m)*rho(j,i,:)
      Sz (m,:) = Sz (m,:) + spinZ (i,j,m)*rho(j,i,:)

     enddo
     enddo
     
     

     


! Projection of spin along the local magnetic field axis
      mod_h = sqrt(hx (m)**2 + hy (m)**2 + hz (m)**2)
      if (mod_h==0) then
      print *, ' '
      print *, 'WARNING!! Magnetic field on spin number:', m
      print *, 'WARNING!! is ZERO! We cannot project along that axis!'
      print *, 'WARNING!! We set the projection to axis X'
      print *, ' '
      Sh (m,:) = Sx (m,:)
      else
      Sh (m,:) = Sx (m,:)*hx (m) + Sy (m,:)*hy (m) + Sz (m,:)*hz (m) 
      Sh (m,:) = Sh (m,:) / mod_h
      endif
     enddo

     

     
     
     
     
      do j = 1, Ntime
      sumrule  = 0
       do i = 1, Ndim
      sumrule  = sumrule  + real (rho(i,i,j))
       enddo
      write (unit_spinoutput,*) Time(j)*time_unit , &
    &  (real(Sx(i,j)), real(Sy(i,j)), real(Sz(i,j)), real(Sh(i,j)), &
    &  i=1, Nm),real(spin2_time(j)),real(sqrt(1+4*(spin2_time(j)))-1)*0.5, sumrule ! nanoseconds
      enddo

      close (unit_spinoutput)
      return

      end subroutine SpinDynamics
end module SpinTools
