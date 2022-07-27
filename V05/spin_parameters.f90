module spin_parameters
use declarations, only: q, qc, ur, ui, zero
use declarations, only: sqr2, sqr3, sqr5, sqr6, sqr7, sqr8
CONTAINS
   ! Pauli matrices initialized in main
   Subroutine Pauli (Nm, N_in, Sx_u, Sy_u, Sz_u)
   implicit none
   integer :: i 
   integer, parameter :: Ndim=6
   integer, intent (in) :: Nm
   integer, intent (in) :: N_in (:)
   complex (qc), allocatable :: Sx_u (:,:,:), Sy_u (:,:,:), Sz_u (:,:,:)
   
  allocate (Sx_u(Nm,Ndim,Ndim),Sy_u(Nm,Ndim,Ndim),Sz_u(Nm,Ndim,Ndim))


! we divide by 2 at the end of the subroutine
   Sx_u = zero
   Sy_u = zero
   Sz_u = zero

! first site
   Sx_u (1,1,2)=ur
   Sx_u (1,2,1)=ur
   Sy_u (1,1,2)=-ui
   Sy_u (1,2,1)=ui
   Sz_u (1,1,1)=ur
   Sz_u (1,2,2)=-ur
! remember dimension is larger because it includes electronic levels
   do i=2,Nm

! spin 1/2
   if (N_in (i) == 2) then
      Sx_u (i,1,2)=ur
      Sx_u (i,2,1)=ur

      Sy_u (i,1,2)=-ui
      Sy_u (i,2,1)=ui

      Sz_u (i,1,1)=ur
      Sz_u (i,2,2)=-ur
!spin 1
  else if (N_in (i) == 3) then
      Sx_u (i,1,2)=sqr2 
      Sx_u (i,2,1)=sqr2
      Sx_u (i,2,3)=sqr2
      Sx_u (i,3,2)=sqr2

      Sy_u (i,1,2)=-ui*sqr2
      Sy_u (i,2,1)=ui*sqr2
      Sy_u (i,2,3)=-ui*sqr2
      Sy_u (i,3,2)=ui*sqr2
   
      Sz_u (i,1,1)=ur*2
      Sz_u (i,3,3)=-ur*2

   else if (N_in (i) == 4) then

      Sx_u (i,1,2)=sqr3
      Sx_u (i,2,1)=sqr3
      Sx_u (i,2,3)=2._q
      Sx_u (i,3,2)=2._q
      Sx_u (i,3,4)=sqr3
      Sx_u (i,4,3)=sqr3

      Sy_u (i,1,2)=-ui*sqr3
      Sy_u (i,2,1)=ui*sqr3
      Sy_u (i,2,3)=-ui*2._q
      Sy_u (i,3,2)=ui*2._q
      Sy_u (i,3,4)=-ui*sqr3
      Sy_u (i,4,3)=ui*sqr3
   
      Sz_u (i,1,1)=3._q
      Sz_u (i,2,2)=1._q
      Sz_u (i,3,3)=-1._q
      Sz_u (i,4,4)=-3._q

   else if (N_in (i) == 5) then

      Sx_u (i,1,2)=2._q
      Sx_u (i,2,1)=2._q
      Sx_u (i,2,3)=sqr6
      Sx_u (i,3,2)=sqr6
      Sx_u (i,3,4)=sqr6
      Sx_u (i,4,3)=sqr6
      Sx_u (i,4,5)=2._q
      Sx_u (i,5,4)=2._q


      Sy_u (i,1,2)=-ui*2._q
      Sy_u (i,2,1)=ui*2._q
      Sy_u (i,2,3)=-ui*sqr6
      Sy_u (i,3,2)=ui*sqr6
      Sy_u (i,3,4)=-ui*sqr6
      Sy_u (i,4,3)=ui*sqr6
      Sy_u (i,4,5)=-ui*2._q
      Sy_u (i,5,4)=ui*2._q
       
      Sz_u (i,1,1)=4._q
      Sz_u (i,2,2)=2._q
      Sz_u (i,4,4)=-2._q
      Sz_u (i,5,5)=-4._q

   else if (N_in (i) == 6) then

      Sx_u (i,1,2)=sqr5
      Sx_u (i,2,1)=sqr5
      Sx_u (i,2,3)=sqr8
      Sx_u (i,3,2)=sqr8
      Sx_u (i,3,4)=3._q
      Sx_u (i,4,3)=3._q
      Sx_u (i,4,5)=sqr8
      Sx_u (i,5,4)=sqr8
      Sx_u (i,5,6)=sqr5
      Sx_u (i,6,5)=sqr5

      Sy_u (i,1,2)=-ui*sqr5
      Sy_u (i,2,1)=ui*sqr5
      Sy_u (i,2,3)=-ui*sqr8
      Sy_u (i,3,2)=ui*sqr8
      Sy_u (i,3,4)=-ui*3._q
      Sy_u (i,4,3)=ui*3._q
      Sy_u (i,4,5)=-ui*sqr8
      Sy_u (i,5,4)=ui*sqr8
      Sy_u (i,5,6)=-ui*sqr5
      Sy_u (i,6,5)=ui*sqr5

      Sz_u (i,1,1)=5._q
      Sz_u (i,2,2)=3._q
      Sz_u (i,3,3)=1._q
      Sz_u (i,4,4)=-1._q
      Sz_u (i,5,5)=-3._q
      Sz_u (i,6,6)=-5._q

   else if (N_in (i) > 6 ) then

      write (*,*) ' '
      write (*,*) ' ERROR:'
      write (*,*) ' We have only considered spins up to 5/2'
      write (*,*) ' please go to subroutine  Pauli '
      write (*,*) ' and write the corresponding Pauli-like matrices!!'
      write (*,*) ' Stop.'      
      write (*,*) ' '
      stop
    
   endif

   enddo

    Sx_u = 0.5_q*Sx_u
    Sy_u = 0.5_q*Sy_u
    Sz_u = 0.5_q*Sz_u


   end Subroutine Pauli
end module  spin_parameters
