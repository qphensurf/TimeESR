module spin_parameters
   
   use declarations, only: q, qc, ur, ui, zero, Spin
   use error_handler, only: errore

contains

   subroutine Pauli( Nm, N_in, Sx_u, Sy_u, Sz_u )
   ! Build pauli matrices. This routine should be initialized in main.
   
      implicit none
       
      ! Arguments
      integer, intent(in) :: Nm
      integer, intent(in) :: N_in(:) ! 2s+1
      complex(qc), allocatable, intent(inout) :: Sx_u(:,:,:), Sy_u(:,:,:), Sz_u(:,:,:)
      
      ! Local variables
      integer :: i, ierr, mdim, ndim, imol, im, jm 
      real(q) :: sval, mval
      complex(qc), allocatable :: Sp_u(:,:,:), Sm_u(:,:,:)
      
      ! Maximum dimension of spin matrices (must be at least as large as spin 1/2 space + 2 electronic levels)
      ndim = maxval(N_in)
      
      ! Safe allocate matrices
      if(allocated(Sp_u)) deallocate(Sp_u)
      if(allocated(Sm_u)) deallocate(Sm_u)
      if(allocated(Sx_u)) deallocate(Sx_u)
      if(allocated(Sy_u)) deallocate(Sy_u)
      if(allocated(Sz_u)) deallocate(Sz_u)
      
      allocate(Sp_u(Nm,ndim,ndim),stat=ierr)
      if(ierr .ne. 0) call errore( 'Pauli', 'memory allocation failed on Sp_u', ierr)
      allocate(Sm_u(Nm,ndim,ndim),stat=ierr)
      if(ierr .ne. 0) call errore( 'Pauli', 'memory allocation failed on Sm_u', ierr)
      allocate(Sx_u(Nm,ndim,ndim),stat=ierr)
      if(ierr .ne. 0) call errore( 'Pauli', 'memory allocation failed on Sx_u', ierr)
      allocate(Sy_u(Nm,ndim,ndim),stat=ierr)
      if(ierr .ne. 0) call errore( 'Pauli', 'memory allocation failed on Sy_u', ierr)
      allocate(Sz_u(Nm,ndim,ndim),stat=ierr)
      if(ierr .ne. 0) call errore( 'Pauli', 'memory allocation failed on Sz_u', ierr)
      
      Sp_u = zero
      Sm_u = zero
      Sx_u = zero
      Sy_u = zero
      Sz_u = zero
      
      ! Fill matrices for each spin object
      do imol = 1, Nm
         
         ! Spin value of object
         sval = Spin(imol)
         mdim = int((2._q * sval) + 1._q)
         
         do im = 1, mdim-1 ! row
            
            mval = zero
            mval = sval - REAL(im,q)
            jm = im + 1 ! column
            
            ! Fill out the plus matrix according to
            ! Splus[1:] = [ 0 \sqrt(s(s+1) - s(s+1)) 0 ...], Splus[2:] = [ 0 0 \sqrt(s(s+1) - (s-1)((s-1)+1)) 0 ...], Splus[3:] = [ 0 0 0 \sqrt(s(s+1) - m(m+1)) ...]
            Sp_u(imol,im,jm) = SQRT(( sval * ( sval + 1._q )) - ( mval * ( mval + 1._q ) ))
            
            ! Sminus is transpose of Splus
            Sm_u(imol,jm,im) = Sp_u(imol,im,jm)
            
            ! Sz is easy
            Sz_u(imol,im,im) = mval + ur
         
         end do
         
         Sz_u(imol,mdim,mdim) = -ur * sval
         
      end do
      
      ! Factors of 1/2 where appropriate
      Sx_u = 0.5_q * ( Sp_u + Sm_u )
      Sy_u = -ui * 0.5_q * ( Sp_u - Sm_u )
      
      ! Deallocate 
      deallocate( Sp_u )
      deallocate( Sm_u )

#ifdef __DEBUG
      do imol = 1, Nm
         sval = Spin(imol)
         mdim = int((2._q * sval) + 1._q)
         write(*,'(/A,I1,A,F3.1)') "Spin of mol ", imol,": ",sval
         write(*,'(A,I1)') "Sx_u for mol ", imol
         do im = 1, mdim
            do jm = 1, mdim
               write(*,'(2I5,2F20.10)') im, jm, REAL(Sx_u(imol,im,jm),q), IMAG(Sx_u(imol,im,jm))
            end do
         end do
         
         write(*,'(A,I1)') "Sy_u for mol ", imol
         do im = 1, mdim
            do jm = 1, mdim
               write(*,'(2I5,2F20.10)') im, jm, REAL(Sy_u(imol,im,jm),q), IMAG(Sy_u(imol,im,jm))
            end do
         end do
         
         write(*,'(A,I1)') "Sz_u for mol ", imol
         do im = 1, mdim
            do jm = 1, mdim
               write(*,'(2I5,2F20.10)') im, jm, REAL(Sz_u(imol,im,jm),q), IMAG(Sz_u(imol,im,jm))
            end do
         end do
      end do
#endif

   end subroutine Pauli
   
end module  spin_parameters
