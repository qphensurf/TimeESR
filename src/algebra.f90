module algebra
! contains interfaces to call LAPACK
! for inversion and diagonalization
!
! does not contain matrix operations making use of BLAS 
!   I had the impression that BLAS3 was an overkill for
!   what we do here, work in progress... N.L. 6 May 2013
!
! contains linear interpolation of arrays

Use declarations,only:q,qc

interface initialize_zero
    module procedure initialize_zero_1, initialize_zero_2
end interface
 
interface num_mat_mult
    module procedure num_mat_mult1, num_mat_mult2
end interface 

CONTAINS
        SUBROUTINE INVERSION (A)
        IMPLICIT NONE
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        integer (I4B), save :: M
        integer (I4B), allocatable, save :: IPIV (:)
        integer (I4B), save :: INFO, LWORK
        complex (qc), DIMENSION (:,:) :: A
        complex (qc), allocatable, save :: WORK(:)

         M=size(A,1)

        LWORK=M

        ALLOCATE (IPIV(M))
        ALLOCATE (WORK(LWORK))


        call ZGETRF (M, M, A, M, IPIV, INFO)

        if (INFO /= 0) then

                print *, 'factorization failed, INFO=', INFO
                print *, 'in INVERSION subroutine; STOPPING'
                stop
        endif

        call ZGETRI (M, A, M, IPIV, WORK, LWORK, INFO)

        if (INFO /= 0) then

                print *, 'inversion failed, INFO=', INFO
                print *, 'STOPPING'

                stop

        endif

        DEALLOCATE (IPIV)
        DEALLOCATE (WORK)

        RETURN

        END SUBROUTINE INVERSION
!
! DIAGONALIZATION 
!
        SUBROUTINE DIAGONALIZE(N_points,A,W)
        IMPLICIT NONE
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
 ! To be consistent with the rest of the code:
        complex(qc), DIMENSION(:,:),INTENT(INOUT) :: A
        real(q),    DIMENSION(:),  INTENT(INOUT) :: W
        character ( len = 1 ) :: UPLO, JOB
        integer (I4B) :: N_points
        integer (I4B) :: INFO, LWORK
        REAL (q),    DIMENSION(:), allocatable   :: RWORK
        complex (qc), DIMENSION(:), allocatable   :: WORK

        UPLO='U' !upper triangle of A is stored
        JOB='V' ! vecteurs et valeurs propres

        LWORK=2*N_points
        INFO = 0

         print *, 'N_points:', N_points
        allocate (WORK(LWORK))
        allocate (RWORK(max(1,3*N_points-2)))
       
       call zheev (JOB, UPLO, N_points, A, N_points, &
&           W, WORK, LWORK,RWORK,INFO)

         print *, 'N_points:', N_points

       if ( INFO /= 0 ) then
         print *, 'the diagonalization had problems'
         print *, 'INFO =', INFO
         print *, 'check your input data'
         print *, 'but we continue the calculation.'
       endif
       deallocate(RWORK,WORK)

        END SUBROUTINE DIAGONALIZE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!    Matrix functions intended to improve portability   !
!    and perhaps performance                            !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine matrix_copy (A,B)
        IMPLICIT NONE
! copies A in B
    integer :: i, j, N
    complex (qc) ::  A (:,:)
    complex (qc) :: B (:,:)

    N = size(A,1)
    
      do j = 1, N
    do i = 1, N
      B (i,j) = A (i,j) 
      enddo
    enddo
    end subroutine matrix_copy
!
! add matrices
!
    subroutine matrix_add (A,B,C)
! adds A+B
        IMPLICIT NONE
    integer :: i, j, N
    complex (qc), intent (in) ::  A (:,:), B (:,:)
    complex (qc), intent (out) :: C (:,:)

    N = size(A,1)
    
    do i = 1, N
      do j = 1, N
       C (i,j) = A (i,j) + B (i,j)
      enddo
    enddo

    end subroutine matrix_add
!
! initialize an array to zero
!
    function initialize_zero_1 (N)
        IMPLICIT NONE
    integer :: i, N
    complex (qc) :: initialize_zero_1 (N)
 
    do i = 1, N
           initialize_zero_1 (i) = (0._q,0._q)
    enddo
    end function initialize_zero_1
!
! initialize an array to zero
!
    function initialize_zero_2 (N,M)
        IMPLICIT NONE
    integer :: i, j, N, M
    complex (qc) :: initialize_zero_2 (N,M)
 
      do j = 1, M
    do i = 1, N
           initialize_zero_2 (i,j) = (0._q,0._q)
      enddo
    enddo
    end function initialize_zero_2
!
! multiply scalar by matrix
!
! real scalar
    function num_mat_mult1 (z, A)
        IMPLICIT NONE
    integer :: i, j, N
    real (q) :: z
    complex (qc) :: A (:,:)
!   complex (qc), allocatable :: num_mat_mult1 (:,:)
    complex (qc) :: num_mat_mult1 (size (A, 1), size (A, 1))

    N = size (A, 1)
!   allocate ( num_mat_mult1 (N,N) )
      do j = 1, N
    do i = 1, N
      num_mat_mult1 (i,j) = A (i,j) * z
      enddo
    enddo
    end function 
! complex scalar
    function num_mat_mult2 (z, A)
        IMPLICIT NONE
    integer :: i, j, N
    complex (qc) :: z, A (:,:)
!   complex (qc), allocatable :: num_mat_mult2 (:,:)
    complex (qc) :: num_mat_mult2 (size (A, 1), size (A, 1))

    N = size (A, 1)
!   allocate ( num_mat_mult2 (N,N) )
      do j = 1, N
    do i = 1, N
      num_mat_mult2 (i,j) = A (i,j) * z
      enddo
    enddo
    end function 
    
end module algebra
