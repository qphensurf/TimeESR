module krondprod
use declarations, only: qc, zero
use spin_parameters !contains definition for Sx_u etc
CONTAINS

! does the Kronecker product
! 1 x S x 1
    subroutine Kronecker_Product (N, Nm, N_in, N_block, Ss, Sx_u, Sy_u, Sz_u)
    implicit none
    integer :: i_, N, Nm, i1, i2, i3, i4
    integer :: N_in (:), N_block(:)
    complex (qc):: Ss(:,:,:,:)
    complex (qc), allocatable, intent (out) :: Sx_u (:,:,:), Sy_u (:,:,:), Sz_u (:,:,:)


! call to Pauli matrices
     call Pauli (Nm, N_in, Sx_u, Sy_u, Sz_u)


       do i2 = 1, N_in (1)
        do i3 = 1, N_in (1)

            do i4 = 1, N/N_block(1)

       Ss (1,1,i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  Sx_u (1,i2,i3);

       Ss (1,2,i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  Sy_u (1,i2,i3);

       Ss (1,3,i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  Sz_u (1,i2,i3);

            enddo
         enddo
        enddo


    do i_ = 2, Nm
 
    do i1=1, N_block(i_-1)

        do i2 = 1, N_in (i_)
        do i3 = 1, N_in (i_)
            do i4 = 1, N/N_block(i_)

       Ss (i_, 1, i4 + (i2-1)*(N/N_block(i_)) + (i1-1)*N_in(i_)*(N/N_block(i_)) ,  &
       &  i4 + (i3-1)*(N/N_block(i_)) + (i1-1)*(N_in(i_)*(N/N_block(i_))))  =  &
       &  Sx_u (i_,i2,i3); 

       Ss (i_, 2, i4 + (i2-1)*(N/N_block(i_)) + (i1-1)*N_in(i_)*(N/N_block(i_)) ,  &
       &  i4 + (i3-1)*(N/N_block(i_)) + (i1-1)*(N_in(i_)*(N/N_block(i_))))  =  &
       &  Sy_u (i_,i2,i3); 

       Ss (i_, 3, i4 + (i2-1)*(N/N_block(i_)) + (i1-1)*N_in(i_)*(N/N_block(i_)) ,  &
       &  i4 + (i3-1)*(N/N_block(i_)) + (i1-1)*(N_in(i_)*(N/N_block(i_))))  =  &
       &  Sz_u (i_,i2,i3); 

            enddo
        enddo
        enddo
     enddo
     enddo

    return
    end subroutine Kronecker_Product

! Matrix product

   Subroutine MatProdSpin (N, i, mol1, mol2, Ss, SprodS)
    implicit none
   integer :: i, j1, j2, j3, j4, N
   complex (qc) :: SprodS (:,:,:), Ss (:,:,:,:)
   complex (qc) :: produc
   integer :: mol1(:), mol2(:)

   

   do j1 =1, 3

   do j2=1, N
   do j4=1, N

   produc = zero
   do j3=1, N
   produc = produc + Ss (mol1(i),j1, j2, j3) * Ss (mol2(i), j1, j3, j4)
   enddo
   SprodS (j1,j2,j4) = produc

   enddo
   enddo

   enddo

   return
   end subroutine MatProdSpin

   Subroutine MatSquareSpin (N, Ss, SprodS)
    implicit none
   integer ::  j1, j2, j3, j4, N
   complex (qc) :: SprodS (:,:,:), Ss (:,:,:)
   complex (qc) :: produc



   do j1 =1, 3

   do j2=1, N
   do j4=1, N

   produc = zero

   do j3=1, N
   produc = produc + Ss (j1, j2, j3) * Ss (j1, j3, j4)
   enddo

   SprodS (j1,j2,j4) = produc

   enddo
   enddo

   enddo

   return
   end subroutine MatSquareSpin

   Subroutine MatSpSm (N, Sp, Sm, Sp2, Sm2)
    implicit none
   integer ::  j1, j2, j3, j4, N
   complex (qc) :: Sp (:,:), Sm (:,:)
   complex (qc) :: Sp2 (:,:), Sm2 (:,:)
   complex (qc) :: produc, produs

   do j1=1, N
   do j2=1, N

   produc = zero
   produs = zero

   do j3=1, N
   produc = produc + Sp (j1,j3)*Sp (j3,j2)
   produs = produs + Sm (j1,j3)*Sm (j3,j2)
   enddo

   Sp2 (j1,j2) = produc
   Sm2 (j1,j2) = produs

   enddo
   enddo

   return
   end subroutine MatSpSm

end module krondprod
