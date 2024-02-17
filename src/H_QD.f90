module H_QD
Use declarations
Use OpenFiles
Use algebra
Use KrondProd
implicit none

CONTAINS

   subroutine Hamiltonian (N, lambda, Delta, Eigenvalues, H, Ss, runs, Nm, hx, hy, hz)
   
   USE io, ONLY: print_2darray_complex16
   implicit none   
   integer, intent(out) :: N, Nm
   real (q), allocatable, intent (out) :: hx (:), hy(:), hz(:)
   real (q), allocatable ::  Eigenvalues (:), Delta (:,:)
   complex (qc), intent (out), allocatable :: lambda (:,:,:), H(:,:), Ss (:,:,:,:)
   logical :: runs, prediag_hamiltonian, eigenvectors


!
! Create many body configurations for spin excitation
! or ESR dynamics
!
! Brute force diagonalization of a model Hamiltonian
! that contains
!         spins, anisotropic exchange, local magnetic fields, Stephen Operators,
!         an exchange interaction between electron site and first spin
!         finite Hubbard U
!
!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!
! in this code, sites or spins are considered to be molecules, this explains the notation

! input file must exist otherwise stop and politely ask to fill in the info

inquire (file ='H_QD.input', exist = presence)

if (.not.presence) then
   write (*,*) '**************************************************'
   write (*,*) 'ERROR:'
   write (*,*) 'Please create an H_QD.input file to run this code.'
   write (*,*) 'This file will contain all the system information.'
   write (*,*) '              Stop.'
   write (*,*) '**************************************************'
else
!  write (*,*) '**************************************************'
!  write (*,*) 'INFO:'
!  write (*,*) 'Reading H_QD.input'
!  write (*,*) ' '
   open (unit_input, file='H_QD.input', status='old')
   ! read input file
     read (unit_input,*)
     read (unit_input,*)
     read (unit_input,*)
     read (unit_input,*) Nm  ! Number of molecules (sites or spins), including the electron molecule (site)
                     ! for example, no spin, only one electron level is Nm=1
                     ! one single Ti is Nm=1
                     ! one Ti and one Fe is Nm=2
   
   ! we read info of each site:
   ! the electron level is always the first site
      allocate (Spin (Nm))! allocate spin of each site
      allocate (hx (Nm), hy (Nm), hz (Nm)) ! allocate one local magnetic field per site
      allocate (gx (Nm), gy (Nm), gz (Nm)) ! allocate gyromagnetic factors per site
      allocate (nx (Nm), ny (Nm), nz (Nm)) ! allocate one axis per site
      allocate (B20 (Nm), B22 (Nm), B40 (Nm), B44 (Nm)) ! Stephen Operators per site
         read (unit_input,*)
      do i_m = 1, Nm
         read (unit_input,*) Spin (i_m) ! read spin of each site including electron site

         read (unit_input,*) B20 (i_m), B22 (i_m), B40 (i_m), B44 (i_m) ! read local Stephen coefficients
                                                               ! Stephen coefficient units should be meV
         read (unit_input,*) nx (i_m), ny (i_m), nz (i_m) ! axis of the Stephen operators
                                                
         read (unit_input,*) hx (i_m), hy (i_m), hz (i_m) ! local magnetic field per site
                                                 ! in Teslas
         read (unit_input,*) gx (i_m), gy (i_m), gz (i_m) ! gyromagnetic vector
         read (unit_input,*)
      enddo
   ! we read the connection between sites
   ! in this code the connections are just anisotropic dipolar exchange interactions
     read (unit_input, *) Np ! Number of connected pairs including electronic site
      if (Np /=0) then
         allocate (mol1(Np), mol2(Np)) 
         allocate (Jexch(Np, 3)) ! anisotropic exchange: needs three components
                                 ! in GHz
         do i=1, Np
            read (unit_input,*) mol1 (i), mol2 (i) ! Indices of the exchange-connected molecules
            read (unit_input,*) Jexch(i,1), Jexch(i,2), Jexch(i,3) ! Three-component exchange in GHz
         enddo
      endif
            read (unit_input,*)
   ! finally, information on the electronic level:
     if (int(2*Spin(1)+1) /= 2) then
        write (*, *) ' '
        write (*,*) 'ERROR: '
        write (*,*) ' You are using a spin different from 1/2 for the transport electron!'
        write (*,*) ' Stop. '
        write (*,*) '  '
        stop
     endif
     read (unit_input, *) eps_QD  ! electronic level in meV
     read (unit_input, *) U_Hubbard ! Hubbard U (meV) leading to spin polarization SU(2) symmetry is not broken
     read (unit_input,*)
     read (unit_input,*) Name_output
     read (unit_input,*) Nplot ! Number of states to print in Spin_distribution.dat
     read (unit_input,*) prediag_hamiltonian ! .true. write pre-diagonalized Hamiltonian to PD_HAMIL.dat
     read (unit_input,*) eigenvectors        ! .true. write all eigvenctors to EIGENVECT.dat
   close (unit_input)
   write (*,*) 'INFO:'
   write (*,*) 'Reading Hamiltonian parametres finished '
   write (*,*) '**************************************************'
   write (*,*) ' '
   ! reading is finished
   
endif
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    Transform all quantities into atomic units                                                       !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tesla to a.u.
      hx = gx*hx*BohrMagneton*1000./Hartree; hy = gy*hy*BohrMagneton*1000./Hartree
      hz = gz*hz*BohrMagneton*1000./Hartree
!     print *, hx, hy, hz
! meV to a.u.
      B20 = B20/Hartree; B22=B22/Hartree; B40=B40/Hartree; B44=B44/Hartree
!GHz to a.u.
      if (Np/=0) then
      Jexch=Jexch*GHz/Hartree
      endif
! meV to a.u.
      eps_QD = eps_QD / Hartree; U_Hubbard=U_Hubbard / Hartree

! Check for sanity against the cutoff range of the Fermi integral
   if (Cutoff .lt. (2.*eps_QD+U_Hubbard)) then
        print *, 'WARNING: Cutoff value is smaller than 2eps+U. Fermi integrals may not converge.'
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
! Electronic Hamiltonian: first site, dimension N_in (1) = 4                                          !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate (N_in(Nm))
        N_in(1)=4  !int(2*Spin(1)+1)+1+1=4
        allocate (H_el(N_in(1),N_in(1)))
        H_el = zero
        H_el(1,1) = eps_QD; H_el(2,2) = eps_QD; H_el(3,3)=0._q; H_el(4,4) = 2*eps_QD+U_Hubbard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    Establish basis set                                                                              !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the basis here is a tensorial product of the first state times the second times...
! each state is given by a spin matrix
! the spin matrix is called Ss (:,:,:,:)
!                  the first entry Ss (i,:,:,:) refers to the site
!                  the second entry Ss (:,i,:,:) refers to the x, y, z component of the spin
!                  the third and fourth entry coincide with the Hamiltonian dimensions Hamiltonian (:,:) 

! Total dimension (Hamiltonian entries) :: N
! each spin is of dimension N_in
! up to a spin i_m the total dimension is the product of previous N_in, stored in N_block
      allocate (N_block(Nm))
          N_block (1) = N_in (1)
! Test to remove after debugging
!         print *, 'N_in',1, N_in(1)
!         print *, 'N_block',1, N_block(1)
       do i_m = 2, Nm
          N_in (i_m)=int(2*Spin(i_m)+1)
          N_block (i_m) = N_block (i_m-1) * N_in (i_m)
!         print *, 'N_in',i_m, N_in(i_m)
!         print *, 'N_block',i_m, N_block(i_m)
       enddo
          N = N_block (Nm) ! full dimension
! MEMORY
      allocate (H (N,N)) ! contains the Hamiltonian before diag and the Eigenstates after
      allocate (W (N)) ! Eigenvalues
      allocate (Identity (N,N))
      allocate (Ss (Nm, 4, N, N)) ! The tensorial spin matrix -basis-
      allocate (Sn (3, N, N)) ! rotated spin matrix
      allocate (SProdS (3, N,N)) ! The following are operated spins 
      allocate (SProdS2 (3, N,N)) 
      allocate (Sp (N,N))
      allocate (Sm (N,N))
      allocate (Sp2 (N,N))
      allocate (Sm2 (N,N))
      allocate (Sp4 (N,N))
      allocate (Sm4 (N,N))

! initializations
     H = zero
     Ss = zero
     Sn = zero
     SProdS = zero
     SProdS2 = zero
     Identity = zero
        do i =1,N
           Identity (i,i) = ur
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the tensorial product of spin matrices (aka the basis set):                          !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call Kronecker_Product (N, Nm, N_in, N_block, Ss, Sx_u, Sy_u, Sz_u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the Hamiltonian                                                                      !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     STEPS:
!
!     ZERO: The electronic contribution to the Hamiltonian
!           we use the same algo as for the spin Ss

       do i2 = 1, N_in (1)
        do i3 = 1, N_in (1)

            do i4 = 1, N/N_block(1)

       H (i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  H_el(i2,i3);

            enddo
         enddo
        enddo



!     FIRST: local magnetic fields

      do i_ = 1, Nm ! Loop on molecules or sites
   
      H (:,:) = H (:,:)+Ss (i_, 1, :,:)*hx(i_)+Ss (i_, 2, :,:)*hy(i_)+Ss (i_, 3, :,:)*hz(i_)

!     SECOND: local Stephen operators
!             We need to compute spin matrices up to 4th power.

      ! We define new directions following the preferential axis given in the input file
      ! in this way we can take the internal "z" axis to be always along the preferential axis
      ! First normalize vector input:
          p_mod=sqrt(nx(i_)**2+ny(i_)**2+nz(i_)**2)
          nx (i_) = nx (i_)/p_mod; ny (i_) = ny (i_)/p_mod; nz (i_) = nz (i_)/p_mod
      ! we make vectors px,py,pz,pxx,... perpendicular to nx,ny,nz
       if (nz (i_) /= 0) then
          if (ny (i_) /= 0 ) then
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=nx (i_); pyy=ny (i_);pzz=-(nx(i_)**2+ny (i_)**2)/nz(i_);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else if (nx (i_) /=0) then
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=nx (i_); pyy=ny (i_);pzz=-(nx(i_)**2+ny (i_)**2)/nz(i_);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else
          px = 1.0; py = 0.0; pz=0.0; pxx = 0.0; pyy = 1.0; pzz=0.0
          endif
       else 
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=0; pyy=0;pzz=1._q
       endif


      ! rotated spins, inside the site loop on i_

       Sn (3,:,:)=Ss(i_,1,:,:)*nx (i_)+Ss(i_,2,:,:)*ny (i_)+Ss(i_,3,:,:)*nz (i_)
       Sn (2,:,:)=Ss(i_,1,:,:)*px+Ss(i_,2,:,:)*py+Ss(i_,3,:,:)*pz
       Sn (1,:,:)=Ss(i_,1,:,:)*pxx+Ss(i_,2,:,:)*pyy+Ss(i_,3,:,:)*pzz

      ! We perform matrix multiplication and use the definitions of Stephens Operators
      ! choosing a few ones up to 4th order -See our paper in J. Phys. Chem. A 124, 2318  (2020)
      ! We add over all molecules the contribution of each local anisotropy, this seems
      ! to work, at least for not extremely coupled spins


        call  MatSquareSpin (N, Sn, SprodS)
        call  MatSquareSpin (N, SprodS, SprodS2)

        H (:,:) = H (:,:)+ 3.*B20(i_)*SprodS(3,:,:)  !Longitudinal anisotropy
        H (:,:) = H (:,:)+ B22(i_)*(SprodS(1,:,:)-SprodS(2,:,:)) ! Transversal anisotropy
        H (:,:) = H (:,:)+ B40(i_)*35.*SprodS2(3,:,:) ! Fourth order: B40 ain't pretty
        H (:,:) = H (:,:)- B40(i_)*30.*(Spin(i_)*(Spin(i_)+1._q)-2*Spin(i_))*SprodS(3,:,:)
        H (:,:) = H (:,:)+ B40(i_)*(3*(Spin(i_)*(Spin(i_)+1))**2-6*(Spin(i_)*(Spin(i_)+1)))*Identity (:,:)

        ! change to circular spins: S+ and S-
        Sp (:,:) = Sn (1, :,:)+ui*Sn (2, :,:)
        Sm (:,:) = Sn (1, :,:)-ui*Sn (2, :,:)

        call MatSpSm (N, Sp, Sm, Sp2, Sm2)
        call MatSpSm (N, Sp2, Sm2, Sp4, Sm4)

        H = H - B44(i_)*0.5*ui*(Sp4-Sm4) ! Fourth order: B44, sort of fourth-order longitudinal anisotropy

   
      enddo

!     THIRD: anisotropic intermolecular exchange interactions
!     we keep the original axis, not the one of the anisotropy

      do i = 1, Np ! loop on pairs

!     if (i ==1) then
!     print *, 'computing first pair!!!'
!     else if (i == 2) then
!     print *, 'computing second pair!!!'
!     else
!     print *, 'computing', i,'th pair!!!'
!     endif

      call MatProdSpin (N, i, mol1, mol2, Ss, SProdS) ! mol1 and mol2 contain the indices of the paired molecules

      H (:,:) = H (:,:)+Jexch (i, 1)*SProdS (1,:,:)+Jexch (i, 2)*SProdS (2,:,:)+Jexch (i, 3)*SProdS (3,:,:)

      enddo

     if(prediag_hamiltonian) call print_2darray_complex16( H, 'PD_HAMIL.dat', Hartree/GHz )

     if (runs) then
! DO NOT diagonalize if Name_output exists, and use instead the values
! of a previous run
    inquire (file =Name_output, exist = presence)
    else
        presence = .false. !if runs =.false. do not read
    endif
    if (presence) then
      open (unit_spinoutput, file=Name_output, form='unformatted')
      read (unit_spinoutput) H, W
      close (unit_spinoutput)

           print *, ' '
           print *, 'WARNING!!!!!'
           print *, ' '
           print *, ' Reading the Hamiltonian from a previous run!!!'
           print *, ' '
           print *, 'WARNING!!!!!'
           print *, ' '


       ! check the runs are compatible
         if (N /= size (H,1)) then
           print *, ' '
           print *, Name_output, 'IS NOT COMPATIBLE with H_QD.input'
           print *, 'remove', Name_output, 'from the running folder.'
           print *, 'STOP.'
           print *, ' '
           stop
         endif
    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        DIAGONALIZE                                                                                  !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   print *, ' '
   print *, 'Beginning of diagonalization'
   print *, ' '
         print *, 'size H:', size(H,1), N
      call DIAGONALIZE (N,H,W)
   print *, 'End of diagonalization, these are the 8 first eigenvalues (GHz):'
   write (*, '(8g14.4)') ((W(i)-W(1))*Hartree/GHz, i=1,min(N,8))
   print *, 'End of diagonalization, 8 first eigenvalues (meV):'
   write (*, '(8g14.4)') ((W(i)-W(1))*Hartree, i=1,min(N,8))
   print *, ' '

! Save it, unformatted so it does not need to be recalculated if it 
! exists on the running folder

      open (unit_spinoutput, file=Name_output, form='unformatted')
      write (unit_spinoutput) H, W
      close (unit_spinoutput)

   print *, 'Unformated states and eigenvalues written in file:  ', Name_output
   print *, ' '
! We go ahead and write the ouput 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        Write output in a useful way                                                                 !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print *, 'Eigenvalues in meV are written in Hamiltonian_Eigen.dat'
   print *, 'together with the 8 first complex coefficients of the Eigenstates.'
   print *, ' '
       open (unit_spinoutput, file='Hamiltonian_Eigen.dat')
       do i =1, N
          write (unit_spinoutput,'(3g14.4)') i, W(i)*Hartree, (W(i)-W(1))*Hartree
          write (unit_spinoutput,'(16g14.4)') (real(H(j,i)),aimag(H(j,i)), j= 1, N)
       enddo
       close (unit_spinoutput)

! Plotting the spin distribution for the first Nplot states 


    if (Nplot > N) then
       Nplot = N
    endif

    allocate (spinX(Nplot,Nm),spinY(Nplot,Nm),spinZ(Nplot,Nm))
    allocate (spin2(Nplot))
      spinX = zero
      spinY = zero
      spinZ = zero

      open (unit_spinoutput, file='Spin_distribution.dat')
    do i=1, Nplot
      write (unit_spinoutput,*) 'State=',i
      write (unit_spinoutput,*) '#  Site, Sx, Sy, Sz'

      do j = 1, Nm

      do j1=1, N
      do j2=1, N
      spinX (i,j)=spinX (i,j)+conjg(H(j1,i))*Ss (j, 1, j1, j2)*H(J2, i)
      spinY (i,j)=spinY (i,j)+conjg(H(j1,i))*Ss (j, 2, j1, j2)*H(J2, i)
      spinZ (i,j)=spinZ (i,j)+conjg(H(j1,i))*Ss (j, 3, j1, j2)*H(J2, i)
      enddo
      enddo

      write (unit_spinoutput,'(2x,I2,3x,3g14.4)') j, real(spinX(i,j)), real(spinY(i,j)), real(spinZ(i,j))

      enddo
    enddo
   print *, 'Spins written in file:  Spin_distribution.dat'
   print *, ' '

      spin2 = zero ! spin square

      write (unit_spinoutput,*) ' State ', ' Excitation Energy (GHz) ', ' (meV) ', ' Spin^2 ', ' Spin '

    do i = 1, Nplot

      do j = 1, Nm
       do j4 = 1, Nm
        do j1=1, N
         do j2=1, N
          do j3=1,N
               spin2(i)=spin2(i)+conjg(H(j1,i))*(Ss (j, 1, j1, j3)*Ss (j4, 1, j3, j2)+  &
      &        Ss (j, 2, j1, j3)*Ss (j4, 2, j3, j2)+ Ss (j, 3, j1, j3)*Ss (j4, 3, j3, j2) &
      &        )*H(J2, i)
          enddo
         enddo
        enddo
       enddo
      enddo

      write (unit_spinoutput,'(2x,I2,3x,4g14.4)') i, (W(i)-W(1))*Hartree/GHz, (W(i)-W(1))*Hartree, &
        real(spin2(i)), 0.5*(sqrt(1+4*real(spin2(i)))-1)
    enddo
      close (unit_spinoutput)
      
   if(eigenvectors) then
      print *, 'Complete set of eigenvectors and corresponding eigenvalues (GHz) written to EIGENVECT.dat'
      open (unit_spinoutput, file='EIGENVECT.dat')
      do i =1, N
         write (unit_spinoutput,'(3g20.8)') i, W(i)*Hartree/GHz, (W(i)-W(1))*Hartree/GHz
         write (unit_spinoutput,*) real(H(:,i)), imag(H(:,i))
      enddo
      close (unit_spinoutput)
   end if
   
   print *, 'The Hamiltonian calculation is DONE!'
   write (*,*) '**************************************************'
  endif !end of diverting if Name_output exists

! Calculation of output

    allocate (Eigenvalues (N))
    allocate (Delta (N,N))
    allocate (lambda (N,N,2))

        Eigenvalues = W
    do i = 1, N
    do j = 1, N
        Delta (i,j) = W(i)-W(j)
    enddo
    enddo

     
! Calculation of lambda

    lambda = zero

    do i = 1, N
    do j = 1, N
    do i_sigma= 1, 2 !1 spin down as always in this code

! contribution from |0Xsigma|

      i2 = 3 ! 3 is |0>
      i3 = i_sigma ! 1 is down and 2 is up as corresponds to the basis set

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo

! contribution from |\bar{sigm}aX4| where 4 is the doubly occupied state (singlet)

      i2 = (-1)**(i_sigma+1)+i_sigma! if i_sigma=1 then this yields 2
                            ! if i_sigma=2 then this yields 1
      i3 = 4 ! 4 is the doubly occupied state

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo


    enddo
     ! test lambda
     !write (123,*) i,j, 1, lambda (i,j,1)
     !write (123,*) i,j, 2, lambda (i,j,2)

    enddo
    enddo

    

   end subroutine Hamiltonian
!
! Cutting off dimensions for sizeable calculations
! 
     subroutine redimensioning (Ndim, Delta, bias_R, bias_L, GammaC)
     implicit none
     integer :: Ndim, i, l
     real (q) :: bias_R, bias_L, GammaC
     real (q), dimension (:,:) :: Delta

       l = 0

     do i = 1, Ndim

       if (Delta (i,1) <= abs(bias_R-bias_L)+10*GammaC) then
           l = l+1
       endif

     enddo

       Ndim = l

     return
     end subroutine redimensioning


end module H_QD
