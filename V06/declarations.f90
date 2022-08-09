module declarations

  implicit none

! PARAMETERS

  integer, parameter :: kk = SELECTED_INT_KIND (10)
  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: qs = SELECTED_REAL_KIND(5)
  integer, parameter :: qc = SELECTED_REAL_KIND(10)
  real (q), parameter :: sqr2=1.41421356237309504880_q
  real (q), parameter :: sqr3=1.73205080756887729352_q
  real (q), parameter :: sqr5=2.23606797749978969640_q
  real (q), parameter :: sqr6=2.44948974278317809819_q
  real (q), parameter :: sqr7=2.64575131106459059050_q
  real (q), parameter :: sqr8=2.82842712474619009760_q
  real (q), parameter :: pi_d = 3.14159265358979323846_q
  real (q), parameter :: sqrtpi_d = 1.77245385090551602729_q
  real (q), parameter :: Hartree = 27211.6_q ! meV
  real (q), parameter :: BohrMagneton=5.7883818066E-5_q !eV/T
  real (q), parameter :: GHz = 4.135665E-3_q ! meV
  real (q), parameter :: time_unit = 2.4188843266E-8_q ! nanoseconds
  complex (qc), parameter :: zero=(0._q,0._q), ui = (0._q,1._q)
  complex (qc), parameter ::  ur = (1._q, 0._q)

! numbers
  integer :: i_m, Nplot, Ndim, Nfreq, Ntime, N_int
  integer :: i, N, Nm, Np, i_, j, j1, j2, l
  integer :: j3, j4, i1, i2, i3, i4, u, v, i_sigma
  integer :: Electrode
  real (q) :: eps_QD, U_Hubbard, p_mod, stept
  real (q) :: px, py, pz, pxx, pyy, pzz, suma, gammaC
  real (q) :: bias_R, bias_L, Spin_polarization_R, Spin_polarization_L
  real (q) :: Temperature, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, Cutoff
! arrays
  integer, allocatable :: mol1(:), mol2(:), N_in (:), N_block (:)
  integer, allocatable :: t_seq (:)
  real (q), allocatable :: curr (:)
  real (q), allocatable :: W (:), Jexch (:,:), nx(:), ny(:), nz(:)
  real (q), allocatable :: hx (:), hy (:), hz (:)
  real (q), allocatable :: gx (:), gy (:), gz (:)
  real (q), allocatable :: B20 (:), B22 (:), B40 (:), B44 (:)
  real (q), allocatable :: Spin (:), Eigenvalues (:), Delta (:,:)
  real (q), allocatable :: t0 (:), t1 (:), time (:)
  real (q), allocatable :: Amplitude_seq (:,:), Freq_seq (:,:)
  real (q), allocatable :: Phase_seq (:) 
  complex (qc), allocatable :: spin2 (:)
  complex (qc), allocatable :: Ss(:,:,:,:), Sn (:,:,:), H_el(:,:)
  complex (qc), allocatable :: H (:,:), SprodS (:,:,:), SprodS2 (:,:,:)
  complex (qc), allocatable :: Identity (:,:)
  complex (qc), allocatable :: Sp (:,:), Sm (:,:)
  complex (qc), allocatable :: Sp2 (:,:), Sm2 (:,:)
  complex (qc), allocatable :: Sp4 (:,:), Sm4 (:,:)
  complex (qc), allocatable :: Sx_u (:,:,:), Sy_u (:,:,:), Sz_u (:,:,:)
  complex (qc), allocatable :: spinX (:,:), spinY(:,:), spinZ(:,:)
  complex (qc), allocatable :: lambda (:,:,:)
  complex (qc), allocatable :: G (:,:,:,:,:), GC (:,:,:,:,:)
  complex (qc), allocatable :: rho (:,:,:)
  complex (qc), allocatable:: fermiR_a(:,:), fermiL_a(:,:)
  complex (qc), allocatable:: ufermiR_a(:,:), ufermiL_a(:,:)



! logical
  logical ::  presence, runs, population, spindyn, redimension

! character
   character ( len = 100 ) :: Name_output, output_file, output_fourier, output_ESR


end module declarations
