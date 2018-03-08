module global_variables

   real*8, parameter :: pi = 3.141592654
   integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_N_ATOM=90000 ,  MAX_N_ATOM_TYPE=20

   integer, parameter :: max_q_form = 20000
   real*8  :: dq_form
   real*8 , dimension( max_q_form ) :: q_grid_form
   real*8 , dimension( MAX_N_ATOM_TYPE , max_q_form ) :: atomtype_form

   integer :: n_atom_type
   character(MAX_ANAME) , dimension( MAX_N_ATOM_TYPE ) ::  atype_name
   integer , dimension( MAX_N_ATOM ) :: atom_index

   character(3), parameter :: Charge_density_Sq = 'yes'  ! 'yes' for charge density, 'no'  for electron density
   ! if we split S(q) into cation/anion contributions, use these names to
   ! identify ions
   character(5)            :: cation_name, anion_name

   !*** this is for computing total charge density
   real*8, dimension(:), allocatable :: charges

   !************* q grid
   real*8,parameter :: q_spacing=0.01, q_max = 0.8, q_min=0.003 ! q_min should be greater than the smallest value of q used in the atomic form factor grids
   integer, parameter :: n_q = 200

   ! FFT
   integer, parameter :: pme_max_print=60  ! this should be less than or equal to the setting of pme_grid
   integer,parameter :: pme_grid=60
   integer,parameter :: spline_order=6
   integer, parameter :: spline_grid = 100000
   real*8,dimension(spline_grid) :: B6_spline, B5_spline, B4_spline, B3_spline 
   real*8,dimension(:,:,:) , allocatable :: Q_grid
   complex*16,dimension(:,:,:) , allocatable :: SQ, B

   ! threads
   integer, parameter :: n_threads=4

   ! set to yes if averaging over a trajectory
   character(3) :: traj_avg="yes"
   ! if there's an additional command line argument (integer), code will calculate
   ! a partial structure factor for that atom type
   character(3) :: partial_structure_factor
   integer :: partial_structure_factor_index


end module global_variables
