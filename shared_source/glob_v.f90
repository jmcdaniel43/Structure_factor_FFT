module global_variables

   real*8, parameter :: pi = 3.141592654
   integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_N_ATOM=14000 ,  MAX_N_ATOM_TYPE=20

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

   !************* q grid for output structure factor
   real*8,parameter :: q_space=0.001

   ! FFT
   integer, parameter :: pme_max_print=20  ! this should be less than or equal to the setting of pme_grid
   integer,parameter :: pme_grid=20
   integer,parameter :: spline_order=6
   integer, parameter :: spline_grid = 100000
   real*8,dimension(spline_grid) :: B6_spline, B5_spline, B4_spline, B3_spline 
   real*8,dimension(:,:,:) , allocatable :: Q_grid
   complex*16,dimension(:,:,:) , allocatable :: SQ, B

   ! for saving Sq trajectory for correlation functions
   complex*16, dimension(:,:,:), allocatable  :: SQt
   real*8 ,  dimension(:,:,:), allocatable  :: SQ_Ct
   ! this array contains the corresponding wavevector magnitude
   ! for index i_index in SQt(i_index,:,:)
   real*8 , dimension(:), allocatable :: kmag_1D
   ! this stores all the kmag values from all vectors
   real*8 , dimension(:), allocatable :: kmag_1Dall
    ! this is not the most efficient way of doing mapping, but it works.  Each
    ! array element stores integer index
   integer, dimension(pme_grid,pme_grid,pme_grid) :: SQt_map_index
   integer, dimension(3)  :: max_index_list


   ! threads
   integer, parameter :: n_threads=16

   ! if there's an additional command line argument (integer), code will calculate
   ! a partial structure factor for that atom type
   character(3) :: partial_structure_factor
   integer :: partial_structure_factor_index


end module global_variables
