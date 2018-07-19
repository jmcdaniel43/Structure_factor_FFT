module structure_factor

contains

  !*******************************
  !  This computes both number density and charge density structure factors,
  !  and correlation functions thereof, e.g. S(q,t)
  !  Because many things are computed, there are many similar
  !  datastructures.  Our notation for structure factor data structures
  !  is as follows:
  ! 
  !  1) "Q" denotes vector, "q" denotes scalar, e.g. q=|Q|
  !     thus SQ(Q1,Q2,Q3)(at least 3D array), and Sq(q)
  !    
  !  2) "_a" denotes atom type (or ion type dependence) so a structure factor
  !     that depends on two ion types could be called SQ_a_b, or if we've 
  !     condensed the indices into one array, SQ_ab
  !
  !  3) dynamic structure factor is labeled with "t", such as Sqt. The
  !     associated correlation function of this quantity is labeled with "Ct",
  !     such as SqCt
  !
  !  4) we compute both total charge structure factors, and electron number
  !     density structure factors.  The former is denoted with "c", the later
  !     with "n", so e.g. SQc, SQn
  !     this additionally applies to all other data structures used to compute
  !     structure factor, for example q_1rc , q_1dn 
  !*******************************
  subroutine generate_structure_factor( n_atom, dfti_desc,dfti_desc_inv, traj_file )
    use global_variables
    use routines
    use pme_routines
    use MKL_DFTI
    integer, intent(in) :: n_atom
    TYPE(DFTI_DESCRIPTOR), pointer,intent(inout):: dfti_desc,dfti_desc_inv
    character(*), intent(in)  :: traj_file
    complex*16,dimension(:,:,:),allocatable::FQc, FQn
    real*8, dimension(:,:,:),allocatable :: SqCtc_avg, SqCtn_avg
    real*8, dimension(:,:,:,:,:),allocatable :: SQQc_a_b, SQQn_a_b  ! fortran not case sensitive QQ and q2
    real*8, dimension(:,:,:), allocatable    :: Sq2c_a_b, Sq2n_a_b  ! denote same thing but distinct fortran names
    complex*16, dimension(:,:,:,:),allocatable :: SQc_a, SQn_a
    ! time dependent Sq
    complex*16, dimension(:,:,:), allocatable  :: Sqtc_a, Sqtn_a
    ! correlation functions S(q,t)
    real*8 ,  dimension(:,:,:), allocatable  :: SqCtc_aa, SqCtn_aa
    ! Qgrids
    real*8, dimension(:,:,:), allocatable      :: Qc_grid, Qn_grid
    real*8,dimension(:), allocatable:: kmag_1D_avg, kmag_1Dall_avg
    real*8,dimension(:), allocatable:: q_1rc, q_1rn
    complex*16,dimension(:), allocatable::q_1dc, q_1dn
    real*8,dimension(3,3) :: kk, kk_avg, box
    integer :: n, K
    real*8,dimension(:), allocatable :: charge_iontype, atomic_number_iontype
    real*8,dimension(:,:),allocatable :: xyz, xyz_scale
    integer :: n_atom_kind, nkmag, nkmag_all, i_type, status, n_traj, nmax_tcf, i_step, ifile=99, i_atom
    real*8 :: vol

  
    !  if n_atom_type == 3, then we have solvent molecules present.
    !  if n_atom_type == 2, just cations and anions
    Select Case(n_atom_type)
    Case(2)
      write(*,*) ""
      write(*,*) "No solvent detected, Computing Sq's for pure ionic liquid"
      write(*,*) ""    
    Case(3)
      write(*,*) ""
      write(*,*) "Solvent detected, total Sq will contain solvent contribution"
      write(*,*) ""
    End Select

    ! note in this code we have disabled the use of Charge_density_Sq tag.  It
    ! should always be set to 'yes'
    Select Case(Charge_density_Sq)
    Case('yes')
       continue
    case default
       write(*,*) "must have Charge_density_Sq='yes'. Please change setting in glob_v.f90 "
       stop
    End Select

    ! n_threads set in glob_v.f90
    write(*,*) ""
    write(*,*) "attempting to use ", n_threads , " threads in Q_grid calculation"
    write(*,*) ""

    ! Here we approximate the q-dependence of the atomic form factors using in
    ! the number density scattering structure factor (X-ray) .  This is because
    ! if we were using rigorously elemental form factors, we would have to
    ! compute partial structure factors for all elements individually.
    ! for short wavevectors ( < 1-2 inverse angtroms ), it is a good
    ! approximation to assume that fi(q) = qi * G(q), where G(q) is an
    ! element-independent function.  This is what we do here....
    write(*,*) "Using approximate functional form for elemental form factors"
    write(*,*) "this should provide a good approximation for small wavevectors,"
    write(*,*) "e.g. < 1-2 inverse angstroms.  See code for details..."

    n=spline_order
    K=pme_grid
    ! here we average wavevectors as these change slightly during NPT
    kk_avg=0d0    
 
    allocate( Qc_grid(pme_grid,pme_grid,pme_grid), Qn_grid(pme_grid,pme_grid,pme_grid) )
    allocate( FQc(pme_grid,pme_grid,pme_grid), q_1rc(pme_grid**3), q_1dc(pme_grid**3) )
    allocate( FQn(pme_grid,pme_grid,pme_grid), q_1rn(pme_grid**3), q_1dn(pme_grid**3) )

    ! structure factor data structures
    allocate(  SQc_a(n_atom_type,pme_grid,pme_grid,pme_grid), SQQc_a_b(n_atom_type,n_atom_type,pme_grid,pme_grid,pme_grid) )
    allocate(  SQn_a(n_atom_type,pme_grid,pme_grid,pme_grid), SQQn_a_b(n_atom_type,n_atom_type,pme_grid,pme_grid,pme_grid) )
    allocate(  Sq2c_a_b(pme_grid*pme_grid*pme_grid,n_atom_type,n_atom_type), Sq2n_a_b(pme_grid*pme_grid*pme_grid,n_atom_type,n_atom_type) )   

    SQc_a=0d0;SQn_a=0d0;
    SQQc_a_b=0d0;SQQn_a_b=0d0;

    allocate( kmag_1Dall(pme_grid*pme_grid*pme_grid) )

    ! get number of trajectory snapshots
    call get_n_trajectories( n_traj, n_atom, traj_file )

    ! open trajectory file, so as we read it in on the fly
    open( ifile, file=traj_file, status='old' )

    ! create scaled coordinates
    allocate(xyz(n_atom,3), xyz_scale(n_atom,3), charge_iontype(n_atom), atomic_number_iontype(n_atom))

    ! now loop over trajectories

    do i_step =1, n_traj
       ! read coordinates from trajectory file
       call read_trajectory_snapshot( ifile , xyz , box, n_atom )
       call construct_reciprocal_lattice_vector(kk,vol, box)

       ! get average of reciprocal lattice vectors for end
       kk_avg = kk_avg + kk

       if ( i_step == 1 ) then
          ! initialize Sqt, TCF datastructures
          call Sqt_initialize_datastructures(SqCtc_aa, SqCtn_aa, Sqtc_a, Sqtn_a, n_traj, kk, nmax_tcf )
       end if

          do i_type = 1 , n_atom_type  ! n_atom_type=2, cation atoms=1, anion atoms=2

             ! note this only creates coordinates for atomtype "i_type"
             ! charge_iontype stores charges for atoms in xyz_scale array,
             ! with corresponding indices
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K, charge_iontype, atomic_number_iontype)
             ! two different Q_grids, one for electron number density, Qn_grid
             ! one for total charge density, Qc_grid
             call grid_Q(Qn_grid,Qc_grid,xyz_scale,n_atom_kind,K,n,charge_iontype, atomic_number_iontype)
             q_1rn=RESHAPE(Qn_grid, (/K**3/) )
             q_1rc=RESHAPE(Qc_grid, (/K**3/) )
             q_1dn=cmplx(q_1rn,0.,16)
             q_1dc=cmplx(q_1rc,0.,16)
             status=DftiComputeForward(dfti_desc, q_1dn)
             status=DftiComputeForward(dfti_desc, q_1dc)
             FQn=RESHAPE(q_1dn, (/K,K,K/) )
             FQc=RESHAPE(q_1dc, (/K,K,K/) )
             ! structure factor = B * FQ
             SQn_a(i_type,:,:,:) = FQn*B
             SQc_a(i_type,:,:,:) = FQc*B
          enddo

          ! normalize by volume.  Stot=Sa*Sb/Volume, so multiply each  1 / vol**(1/2)
          SQn_a = SQn_a / vol**(0.5d0)
          SQc_a = SQc_a / vol**(0.5d0)

          ! now create all the cross SQ2 structure factors for this snapshot
          call combine_partial_structure_factors( SQQn_a_b , SQn_a )
          call combine_partial_structure_factors( SQQc_a_b , SQc_a )


       write(*,*) "i_step", i_step
       ! save SQc_a, SQn_a for desired kmag wavevectors
       call Sq_TCF_store( i_step, Sqtn_a, SQn_a )
       call Sq_TCF_store( i_step, Sqtc_a, SQc_a )

    enddo

       kk_avg = kk_avg / dble(n_traj)
       SQQn_a_b = SQQn_a_b / dble(n_traj)
       SQQc_a_b = SQQc_a_b / dble(n_traj)

    ! NPT, calculate average reciprocal lattice vectors
    write(*,*) " average reciprocal lattice_vectors"
    write(*,*) kk_avg(1,:)
    write(*,*) kk_avg(2,:)
    write(*,*) kk_avg(3,:)


    ! compute correlation functions for both number and charge
    ! structure factors
    call Sq_TCF_compute( nmax_tcf, kk_avg, SqCtn_aa, Sqtn_a  )
    call Sq_TCF_compute( nmax_tcf, kk_avg, SqCtc_aa, Sqtc_a  )

    ! average over 3D wavectors to get 1D SQ
    allocate( SqCtc_avg(size(SqCtc_aa(:,1,1)),size(SqCtc_aa(1,:,1)),size(SQCtc_aa(1,1,:)) ), SqCtn_avg(size(SqCtn_aa(:,1,1)),size(SqCtn_aa(1,:,1)),size(SQCtn_aa(1,1,:)) )  )
    allocate( kmag_1D_avg(size(kmag_1D)) )

    ! we need to copy kmag_1D, because we sort it twice.  No reason for this
    ! except different data structures
    allocate(  kmag_1Dall_avg(size(kmag_1Dall)) )

    call Sq_collapse_1D( nkmag, SqCtn_avg , SqCtn_aa , kmag_1D_avg, kmag_1D, q_space )
    call Sq_collapse_1D( nkmag, SqCtc_avg , SqCtc_aa , kmag_1D_avg, kmag_1D, q_space )

    call Sq_collapse_SQ2( nkmag_all, Sq2n_a_b , SQQn_a_b , kmag_1Dall_avg , kmag_1Dall , q_space )
    call Sq_collapse_SQ2( nkmag_all, Sq2c_a_b , SQQc_a_b , kmag_1Dall_avg , kmag_1Dall , q_space )

    ! now print, note because we took the FT in reduced coordinates, we need to
    ! convert to the physical wavevectors
  
    ! when we print correlation functions, we will also print power spectrum so
    ! we will need to do Fourier transform.  Re-initialize the FFT descriptors
    ! to 1D, don't need the old Descriptor anymore
    
    call initialize_FFT_1D_temporal( size(SqCtn_avg(1,1,:)), dfti_desc )
 
    ! print number density S(q)'s, output files will be labeled with suffix='n' 
    call Sq_TCF_print( Sq2n_a_b, kmag_1Dall_avg , nkmag_all, SqCtn_avg, kmag_1D_avg, nmax_tcf, nkmag, "n", dfti_desc )
    ! print charge density S(q)'s, output files will be labeled with suffix='c'
    call Sq_TCF_print( Sq2c_a_b, kmag_1Dall_avg , nkmag_all, SqCtc_avg, kmag_1D_avg, nmax_tcf, nkmag, "c", dfti_desc  )


    deallocate( Qc_grid, Qn_grid, FQc, FQn, q_1rc, q_1rn, q_1dc, q_1dn, SQc_a, SQn_a, SQQc_a_b, SQQn_a_b, Sq2c_a_b, Sq2n_a_b, SqCtc_aa, SqCtn_aa, Sqtc_a, Sqtn_a, SqCtc_avg, SqCtn_avg )

    close( ifile )

  end subroutine generate_structure_factor


  !*********************************
  ! this function gives the approximate form-factor
  ! dependence on q.  We assume this functional form
  ! is element-independent, which is not true, but is
  ! a fairly good approximation for small q < 1-2 inverse angstrom
  ! this approximation saves a ton of computer time, as it means
  ! we don't have to compute partial-structure factors as the
  ! Fi(q) dependence of the form factor comes out of the sum over
  ! atomtypes.  We use a functional form of
  ! F(q) = Atomic_number * exp(-a0 *q**a1 ).  We have determined
  ! a0 and a1 from fits to the form factor for carbon.  Thes parameters are:
  !             a0=0.01 nm^(a1) and a1=1.42.
  !            note we have a conversion below for angstroms
  ! we have verified that these parameters
  ! are also reasonable for oxygen and hydrogen
  real*8 function form_factor_approx( kmag )
       real*8, intent(in) :: kmag
       real*8  :: arg
       real*8, parameter :: a0=0.01 , a1=1.42  ! nm units
   ! kmag is in inverse angstroms, so make sure to convert
   ! also, don't need Atomic number prefactor, as this was
   ! included in the Qgrid and Fourier transform
   arg = kmag * 10d0  ! convert angstrom^-1 to nm^-1

   form_factor_approx = exp(-a0*arg**a1)

  end function form_factor_approx


  !*********************************
  ! this subroutine initializes data structures required to
  ! map from 3D vector arrays "SQ" to 1D vector arrays "Sq"
  ! it also initializes the time dependent structure factor
  ! data arrays
  !
  ! we also construct the array kmag_1D which provides
  ! the magnitude of the wavevector mapped to index  i_index in
  ! Sqt(i_index,i_type,i_step)
  !*********************************
  subroutine Sqt_initialize_datastructures(SqCtc_aa, SqCtn_aa, Sqtc_a, Sqtn_a,n_traj, kk, nmax_tcf )
    use global_variables
    complex*16, allocatable, dimension(:,:,:), intent(inout) :: Sqtc_a, Sqtn_a
    real*8, allocatable, dimension(:,:,:), intent(inout) :: SqCtc_aa, SqCtn_aa
    real*8,dimension(:,:), intent(in) :: kk
    integer, intent(in)               :: n_traj
    integer, intent(out)              :: nmax_tcf

    integer :: i_k , j_k , l_k, i_index, i_index_tot, i_type
    real*8                 :: k_mag, k_vec(3)


    ! use nmax_tcf_frac times the trajectory size for tau
    ! nmax_tcf_frac is read in as input parameter
    nmax_tcf = int( (n_traj-1) * nmax_tcf_frac)


       ! initialize
       max_index_list=0
       i_index=1
       i_index_tot=1  ! this is just used for kmag_1Dall
       do i_k=1,pme_grid
          do j_k=1,pme_grid
             do l_k=1,pme_grid

                ! convert wavevector, note the reciprocal lattice vectors kk
                ! don't have the 2*pi factor
                k_vec(:) = 2 * pi * ( dble(i_k-1) * kk(1,:) +  dble(j_k-1) * kk(2,:) +  dble(l_k-1) * kk(3,:)  )
                k_mag = sqrt(dot_product(k_vec,k_vec))

                ! this just stores all kvec magnitudes
                kmag_1Dall(i_index_tot) = k_mag
                i_index_tot = i_index_tot + 1

                ! see if we're using this wavevector
                if ( k_mag < k_mag_max ) then
                   SQq_map_index(i_k,j_k,l_k) = i_index
                   i_index = i_index + 1
                   if ( i_k > max_index_list(1) ) then
                      max_index_list(1) = i_k
                   endif
                   if ( j_k > max_index_list(2) ) then
                      max_index_list(2) = j_k
                   endif
                   if ( l_k > max_index_list(3) ) then
                      max_index_list(3) = l_k
                   endif
                else
                   SQq_map_index(i_k,j_k,l_k) = 0
                end if

             enddo
          enddo
       enddo

       ! now allocate, in Fortran column major, second array dimension is number
       ! of atom types, should be 2
       i_index = i_index -1
       write(*,*) "will compute correlation function for  ", i_index, "wavevectors"
       allocate( Sqtc_a(i_index, n_atom_type, n_traj), Sqtn_a(i_index, n_atom_type, n_traj)  )
       ! here we assume we are computing correlation function for cation-cation,
       ! anion-anion, and cation-anion, and total, hence four indices
       allocate( SqCtc_aa(i_index, 4, nmax_tcf), SqCtn_aa(i_index, 4, nmax_tcf), kmag_1D(i_index)  )

  end subroutine  Sqt_initialize_datastructures

  


  !*********************************
  ! this subroutine stores structure factors during the course of the simulation
  !*********************************
  subroutine Sq_TCF_store( i_step, Sqt_a, SQ_a )
    use global_variables
    integer, intent(in) :: i_step
    complex*16, dimension(:,:,:), intent(inout) :: Sqt_a
    complex*16, dimension(:,:,:,:), intent(in) :: SQ_a 

    integer :: i_k , j_k , l_k, i_index, i_type

       ! storing Sqt from trajectory
       do i_k=1,max_index_list(1)
          do j_k=1,max_index_list(2)
             do l_k=1,max_index_list(3)
                i_index =  SQq_map_index(i_k,j_k,l_k)
                if ( i_index > 0 ) then
                   do i_type =1, n_atom_type
                      Sqt_a(i_index,i_type,i_step) =  SQ_a(i_type,i_k,j_k,l_k)
                   enddo
                endif
             enddo
          enddo
       enddo

  end subroutine Sq_TCF_store


  !*********************************
  ! this subroutine computes correlation functions of the Structure factors
  ! the structure factors are stored in global variables
  !*********************************
  subroutine Sq_TCF_compute( nmax_tcf, kk_avg, SqCt_aa, Sqt_a  )
    use global_variables
    complex*16, dimension(:,:,:), intent(inout)  :: Sqt_a
    real*8 ,  dimension(:,:,:), intent(inout)  :: SqCt_aa
    integer, intent(in) :: nmax_tcf
    real*8,dimension(:,:), intent(in) :: kk_avg

    integer :: i_k , j_k , l_k, i_index, i_type, i_t0, i_t1, i_tau, i_max
    integer, dimension(:), allocatable :: n_avg
    real*8                 :: norm,k_mag, k_vec(3), SQlocal(4)
    complex*16             :: Sqtotalt0, Sqtotalt1

    ! compute correlation functions
    i_max = size(Sqt_a(1,1,:))-1

    allocate(n_avg(i_max))
    SqCt_aa=0d0
    n_avg=0
    do i_t0 = 1, i_max
       ! use 1/4 the trajectory size for tau for statistics reasons
       do i_tau = 1 , nmax_tcf
          i_t1 = i_t0 + i_tau - 1
          if ( i_t1 <= i_max) then
             n_avg(i_tau) = n_avg(i_tau) + 1
             do i_index=1, size(Sqt_a(:,1,1))
                ! cation-cation
                SqCt_aa(i_index,1,i_tau) = SqCt_aa(i_index,1,i_tau) + 2*dble(Sqt_a(i_index,1,i_t1) * conjg(Sqt_a(i_index,1,i_t0)))
                ! anion-anion
                SqCt_aa(i_index,2,i_tau) = SqCt_aa(i_index,2,i_tau) + 2*dble(Sqt_a(i_index,2,i_t1) * conjg(Sqt_a(i_index,2,i_t0)))
                ! cation-anion
                SqCt_aa(i_index,3,i_tau) = SqCt_aa(i_index,3,i_tau) + dble(Sqt_a(i_index,2,i_t1) * conjg(Sqt_a(i_index,1,i_t0)))
                SqCt_aa(i_index,3,i_tau) = SqCt_aa(i_index,3,i_tau) + dble(Sqt_a(i_index,1,i_t1) * conjg(Sqt_a(i_index,2,i_t0)))
                ! total--this depends on whether there is solvent or not
                Select Case(n_atom_type)
                Case(2)
                  ! no solvent
                  Sqtotalt0 = Sqt_a(i_index,1,i_t0) + Sqt_a(i_index,2,i_t0)  
                  Sqtotalt1 = Sqt_a(i_index,1,i_t1) + Sqt_a(i_index,2,i_t1)
                  SqCt_aa(i_index,4,i_tau) = SqCt_aa(i_index,4,i_tau) + 2*dble(Sqtotalt1 * conjg(Sqtotalt0) )
                Case(3)
                  ! solvent present, total structure factor is stored in 3rd index
                  SqCt_aa(i_index,4,i_tau) = SqCt_aa(i_index,4,i_tau) + 2*dble(Sqt_a(i_index,3,i_t1) * conjg(Sqt_a(i_index,3,i_t0)))                  
                End Select
             enddo
          endif
       enddo
    enddo

    ! now average over trajectory
    do i_tau=1, nmax_tcf
       SqCt_aa(:,:,i_tau) = SqCt_aa(:,:,i_tau) / dble(n_avg(i_tau))
    enddo

    ! now normalize by SqCt(t=0)
    do i_index=1, size(SqCt_aa(:,1,1))
        do i_type=1,4
            norm = SqCt_aa(i_index,i_type,1)
           do i_tau=1, nmax_tcf
                SqCt_aa(i_index,i_type,i_tau) = SqCt_aa(i_index,i_type,i_tau) / norm
           enddo
        enddo 
    enddo

    ! now loop through k indices and fill in kmag_1D array from kk_avg
    do i_k=1,max_index_list(1)
       do j_k=1,max_index_list(2)
          do l_k=1,max_index_list(3)
             i_index =  SQq_map_index(i_k,j_k,l_k)
             if ( i_index > 0 ) then
                k_vec(:) = 2 * pi * ( dble(i_k-1) * kk_avg(1,:) +  dble(j_k-1) * kk_avg(2,:) +  dble(l_k-1) * kk_avg(3,:)  )
                k_mag = sqrt(dot_product(k_vec,k_vec))
                kmag_1D(i_index) = k_mag
             endif
          enddo
       enddo
    enddo


  end subroutine Sq_TCF_compute


  !******************************
  ! IMPORTANT Note:  Input arrays Sq_in, kmag_in will be
  ! changed because they are being sorted
  !
  ! this subroutine maps Sq as a function of 3D kvector
  ! to Sq as a function of kmag, by averaging over orientation
  ! Sq_in and kmag_in are 1D arrays for the kvector index,
  ! where kmag_in(index) stores the kvector magnitude.  Note
  ! there are duplicate magnitudes since we have stored
  ! all the 3D vector info in these arrays
  !******************************
  subroutine Sq_collapse_1D( nkmag, Sq_avg , Sq_in, kmag_avg, kmag_in, q_space  )
    use routines
    integer , intent(out)               :: nkmag
    real*8,dimension(:,:,:), intent(out):: Sq_avg   
    real*8,dimension(:),     intent(out):: kmag_avg
    real*8,dimension(:,:,:), intent(inout) :: Sq_in
    real*8,dimension(:), intent(in) :: kmag_in
    real*8, intent(in)   :: q_space

    real*8, dimension(:), allocatable :: kmag_loc
    integer :: i, j, k, n_avg
    integer :: index_store, index_current, index_local
    real*8  :: kmag_current, kmag_new  

    allocate(kmag_loc(size(kmag_in)) )
    kmag_loc = kmag_in
    ! first, sort the arrays by wavevector magnitude
    call Sort(kmag_loc, Sq_in) 

    ! now average use q_space resolution, if wavevectors are
    ! closer in magnitude, then average


    index_store=1
    index_current=1
    do i=1, size(kmag_loc)-1
       ! look for similar kvectors and average
       if ( index_current >= size(kmag_loc)-1 ) exit  ! reached the end
       n_avg=1
       ! store values for this wavevector, will add to this if averaging
       kmag_current = kmag_loc(index_current)
       kmag_avg(index_store) = kmag_current
       Sq_avg(index_store,:,:) = Sq_in(index_current,:,:)
       ! now look for close value wavevectors
       do j=1, size(kmag_loc)
          index_local = index_current + j
          if ( index_local > size(kmag_loc)) exit
          ! see if differences in kvecs are less than resolution
          if ( q_space > abs( kmag_current - kmag_loc(index_local) )  )  then
             ! add to average
             n_avg = n_avg + 1
             kmag_avg(index_store) = kmag_avg(index_store) + kmag_loc(index_local)
             Sq_avg(index_store,:,:)  =  Sq_avg(index_store,:,:) + Sq_in(index_local,:,:)
          else
             exit
          endif
       enddo
       ! now average for this wavevector and increment indices
       kmag_avg(index_store) = kmag_avg(index_store) / dble(n_avg)
       Sq_avg(index_store,:,:) =  Sq_avg(index_store,:,:) / dble(n_avg)
       index_current = index_current + n_avg
       index_store = index_store + 1
    enddo

    nkmag = index_store -1

  end subroutine Sq_collapse_1D

  !******
  ! this subroutine performs the same kvec averaging
  ! as the above Sq_collapse_Ct subroutine, but needs
  ! to rearrange data structures first, and then calls the routine
  !*************
  subroutine Sq_collapse_SQ2( nkmag_all, Sq2_a_b , SQ2, kmag_1Dall_avg , kmag_1Dall , q_space )
    integer, intent(out)                 :: nkmag_all
    real*8, dimension(:,:,:), intent(out) :: Sq2_a_b
    real*8, dimension(:,:,:,:,:), intent(in) :: SQ2
    real*8, dimension(:), intent(out) :: kmag_1Dall_avg
    real*8, dimension(:), intent(in)  :: kmag_1Dall
    real*8, intent(in)                :: q_space

    real*8, dimension(:,:,:), allocatable :: SQ2_temp
    integer :: i_k, j_k, l_k, i_index_tot, K

    K = size(SQ2(1,1,1,1,:))

    ! this is temporary array for data structure compatability
    allocate(SQ2_temp(K*K*K,size(SQ2(:,1,1,1,1)),size(SQ2(1,:,1,1,1)) )  )

    i_index_tot=1
    do i_k=1,K
       do j_k=1,K
          do l_k=1,K
             SQ2_temp(i_index_tot,:,:) = SQ2(:,:,i_k,j_k,l_k)
             i_index_tot = i_index_tot + 1
          enddo
       enddo
    enddo

    call Sq_collapse_1D( nkmag_all, Sq2_a_b , SQ2_temp, kmag_1Dall_avg, kmag_1Dall, q_space  )

    deallocate(SQ2_temp)

  end subroutine Sq_collapse_SQ2




  !*********************************
  ! this subroutine prints correlation functions of the Structure factors
  ! see definites of structure factor data types in the main subroutine in this
  ! module
  !
  ! we use this subroutine to print both number and charge structure factors
  ! input suffix "n" or "c" will signify what we are printing, and will tag
  ! this suffix onto the file names accordingly
  !
  !  Note that if we're printing number density structure factors,
  !  S(q) will be multipled by generic form factor function F(q)
  !  which is approximation of element-independent function form
  !  see function 'form_factor_approx'
  !
  !  Sq2n_a_b and Sq2c_a_b are "number" and "charge" square structure factors as
  !  a function of (1D) wavevector magnitude and iontypes "a" and "b"
  !
  !  SqCtn_avg, and SqCtc_avg are S(q,t) TCFs for 1D S(q) for cation-cation,
  !  anion-anion, and cation-anion cross structure factors
  !  
  !*********************************
  subroutine Sq_TCF_print( Sq2_a_b, kmag_1Dall_avg , nkmag_all, SqCt_avg, kmag_1D_avg, nmax_tcf, nkmag, suffix, dfti_desc )
    use pme_routines
    use MKL_DFTI
    use global_variables
    real*8,dimension(:,:,:), intent(in) :: Sq2_a_b, SqCt_avg
    real*8,dimension(:) ,    intent(in) :: kmag_1Dall_avg, kmag_1D_avg
    integer, intent(in)                 :: nmax_tcf, nkmag, nkmag_all
    character(1), intent(in)            :: suffix   ! this should be "n" or "c" corresponding to number or charge S(q)
    TYPE(DFTI_DESCRIPTOR), pointer,intent(inout):: dfti_desc

    real*8,dimension(:), allocatable  :: Comega
    integer        :: i_k, i_type, i_t
    real*8         :: time, omega, norm,Sq_print
    character(100) :: filecatCt_out, fileanCt_out, filecatanCt_out, filetotCt_out
    character(100) :: filecatCw_out, fileanCw_out, filecatanCw_out, filetotCw_out
    character(100) :: filecatSq_out, fileanSq_out, filecatanSq_out, filetotSq_out
    character(100) :: filecatSqq2_out, fileanSqq2_out, filecatanSqq2_out, filetotSqq2_out
   
    integer  :: ifile1, ifile2, ifile3, ifile4, ifile5 , ifile6, ifile7, ifile8, ifile9, ifile10, ifile11, ifile12, ifile13, ifile14, ifile15, ifile16

    ! this is for printing correlation function
    !integer, parameter :: nk_print=40

    ! these are output file names
    filecatCt_out = 'SqCt_cation' // suffix // '.xvg'
    fileanCt_out = 'SqCt_anion' // suffix // '.xvg'
    filecatanCt_out = 'SqCt_cat_an' // suffix // '.xvg'
    filetotCt_out = 'SqCt_total' // suffix // '.xvg'
    filecatCw_out = 'SqComega_cation' // suffix // '.xvg'
    fileanCw_out = 'SqComega_anion' // suffix // '.xvg'
    filecatanCw_out = 'SqComega_cat_an' // suffix // '.xvg'
    filetotCw_out = 'SqComega_total' // suffix // '.xvg'
    filecatSq_out = 'Sqavg_cation' // suffix // '.xvg'
    fileanSq_out = 'Sqavg_anion' // suffix // '.xvg'
    filecatanSq_out = 'Sqavg_cat_an' // suffix // '.xvg'
    filetotSq_out = 'Sqavg_total' // suffix // '.xvg'
    filecatSqq2_out = 'Sqq2avg_cation' // suffix // '.xvg'
    fileanSqq2_out = 'Sqq2avg_anion' // suffix // '.xvg'
    filecatanSqq2_out = 'Sqq2avg_cat_an' // suffix // '.xvg'
    filetotSqq2_out = 'Sqq2avg_total' // suffix // '.xvg'

    ifile1=86
    ifile2=87
    ifile3=88
    ifile4=89
    ifile5=90
    ifile6=91
    ifile7=92
    ifile8=93
    ifile9=94
    ifile10=95
    ifile11=96
    ifile12=97
    ifile13=98
    ifile14=99
    ifile15=100
    ifile16=101

    open( ifile1, file=filecatCt_out, status='new' )
    open( ifile2, file=fileanCt_out, status='new' )
    open( ifile3, file=filecatanCt_out, status='new' )
    open( ifile4, file=filetotCt_out, status='new' )
    open( ifile5, file=filecatSq_out, status='new' )
    open( ifile6, file=fileanSq_out, status='new' )
    open( ifile7, file=filecatanSq_out, status='new' )
    open( ifile8, file=filetotSq_out, status='new' )
    open( ifile9, file=filecatSqq2_out, status='new' )
    open( ifile10, file=fileanSqq2_out, status='new' )
    open( ifile11, file=filecatanSqq2_out, status='new' )
    open( ifile12, file=filetotSqq2_out, status='new' )
    open( ifile13, file=filecatCw_out, status='new' )
    open( ifile14, file=fileanCw_out, status='new' )
    open( ifile15, file=filecatanCw_out, status='new' )
    open( ifile16, file=filetotCw_out, status='new' )


    ! first print average Sq, Sq/q2
    ! skip first wavevector which is mag 0
    do i_k=2, nkmag_all 

       ! if number density structure factor, multiply by approximate f(q)
       Select Case(suffix)
       Case('n')
           ! normalize by 1/sum_i(fq_i**2), thus form factors cancel out in our
           ! approximation of same q dependence for f(q) of each element
           !norm = form_factor_approx( kmag_1Dall_avg(i_k) )
           norm = 1d0 / dble(fq_norm)
       Case default
           ! convert from q^2/A to kJ/mol
           norm = 2625.4996 / 1.88973
           ! normalize for sum-rule
           ! for non-polarizable simulations, we have 4*pi/kT * Sq/q**2 = 1
           norm = norm * 4d0 * pi / ( 0.008314462d0 * temperature )
       End Select

       ! cation
       Sq_print = norm*Sq2_a_b(i_k,1,1)
       write(ifile5,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       Sq_print = norm*Sq2_a_b(i_k,1,1) / kmag_1Dall_avg(i_k)**2
       write(ifile9,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! anion
       Sq_print = norm*Sq2_a_b(i_k,2,2)
       write(ifile6,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       Sq_print = norm*Sq2_a_b(i_k,2,2) / kmag_1Dall_avg(i_k)**2
       write(ifile10,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! cation-anion
       Sq_print = norm*Sq2_a_b(i_k,1,2)
       write(ifile7,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       Sq_print = norm*Sq2_a_b(i_k,1,2) / kmag_1Dall_avg(i_k)**2
       write(ifile11,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! total
       ! depends on whether we have solvent or not...
       Select Case(n_atom_type)
         Case(2)
         ! no solvent
         Sq_print = norm*(Sq2_a_b(i_k,1,1) + Sq2_a_b(i_k,2,2) + 2 * Sq2_a_b(i_k,1,2))
         Case(3)
         ! solvent presenti
       Sq_print = norm*Sq2_a_b(i_k,3,3)
       End Select
       write(ifile8,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       Sq_print = Sq_print / kmag_1Dall_avg(i_k)**2
       write(ifile12,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
    enddo

    ! this is for FT
    allocate( Comega(size(SqCt_avg(1,1,:))) )

    ! now print correlation function, all wavevectors to each file
    ! don't need norm prefactor here, as we normalize to Ct(0)

!    if ( nk_print > nkmag ) then
!       write(*,*) "nk_print > nkmag , please increase setting of nkmag "
!       write(*,*) "ACF calculation"
!       stop
!    end if

    do i_k=2, nkmag   ! this would print a lot of correlation functions
    !   dt is in ps
    !do i_k=2, nk_print
       write(ifile1,*) ""
       write(ifile1,'(A30,F14.4)') "# C(t) (ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile2,*) ""
       write(ifile2,'(A30,F14.4)') "# C(t) (ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile3,*) ""
       write(ifile3,'(A30,F14.4)') "# C(t) (ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile4,*) ""
       write(ifile4,'(A30,F14.4)') "# C(t) (ps) for wavevector" , kmag_1D_avg(i_k)

       write(ifile13,*) ""
       write(ifile13,'(A30,F14.4)') "# C(omega) (rad/ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile14,*) ""
       write(ifile14,'(A30,F14.4)') "# C(omega) (rad/ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile15,*) ""
       write(ifile15,'(A30,F14.4)') "# C(omega) (rad/ps) for wavevector" , kmag_1D_avg(i_k)
       write(ifile16,*) ""
       write(ifile16,'(A30,F14.4)') "# C(omega) (rad/ps) for wavevector" , kmag_1D_avg(i_k)

       ! cation C(t)
       do i_t=1,nmax_tcf
          time=(i_t-1)*time_step
          write(ifile1,'(F14.1, E20.6)') time , SqCt_avg(i_k,1,i_t) 
       enddo
       ! cation C(w)
       Comega(:) = SqCt_avg(i_k,1,:)
       call compute_FFT_1D_temporal( Comega , dfti_desc )
       do i_t=1,nmax_tcf
             omega = 2d0 * pi * dble(i_t - 1) / (dble(nmax_tcf) * time_step ) 
             write(ifile13,'(E14.6, E20.6)') omega , Comega(i_t)
       enddo

       ! anion C(t)
       do i_t=1,nmax_tcf
          time=(i_t-1)*time_step
          write(ifile2,'(F14.1, E20.6)') time , SqCt_avg(i_k,2,i_t)
       enddo
       ! anion C(w)
       Comega(:) = SqCt_avg(i_k,2,:)
       call compute_FFT_1D_temporal( Comega , dfti_desc )
       do i_t=1,nmax_tcf
             omega = 2d0 * pi * dble(i_t - 1) / (dble(nmax_tcf) * time_step )
             write(ifile14,'(E14.6, E20.6)') omega , Comega(i_t)
       enddo

       ! cation-anion C(t)
       do i_t=1,nmax_tcf
          time=(i_t-1)*time_step
          write(ifile3,'(F14.1, E20.6)') time , SqCt_avg(i_k,3,i_t)
       enddo
       ! cation-anion C(w)
       Comega(:) = SqCt_avg(i_k,3,:)
       call compute_FFT_1D_temporal( Comega , dfti_desc )
       do i_t=1,nmax_tcf
             omega = 2d0 * pi * dble(i_t - 1) / (dble(nmax_tcf) * time_step )
             write(ifile15,'(E14.6, E20.6)') omega , Comega(i_t)
       enddo


       ! total C(t)
       do i_t=1,nmax_tcf
          time=(i_t-1)*time_step
          write(ifile4,'(F14.1, E20.6)') time , SqCt_avg(i_k,4,i_t)
       enddo
       ! total C(w)
       Comega(:) = SqCt_avg(i_k,4,:)
       call compute_FFT_1D_temporal( Comega , dfti_desc )
       do i_t=1,nmax_tcf
             omega = 2d0 * pi * dble(i_t - 1) / (dble(nmax_tcf) * time_step )
             write(ifile16,'(E14.6, E20.6)') omega , Comega(i_t)
       enddo


    enddo

    close(ifile1)
    close(ifile2)
    close(ifile3)
    close(ifile4)
    close(ifile5)
    close(ifile6)
    close(ifile7)
    close(ifile8)
    close(ifile9)
    close(ifile10)
    close(ifile11)
    close(ifile12)
    close(ifile13)
    close(ifile14)
    close(ifile15)
    close(ifile16)


  end subroutine Sq_TCF_print


  subroutine  combine_partial_structure_factors( SQ2 , SQ )
    use global_variables
    real*8,dimension(:,:,:,:,:),intent(inout) :: SQ2
    complex*16,dimension(:,:,:,:),intent(in) :: SQ

    integer :: i_type, j_type

    do i_type=1,n_atom_type
       do j_type = i_type, n_atom_type

          if ( i_type == j_type ) then
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + dble(SQ(i_type,:,:,:))**2+aimag(SQ(i_type,:,:,:))**2
          else
             ! here we have unlike types, so we want SQi(c.c.) * SQj + SQi * SQj(c.c.) , where
             ! (c.c.) is complex conjugate, which equals 2 * real ( SQi * SQj(c.c.) )
             !SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + 2d0 * dble( SQ(i_type,:,:,:) * conjg((SQ(j_type,:,:,:))) )
             ! multiply by 0.5 for double counting
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + 0.5d0 *  2d0 * dble( SQ(i_type,:,:,:) * conjg((SQ(j_type,:,:,:))) )
          end if
       enddo
    enddo

  end subroutine combine_partial_structure_factors



end module structure_factor
