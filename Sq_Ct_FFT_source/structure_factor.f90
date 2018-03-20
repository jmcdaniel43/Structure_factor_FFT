module structure_factor

contains

  !*******************************
  !  This computes both number density and charge density structure factors,
  !  and correlation functions thereof, e.g. S(q,t)
  !  Because many things are computed, there are many similar
  !  datastructures that could create confusion.
  !  Here we outline the data structures
  ! 
  !  SQ(:,:,:) is a structure factor for a specific snapshot, for a specific
  !  atom type (cations, anions, other atom type, etc.)
  ! 
  !  SQ_store(:,:,:,:) stores the structure factors for all atomtypes
  !  e.g. SQ_store(i_type,:,:,:), for a snapshot
  !
  !  SQ2_store keeps a running sum of atomtype-atomtype scattering structure
  !  factors throughout the trajectory, e.g. SQ_store(i_type,j_type,:,:,:) +=
  !  SQ_store(i_type,:,:,:) *complexconjugate SQ_store(j_type,:,:,:)
  !
  !  SQt stores Structure factor trajectory for select wavevectors,
  !  see Sq_TCF subroutine (this array is stored as global variable)
  !*******************************
  subroutine generate_structure_factor( n_atom, dfti_desc,dfti_desc_inv, traj_file )
    use global_variables
    use routines
    use pme_routines
    use MKL_DFTI
    integer, intent(in) :: n_atom
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    character(*), intent(in)  :: traj_file

    complex*16,dimension(:,:,:),allocatable::FQ
    real*8, dimension(:,:,:),allocatable :: SQ2, SQ_Ct_avg
    real*8, dimension(:,:,:,:,:),allocatable :: SQ2_store
    real*8, dimension(:,:,:), allocatable    :: SQ2_avg
    complex*16, dimension(:,:,:,:),allocatable :: SQ_store
    real*8,dimension(:), allocatable::q_1r, kmag_1D_avg, kmag_1Dall_avg
    complex*16,dimension(:), allocatable::q_1d
    real*8,dimension(3,3) :: kk, kk_avg, box
    integer :: n, K
    real*8,dimension(:), allocatable :: charge_iontype
    real*8,dimension(:,:),allocatable :: xyz, xyz_scale
    integer :: n_atom_kind, nkmag, nkmag_all, i_type, status, n_traj, nmax_tcf, i_step, ifile=99, i_atom

    n=spline_order
    K=pme_grid
    ! here we average wavevectors as these change slightly during NPT
    kk_avg=0d0    

    allocate( FQ(pme_grid,pme_grid,pme_grid), SQ2(pme_grid,pme_grid,pme_grid), q_1r(pme_grid**3), q_1d(pme_grid**3) )
    SQ2=0d0

    ! we need these arrays if doing a full scattering structure factor and also
    ! if doing an electron density structure factor.  The only time we don't
    ! need the arrays is a partial structure factors, so might as well allocate
    ! them no matter what
    allocate(  SQ_store(n_atom_type,pme_grid,pme_grid,pme_grid), SQ2_store(n_atom_type,n_atom_type,pme_grid,pme_grid,pme_grid) )
    SQ_store=0d0
    SQ2_store=0d0

    allocate( kmag_1Dall(pme_grid*pme_grid*pme_grid) )
    ! this will store the radial average over wavevectors of SQ2_store
    allocate(  SQ2_avg(pme_grid*pme_grid*pme_grid,n_atom_type,n_atom_type) )


    ! get number of trajectory snapshots
    call get_n_trajectories( n_traj, n_atom, traj_file )
    ! open trajectory file, so as we read it in on the fly
    open( ifile, file=traj_file, status='old' )

    ! create scaled coordinates
    allocate(xyz(n_atom,3), xyz_scale(n_atom,3), charge_iontype(n_atom))

    ! now loop over trajectories

    do i_step =1, n_traj
       ! read coordinates from trajectory file
       call read_trajectory_snapshot( ifile , xyz , box, n_atom )
       call construct_reciprocal_lattice_vector(kk, box)

       ! get average of reciprocal lattice vectors for end
       kk_avg = kk_avg + kk

       if ( i_step == 1 ) then
          ! initialize Sq_TCF with initial kvectors, n_traj
          call Sq_TCF_store( i_step, SQ_store, kk, n_traj )
       end if

       !*************** Decide whether we're doing total charge structure factor
       ! or X-ray scattering, electron density structure factor
       Select Case(Charge_density_Sq)
       Case('yes')
          do i_type = 1 , n_atom_type  ! n_atom_type=2, cation atoms=1, anion atoms=2

             ! note this only creates coordinates for atomtype "i_type"
             ! charge_iontype stores charges for atoms in xyz_scale array,
             ! with corresponding indices
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K, charge_iontype)
             call grid_Q(Q_grid,xyz_scale,n_atom_kind,K,n, charge_iontype)
             q_1r=RESHAPE(Q_grid, (/K**3/) )
             q_1d=cmplx(q_1r,0.,16)
             status=DftiComputeForward(dfti_desc, q_1d)
             FQ=RESHAPE(q_1d, (/K,K,K/) )
             ! structure factor = B * FQ
             SQ = FQ*B
             SQ_store(i_type,:,:,:) = SQ(:,:,:)
          enddo

          ! now create all the cross SQ2 structure factors for this snapshot
          call combine_partial_structure_factors( SQ2_store , SQ_store )

       Case default

          Select Case(partial_structure_factor)
          Case("no")
             !****************** computing full structure factor
             ! loop over  atom types for structure factor
             do i_type = 1 , n_atom_type
                ! note this only creates coordinates for atomtype "i_type"
                call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
                call grid_Q(Q_grid,xyz_scale,n_atom_kind,K,n)
                q_1r=RESHAPE(Q_grid, (/K**3/) )
                q_1d=cmplx(q_1r,0.,16)
                status=DftiComputeForward(dfti_desc, q_1d)
                FQ=RESHAPE(q_1d, (/K,K,K/) )
                ! structure factor = B * FQ
                SQ = FQ*B
                SQ_store(i_type,:,:,:) = SQ(:,:,:)
             enddo

             ! now create all the cross SQ2 structure factors for this snapshot
             call combine_partial_structure_factors( SQ2_store , SQ_store )

          Case("yes")
             !********************* computing partial structure factor
             i_type = partial_structure_factor_index
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
             call grid_Q(Q_grid,xyz_scale,n_atom_kind,K,n)
             q_1r=RESHAPE(Q_grid, (/K**3/) )
             q_1d=cmplx(q_1r,0.,16)
             status=DftiComputeForward(dfti_desc, q_1d)
             FQ=RESHAPE(q_1d, (/K,K,K/) )
             ! structure factor = B * FQ
             SQ = FQ*B
             SQ2 = SQ2 + dble(SQ)**2+aimag(SQ)**2
          End Select


       End Select

       write(*,*) "i_step", i_step
       ! save SQ_store at this snapshot to compute TCF
       call Sq_TCF_store( i_step, SQ_store, kk )

    enddo


    Select Case(Charge_density_Sq)
    Case('yes')
       kk_avg = kk_avg / dble(n_traj)
       SQ2_store = SQ2_store / n_traj
    Case('no')
       ! now average
       Select Case(partial_structure_factor)
       Case("no")
          kk_avg = kk_avg / dble(n_traj)
          SQ2_store = SQ2_store / n_traj
          ! now add partial structure factors
          call add_partial_structure_factors( SQ2 , SQ2_store, kk )
       Case("yes")
          kk_avg = kk_avg / dble(n_traj)
          SQ2 = SQ2 / n_traj
       end Select
    End Select

    ! NPT, calculate average reciprocal lattice vectors
    write(*,*) " average reciprocal lattice_vectors"
    write(*,*) kk_avg(1,:)
    write(*,*) kk_avg(2,:)
    write(*,*) kk_avg(3,:)


    ! compute correlation functions
    call Sq_TCF_compute( nmax_tcf, kk_avg  )

    ! average over 3D wavectors to get 1D SQ
    allocate( SQ_Ct_avg(size(SQ_Ct(:,1,1)),size(SQ_Ct(1,:,1)),size(SQ_Ct(1,1,:)) ), kmag_1D_avg(size(kmag_1D)) )
    ! we need to copy kmag_1D, because we sort it twice.  No reason for this
    ! except different data structures
    allocate(  kmag_1Dall_avg(size(kmag_1Dall)) )
    call Sq_collapse_1D( nkmag, SQ_Ct_avg , SQ_Ct , kmag_1D_avg, kmag_1D, q_space )
    call Sq_collapse_SQ2( nkmag_all, SQ2_avg , SQ2_store , kmag_1Dall_avg , kmag_1Dall , q_space )

    ! now print, note because we took the FT in reduced coordinates, we need to
    ! convert to the physical wavevectors
    Select Case(Charge_density_Sq)
    Case('yes')
       call Sq_TCF_print( SQ2_avg, kmag_1Dall_avg , nkmag_all, SQ_Ct_avg, kmag_1D_avg, nmax_tcf, nkmag  )
    case default 
       write(*,*) " k vec ,  SQ^2 "
       call print_SQ( SQ2 , kk_avg, pme_max_print )
    End Select

    deallocate( FQ, SQ2, q_1r, q_1d )

    Select Case(partial_structure_factor)
    Case("no")
       deallocate( SQ_store, SQ2_store)
    End Select

    close( ifile )

  end subroutine generate_structure_factor



  subroutine add_partial_structure_factors( SQ2 , SQ2_store,kk )
    use global_variables
    real*8,dimension(:,:,:), intent(out) :: SQ2
    real*8,dimension(:,:,:,:,:), intent(in) :: SQ2_store
    real*8,dimension(3,3) , intent(in) :: kk

    integer :: i_type, j_type, i_k, j_k, l_k
    real*8  :: fi, fj, k_vec(3), k_mag

    SQ2=0d0

    ! here we add only those structure factors that will be printed, to save time
    do l_k = 1, pme_max_print
       do j_k = 1, pme_max_print
          do i_k = 1, pme_max_print

             ! convert wavevector, note the reciprocal lattice vectors kk don't have the 2*pi factor
             k_vec(:) = 2 * pi * ( dble(i_k-1) * kk(1,:) +  dble(j_k-1) * kk(2,:) +  dble(l_k-1) * kk(3,:)  )
             k_mag = sqrt(dot_product(k_vec,k_vec))

             do i_type=1,n_atom_type
                do j_type = i_type, n_atom_type

                   fi = get_form_fac( i_type, k_mag )
                   fj = get_form_fac( j_type, k_mag )

                   SQ2(i_k,j_k,l_k) = SQ2(i_k,j_k,l_k) + fi * fj * SQ2_store(i_type,j_type,i_k,j_k,l_k)
		enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine add_partial_structure_factors


  !*********************************
  ! this subroutine stores structure factors during the course of the simulation
  ! and finally computes correlation functions of these quantities
  !
  ! when we initialize, we construct the array kmag_1D which provides
  ! the magnitude of the wavevector mapped to index  i_index in
  ! SQt(i_index,i_type,i_step)
  !*********************************
  subroutine Sq_TCF_store( i_step, SQ_store, kk, n_traj )
    use global_variables
    integer, intent(in) :: i_step
    complex*16, dimension(:,:,:,:), intent(in) :: SQ_store 
    real*8,dimension(:,:), intent(in) :: kk
    integer, intent(in), optional     :: n_traj

    integer :: i_k , j_k , l_k, i_index, i_index_tot, i_type, i_t0, i_t1, i_tau, i_max
    real*8, parameter :: k_mag_max = 1.0  ! maximum wavevector for computing correlation function in inverse angstrom
    real*8                 :: k_mag, k_vec(3)

    if ( present(n_traj) ) then
       write(*,*) "initializing SQt array mapping"
       ! initialize
       max_index_list=0
       i_index=1
       i_index_tot=1  ! this is just used for kmag_1Dall
       do i_k=1,pme_grid
          do j_k=1,pme_grid
             do l_k=1,pme_grid

                ! convert wavevector, note the reciprocal lattice vectors kk don't
                ! have the 2*pi factor
                k_vec(:) = 2 * pi * ( dble(i_k-1) * kk(1,:) +  dble(j_k-1) * kk(2,:) +  dble(l_k-1) * kk(3,:)  )
                k_mag = sqrt(dot_product(k_vec,k_vec))

                ! this just stores all kvec magnitudes
                kmag_1Dall(i_index_tot) = k_mag
                i_index_tot = i_index_tot + 1

                ! see if we're using this wavevector
                if ( k_mag < k_mag_max ) then
                   SQt_map_index(i_k,j_k,l_k) = i_index
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
                   SQt_map_index(i_k,j_k,l_k) = 0
                end if

             enddo
          enddo
       enddo

       ! now allocate, in Fortran column major, second array dimension is number of
       ! atom types, should be 2
       i_index = i_index -1
       write(*,*) "will compute correlation function for  ", i_index, " wavevectors"
       allocate( SQt(i_index, n_atom_type, n_traj)  )
       ! here we assume we are computing correlation function for cation-cation,
       ! anion-anion, and cation-anion, hence three indices
       allocate( SQ_Ct(i_index, 3, n_traj), kmag_1D(i_index)  )   

    else

       ! storing Sqt from trajectory
       do i_k=1,max_index_list(1)
          do j_k=1,max_index_list(2)
             do l_k=1,max_index_list(3)
                i_index =  SQt_map_index(i_k,j_k,l_k)
                if ( i_index > 0 ) then
                   do i_type =1, n_atom_type
                      SQt(i_index,i_type,i_step) =  SQ_store(i_type,i_k,j_k,l_k)
                   enddo
                endif
             enddo
          enddo
       enddo


       !  if ( i_step == size(SQt(1,1,:)) ) then
       !    write(*,*) max_index_list
       !    i_index =  SQt_map_index(1,1,1)
       !    do i_k = 1, 2
       !       write(*,*) SQt(i_index,i_k, i_step)
       !       write(*,*) SQt(i_index,i_k, i_step-1)
       !    enddo
       !   stop
       !  end if


    end if

  end subroutine Sq_TCF_store


  !*********************************
  ! this subroutine computes correlation functions of the Structure factors
  ! the structure factors are stored in global variables
  !*********************************
  subroutine Sq_TCF_compute( nmax_tcf, kk_avg  )
    use global_variables
    integer, intent(out) :: nmax_tcf
    real*8,dimension(:,:), intent(in) :: kk_avg

    integer :: i_k , j_k , l_k, i_index, i_type, i_t0, i_t1, i_tau, i_max
    integer, dimension(:), allocatable :: n_avg
    real*8                 :: k_mag, k_vec(3), SQlocal(4)

    ! compute correlation functions
    i_max = size(SQt(1,1,:))-1
    ! use 1/4 the trajectory size for tau for statistics reasons
    nmax_tcf = i_max/4

    allocate(n_avg(i_max))
    SQ_Ct=0d0
    n_avg=0
    do i_t0 = 1, i_max
       ! use 1/4 the trajectory size for tau for statistics reasons
       do i_tau = 1 , nmax_tcf
          i_t1 = i_t0 + i_tau
          if ( i_t1 <= i_max) then
             n_avg(i_tau) = n_avg(i_tau) + 1
             do i_index=1, size(SQt(:,1,1))
                ! cation-cation
                SQ_Ct(i_index,1,i_tau) = SQ_Ct(i_index,1,i_tau) + 2*dble(SQt(i_index,1,i_t1) * conjg(SQt(i_index,1,i_t0)))
                ! anion-anion
                SQ_Ct(i_index,2,i_tau) = SQ_Ct(i_index,2,i_tau) + 2*dble(SQt(i_index,2,i_t1) * conjg(SQt(i_index,2,i_t0)))
                ! cation-anion
                SQ_Ct(i_index,3,i_tau) = SQ_Ct(i_index,3,i_tau) + dble(SQt(i_index,2,i_t1) * conjg(SQt(i_index,1,i_t0)))
                SQ_Ct(i_index,3,i_tau) = SQ_Ct(i_index,3,i_tau) + dble(SQt(i_index,1,i_t1) * conjg(SQt(i_index,2,i_t0)))
             enddo
          endif
       enddo
    enddo

    ! now average over trajectory
    do i_tau=1, nmax_tcf
       SQ_Ct(:,:,i_tau) = SQ_Ct(:,:,i_tau) / dble(n_avg(i_tau))
    enddo


    ! now loop through k indices and fill in kmag_1D array from kk_avg
    do i_k=1,max_index_list(1)
       do j_k=1,max_index_list(2)
          do l_k=1,max_index_list(3)
             i_index =  SQt_map_index(i_k,j_k,l_k)
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
    real*8,dimension(:), intent(inout) :: kmag_in
    real*8, intent(in)   :: q_space

    integer :: i, j, k, n_avg
    integer :: index_store, index_current, index_local
    real*8  :: kmag_current, kmag_new  

    ! first, sort the arrays by wavevector magnitude
    call Sort(kmag_in, Sq_in) 

    ! now average use q_space resolution, if wavevectors are
    ! closer in magnitude, then average


    index_store=1
    index_current=1
    do i=1, size(kmag_in)-1
       ! look for similar kvectors and average
       if ( index_current >= size(kmag_in)-1 ) exit  ! reached the end
       n_avg=1
       ! store values for this wavevector, will add to this if averaging
       kmag_current = kmag_in(index_current)
       kmag_avg(index_store) = kmag_current
       Sq_avg(index_store,:,:) = Sq_in(index_current,:,:)
       ! now look for close value wavevectors
       do j=1, size(kmag_in)
          index_local = index_current + j
          if ( index_local > size(kmag_in)) exit
          ! see if differences in kvecs are less than resolution
          if ( q_space > abs( kmag_current - kmag_in(index_local) )  )  then
             ! add to average
             n_avg = n_avg + 1
             kmag_avg(index_store) = kmag_avg(index_store) + kmag_in(index_local)
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
  subroutine Sq_collapse_SQ2( nkmag_all, SQ2_avg , SQ2_store , kmag_1Dall_avg , kmag_1Dall , q_space )
    integer, intent(out)                 :: nkmag_all
    real*8, dimension(:,:,:), intent(out) :: SQ2_avg
    real*8, dimension(:,:,:,:,:), intent(in) :: SQ2_store
    real*8, dimension(:), intent(out) :: kmag_1Dall_avg
    real*8, dimension(:), intent(inout)  :: kmag_1Dall
    real*8, intent(in)                :: q_space

    real*8, dimension(:,:,:), allocatable :: SQ2_temp
    integer :: i_k, j_k, l_k, i_index_tot, K

    K = size(SQ2_store(1,1,1,1,:))

    ! this is temporary array for data structure compatability
    allocate(SQ2_temp(K*K*K,size(SQ2_store(:,1,1,1,1)),size(SQ2_store(1,:,1,1,1)) )  )

    i_index_tot=1
    do i_k=1,K
       do j_k=1,K
          do l_k=1,K
             SQ2_temp(i_index_tot,:,:) = SQ2_store(:,:,i_k,j_k,l_k)
             i_index_tot = i_index_tot + 1
          enddo
       enddo
    enddo

    call Sq_collapse_1D( nkmag_all, SQ2_avg , SQ2_temp, kmag_1Dall_avg, kmag_1Dall, q_space  )

    deallocate(SQ2_temp)

  end subroutine Sq_collapse_SQ2




  !*********************************
  ! this subroutine prints correlation functions of the Structure factors
  !*********************************
  subroutine Sq_TCF_print( SQ2_avg, kmag_1Dall_avg , nkmag_all, SQ_Ct_avg, kmag_1D_avg, nmax_tcf, nkmag  )
    use global_variables
    real*8,dimension(:,:,:), intent(in) :: SQ2_avg, SQ_Ct_avg
    real*8,dimension(:) ,    intent(in) :: kmag_1Dall_avg, kmag_1D_avg
    integer, intent(in)                 :: nmax_tcf, nkmag, nkmag_all

    integer        :: i_k, i_type, i_t
    real*8         :: Sq_print
    character(100) :: filecatCt_out, fileanCt_out, filecatanCt_out, filetotCt_out
    character(100) :: filecatSq_out, fileanSq_out, filecatanSq_out, filetotSq_out
    character(100) :: filecatSqq2_out, fileanSqq2_out, filecatanSqq2_out, filetotSqq2_out
    integer  :: ifile1, ifile2, ifile3, ifile4, ifile5 , ifile6, ifile7, ifile8, ifile9, ifile10, ifile11, ifile12

    ! this is for printing correlation function
    integer, parameter :: nk_print=15

    ! these are output file names
    filecatCt_out = 'SqCt_cation.xvg'
    fileanCt_out = 'SqCt_anion.xvg'
    filecatanCt_out = 'SqCt_cat_an.xvg'
    filetotCt_out = 'SqCt_total.xvg'
    filecatSq_out = 'Sqavg_cation.xvg'
    fileanSq_out = 'Sqavg_anion.xvg'
    filecatanSq_out = 'Sqavg_cat_an.xvg'
    filetotSq_out = 'Sqavg_total.xvg'
    filecatSqq2_out = 'Sqq2avg_cation.xvg'
    fileanSqq2_out = 'Sqq2avg_anion.xvg'
    filecatanSqq2_out = 'Sqq2avg_cat_an.xvg'
    filetotSqq2_out = 'Sqq2avg_total.xvg'

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



    ! first print average Sq, Sq/q2
    ! skip first wavevector which is mag 0
    do i_k=2, nkmag_all 
       ! cation
       write(ifile5,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), SQ2_avg(i_k,1,1)
       Sq_print = SQ2_avg(i_k,1,1) / kmag_1Dall_avg(i_k)**2
       write(ifile9,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! anion
       write(ifile6,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), SQ2_avg(i_k,2,2)
       Sq_print = SQ2_avg(i_k,2,2) / kmag_1Dall_avg(i_k)**2
       write(ifile10,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! cation-anion
       write(ifile7,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), SQ2_avg(i_k,1,2)
       Sq_print = SQ2_avg(i_k,1,2) / kmag_1Dall_avg(i_k)**2
       write(ifile11,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       ! total
       Sq_print = SQ2_avg(i_k,1,1) + SQ2_avg(i_k,2,2) + 2 * SQ2_avg(i_k,1,2)
       write(ifile8,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
       Sq_print = Sq_print / kmag_1Dall_avg(i_k)
       write(ifile12,'(F14.6, E20.6)') kmag_1Dall_avg(i_k), Sq_print
    enddo


    ! now print correlation function, all wavevectors to each file
    !   do i_k=2, nkmag   ! this would print a lot of correlation functions
    do i_k=2, nk_print
       write(ifile1,*) ""
       write(ifile1,*) "# C(t) for wavevector" , kmag_1D_avg(i_k)
       write(ifile2,*) ""
       write(ifile2,*) "# C(t) for wavevector" , kmag_1D_avg(i_k)
       write(ifile3,*) ""
       write(ifile3,*) "# C(t) for wavevector" , kmag_1D_avg(i_k)
       write(ifile4,*) ""
       write(ifile4,*) "# C(t) for wavevector" , kmag_1D_avg(i_k)
       ! cation C(t)
       do i_t=1,nmax_tcf
          write(ifile1,'(I8, E20.6)') i_t , SQ_Ct_avg(i_k,1,i_t) 
       enddo
       ! anion C(t)
       do i_t=1,nmax_tcf
          write(ifile2,'(I8, E20.6)') i_t , SQ_Ct_avg(i_k,2,i_t)
       enddo
       ! cation-anion C(t)
       do i_t=1,nmax_tcf
          write(ifile3,'(I8, E20.6)') i_t , SQ_Ct_avg(i_k,3,i_t)
       enddo
       ! total C(t)
       do i_t=1,nmax_tcf
          Sq_print = SQ_Ct_avg(i_k,1,i_t) + SQ_Ct_avg(i_k,2,i_t) + 2 * SQ_Ct_avg(i_k,3,i_t)
          write(ifile4,'(I8, E20.6)') i_t , Sq_print
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

  end subroutine Sq_TCF_print


  subroutine  combine_partial_structure_factors( SQ2 , SQ_store )
    use global_variables
    real*8,dimension(:,:,:,:,:),intent(inout) :: SQ2
    complex*16,dimension(:,:,:,:),intent(in) :: SQ_store

    integer :: i_type, j_type

    do i_type=1,n_atom_type
       do j_type = i_type, n_atom_type

          if ( i_type == j_type ) then
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + dble(SQ_store(i_type,:,:,:))**2+aimag(SQ_store(i_type,:,:,:))**2
          else
             ! here we have unlike types, so we want SQi(c.c.) * SQj + SQi * SQj(c.c.) , where
             ! (c.c.) is complex conjugate, which equals 2 * real ( SQi * SQj(c.c.) )
             !SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + 2d0 * dble( SQ_store(i_type,:,:,:) * conjg((SQ_store(j_type,:,:,:))) )
             ! multiply by 0.5 for double counting
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + 0.5d0 *  2d0 * dble( SQ_store(i_type,:,:,:) * conjg((SQ_store(j_type,:,:,:))) )
          end if
       enddo
    enddo

  end subroutine combine_partial_structure_factors



  subroutine print_SQ( SQ2 , kk, K )
    use global_variables
    real*8,dimension(:,:,:),intent(in) :: SQ2
    real*8,dimension(3,3),intent(in) :: kk
    integer, intent(in) :: K

    integer :: i, j , l , n
    real*8, dimension(3) :: k_vec

    do i=1, K
       do j=1,K
          do l=1,K
             ! convert wavevector, note the reciprocal lattice vectors kk don't have the 2*pi factor
             k_vec(:) = 2 * pi * ( dble(i-1) * kk(1,:) +  dble(j-1) * kk(2,:) +  dble(l-1) * kk(3,:)  )
             Select Case(Charge_density_Sq)
             Case('yes')
                write(*,'(3F14.6, E20.6)') k_vec, SQ2(i,j,l)
             case default
                write(*,'(3F14.6, F20.6)') k_vec, SQ2(i,j,l)
             end select
          enddo
       enddo
    enddo

  end subroutine print_SQ


  ! *************************
  ! if q_mag is greater than the longest qvec for which we
  ! have the form factor, then return 0
  !*************************
  real*8 function get_form_fac( i_index, q_mag )
    use global_variables
    integer, intent(in) :: i_index
    real*8, intent(in) :: q_mag

    integer :: i, index, flag

    !    flag=0
    !    do i=1, max_q_form-1
    !       if ( ( q_grid_form(i) < q_mag ) .and. ( q_mag < q_grid_form(i+1)  ) ) then
    !          index = i
    !          flag=1
    !          exit
    !       endif
    !    enddo

    !    if ( flag == 0 ) then
    !       get_form_fac=0d0
!!$       write(*,*) "couldn't find atomic form factor for q value ", q_mag
!!$       stop
    !   else
    index = ceiling( q_mag / dq_form )
    if ( index < max_q_form ) then
       get_form_fac = atomtype_form( i_index, index )
    else
       get_form_fac = 0d0
    endif
    !   endif

  end function get_form_fac




end module structure_factor
