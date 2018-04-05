module structure_factor

contains

  !*******************************
  ! this subroutine computes either number density or charge density structure factors
  ! based on the input and parameter settings.
  ! PME routines are setup to return both number density Qn, and charge density Qc,
  ! interpolated grids to save time.  We may or may not be using both
  ! structure factors will be returned in 3D , S(qvec), which is the most
  ! general data for non isotropic systems
  !*******************************
  subroutine generate_structure_factor( n_atom, dfti_desc,dfti_desc_inv, traj_file )
    use global_variables
    use routines
    use pme_routines
    use MKL_DFTI
    integer, intent(in) :: n_atom
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    character(*), intent(in) :: traj_file

    complex*16,dimension(:,:,:),allocatable::FQ, SQ
    real*8, dimension(:,:,:),allocatable :: Qn, Qc, SQ2
    real*8, dimension(:,:,:,:,:),allocatable :: SQ2_store
    complex*16, dimension(:,:,:,:),allocatable :: SQ_store
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8,dimension(3,3) :: kk, kk_avg, box
    real*8  ::vol, dt
    integer :: n, K
    real*8,dimension(:), allocatable :: charge_iontype, atomic_number_iontype
    real*8,dimension(:,:),allocatable :: xyz, xyz_scale
    integer :: n_atom_kind, i_type, status, n_traj, i_step, ifile=99, i_atom

    n=spline_order
    K=pme_grid
    ! here we average wavevectors as these change slightly during NPT
    kk_avg=0d0    

    allocate( Qn(pme_grid,pme_grid,pme_grid),Qc(pme_grid,pme_grid,pme_grid),FQ(pme_grid,pme_grid,pme_grid),SQ(pme_grid,pme_grid,pme_grid),SQ2(pme_grid,pme_grid,pme_grid), q_1r(pme_grid**3), q_1d(pme_grid**3) )
    SQ2=0d0


    ! we need these arrays if doing a full scattering structure factor and also
    ! if doing an electron density structure factor.  The only time we don't
    ! need the arrays is a partial structure factors, so might as well allocate
    ! them no matter what
    allocate(  SQ_store(n_atom_type,pme_grid,pme_grid,pme_grid), SQ2_store(n_atom_type,n_atom_type,pme_grid,pme_grid,pme_grid) )
    SQ_store=0d0
    SQ2_store=0d0

    ! get number of trajectory snapshots
    call get_n_trajectories( n_traj, dt, n_atom, traj_file )
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


       !*************** Decide whether we're doing total charge structure factor
       ! or X-ray scattering, electron density structure factor
       Select Case(Charge_density_Sq)
       Case('yes')
          do i_type = 1 , n_atom_type  ! n_atom_type=2, cation atoms=1, anion atoms=2
             ! note this only creates coordinates for atomtype "i_type"
             ! charge_iontype stores charges for atoms in xyz_scale array,
             ! with corresponding indices
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K, charge_iontype,atomic_number_iontype)
             ! returns both Qn number density and Qc charge density grids.
             ! here, we are only using Qc
             call grid_Q(Qn,Qc,xyz_scale,n_atom_kind,K,n, charge_iontype,atomic_number_iontype)
             q_1r=RESHAPE(Qc, (/K**3/) )
             q_1d=cmplx(q_1r,0.,16)
             status=DftiComputeForward(dfti_desc, q_1d)
             FQ=RESHAPE(q_1d, (/K,K,K/) )
             ! structure factor = B * FQ
             SQ_store(i_type,:,:,:) = FQ*B
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
                ! returns both Qn number density and Qc charge density grids.
                ! here, we are only using Qn
                ! input arrays charge_iontype and atomic_number_iontype are
                ! allocated but uninitialized.  This is fine, as they will not
                ! be used here, just token inputs
                call grid_Q(Qn,Qc,xyz_scale,n_atom_kind,K,n,charge_iontype,atomic_number_iontype)
                q_1r=RESHAPE(Qn, (/K**3/) )
                q_1d=cmplx(q_1r,0.,16)
                status=DftiComputeForward(dfti_desc, q_1d)
                FQ=RESHAPE(q_1d, (/K,K,K/) )
                ! structure factor = B * FQ
                SQ_store(i_type,:,:,:) = FQ*B
             enddo

             ! now create all the cross SQ2 structure factors for this snapshot
             call combine_partial_structure_factors( SQ2_store , SQ_store )

          Case("yes")
             !********************* computing partial structure factor
             i_type = partial_structure_factor_index
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
             ! returns both Qn number density and Qc charge density grids.
             ! here, we are only using Qn
             ! input arrays charge_iontype and atomic_number_iontype are
             ! allocated but uninitialized.  This is fine, as they will not
             ! be used here, just token inputs
             call grid_Q(Qn,Qc,xyz_scale,n_atom_kind,K,n,charge_iontype,atomic_number_iontype)
             q_1r=RESHAPE(Qn, (/K**3/) )
             q_1d=cmplx(q_1r,0.,16)
             status=DftiComputeForward(dfti_desc, q_1d)
             FQ=RESHAPE(q_1d, (/K,K,K/) )
             ! structure factor = B * FQ
             SQ = FQ*B
             SQ2 = SQ2 + dble(SQ)**2+aimag(SQ)**2
          End Select

       End Select

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


    ! now print, note because we took the FT in reduced coordinates, we need to
    ! convert to the physical wavevectors
    Select Case(Charge_density_Sq)
    Case('yes')
       write(*,*) " k vec ,  SQ cation * SQ cation "
       SQ2 = SQ2_store(1,1,:,:,:)
       call print_SQ( SQ2 , kk_avg, pme_max_print )   
       write(*,*) " k vec, SQ anion * SQ anion "
       SQ2 = SQ2_store(2,2,:,:,:)
       call print_SQ( SQ2 , kk_avg, pme_max_print )
       write(*,*) " k vec, SQ cation * SQ anion "
       SQ2 = SQ2_store(1,2,:,:,:)
       call print_SQ( SQ2 , kk_avg, pme_max_print )

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
    ! could be 0 if q_mag is zero.  We don't care about scattering at k=0
    ! anyway, so might as well set to zero
    if (index == 0 ) then
       get_form_fac = 0d0
    elseif ( index < max_q_form ) then
       get_form_fac = atomtype_form( i_index, index )
    else
       get_form_fac = 0d0
    endif
    !   endif

  end function get_form_fac




end module structure_factor
