program compute_structure_factor
  use routines
  use global_variables
  use structure_factor
  use MKL_DFTI
  use pme_routines

  character(MAX_FN) :: ifile_conf, traj_file, charge_file, form_file_dir
  character(MAX_ANAME), dimension(MAX_N_ATOM) :: alist
  integer      :: n_atom
  real*8, dimension(MAX_N_ATOM,3) :: xyz
  real*8, dimension(3,3) :: box, kk
  TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc,dfti_desc_inv

  integer :: i_atom, i_type
  !************************************
  ! IMPORTANT: this determines whether we're computing
  ! the CHARGE DENSITY structure factor (PME), or
  ! ELECTRON DENSITY structure factor (X-Ray) scattering.
  ! the first is used to evaluate electrolyte sum-rules, and
  ! requires input charges, while the second requires input form-factors
  !*************************************
  Select Case(Charge_density_Sq)
  Case("yes")
      write(*,*) "Computing S(q) of TOTAL charge density.  This requires input"
      write(*,*) "charges, and can be used to evaluate sum rules  "
  Case default
      write(*,*) "Computing S(q) of electron density, using input form factors"
      write(*,*) "this is the observable in X-Ray diffraction "   
  End Select

  ! charge_file is only present if we're computing charge density structure factor
  ! form_file_dir is directory with scattering form functions are located, and
  ! is only input if we're computing number density structure factor

  call sort_input_files( ifile_conf, traj_file, charge_file , form_file_dir )
  call read_gro( ifile_conf, n_atom, xyz, alist, box )

  Select Case(Charge_density_Sq)
  Case("yes")
     allocate(charges(n_atom))
     call read_charges( charge_file, charges , n_atom ) 
  Case default
     call create_atom_index( n_atom, alist )
     call get_atomic_form_factor( form_file_dir )
  End Select

  ! initialize bspline interpolation and FFT
  call initialize_spline_FFT(dfti_desc,dfti_desc_inv)

  ! compute structure factors
  call generate_structure_factor( n_atom, dfti_desc,dfti_desc_inv, traj_file )


end program compute_structure_factor
