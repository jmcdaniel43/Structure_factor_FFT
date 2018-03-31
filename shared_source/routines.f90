module routines

contains


  !*******************************************
  ! this subroutine decides which command line arguments
  ! correspond to which input files, depending on parameter settings
  !******************************************
  subroutine sort_input_files( ifile_conf, traj_file, charge_file, form_file_dir )
    use global_variables
    character(*), intent(out) :: ifile_conf, traj_file, charge_file, form_file_dir
    character(5) :: temp
    integer :: n_files

    n_files = IARGC ()

    ! first file is always .gro file
    call getarg( 1, ifile_conf )

    ! See if we're doing charge density or electron density structure factor
    Select Case(Charge_density_Sq)
    Case('yes')
       call getarg( 2, traj_file )
       call getarg( 3, charge_file )
       ! cation/anion S(q) decomposition
       call getarg( 4, cation_name)
       call getarg( 5, anion_name)
       call trim_head( cation_name )
       call trim_head( anion_name )

       write(*,*) "Identifying  ", cation_name , " as cations"
       write(*,*) "Identifying  ", anion_name ,  " as anions"
       write(*,*) "Will decompose S(q) into separate ion contributions"
    Case default
       ! this is number density structure factor for X-Ray scattering
       ! here trajectory file should be present
       call getarg( 2, traj_file )
       call getarg( 3, form_file_dir )
       ! if an additional argument, compute partial structure factor for requested atom index
       if ( n_files == 4 ) then
          partial_structure_factor="yes"
          call getarg(4, temp )
          read(temp,'(I5)') partial_structure_factor_index
       else
          partial_structure_factor="no"
       end if

    End Select


  end subroutine sort_input_files




  !***********************************************************************
  ! this subroutine reads a .gro file and stores the atomic names and coordinates
  !***********************************************************************
  subroutine read_gro( ifn, n_atom, xyz, alist, box )
    use global_variables
    implicit none
    character(MAX_FN), intent(in) :: ifn
    character(MAX_ANAME), intent(out), dimension(MAX_N_ATOM) :: alist
    integer, intent(out) :: n_atom
    real*8, intent(out), dimension(MAX_N_ATOM,3) :: xyz
    real*8, intent(out), dimension(3,3) :: box
    integer :: ifile, i_atom, i_mole, index, ion_type, nargs, inputstatus
    character(40),dimension(9)  :: args
    character(MAX_ANAME) :: aname, mname
    character(200)::line
    real*8 :: vel(3), tmp1, tmp2, tmp3
    integer :: ios

    ifile=99
    open( ifile, file=ifn, status='old' )
    read( ifile, * ) line
    read( ifile, '(I)' ) n_atom

    ! allocate atomic_number array if needed. this is global variable
    Select Case(Charge_density_Sq)
    Case('yes')
        allocate(atomic_number(n_atom))
    end Select

    do i_atom = 1 , n_atom
       if ( i_atom > MAX_N_ATOM ) then
          stop " please increase setting of MAX_N_ATOM"
       endif
       ! here we assume trjconv was used, which makes a gro file with higher precision F9.4 coordinates
       !read( ifile, '(I5,2A5,I5,3F9.4)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3)
       read( ifile, '(I5,2A5,I5,3F8.3)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3)
       ! this format may also be possible
       !         read( ifile, '(I5,2A5,I5,3F10.5)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3)

       Select Case(Charge_density_Sq)
       Case('yes')
          call trim_head( mname )
          call match_ion_type( mname , ion_type )
          atom_index( i_atom ) = ion_type
          n_atom_type = 2  ! either a cation atom, or anion atom
          call trim_head( aname )
          call fill_atomic_number( i_atom, aname )
       End Select
       ! convert nm to angstoms
       xyz(i_atom,:) = xyz(i_atom,:) * 10d0
       call trim_end( aname )
       alist( i_atom ) = aname
    end do

    box=0d0

    ! now box
    Read(ifile,'(A)',Iostat=inputstatus) line
    ! if its orthogonal, 3 arguments, if not 9 arguments
    call parse(line," ",args,nargs)

    Select Case(nargs)
    Case(3)
       read(args(1),*) box(1,1)
       read(args(2),*) box(2,2)
       read(args(3),*) box(3,3)
    Case(9)
       read(args(1),*) box(1,1)
       read(args(2),*) box(2,2)
       read(args(3),*) box(3,3)
       read(args(4),*) box(1,2)
       read(args(5),*) box(1,3)
       read(args(6),*) box(2,1)
       read(args(7),*) box(2,3)
       read(args(8),*) box(3,1)
       read(args(9),*) box(3,2)
    case default
       stop "error reading box in read_trajectory_snapshot subroutine"
    End Select


    ! convert to angstroms
    box(:,:) = box(:,:) * 10d0


    close( ifile )
  end subroutine read_gro




  !***********************************
  ! returns "1" for cations, "2" for anions
  subroutine match_ion_type( mname , ion_type )
    use global_variables    
    character(*),intent(in) :: mname
    integer, intent(out)  :: ion_type

    if ( mname .eq. cation_name) then
       ion_type = 1
    elseif( mname .eq. anion_name) then
       ion_type = 2
    else
       write(*,*) "can't recognize ion type ", mname
       stop
    end if

  end subroutine match_ion_type





  !********************************************
  ! this just reads a list of charges for all atoms in the simulation
  ! we don't use atomtypes here.
  !********************************************
  subroutine read_charges( charge_file, charges , n_atom )
    character(*), intent(in) :: charge_file
    real*8, dimension(:), intent(out) :: charges     
    integer, intent(in)      :: n_atom
    integer  :: ifile, i_atom

    ifile=99
    open( ifile, file=charge_file, status='old' )

    do i_atom = 1 , n_atom
       read( ifile, * ) charges(i_atom)
    enddo

    close( ifile )
  end subroutine read_charges


  !*******************************************
  ! this subroutine reads the number of trajectories from a .gro file
  !*******************************************
  subroutine get_n_trajectories( n_traj, dt, n_atom, traj_file )
    use global_variables
    integer, intent(out) :: n_traj
    real*8, intent(out)  :: dt
    character(*), intent(in) :: traj_file
    integer, intent(in) :: n_atom

    character(200) :: line
    integer :: ifile=99, inputstatus, i_line, i_atom, nargs
    real*8  :: t0, t1
    character(80),dimension(10)  :: args
    character(1) :: delim='='
    write(*,*) " figuring out number of snapshots from trajectory ...."
    open( ifile, file=traj_file, status='old' )

    n_traj = 0
    do
       Read(ifile,'(A)',Iostat=inputstatus) line
       ! get dt from 1st and 2nd frames
       if ( n_traj == 0 ) then
          ! parse on '=', could be different number of words here...
          call parse(line,delim,args,nargs)
          ! make sure we're reading this correctly
          if ( nargs /= 2 ) then
             write(*,*) "can't read timestamp in 1st frame of trajectory"
             stop
          endif
          read(args(2),*) t0 
       elseif ( n_traj == 1 ) then
          call parse(line,delim,args,nargs) 
          read(args(2),*) t1     
       endif
       If(inputstatus < 0) Exit
       n_traj = n_traj + 1
       Read(ifile,'(A)',Iostat=inputstatus) line
       do i_atom=1,n_atom
          Read(ifile,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
       enddo
       Read(ifile,'(A)',Iostat=inputstatus) line
    enddo

    dt = t1-t0
    write(*,*) " number of snapshots : ", n_traj
    write(*,*) " dt between snapshots ", dt , " ps"
    close( ifile )

  end subroutine get_n_trajectories



  !*******************************************
  ! this subroutine reads xyz coordinates and box vectors for a
  ! snapshot from a .gro file
  !
  ! here we read n_strides snapshots at once, and store in
  ! memory.  This prevents us from too much disk i/o 
  !*******************************************
  subroutine read_trajectory_snapshot( ifile, xyz, box, n_atom )
    use global_variables
    integer,intent(in) :: ifile
    integer, intent(in) :: n_atom
    real*8,dimension(:,:),intent(out) :: box
    real*8,dimension(:,:),intent(out) :: xyz
    character(20) :: junk1, junk2, junk3, junk4

    character(200) :: line
    integer :: inputstatus, i_line, i_atom, i_mole, index, nargs
    character(40),dimension(9)  :: args
    character(5) :: mname, aname
    !**************** these data structures are for reading in n_strides of data
    integer           :: i_read
    integer,parameter :: n_stride=20
    integer, save     :: i_stride=1, initialize=0
    real*8,dimension(:,:,:), allocatable,save :: box_stride, xyz_stride

    ! if initialize, allocate data structures
    if ( initialize == 0 ) then
       allocate(box_stride(size(box(:,1)),size(box(1,:)),n_stride))
       allocate(xyz_stride(size(xyz(:,1)),size(xyz(1,:)),n_stride))
       initialize=1
    end if


    ! if i_stride ==1, that means we need to read n_stride new trajectories

    if ( i_stride == 1 ) then
       !zero box in case orthorhombic
       box_stride(:,:,:)=0d0

       do i_read = 1 , n_stride

          Read(ifile,'(A)',Iostat=inputstatus) line
          if ( inputstatus < 0 ) then
             ! end of file.  This is probably ok, as we don't have an integer *
             ! n_stride number of trajectories and so this will naturally happen
             exit
          end if

          Read(ifile,'(A)',Iostat=inputstatus) line
          if ( inputstatus < 0 ) then
             write(*,*) "error reading trajectory file"
             stop
          endif

          do i_atom = 1 , n_atom
             read( ifile, '(I5,2A5,I5,3F8.3)' ), i_mole, mname, aname, index, xyz_stride(i_atom,1,i_read), xyz_stride(i_atom,2,i_read), xyz_stride(i_atom,3,i_read)
          enddo

          ! now box
          Read(ifile,'(A)',Iostat=inputstatus) line
          ! if its orthogonal, 3 arguments, if not 9 arguments
          call parse(line," ",args,nargs)

          Select Case(nargs)
          Case(3)
             read(args(1),*) box_stride(1,1,i_read)
             read(args(2),*) box_stride(2,2,i_read)
             read(args(3),*) box_stride(3,3,i_read)
          Case(9)
             read(args(1),*) box_stride(1,1,i_read)
             read(args(2),*) box_stride(2,2,i_read)
             read(args(3),*) box_stride(3,3,i_read)
             read(args(4),*) box_stride(1,2,i_read)
             read(args(5),*) box_stride(1,3,i_read)
             read(args(6),*) box_stride(2,1,i_read)
             read(args(7),*) box_stride(2,3,i_read)
             read(args(8),*) box_stride(3,1,i_read)
             read(args(9),*) box_stride(3,2,i_read)
          case default
             stop "error reading box in read_trajectory_snapshot subroutine"
          End Select

       enddo

       ! convert to angstroms
       xyz_stride = xyz_stride * 10d0
       box_stride = box_stride * 10d0

       ! now assign to output data structures
       xyz(:,:) = xyz_stride(:,:,i_stride)
       box(:,:) = box_stride(:,:,i_stride)


    else

       ! *****************************************
       ! here we already have read in this data, copy it 
       ! from memory
       ! *******************************************
       xyz(:,:) = xyz_stride(:,:,i_stride)
       box(:,:) = box_stride(:,:,i_stride)

    endif

    !************ advance the strid counter
    i_stride = i_stride + 1
    if ( i_stride > n_stride ) then
       i_stride = 1
    endif


  end subroutine read_trajectory_snapshot




  subroutine create_atom_index( n_atom, alist )
    use global_variables
    character(MAX_ANAME), intent(in), dimension(MAX_N_ATOM) :: alist
    integer, intent(in) :: n_atom   

    integer :: i_atom, i_type, flag_new

    ! get first atom type
    n_atom_type = 1

    atom_index(1) = 1
    atype_name(1) = alist(1)

    do i_atom = 2, n_atom
       ! see if this a new atom type
       flag_new=1
       do i_type = 1, n_atom_type
          if ( alist( i_atom ) .eq. atype_name(i_type ) ) then
             flag_new=0
             atom_index(i_atom ) = i_type
             exit
          end if
       enddo

       ! if new atom type
       if ( flag_new == 1 ) then
          n_atom_type = n_atom_type + 1
          atype_name(n_atom_type) = alist( i_atom )
          atom_index( i_atom ) = n_atom_type
       end if
    enddo

  end subroutine create_atom_index



  !*********************************
  ! this subroutine gets the atomic form factors
  ! for every atomtype in the system
  !*********************************

  subroutine get_atomic_form_factor( form_file_dir )
    use global_variables
    character(*), intent(inout) :: form_file_dir  

    character(MAX_FN) :: form_file
    character(4) :: prefix, suffix
    character(MAX_ANAME) :: aname

    integer :: i, i_type, ifile=99

    prefix="AFF_"
    suffix=".out"

    ! this is directory with form factor files
    call trim_head(form_file_dir)

    ! loop over atom types
    do i_type = 1, n_atom_type
       aname = atype_name(i_type)
       call trim_head( aname )   

       form_file = trim(form_file_dir)//prefix//trim(aname)//suffix

       write(*,*) "getting form factor for atom type ", aname , "  ....."

       open( ifile, file=form_file, status='old' )

       do i=1, max_q_form
          read(ifile,*) q_grid_form(i) , atomtype_form(i_type, i )

          ! convert from nm^-1, to angstrom^-1
          q_grid_form(i) = q_grid_form(i) / 10d0

       enddo
       dq_form = q_grid_form(2) - q_grid_form(1)

       close(ifile)
    enddo


  end subroutine get_atomic_form_factor


  !****************************************
  ! this subroutine fills in atomic numbers to be
  ! used in number structure factor calculation
  ! this is based on atom name in gromacs, so 
  ! it's imperfect, as atom names don't always reflect elements...
  !****************************************
  subroutine fill_atomic_number( i_atom, aname )
    use global_variables
    integer, intent(in) :: i_atom
    character(*), intent(in) :: aname    

    ! guess element based on atom name
    ! two letter dictionaries first...
    if ( ( aname(1:2) .eq. "Sh" ) .or. ( aname(1:2) .eq. "SH" ) ) then
       atomic_number(i_atom) = 0.0  ! drude oscillator
    elseif ( ( aname(1:2) .eq. "Cl" ) .or. ( aname(1:2) .eq. "CL" ) ) then
       atomic_number(i_atom) = 17.0
    elseif ( ( aname(1:2) .eq. "Na" ) .or. ( aname(1:2) .eq. "NA" ) ) then
       atomic_number(i_atom) = 11.0
    ! now one letter if couldn't match two letters
    elseif ( aname(1:1) .eq. "S" ) then
       atomic_number(i_atom) = 16.0
    elseif ( aname(1:1) .eq. "P" ) then
       atomic_number(i_atom) = 15.0
    elseif ( aname(1:1) .eq. "F" ) then
       atomic_number(i_atom) = 9.0
    elseif ( aname(1:1) .eq. "O" ) then
       atomic_number(i_atom) = 8.0
    elseif ( aname(1:1) .eq. "N" ) then
       atomic_number(i_atom) = 7.0
    elseif ( aname(1:1) .eq. "C" ) then
       atomic_number(i_atom) = 6.0
    elseif ( aname(1:1) .eq. "B" ) then
       atomic_number(i_atom) = 5.0
    elseif ( aname(1:1) .eq. "H" ) then
       atomic_number(i_atom) = 1.0
    else
       write(*,*) "can't determine element for atomtype ", aname
       stop
    endif

  end subroutine fill_atomic_number



  !********************************************
  ! this subroutine creates direct coordinates, scaled
  ! by input integer "K" (pme_grid), using the
  ! reciprocal lattice vectors
  ! 
  ! note, only coordinates are returned for atoms of type "i_type"
  !********************************************
  subroutine create_scaled_direct_coordinates(i_type,xyz_scale, xyz, n_atom, n_atom_kind, kk, K,charge_iontype)
    use global_variables
    integer  , intent(in)  :: i_type
    real*8,dimension(:,:),intent(out) :: xyz_scale
    real*8,dimension(:,:),intent(in) :: xyz
    integer, intent(in) :: n_atom
    integer, intent(out) :: n_atom_kind
    real*8,dimension(:,:),intent(in) :: kk
    integer, intent(in) :: K
    real*8,dimension(:),intent(out),optional:: charge_iontype

    integer :: i_atom,j,l, index
    real*8,parameter :: small=1D-6

    n_atom_kind=0
    do j=1,n_atom
       index = atom_index(j)
       ! if desired atom type
       if ( index == i_type ) then
          n_atom_kind = n_atom_kind + 1
          i_atom = n_atom_kind
          ! if outputting charge_iontype
          if( present(charge_iontype) ) then
             charge_iontype(i_atom) = charges(j)  ! charges is global array
          endif
          do l=1,3
             xyz_scale(i_atom,l)=dble(K)*dot_product(kk(l,:),xyz(j,:))
             ! if atoms are out of grid, shift them back in
             if (xyz_scale(i_atom,l)<0.) then
                xyz_scale(i_atom,l)=xyz_scale(i_atom,l)+dble(K)
             else if(xyz_scale(i_atom,l)>= dble(K)) then
                xyz_scale(i_atom,l)=xyz_scale(i_atom,l)-dble(K)
             endif
             ! make sure scaled coordinates are not numerically equal to zero, otherwise this will screw up Q grid routine
             if ( abs(xyz_scale(i_atom,l)) < small ) then
                xyz_scale(i_atom,l) = small
             end if
          enddo
       end if
    enddo

  end subroutine create_scaled_direct_coordinates



  !******************************************
  ! reciprocal lattice vector.  This is essentially the same
  ! as initialize_non_orth_transform subroutine, but we keep both
  ! in for compatibility with older code
  !******************************************
  subroutine construct_reciprocal_lattice_vector(kk,vol,box)
    real*8,dimension(:,:),intent(out) :: kk
    real*8, intent(out)   :: vol
    real*8,dimension(:,:),intent(in) :: box

    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3)

    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

  end subroutine construct_reciprocal_lattice_vector



  subroutine crossproduct( a,b,ans )
    implicit none
    real*8, intent(in) :: a(3),b(3)
    real*8, intent(out) :: ans(3)
    ans(1) = a(2)*b(3)-a(3)*b(2)
    ans(2) = -a(1)*b(3)+a(3)*b(1)
    ans(3) = a(1)*b(2)-a(2)*b(1)
  end subroutine crossproduct

  real function volume( a, b, c )
    implicit none
    real*8, intent(in) :: a(3), b(3), c(3)
    volume = a(1) * (b(2)*c(3)-b(3)*c(2)) - a(2) * (b(1)*c(3)-b(3)*c(1)) + a(3) * (b(1)*c(2)-b(2)*c(1))
    volume = abs(volume)
  end function volume


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_end( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space=0
    flag = 0
    do i = n, 1, -1
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i+n_space:i+n_space) = aname(i:i)
       else
          aname(i+n_space:i+n_space) = aname(i:i)
       end if
    end do
    do i = 1, n_space
       aname(i:i) = ' '
    end do
  end subroutine trim_end


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_head( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space = 0
    flag = 0
    do i = 1, n
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i-n_space:i-n_space) = aname(i:i)
       else 
          aname(i-n_space:i-n_space) = aname(i:i)
       end if
    end do
    do i = n-n_space+1, n
       aname(i:i) = ' ' 
    end do
  end subroutine trim_head


  !**************************************************
  !   these are routines for sorting arrays in terms of
  !   increasing element value.  Here we sort 2 arrays
  !   simultaneously based on the first array values
  !**************************************************


  ! --------------------------------------------------------------------
  ! This simultaneously sorts x,y arrays based on corresponding elements of x
  ! right now we are inputing a 3D array for y, and are sorting all elements
  ! based on first dimension
  ! --------------------------------------------------------------------

  SUBROUTINE  Sort(x, y)
    IMPLICIT  NONE
    Real*8, DIMENSION(:), INTENT(INOUT) :: x
    real*8, dimension(:,:,:), intent(inout) :: y
    INTEGER                               :: i, max_step 
    INTEGER                               :: Location

    max_step = size(x)-1              ! except for the last
    DO i = 1, max_step
       Location = FindMinimum(x, i)     ! find min from this to last
       CALL  Swap(x(i), x(Location))  ! swap this and the minimum
       CALL  Swap2D(y(i,:,:), y(Location,:,:))
    END DO
  END SUBROUTINE  Sort



  ! --------------------------------------------------------------------
  ! INTEGER FUNCTION  FindMinimum():
  !    This function returns the location of the minimum in the section
  ! between Start and End.
  ! --------------------------------------------------------------------

  INTEGER FUNCTION  FindMinimum(x, Start)
    IMPLICIT  NONE
    Real*8, DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: Start
    Real*8                             :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i

    Minimum  = x(Start)               ! assume the first is the min
    Location = Start                  ! record its position
    DO i = Start+1, size(x)               ! start with next elements
       IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
          Minimum  = x(i)             !      Yes, a new minimum found
          Location = i                !      record its position
       END IF
    END DO
    FindMinimum = Location            ! return the position
  END FUNCTION  FindMinimum

  ! --------------------------------------------------------------------
  ! SUBROUTINE  Swap():
  !    This subroutine swaps the values of its two formal arguments.
  ! --------------------------------------------------------------------

  SUBROUTINE  Swap(a, b)
    IMPLICIT  NONE
    real*8 , INTENT(INOUT) :: a, b
    real*8                 :: Temp

    Temp = a
    a    = b
    b    = Temp
  END SUBROUTINE  Swap



  SUBROUTINE  Swap2D(a, b)
    IMPLICIT  NONE
    real*8 , dimension(:,:), INTENT(INOUT) :: a, b
    real*8 , dimension(size(a(:,1)),size(a(1,:))):: Temp

    Temp = a
    a    = b
    b    = Temp
  END SUBROUTINE  Swap2D



  !******************************************************
  !        Here are some string manipulation routines
  !        written by Dr. George Benthien and taken from
  !        http://www.gbenthien.net/strings/index.html
  !******************************************************

  subroutine parse(str,delims,args,nargs)

    ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    ! the delimiters contained in the string 'delims'. Preceding a delimiter in
    ! 'str' by a backslash (\) makes this particular instance not a delimiter.
    ! The integer output variable nargs contains the number of arguments found.

    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:) :: args
    integer,intent(out) ::nargs
    integer :: i,k,na,lenstr


    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
       args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    k=0

    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       call split(str,delims,args(nargs))
       call removebksl(args(nargs))
    end do
    str=strsav

  end subroutine parse



  subroutine compact(str)

    ! Converts multiple spaces and tabs to single spaces; deletes control
    ! characters;
    ! removes initial spaces.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str)):: outstr
    integer :: i,k,ich,isp,lenstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    isp=0
    k=0

    do i=1,lenstr
       ch=str(i:i)
       ich=iachar(ch)

       select case(ich)

       case(9,32)     ! space or tab character
          if(isp==0) then
             k=k+1
             outstr(k:k)=' '
          end if
          isp=1

        case(33:)      ! not a space, quote, or control character
           k=k+1
           outstr(k:k)=ch
           isp=0

       end select

    end do

    str=adjustl(outstr)

  end subroutine compact


  subroutine split(str,delims,before,sep)

    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the
    ! found delimiter. A delimiter in 'str' is treated like an ordinary
    ! character if it is preceded by a backslash (\). If the backslash
    ! character is desired in 'str', then precede it with another backslash.

    character(len=*) :: str,delims,before
    character,optional :: sep
    logical :: pres
    character :: ch,cha
    integer:: i,k,lenstr,ibsl,ipos,iposa

    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if(lenstr == 0) return        ! string str is empty
    k=0
    ibsl=0                        ! backslash initially inactive
    before=' '
    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          before(k:k)=ch
          ibsl=0
          cycle
       end if
       if(ch == '\') then          ! backslash with backslash inactive
          k=k+1
          before(k:k)=ch
          ibsl=1
          cycle
       end if
       ipos=index(delims,ch)
       if(ipos == 0) then          ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
       end if
       if(ch /= ' ') then          ! character is a delimiter that is not a  space
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
       cha=str(i+1:i+1)            ! character is a space delimiter
       iposa=index(delims,cha)
       if(iposa > 0) then          ! next character is a delimiter
          str=str(i+2:)
          if(pres) sep=cha
          exit
       else
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
    end do
    if(i >= lenstr) str=''
    str=adjustl(str)              ! remove initial spaces
    return

  end subroutine split

  !**********************************************************************

  subroutine removebksl(str)

    ! Removes backslash (\) characters. Double backslashes (\\) are replaced
    ! by a single backslash.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr
    integer :: i,k,ibsl,lenstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0
    ibsl=0                        ! backslash initially inactive

    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          outstr(k:k)=ch
          ibsl=0
          cycle
       end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************



end module routines
