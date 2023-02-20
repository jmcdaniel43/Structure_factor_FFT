!  openmmdcd Fortran Interface openmmdcd Example Program with Wrapper
!  2017 (c) Chang Yun Son <cson@chem.wisc.edu>
!
 
module traj
!  use xtc, only: xtcfile
  use dcd, only: openmmdcdfile

  implicit none
  private
  public :: trajfile

!  type(xtcfile) :: trajxtc
  type(openmmdcdfile) :: openmmdcd

  type trajfile
    integer :: STAT, STEP, NATOMS, NFRAMES, INTERVAL
    real :: time, prec, dt
    logical           :: bBoxFlag
    logical           :: btrajxtcType, bOpenmmDcdType
    real*8,dimension(3,3) :: box
    real,dimension(:,:),allocatable :: pos
    
    contains
      procedure :: init => traj_init
      procedure :: read => traj_read
      procedure :: count => traj_count
      procedure :: close => traj_close
  end type trajfile

contains
  subroutine traj_init(this, fname)
    class(trajfile), intent(inout) :: this
    CHARACTER(len=256), intent(in) :: fname
    CHARACTER(len=256)             :: filename
    CHARACTER(len=3)               :: extension
    integer                        :: ppos, isize

    this % btrajxtcType = .false.
    this % bOpenmmDcdType = .false.

    ppos = scan(trim(fname),".",BACK=.true.)
    if (ppos > 0) extension = trim(fname(ppos+1:))
    !if (extension .eq. 'xtc') then
    !  write(*,*) 'Input xtc file : ',trim(fname)
    !  this % btrajxtcType = .true.
    !  call trajxtc % init(fname)
    !  isize = size(trajxtc % pos(1,:))
    !  this % NATOMS = isize
    !  if(.not. allocated(this % pos)) then
    !    allocate(this % pos(3,isize))
    !  endif
    !else if (extension .eq. 'dcd') then
    if (extension .eq. 'dcd') then
      write(*,*) 'Input openmmdcd file : ',trim(fname)
      this % bOpenmmDcdType = .true.
      call openmmdcd % init(fname)
      isize = size(openmmdcd % pos(1,:))
      this % NATOMS = isize
      if(.not. allocated(this % pos)) then
        allocate(this % pos(3,isize))
      endif
      this % NFRAMES = openmmdcd % NFRAMES
    else
      write(*,*) 'Unknown file extension in file ',trim(fname)
      write(*,*) 'Only xtc or openmmdcd file formats are supported'
      stop
    endif
  end subroutine traj_init

  subroutine traj_read(this)
    class(trajfile), intent(inout) :: this
    real*8,dimension(3) :: boxlen, boxang
    real*8 :: sinGamma, cx, cy, cz

    !if(this % btrajxtcType) then
    !  call trajxtc % read
    !  this % NATOMS = trajxtc % NATOMS
    !  this % time   = trajxtc % time
    !  this % STEP   = trajxtc % STEP
    !  this % prec   = trajxtc % prec
    !  this % box    = trajxtc % box
    !  this % pos    = trajxtc % pos
    !  this % STAT   = trajxtc % STAT
    !else if (this % bOpenmmDcdType) then
    if (this % bOpenmmDcdType) then
      call openmmdcd % read
      this % NATOMS = openmmdcd % NATOMS
      this % time   = openmmdcd % time
      this % STEP   = openmmdcd % STEP
      this % prec   = openmmdcd % prec
      this % pos    = openmmdcd % pos / 10.d0
      this % STAT   = openmmdcd % STAT
      boxlen = openmmdcd % box(:,1) / 10.d0
      boxang = openmmdcd % box(:,2)
      ! need to change -> properly reshape triclinic box dimensions into box vectors
      sinGamma = dsin(dacos(boxang(1)))
      this % box(1,:) = (/ boxlen(1), 0.d0,0.d0 /)
      this % box(2,:) = (/ boxang(1), sinGamma,0.d0 /) * boxlen(2)
      cx = boxlen(3) * boxang(2)
      cy = boxlen(3) * ( boxang(3) - boxang(2)*boxang(1))/sinGamma
      cz = dsqrt(boxlen(3)*boxlen(3) + cx*cx + cy*cy)
      this % box(3,:) = (/ cx, cy, cz /)
    endif
  end subroutine traj_read

  subroutine traj_count(this, fname, icnt)
    class(trajfile), intent(inout) :: this
    CHARACTER(len=256), intent(in) :: fname
    CHARACTER(len=256)             :: filename
    integer, intent(out) :: icnt
    real*8,dimension(3) :: boxlen, boxang
    real*8 :: sinGamma, cx, cy, cz

    !if(this % btrajxtcType) then
    !  call trajxtc % count(fname, icnt)
      !icnt = this % NFRAMES
    !else if (this % bOpenmmDcdType) then
    if (this % bOpenmmDcdType) then
      icnt = this % NFRAMES
    endif
  end subroutine traj_count

  subroutine traj_close(this)
    class(trajfile), intent(inout) :: this

    !if(this % btrajxtcType) then
    !  call trajxtc % close
    !else if (this % bOpenmmDcdType) then
    if (this % bOpenmmDcdType) then
      call openmmdcd % close
    endif
  end subroutine traj_close

end module traj
