!  DCD Fortran Interface DCD Example Program with Wrapper
!  2017 (c) Chang Yun Son <cson@chem.wisc.edu>
!
 
module dcd
  implicit none
  private
  public :: openmmdcdfile

  integer :: fp = 10
  real    :: prec = 1e-4

  type openmmdcdfile
    integer :: STAT, STEP, NATOMS, NFRAMES, INTERVAL, FIRSTSTEP
    real :: time, prec, dt, dtframe
    logical           :: bBoxFlag
    real*8,dimension(3,2) :: box
    real,dimension(:,:),allocatable :: pos
    
    contains
      procedure :: init => dcd_init
      procedure :: read => dcd_read
      procedure :: close => dcd_close
  end type openmmdcdfile
contains
  subroutine dcd_init(this, fname)
    class(openmmdcdfile), intent(inout) :: this
    CHARACTER(len=256), intent(in) :: fname
    integer           :: istat,iBoxFlag
    integer           :: dummyi,i,testi
    integer           :: dummyiar(8)
    real              :: dummyr
    character*4       :: dummyc
    character*80      :: dummystr,dummystr2

    write(*,*) 'start reading input file ',trim(fname)
    open(fp,file=trim(fname),status='old',form='unformatted')
    read(fp) dummyc, this % NFRAMES, this % FIRSTSTEP, this % INTERVAL, (dummyi,i=1,6), this % dt, &
            iBoxFlag, (dummyiar(i),i=1,8),testi

    ! convert the dcd time unit to the regular fs unit
    this % dt = this % dt * 0.04888821

    this % STEP = this % FIRSTSTEP
    this % dtframe = this % INTERVAL * this % dt
    this % time = this % FIRSTSTEP * this % dt

    write(*,*) 'dummyc NFRAMES FIRSTSTEP STEP INTERVAL dt',dummyc, this % NFRAMES, &
                this % FIRSTSTEP, this % STEP, this % INTERVAL, this % dt
    !read(fp) iBoxFlag, (dummyi,i=1,8),testi
    !read(fp) iBoxFlag, (dummyiar(i),i=1,8),testi
    write(*,*) 'iBoxFlag',iBoxFlag, (dummyiar(i),i=1,8),testi
    read(fp) dummyi,dummystr,dummystr2
    write(*,*) 'dummystr',dummyi,dummystr,dummystr2
    read(fp) this % NATOMS
    write(*,*) 'NATOMS', this % NATOMS
    
    if(iBoxFlag .ne. 0) this % bBoxFlag = .True.
    if(.not. allocated(this % pos)) then
        allocate(this % pos(3,this % NATOMS))
        write(*,*) 'position array is allocated '
    endif
  end subroutine dcd_init

  subroutine dcd_read(this)
    class(openmmdcdfile), intent(inout) :: this
    integer           :: istat
    integer           :: dummyi, ilength
    real              :: dummyr
    integer           :: natom,i,j

    natom = this % NATOMS
    this % STEP = this % STEP + this % INTERVAL
    this % time = this % time +  this % dtframe

    read(fp,IOSTAT=istat) (this % box(1,j), j=1,2),(this % box(2,j),j=1,2), & 
            (this % box(3,j),j=2,1,-1)
    !write(*,*) 'Box length',(this % box(i,1), i=1,3),(this % box(i,2),i=1,3)
    read(fp,IOSTAT=istat) (this % pos(1,j),j=1,natom)
    read(fp,IOSTAT=istat) (this % pos(2,j),j=1,natom)
    read(fp,IOSTAT=istat) (this % pos(3,j),j=1,natom)
    this % STAT = istat

  end subroutine dcd_read

  subroutine dcd_close(this)
    class(openmmdcdfile), intent(inout) :: this
    close(fp)
  end subroutine dcd_close

end module dcd
