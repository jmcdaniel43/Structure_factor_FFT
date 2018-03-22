module pme_routines
  use routines
  implicit none

  !*******************************************************************
  ! PARTICLE MESH EWALD SUBROUTINES
  !
  ! This module contains PME subroutines for energy and force
  ! pme routines use discrete fourier transforms in MKL library
  !
  !  reference for the pme algorithm is
  !  Essmann et al , J. Chem. Phys. 1995, 103, 8577-8593
  !*******************************************************************


contains


  subroutine initialize_spline_FFT(dfti_desc,dfti_desc_inv)
    use MKL_DFTI
    use global_variables
    TYPE(DFTI_DESCRIPTOR), pointer,intent(out):: dfti_desc,dfti_desc_inv

    integer:: i, length(3), status

    ! allocate arrays
    allocate(B(pme_grid,pme_grid,pme_grid))

    ! set up fourier transform descriptors
    length=pme_grid

    status=DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status=DftiCommitDescriptor(dfti_desc)
    ! don't scale back transform because pick up factor of K^3 from convolution
    status=DftiCreateDescriptor(dfti_desc_inv, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status = DftiCommitDescriptor(dfti_desc_inv)

    ! compute B array
    call B_array(B, pme_grid,spline_order)

!!!!!!!!!!!!!!!!grid B_splines
    if (spline_order .eq. 6) then
       do i=1,spline_grid
          B6_spline(i)=B_spline(6./dble(spline_grid)*dble(i),6)
          B5_spline(i)=B_spline(5./dble(spline_grid)*dble(i),5)
       enddo
    else if (spline_order .eq. 4) then
       do i=1,spline_grid
          B4_spline(i)=B_spline(4./dble(spline_grid)*dble(i),4)
          B3_spline(i)=B_spline(3./dble(spline_grid)*dble(i),3)
       enddo
    else
       stop "requested spline order not implemented"
    endif


  end subroutine initialize_spline_FFT




  !*****************************************************************
  ! This subroutine interpolates charges onto Q grid to be used in pme reciprocal space
  ! routines. 
  !
  ! Qn grid is used to grid number density for Sn(q), and Qc grid is used to grid
  ! charge density for Sc(q).
  !
  ! If 'charge_iontype' is present,
  ! then we are creating two Qgrids, one with charge density, and
  ! one with electron number density (X-ray structure)
  !*****************************************************************
  subroutine grid_Q(Qn,Qc, xyz,n_atom,K,n,charge_iontype)
    use global_variables
    use omp_lib
    integer, intent(in) :: n_atom
    real*8, intent(in), dimension(:,:) :: xyz
    integer,intent(in)::K,n
    real*8,dimension(:),intent(in),optional::charge_iontype
    real*8,dimension(:,:,:),intent(out)::Qn,Qc
    integer::i,j,k1,k2,k3,n1,n2,n3,nn1,nn2,nn3,nearpt(3),splindex(3)
    real*8::sum1,sum2
    real*8,dimension(3)::u,arg


    Qn=0D0
    Qc=0d0
 
    ! parameter spline_grid undeclared, but ok
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_atom,xyz,charge_iontype,atomic_number,n,B6_spline,B4_spline,K) REDUCTION(+:Qn,Qc)
    !$OMP DO SCHEDULE(dynamic, n_threads)
    !$
    do j=1,n_atom
       u=xyz(j,:)
       nearpt=floor(u)
       ! only need to go to k=0,n-1, for k=n, arg > n, so don't consider this
       do k1=0,n-1
          n1=nearpt(1)-k1
          arg(1)=u(1)-dble(n1);
          ! shift index of array storage if < 0
          if(n1<0) then
             n1=n1+K
          endif
          do k2=0,n-1
             n2=nearpt(2)-k2
             arg(2)=u(2)-dble(n2)
             ! shift index of array storage if < 0
             if(n2<0) then
                n2=n2+K
             endif
             do k3=0,n-1
                n3=nearpt(3)-k3
                arg(3)=u(3)-dble(n3);
                ! shift index of array storage if < 0
                if(n3<0) then
                   n3=n3+K
                endif

                sum1=0d0
                sum2=0d0
                splindex = ceiling(arg/6.D0*dble(spline_grid))
               
                ! note 0<arg<n , so arg should always be within bounds of gridded spline
                ! use Case statement to see if we want total charge density
                if(spline_order .eq.6) then
                   Select Case(Charge_density_Sq)
                   Case('yes')
                       ! number density (sum1) and total charge (sum2) structure factors
                       ! Qn, we multiply by atomic_number, as we are using a
                       ! general approximation for atomic form factors to save time
                       sum1=atomic_number(j)*B6_spline(splindex(1))*B6_spline(splindex(2))*B6_spline(splindex(3))
                       sum2=charge_iontype(j)*B6_spline(splindex(1))*B6_spline(splindex(2))*B6_spline(splindex(3))
                   Case default
                       ! we will mutiply S(q) by form factor F(q), so we don't multiply by atomic_number here
                       sum1=B6_spline(splindex(1))*B6_spline(splindex(2))*B6_spline(splindex(3))
                   End Select
                else
                   Select Case(Charge_density_Sq)
                   Case('yes')
                       ! number density (sum1) and total charge (sum2) structure factors
                       ! Qn, we multiply by atomic_number, as we are using a
                       ! general approximation for atomic form factors to save time
                       sum1=atomic_number(j)*B4_spline(splindex(1))*B4_spline(splindex(2))*B4_spline(splindex(3))   
                       sum2=charge_iontype(j)*B4_spline(splindex(1))*B4_spline(splindex(2))*B4_spline(splindex(3))
                   Case default
                       ! we will mutiply S(q) by form factor F(q), so we don't multiply by atomic_number here
                       sum1=B4_spline(splindex(1))*B4_spline(splindex(2))*B4_spline(splindex(3))
                   End Select
                endif

                Qn(n1+1,n2+1,n3+1)=Qn(n1+1,n2+1,n3+1)+sum1
                Qc(n1+1,n2+1,n3+1)=Qc(n1+1,n2+1,n3+1)+sum2

             enddo
          enddo
       enddo
    enddo

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine grid_Q





  !***********************************************
  ! this function calculates B_splines which are used in pme as interpolating
  ! functions.  B_splines are calculated recursively, and therefore it's a good idea
  ! to grid them
  !************************************************
  real*8 function B_spline(u,n)
    real*8,intent(in)::u
    integer,intent(in)::n
    integer::i,j
    real,dimension(n-1,n-1)::mn
    real*8::ui

!!!!!!!! define m2 for n-1 values
    do i=1,n-1
       ui=u-dble(i-1)
       if((ui<0.).or.(ui>2.)) then
          mn(1,i)=0.D0
       else
          mn(1,i)=1.D0-abs(ui-1.D0)
       endif
    enddo

!!!!!!!!!!! define mj recursively for n-1-(j-1) values

    do j=2,n-1
       do i=1,n-j
          ui=u-dble(i-1)
          mn(j,i)=(ui/dble(j))*mn(j-1,i)+((dble(j+1)-ui)/dble(j))*mn(j-1,i+1)
       enddo
    enddo

    B_spline=mn(n-1,1)

  end function B_spline


  !***********************************************************
  ! This is a routine for reciprocal space pme calculation
  !***********************************************************
  subroutine B_array(B,K,n)
    complex*16,dimension(:,:,:),intent(out)::B
    integer,intent(in)::K,n
    real*8, parameter :: pi=3.14159265
    integer::i,j,l,m1,m2,m3
    do i=0,K-1
       do j=0,K-1
          do l=0,K-1
             B(i+1,j+1,l+1)=bm(i,n,K)*bm(j,n,K)*bm(l,n,K)
          enddo
       enddo
    enddo

  end subroutine B_array


  !******************************************************
  ! this is needed in reciprocal space pme calculation
  !******************************************************
  function bm(m,n,K)
    use global_variables
    complex*16::bm
    integer,intent(in)::m,n,K
    integer::i
    complex*16::sum
    real*8::tmp

    sum=0.D0
    do i=0,n-2
       tmp=2.D0*pi*dble(m*i)/dble(K)
       sum=sum+B_spline(dble(i+1),n)*cmplx(cos(tmp),sin(tmp))
!!$     sum=sum+B6_spline(dble(i+1)/6.*dble(spline_grid))*cmplx(cos(tmp),sin(tmp))
    enddo

    tmp=2d0*pi*dble(n-1)*dble(m)/dble(K)
    bm= cmplx(cos(tmp),sin(tmp))/sum

  end function bm




end module pme_routines
