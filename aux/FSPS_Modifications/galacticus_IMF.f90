module Galacticus_IMF_Module
  USE sps_vars
  private
  public :: Get_Galacticus_IMF, Set_Galacticus_IMF_Limits

  logical :: moduleInitialized=.false.

  integer :: nTable
  real(SP), allocatable, dimension(:) :: massTable,imfTable,y2

  real(SP) :: massLow,massHigh

contains

  subroutine Read_Galacticus_IMF()
    implicit none
    integer :: i,ierr,iunit
    real(SP) :: dy
    
    if (.not.moduleInitialized) then
       ! Read the file here.
       iunit=448
       nTable=0
       open(unit=iunit,file=trim(imfFileName),status="old",form="formatted",iostat=ierr)
       do while (ierr.eq.0)
          read (iunit,*,iostat=ierr)
          if (ierr.eq.0) nTable=nTable+1
       end do
       close(iunit)
       allocate(massTable(nTable))
       allocate(imfTable(nTable))
       allocate(y2(nTable+1))
       open(unit=iunit,file=trim(imfFileName),status="old",form="formatted",iostat=ierr)
       do i=1,nTable
          read (iunit,*,iostat=ierr) massTable(i),imfTable(i)
       end do
       close(iunit)
       massLow=massTable(1)
       massHigh=massTable(nTable)
       call gspline(massTable,imfTable,y2)
       moduleInitialized=.true.
    end if
    return
  end subroutine Read_Galacticus_IMF

  subroutine Set_Galacticus_IMF_Limits()
    implicit none
    
     call Read_Galacticus_IMF()
     imf_lower_limit=massLow
     imf_upper_limit=massHigh
     return
  end subroutine Set_Galacticus_IMF_Limits
  
  subroutine Get_Galacticus_IMF(mass,imf)
    implicit none
    real(SP), intent(in), dimension(:) :: mass
    real(SP), intent(out), dimension(:) :: imf
    integer :: i

    call Read_Galacticus_IMF()
    
    ! Do polynomial interpolation to get the IMF at the requested masses.
    do i=1,size(mass)
       if (mass(i).lt.massLow.or.mass(i).gt.massHigh) then
          imf(i)=0.0
       else
          imf(i)=gsplint(massTable,imfTable,y2,mass(i))
       end if
    end do

    return
  end subroutine Get_Galacticus_IMF

  subroutine gspline(x,y,z)
    implicit none
    real(SP), intent(in   ), dimension(0:       ) :: x,y
    real(SP), intent(inout), dimension(0:size(x)) :: z
    real(SP),                dimension(0:size(x)) :: h,b,u,v
    integer                                   :: i,n

    ! Compute number of points.
    n=size(x)
    ! Compute h and b.
    do i=0,n-2
       h(i)= x(i+1)-x(i)
       b(i)=(y(i+1)-y(i))/h(i)
    end do
    h(n-1)=h(n-2)
    b(n-1)=b(n-2)
    ! Gaussian elimination.
    u(1)=2.0*(h(0)+h(1))
    v(1)=6.0*(b(1)-b(0))
    do i=2,n-1
       u(i)=2.0*(h(i-1)+h(i))-h(i-1)**2/u(i-1)
       v(i)=6.0*(b(i)-b(i-1))-h(i-1)*v(i-1)/u(i-1)
    end do
    ! Back-substitution.
    z(n)=0.0
    do i=n-1,1,-1
       z(i)=(v(i)-h(i)*z(i+1))/u(i)
    end do
    z(0)=0.0
    return
  end subroutine gspline
  
  function gsplint(x,y,z,xx)
    implicit none
    real(SP)                            :: gsplint
    real(SP), intent(in), dimension(0:) :: x,y,z
    real(SP), intent(in)                :: xx
    real(SP)                            :: a,b,c,d,h
    integer                             :: i,n

    n=size(x)
    do i=0,n-1
       if (xx <= x(i+1)) then
          exit
       end if
    end do
    h=x(i+1)-x(i)
    a=y(i)
    b=-h*z(i+1)/6.0-h*z(i)/3.0+(y(i+1)-y(i))/h
    c=z(i)/2.0
    d=(z(i+1)-z(i))/6.0/h
    gsplint=a+(xx-x(i))*(b+(xx-x(i))*(c+(xx-x(i))*d))
    return
  end function gsplint

end module Galacticus_IMF_Module
