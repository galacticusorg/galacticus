module Galacticus_IMF_Module
  USE nrtype; use nr
  private
  public :: Get_Galacticus_IMF

  logical :: moduleInitialized=.false.

  integer :: nTable
  real(SP), allocatable, dimension(:) :: massTable,imfTable,y2

  real(SP) :: massLow,massHigh

contains

  subroutine Get_Galacticus_IMF(mass,imf)
    implicit none
    real(SP), intent(in), dimension(:) :: mass
    real(SP), intent(out), dimension(:) :: imf
    integer :: i,ierr,iunit
    real(SP) :: dy,lowLim,upLim

    if (.not.moduleInitialized) then
       ! Read the file here.
       iunit=448
       nTable=0
       open(unit=iunit,file="galacticus.imf",status="old",form="formatted",iostat=ierr)
       do while (ierr.eq.0)
          read (iunit,*,iostat=ierr)
          if (ierr.eq.0) nTable=nTable+1
       end do
       close(iunit)
       allocate(massTable(nTable))
       allocate(imfTable(nTable))
       allocate(y2(nTable))
       open(unit=iunit,file="galacticus.imf",status="old",form="formatted",iostat=ierr)
       do i=1,nTable
          read (iunit,*,iostat=ierr) massTable(i),imfTable(i)
       end do
       close(iunit)
       massLow=massTable(1)
       massHigh=massTable(nTable)
       lowLim=2.0e30
       upLim=2.0e30
       call spline(massTable,imfTable,lowLim,upLim,y2)
       moduleInitialized=.true.
    end if

    ! Do polynomial interpolation to get the IMF at the requested masses.
    do i=1,size(mass)
       if (mass(i).lt.massLow.or.mass(i).gt.massHigh) then
          imf(i)=0.0
       else
          imf(i)=splint(massTable,imfTable,y2,mass(i))
       end if
    end do

    return
  end subroutine Get_Galacticus_IMF

end module Galacticus_IMF_Module
