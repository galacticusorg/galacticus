!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.






!% Contains a module which implements calculations of virial overdensity using the fitting function of
!% \cite{bryan_statistical_1998}.

module Virial_Densities_Bryan_Norman
  private
  public :: Virial_Density_Bryan_Norman_Initialize, Virial_Density_Bryan_Norman_State_Store,&
       & Virial_Density_Bryan_Norman_State_Retrieve

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMinimum=1.0d0, deltaTableTimeMaximum=20.0d0
  integer,          parameter :: deltaTableNPointsPerDecade=100

  ! Labels for different fitting function types.
  integer                     :: fitType
  integer,          parameter :: fitTypeZeroLambda=0, fitTypeFlatUniverse=1

contains

  !# <virialDensityContrastMethod>
  !#  <unitName>Virial_Density_Bryan_Norman_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Virial_Density_Bryan_Norman_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the \cite{bryan_statistical_1998} fitting function module.
    use Input_Parameters
    use ISO_Varying_String
    use Numerical_Comparison
    use Galacticus_Error
    use Cosmological_Parameters
   implicit none
    type(varying_string),          intent(in)    :: virialDensityContrastMethod
    procedure(),          pointer, intent(inout) :: Virial_Density_Contrast_Tabulate
    
    if (virialDensityContrastMethod.eq.'Bryan + Norman') then
       Virial_Density_Contrast_Tabulate => Virial_Density_Bryan_Norman
       ! Check that fitting formulae are applicable to this cosmology.
       if (Omega_DE() == 0.0d0) then
          fitType=fitTypeZeroLambda
       else if (.not.Values_Differ(Omega_0()+Omega_DE(),1.0d0,absTol=1.0d-6)) then
          fitType=fitTypeFlatUniverse
       else
          call Galacticus_Error_Report('Virial_Density_Bryan_Norman_Initialize','no fitting formula available for this cosmology')
       end if
    end if
    return
  end subroutine Virial_Density_Bryan_Norman_Initialize

  subroutine Virial_Density_Bryan_Norman(time,deltaTableNumberPoints,deltaTableTime,deltaTableDelta)
    !% Tabulate the virial density contrast for the \cite{bryan_statistical_1998} fitting function module.
    use Memory_Management
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Numerical_Ranges
    implicit none
    double precision, intent(in)                               :: time
    integer,          intent(out)                              :: deltaTableNumberPoints
    double precision, intent(inout), allocatable, dimension(:) :: deltaTableTime,deltaTableDelta
    integer                                                    :: iTime
    double precision                                           :: x

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(dlog10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    
    ! Deallocate arrays if currently allocated.
    if (allocated(deltaTableTime )) call Dealloc_Array(deltaTableTime )
    if (allocated(deltaTableDelta)) call Dealloc_Array(deltaTableDelta)
    ! Allocate the arrays to current required size.
    call Alloc_Array(deltaTableTime ,deltaTableNumberPoints,'deltaTableTime' )
    call Alloc_Array(deltaTableDelta,deltaTableNumberPoints,'deltaTableDelta')
    
    ! Create set of grid points in time variable.
    deltaTableTime=Make_Range(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints,rangeTypeLogarithmic)
    
    ! Evaluate the fitting formulae of Bryan & Norman at each time to get the density contrast.
    do iTime=1,deltaTableNumberPoints
       x=Omega_Matter(deltaTableTime(iTime))-1.0d0
       select case (fitType)
       case (fitTypeZeroLambda)
          deltaTableDelta(iTime)=(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)/Omega_Matter(deltaTableTime(iTime))
       case (fitTypeFlatUniverse)
          deltaTableDelta(iTime)=(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/Omega_Matter(deltaTableTime(iTime))
       end select
    end do
    return
  end subroutine Virial_Density_Bryan_Norman

end module Virial_Densities_Bryan_Norman
