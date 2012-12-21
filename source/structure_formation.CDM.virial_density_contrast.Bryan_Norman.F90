!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  implicit none
  private
  public :: Virial_Density_Bryan_Norman_Initialize

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
    use ISO_Varying_String
    use Numerical_Comparison
    use Galacticus_Error
    use Cosmological_Parameters
   implicit none
    type(varying_string),          intent(in)    :: virialDensityContrastMethod
    procedure(),          pointer, intent(inout) :: Virial_Density_Contrast_Tabulate
    
    if (virialDensityContrastMethod == 'Bryan-Norman1998') then
       Virial_Density_Contrast_Tabulate => Virial_Density_Bryan_Norman
       ! Check that fitting formulae are applicable to this cosmology.
       if (Omega_DE() == 0.0d0) then
          fitType=fitTypeZeroLambda
       else if (.not.Values_Differ(Omega_Matter()+Omega_DE(),1.0d0,absTol=1.0d-6)) then
          fitType=fitTypeFlatUniverse
       else
          call Galacticus_Error_Report('Virial_Density_Bryan_Norman_Initialize','no fitting formula available for this cosmology')
       end if
    end if
    return
  end subroutine Virial_Density_Bryan_Norman_Initialize

  subroutine Virial_Density_Bryan_Norman(time,deltaVirialTable)
    !% Tabulate the virial density contrast for the \cite{bryan_statistical_1998} fitting function module.
    use Memory_Management
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Numerical_Ranges
    use Tables
    implicit none
    double precision       , intent(in   )              :: time
    class           (table), intent(inout), allocatable :: deltaVirialTable
    integer                                             :: iTime,deltaTableNumberPoints
    double precision                                    :: x

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
        ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(dlog10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    ! Deallocate table if currently allocated.
    if (allocated(deltaVirialTable)) then
       call deltaVirialTable%destroy()
       deallocate(deltaVirialTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaVirialTable)
    select type (deltaVirialTable)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call deltaVirialTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)    
       ! Evaluate the fitting formulae of Bryan & Norman at each time to get the density contrast.
       do iTime=1,deltaTableNumberPoints
          x=Omega_Matter_Total(deltaVirialTable%x(iTime))-1.0d0
          select case (fitType)
          case (fitTypeZeroLambda)
             call deltaVirialTable%populate((18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)/Omega_Matter_Total(deltaVirialTable%x(iTime)),iTime)
          case (fitTypeFlatUniverse)
             call deltaVirialTable%populate((18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/Omega_Matter_Total(deltaVirialTable%x(iTime)),iTime)
          end select
       end do
    end select
    return
  end subroutine Virial_Density_Bryan_Norman

end module Virial_Densities_Bryan_Norman
