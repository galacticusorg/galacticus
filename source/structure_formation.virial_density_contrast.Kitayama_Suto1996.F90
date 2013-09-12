!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
!% \cite{kitayama_semianalytic_1996}.

module Virial_Densities_Kitayama_Suto1996
  implicit none
  private
  public :: Virial_Density_Kitayama_Suto1996_Initialize

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMaximum     =20.0d0, deltaTableTimeMinimum=1.0d0
  integer         , parameter :: deltaTableNPointsPerDecade=100

contains

  !# <virialDensityContrastMethod>
  !#  <unitName>Virial_Density_Kitayama_Suto1996_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Virial_Density_Kitayama_Suto1996_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the \cite{kitayama_semianalytic_1996} fitting function module.
    use ISO_Varying_String
    use Numerical_Comparison
    use Galacticus_Error
    use Cosmology_Parameters
   implicit none
    type     (varying_string                  ), intent(in   )          :: virialDensityContrastMethod
    procedure(Virial_Density_Kitayama_Suto1996), intent(inout), pointer :: Virial_Density_Contrast_Tabulate
    class    (cosmologyParametersClass        )               , pointer :: thisCosmologyParameters

    if (virialDensityContrastMethod == 'Kitayama-Suto1996') then
       Virial_Density_Contrast_Tabulate => Virial_Density_Kitayama_Suto1996
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       ! Check that fitting formula is applicable to this cosmology.
       if (Values_Differ(thisCosmologyParameters%OmegaMatter()+thisCosmologyParameters%OmegaDarkEnergy(),1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Virial_Density_Kitayama_Suto1996_Initialize','no fitting formula available for this cosmology')
    end if
    return
  end subroutine Virial_Density_Kitayama_Suto1996_Initialize

  subroutine Virial_Density_Kitayama_Suto1996(time,deltaVirialTable)
    !% Tabulate the virial density contrast for the \cite{kitayama_semianalytic_1996} fitting function module.
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Tables
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaVirialTable
    class           (cosmologyFunctionsClass), pointer                    :: cosmologyFunctionsDefault
    integer                                                               :: deltaTableNumberPoints   , iTime
    double precision                                                      :: omegaf

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)

    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    ! Deallocate table if currently allocated.
    if (allocated(deltaVirialTable)) then
       call deltaVirialTable%destroy()
       deallocate(deltaVirialTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaVirialTable)
    select type (deltaVirialTable)
    type is (table1DLogarithmicLinear)
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       ! Create the table.
       call deltaVirialTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Evaluate the fitting formulae of Kitayama & Suto at each time to get the density contrast.
       do iTime=1,deltaTableNumberPoints
          omegaf=1.0d0/cosmologyFunctionsDefault%omegaMatterEpochal(deltaVirialTable%x(iTime))-1.0d0
          call deltaVirialTable%populate(18.0d0*Pi**2*(1.0d0+0.4093d0*omegaf**0.9052d0),iTime)
       end do
    end select
    return
  end subroutine Virial_Density_Kitayama_Suto1996

end module Virial_Densities_Kitayama_Suto1996
