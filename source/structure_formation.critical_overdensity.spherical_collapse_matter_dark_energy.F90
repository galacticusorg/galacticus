!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of critical overdensity for collapse based on spherical collapse in a
  !% matter plus dark energy universe.

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseMatterDE" defaultThreadPrivate="yes">
  !#  <description>Critical overdensity for collapse based on the spherical collapse in a matter plus dark energy universe.</description>
  !# </criticalOverdensity>
  use Tables

  type, extends(criticalOverdensitySphericalCollapseMatterLambda) :: criticalOverdensitySphericalCollapseMatterDE
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark energy universe.
     private
   contains
     final     ::                 sphericalCollapseMatterDEDestructor
     procedure :: retabulate   => sphericalCollapseMatterDERetabulate
  end type criticalOverdensitySphericalCollapseMatterDE

  interface criticalOverdensitySphericalCollapseMatterDE
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity for collapse class.
     module procedure sphericalCollapseMatterDEConstructorParameters
     module procedure sphericalCollapseMatterDEConstructorInternal
  end interface criticalOverdensitySphericalCollapseMatterDE

contains

  function sphericalCollapseMatterDEConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(criticalOverdensitySphericalCollapseMatterDE)                :: sphericalCollapseMatterDEConstructorParameters
    type(inputParameters                             ), intent(in   ) :: parameters

    sphericalCollapseMatterDEConstructorParameters%tableInitialized=.false.
    return
  end function sphericalCollapseMatterDEConstructorParameters

  function sphericalCollapseMatterDEConstructorInternal()
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity class.
    implicit none
    type(criticalOverdensitySphericalCollapseMatterDE) :: sphericalCollapseMatterDEConstructorInternal

    sphericalCollapseMatterDEConstructorInternal%tableInitialized=.false.
    return
  end function sphericalCollapseMatterDEConstructorInternal

  subroutine sphericalCollapseMatterDEDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity for collapse class.
    implicit none
    type (criticalOverdensitySphericalCollapseMatterDE), intent(inout) :: self
    
    if (self%tableInitialized) then
       call self%overdensityCritical%destroy()
       deallocate(self%overdensityCritical)
    end if
    return
  end subroutine sphericalCollapseMatterDEDestructor

  subroutine sphericalCollapseMatterDERetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterDE), intent(inout) :: self
    double precision                                              , intent(in   ) :: time
    logical                                                                       :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Dark_Energy_Critical_Overdensity_Tabulate(time,self%overdensityCritical)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%overdensityCritical%x(+1)
       self%tableTimeMaximum=self%overdensityCritical%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterDERetabulate
