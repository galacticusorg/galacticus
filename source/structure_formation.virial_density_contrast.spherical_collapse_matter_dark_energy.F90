!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo virial density contrasts based on spherical collapse in a matter plus dark energy universe.

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseMatterDE">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus dark eneryg universe.</description>
  !# </virialDensityContrast>
  use Tables

  type, extends(virialDensityContrastSphericalCollapseMatterLambda) :: virialDensityContrastSphericalCollapseMatterDE
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark eneryg universe.
     private
   contains
     procedure :: retabulate => sphericalCollapseMatterDERetabulate
  end type virialDensityContrastSphericalCollapseMatterDE

  interface virialDensityContrastSphericalCollapseMatterDE
     !% Constructors for the {\tt sphericalCollapseMatterDE} dark matter halo virial density contrast class.
     module procedure sphericalCollapseMatterDEDefaultConstructor
  end interface virialDensityContrastSphericalCollapseMatterDE

contains

  function sphericalCollapseMatterDEDefaultConstructor()
    !% Default constructor for the {\tt sphericalCollapseMatterDE} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type (virialDensityContrastSphericalCollapseMatterDE), target  :: sphericalCollapseMatterDEDefaultConstructor
    
    sphericalCollapseMatterDEDefaultConstructor%tableInitialized=.false.
    return
  end function sphericalCollapseMatterDEDefaultConstructor

  subroutine sphericalCollapseMatterDERetabulate(self,time)
    !% Recompute the look-up tables for virial density contrast.
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterDE), intent(inout) :: self
    double precision                                                        , intent(in   ) :: time
    logical                                                                                 :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Dark_Energy_Virial_Density_Contrast_Tabulate(time,self%deltaVirial)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%deltaVirial%x(+1)
       self%tableTimeMaximum=self%deltaVirial%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterDERetabulate
