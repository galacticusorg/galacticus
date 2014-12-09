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
     logical                                :: turnaroundInitialized
     double precision                       :: turnaroundTimeMinimum, turnaroundTimeMaximum
     class           (table1D), allocatable :: turnaround
   contains
     procedure :: retabulate                => sphericalCollapseMatterDERetabulate
     procedure :: turnAroundOverVirialRadii => sphericalCollapseMatterDETurnAroundOverVirialRadii
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
    
    sphericalCollapseMatterDEDefaultConstructor%tableInitialized     =.false.
    sphericalCollapseMatterDEDefaultConstructor%turnaroundInitialized=.false.
    return
  end function sphericalCollapseMatterDEDefaultConstructor

  subroutine sphericalCollapseMatterDERetabulate(self,time)
    !% Recompute the look-up tables for virial density contrast.
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterDE), intent(inout) :: self
    double precision                                                , intent(in   ) :: time
    logical                                                                         :: remakeTable

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

  double precision function sphericalCollapseMatterDETurnAroundOverVirialRadii(self,time,expansionFactor,collapsing)
    !% Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    !% constant universe.
    use Galacticus_Error
    use Cosmology_Functions
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterDE), intent(inout)           :: self
    double precision                                                , intent(in   ), optional :: time               , expansionFactor
    logical                                                         , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass                       ), pointer                 :: cosmologyFunctions_
    logical                                                                                   :: remakeTable        , collapsingActual
    double precision                                                                          :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('sphericalCollapseMatterDETurnAroundOverVirialRadii','only one argument can be specified')
       else
          timeActual=time
       end if
    else
       if (present(expansionFactor)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          ! Get the default cosmology functions object.
          cosmologyFunctions_ => cosmologyFunctions()
          timeActual=cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Galacticus_Error_Report('sphericalCollapseMatterDEDETurnAroundOverVirialRadii','at least one argument must be given')
       end if
    end if
    ! Check if we need to recompute our table.
    if (self%turnaroundInitialized) then
       remakeTable=(timeActual < self%turnaroundTimeMinimum .or. timeActual > self%turnaroundTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Dark_Energy_Turnaround_Radius_Tabulate(timeActual,self%turnaround)
       self%turnaroundInitialized=.true.
       self%turnaroundTimeMinimum=self%turnaround%x(+1)
       self%turnaroundTimeMaximum=self%turnaround%x(-1)
    end if
    ! Interpolate to get the ratio.
    sphericalCollapseMatterDETurnAroundOverVirialRadii=self%turnaround%interpolate(timeActual)
    return
  end function sphericalCollapseMatterDETurnAroundOverVirialRadii
