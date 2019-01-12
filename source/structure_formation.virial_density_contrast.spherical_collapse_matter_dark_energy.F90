!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  use Tables

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseMatterDE">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus dark eneryg universe.</description>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastSphericalCollapseMatterLambda) :: virialDensityContrastSphericalCollapseMatterDE
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark eneryg universe.
     private
     logical                                :: turnaroundInitialized
     double precision                       :: turnaroundTimeMinimum, turnaroundTimeMaximum
     class           (table1D), allocatable :: turnaround
     integer                                :: energyFixedAt
   contains
     final     ::                              sphericalCollapseMatterDEDestructor
     procedure :: retabulate                => sphericalCollapseMatterDERetabulate
     procedure :: turnAroundOverVirialRadii => sphericalCollapseMatterDETurnAroundOverVirialRadii
  end type virialDensityContrastSphericalCollapseMatterDE

  interface virialDensityContrastSphericalCollapseMatterDE
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterDE} dark matter halo virial density contrast class.
     module procedure sphericalCollapseMatterDEConstructorParameters
     module procedure sphericalCollapseMatterDEConstructorInternal
  end interface virialDensityContrastSphericalCollapseMatterDE

contains

  function sphericalCollapseMatterDEConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} dark matter halo virial density contrast class that takes a parameter set as input.
    use Input_Parameters
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    type (virialDensityContrastSphericalCollapseMatterDE)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    type (varying_string                                )                :: energyFixedAt
    
    !# <inputParameter>
    !#   <name>energyFixedAt</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('turnaround')</defaultValue>
    !#   <description>Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
    !#     ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
    !#     \citealt{percival_cosmological_2005}; \S8.)</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_" source="parameters"/>
    self=virialDensityContrastSphericalCollapseMatterDE(enumerationDarkEnergySphericalCollapseEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function sphericalCollapseMatterDEConstructorParameters

  function sphericalCollapseMatterDEConstructorInternal(energyFixedAt,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} dark matter halo virial density contrast class.
    implicit none
    type   (virialDensityContrastSphericalCollapseMatterDE)                        :: self
    integer                                                , intent(in   )         :: energyFixedAt
    class  (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="energyFixedAt, *cosmologyFunctions_"/>

    self%tableInitialized     =.false.
    self%turnaroundInitialized=.false.
    return
  end function sphericalCollapseMatterDEConstructorInternal

  subroutine sphericalCollapseMatterDEDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseMatterDE), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    if (self%turnaroundInitialized) then
       call self%turnaround%destroy()
       deallocate(self%turnaround)
    end if
    return
  end subroutine sphericalCollapseMatterDEDestructor

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
       call Spherical_Collapse_Dark_Energy_Virial_Density_Contrast_Tabulate(time,self%energyFixedAt,self%deltaVirial,self%cosmologyFunctions_)
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
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterDE), intent(inout)           :: self
    double precision                                                , intent(in   ), optional :: time       , expansionFactor
    logical                                                         , intent(in   ), optional :: collapsing
    logical                                                                                   :: remakeTable, collapsingActual
    double precision                                                                          :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('only one argument can be specified'//{introspection:location})
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
          timeActual=self%cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Galacticus_Error_Report('at least one argument must be given'//{introspection:location})
       end if
    end if
    ! Check if we need to recompute our table.
    if (self%turnaroundInitialized) then
       remakeTable=(timeActual < self%turnaroundTimeMinimum .or. timeActual > self%turnaroundTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Dark_Energy_Turnaround_Radius_Tabulate(timeActual,self%energyFixedAt,self%turnaround,self%cosmologyFunctions_)
       self%turnaroundInitialized=.true.
       self%turnaroundTimeMinimum=self%turnaround%x(+1)
       self%turnaroundTimeMaximum=self%turnaround%x(-1)
    end if
    ! Interpolate to get the ratio.
    sphericalCollapseMatterDETurnAroundOverVirialRadii=self%turnaround%interpolate(timeActual)
    return
  end function sphericalCollapseMatterDETurnAroundOverVirialRadii
