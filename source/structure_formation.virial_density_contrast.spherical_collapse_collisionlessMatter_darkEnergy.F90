!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus dark eneryg universe.</description>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt) :: virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark eneryg universe.
     private
     integer :: energyFixedAt
   contains
  end type virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy

  interface virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrDrkEnrgy} dark matter halo virial density contrast class.
     module procedure sphericalCollapseCllsnlssMttrDrkEnrgyConstructorParameters
     module procedure sphericalCollapseCllsnlssMttrDrkEnrgyConstructorInternal
  end interface virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy

contains

  function sphericalCollapseCllsnlssMttrDrkEnrgyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrDrkEnrgy} dark matter halo virial density contrast class that takes a parameter set as input.
    use :: Input_Parameters          , only : inputParameter                                , inputParameters
    use :: Spherical_Collapse_Solvers, only : enumerationCllsnlssMttrDarkEnergyFixedAtEncode
    implicit none
    type   (virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy)                :: self
    type   (inputParameters                                           ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                                   ), pointer       :: cosmologyFunctions_
    type   (varying_string                                            )                :: energyFixedAt
    logical                                                                            :: tableStore

    !# <inputParameter>
    !#   <name>tableStore</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
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
    self=virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy(tableStore,enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function sphericalCollapseCllsnlssMttrDrkEnrgyConstructorParameters

  function sphericalCollapseCllsnlssMttrDrkEnrgyConstructorInternal(tableStore,energyFixedAt,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrDrkEnrgy} dark matter halo virial density contrast class.
    use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrDarkEnergy
    implicit none
    type   (virialDensityContrastSphericalCollapseCllsnlssMttrDrkEnrgy)                        :: self
    integer                                                            , intent(in   )         :: energyFixedAt
    logical                                                            , intent(in   )         :: tableStore
    class  (cosmologyFunctionsClass                                   ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="tableStore, energyFixedAt, *cosmologyFunctions_"/>

    self%tableInitialized     =.false.
    self%turnaroundInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrDarkEnergy :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrDarkEnergy)
       !# <referenceConstruct isResult="yes" owner="self" nameAssociated="sphericalCollapseSolver_" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrDarkEnergy(self%energyFixedAt,self%cosmologyFunctions_)"/>
    end select
    return
  end function sphericalCollapseCllsnlssMttrDrkEnrgyConstructorInternal
