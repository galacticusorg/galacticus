!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  An implementation of dark matter halo virial density contrasts based on spherical collapse in a matter plus dark energy universe.
  !!}

  use :: Spherical_Collapse_Solvers, only : enumerationCllsnlssMttrDarkEnergyFixedAtType

  !![
  <virialDensityContrast name="virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy">
   <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus dark energy universe.</description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt) :: virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy
     !!{
     A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark energy universe.
     !!}
     private
     type(enumerationCllsnlssMttrDarkEnergyFixedAtType) :: energyFixedAt
   contains
  end type virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy

  interface virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy
     !!{
     Constructors for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy} dark matter halo virial density contrast class.
     !!}
     module procedure sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters
     module procedure sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal
  end interface virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy

contains

  function sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy} dark matter halo virial density contrast class that takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                                , inputParameters
    use :: Spherical_Collapse_Solvers, only : enumerationCllsnlssMttrDarkEnergyFixedAtEncode
    implicit none
    type   (virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy)                :: self
    type   (inputParameters                                          ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                                  ), pointer       :: cosmologyFunctions_
    type   (varying_string                                           )                :: energyFixedAt
    logical                                                                           :: tableStore

    !![
    <inputParameter>
      <name>tableStore</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    </inputParameter>
    <inputParameter>
      <name>energyFixedAt</name>
      <defaultValue>var_str('turnaround')</defaultValue>
      <description>Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
        ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
        \citealt{percival_cosmological_2005}; \S8.)</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy(tableStore,enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters

  function sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal(tableStore,energyFixedAt,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy} dark matter halo virial density contrast class.
    !!}
    use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrDarkEnergy
    implicit none
    type   (virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy)                        :: self
    type   (enumerationCllsnlssMttrDarkEnergyFixedAtType             ), intent(in   )         :: energyFixedAt
    logical                                                           , intent(in   )         :: tableStore
    class  (cosmologyFunctionsClass                                  ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="tableStore, energyFixedAt, *cosmologyFunctions_"/>
    !!]

    self%tableInitialized     =.false.
    self%turnaroundInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrDarkEnergy :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrDarkEnergy)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="sphericalCollapseSolver_" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrDarkEnergy(self%energyFixedAt,self%cosmologyFunctions_)"/>
       !!]
    end select
    return
  end function sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal
