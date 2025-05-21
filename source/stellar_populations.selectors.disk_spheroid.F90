!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a stellar population selector class which returns a different population for disks and spheroids.
  !!}

  !![
  <stellarPopulationSelector name="stellarPopulationSelectorDiskSpheroid">
   <description>
    A stellar population selector class which selects a different population for disks and spheroids, irrespective of other
    physical conditions. The populations to use are specified by the {\normalfont \ttfamily [stellarPopulationDisk]} and
    {\normalfont \ttfamily [stellarPopulationSpheroid]} parameters.
   </description>
  </stellarPopulationSelector>
  !!]
  type, extends(stellarPopulationSelectorClass) :: stellarPopulationSelectorDiskSpheroid
     !!{
     A stellar population selector class which returns a different population for disks and spheroids.
     !!}
     private
     class(stellarPopulationClass), pointer :: stellarPopulationDisk_ => null(), stellarPopulationSpheroid_ => null()
   contains
     final     ::                                 diskSpheroidDestructor
     procedure :: select                       => diskSpheroidSelect
     procedure :: isStarFormationRateDependent => diskSpheroidIsStarFormationRateDependent
  end type stellarPopulationSelectorDiskSpheroid

  interface stellarPopulationSelectorDiskSpheroid
     !!{
     Constructors for the \refClass{stellarPopulationSelectorDiskSpheroid} stellar population selector class.
     !!}
     module procedure diskSpheroidConstructorParameters
     module procedure diskSpheroidConstructorInternal
  end interface stellarPopulationSelectorDiskSpheroid

contains

  function diskSpheroidConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSelectorDiskSpheroid} stellar population class which takes a parameter list as input.
    !!}
    use :: Input_Parameters   , only : inputParameter   , inputParameters
    use :: Stellar_Populations, only : stellarPopulation, stellarPopulationClass
    implicit none
    type (stellarPopulationSelectorDiskSpheroid)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(stellarPopulationClass               ), pointer       :: stellarPopulationDisk_, stellarPopulationSpheroid_

    !![
    <objectBuilder class="stellarPopulation" parameterName="stellarPopulationDisk"     name="stellarPopulationDisk_"     source="parameters"/>
    <objectBuilder class="stellarPopulation" parameterName="stellarPopulationSpheroid" name="stellarPopulationSpheroid_" source="parameters"/>
    !!]
    self=stellarPopulationSelectorDiskSpheroid(stellarPopulationDisk_,stellarPopulationSpheroid_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarPopulationDisk_"    />
    <objectDestructor name="stellarPopulationSpheroid_"/>
    !!]
    return
  end function diskSpheroidConstructorParameters

  function diskSpheroidConstructorInternal(stellarPopulationDisk_,stellarPopulationSpheroid_) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationSelectorDiskSpheroid} stellar population selector class.
    !!}
    implicit none
    type (stellarPopulationSelectorDiskSpheroid)                        :: self
    class(stellarPopulationClass               ), intent(in   ), target :: stellarPopulationDisk_, stellarPopulationSpheroid_
    !![
    <constructorAssign variables="*stellarPopulationDisk_, *stellarPopulationSpheroid_"/>
    !!]

    return
  end function diskSpheroidConstructorInternal

  subroutine diskSpheroidDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationSelectorDiskSpheroid} stellar population selector class.
    !!}
    implicit none
    type(stellarPopulationSelectorDiskSpheroid), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulationDisk_"    />
    <objectDestructor name="self%stellarPopulationSpheroid_"/>
    !!]
    return
  end subroutine diskSpheroidDestructor

  function diskSpheroidSelect(self,rateStarFormation,abundances_,component)
    !!{
    Return a diskSpheroid stellar population.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponent, nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class           (stellarPopulationClass               ), pointer       :: diskSpheroidSelect
    class           (stellarPopulationSelectorDiskSpheroid), intent(inout) :: self
    double precision                                       , intent(in   ) :: rateStarFormation
    type            (abundances                           ), intent(in   ) :: abundances_
    class           (nodeComponent                        ), intent(in   ) :: component
    !$GLC attributes unused :: rateStarFormation, abundances_

    select type (component)
    class is (nodeComponentDisk    )
       diskSpheroidSelect => self%stellarPopulationDisk_
    class is (nodeComponentSpheroid)
       diskSpheroidSelect => self%stellarPopulationSpheroid_
    class default
       call Error_Report('only disk and spheroid components are supported by this class'//{introspection:location})
    end select
    return
  end function diskSpheroidSelect

  logical function diskSpheroidIsStarFormationRateDependent(self)
    !!{
    Return false indicating that stellar population selection is not dependent on star formation rate.
    !!}
    implicit none
    class(stellarPopulationSelectorDiskSpheroid), intent(inout) :: self
    !$GLC attributes unused :: self

    diskSpheroidIsStarFormationRateDependent=.false.
    return
  end function diskSpheroidIsStarFormationRateDependent

