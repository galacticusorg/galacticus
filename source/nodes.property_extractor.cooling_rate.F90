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
Implements a cooling rate property extractor class.
!!}

  use :: Cooling_Rates, only : coolingRate, coolingRateClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRateCooling">
   <description>
   A cooling rate property extractor class. Extracts the rate at which gas is cooling from the halo (assuming no sources of
   heating) in $M_\odot$ Gyr$^{-1}$.
  </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRateCooling
     !!{
     A rateCooling property extractor class.
     !!}
     private
     class(coolingRateClass), pointer :: coolingRate_ => null()
   contains
     final     ::                rateCoolingDestructor
     procedure :: extract     => rateCoolingExtract
     procedure :: name        => rateCoolingName
     procedure :: description => rateCoolingDescription
     procedure :: unitsInSI   => rateCoolingUnitsInSI
  end type nodePropertyExtractorRateCooling

  interface nodePropertyExtractorRateCooling
     !!{
     Constructors for the \refClass{nodePropertyExtractorRateCooling} output analysis class.
     !!}
     module procedure rateCoolingConstructorParameters
     module procedure rateCoolingConstructorInternal
  end interface nodePropertyExtractorRateCooling

contains

  function rateCoolingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRateCooling} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRateCooling)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(coolingRateClass                ), pointer       :: coolingRate_

    !![
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRateCooling(coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function rateCoolingConstructorParameters

  function rateCoolingConstructorInternal(coolingRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorRateCooling} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRateCooling)                        :: self
    class(coolingRateClass                ), intent(in   ), target :: coolingRate_
    !![
    <constructorAssign variables="*coolingRate_"/>
    !!]

    return
  end function rateCoolingConstructorInternal

  subroutine rateCoolingDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorRateCooling} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRateCooling), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine rateCoolingDestructor

  double precision function rateCoolingExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRateCooling), intent(inout), target   :: self
    type (treeNode                        ), intent(inout), target   :: node
    type (multiCounter                    ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    rateCoolingExtract=self%coolingRate_%rate(node)
    return
  end function rateCoolingExtract

  function rateCoolingName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string                  )                :: rateCoolingName
    class(nodePropertyExtractorRateCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    rateCoolingName=var_str('hotHaloRateCooling')
    return
  end function rateCoolingName

  function rateCoolingDescription(self)
    !!{
    Return a description of the rateCooling property.
    !!}
    implicit none
    type (varying_string                  )                :: rateCoolingDescription
    class(nodePropertyExtractorRateCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    rateCoolingDescription=var_str('The rate of mass cooling in the hot halo [M☉/Gyr].')
    return
  end function rateCoolingDescription

  double precision function rateCoolingUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar
    implicit none
    class(nodePropertyExtractorRateCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    rateCoolingUnitsInSI=massSolar/gigaYear
    return
  end function rateCoolingUnitsInSI


