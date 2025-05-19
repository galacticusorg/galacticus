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
Implements a radiusCooling property extractor class.
!!}

  use :: Cooling_Radii, only : coolingRadius, coolingRadiusClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusCooling">
   <description>
    A cooling radius property extractor class. Extracts the characteristic cooling radius in the halo in Mpc.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusCooling
     !!{
     A cooling radius property extractor class.
     !!}
     private
     class(coolingRadiusClass), pointer :: coolingRadius_ => null()
   contains
     final     ::                radiusCoolingDestructor
     procedure :: extract     => radiusCoolingExtract
     procedure :: name        => radiusCoolingName
     procedure :: description => radiusCoolingDescription
     procedure :: unitsInSI   => radiusCoolingUnitsInSI
  end type nodePropertyExtractorRadiusCooling

  interface nodePropertyExtractorRadiusCooling
     !!{
     Constructors for the {\normalfont \ttfamily radiusCooling} output analysis class.
     !!}
     module procedure radiusCoolingConstructorParameters
     module procedure radiusCoolingConstructorInternal
  end interface nodePropertyExtractorRadiusCooling

contains

  function radiusCoolingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusCooling} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiusCooling)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(coolingRadiusClass                ), pointer       :: coolingRadius_

    !![
    <objectBuilder class="coolingRadius" name="coolingRadius_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusCooling(coolingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRadius_"/>
    !!]
    return
  end function radiusCoolingConstructorParameters

  function radiusCoolingConstructorInternal(coolingRadius_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiusCooling} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusCooling)                        :: self
    class(coolingRadiusClass                ), intent(in   ), target :: coolingRadius_
    !![
    <constructorAssign variables="*coolingRadius_"/>
    !!]

    return
  end function radiusCoolingConstructorInternal

  subroutine radiusCoolingDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiusCooling} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusCooling), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRadius_"/>
    !!]
    return
  end subroutine radiusCoolingDestructor

  double precision function radiusCoolingExtract(self,node,instance)
    !!{
    Implement a cooling radius property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorRadiusCooling), intent(inout), target   :: self
    type (treeNode                          ), intent(inout), target   :: node
    type (multiCounter                      ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    radiusCoolingExtract=self%coolingRadius_%radius(node)
    return
  end function radiusCoolingExtract

  function radiusCoolingName(self)
    !!{
    Return the name of the cooling radius property.
    !!}
    implicit none
    type (varying_string                    )                :: radiusCoolingName
    class(nodePropertyExtractorRadiusCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusCoolingName=var_str('hotHaloRadiusCooling')
    return
  end function radiusCoolingName

  function radiusCoolingDescription(self)
    !!{
    Return a description of the cooling radius property.
    !!}
    implicit none
    type (varying_string                    )                :: radiusCoolingDescription
    class(nodePropertyExtractorRadiusCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusCoolingDescription=var_str('Cooling radius in the hot halo [Mpc].')
    return
  end function radiusCoolingDescription

  double precision function radiusCoolingUnitsInSI(self)
    !!{
    Return the units of the cooling radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusCooling), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusCoolingUnitsInSI=megaParsec
    return
  end function radiusCoolingUnitsInSI


