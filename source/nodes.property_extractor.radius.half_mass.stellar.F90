!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a half-stellar mass radius output analysis property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusHalfMassStellar">
   <description>A half-(stellar) mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusHalfMassStellar
     !!{
     A half-(stellar) mass property extractor output analysis class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                radiusHalfMassStellarDestructor
     procedure :: extract     => radiusHalfMassStellarExtract
     procedure :: name        => radiusHalfMassStellarName
     procedure :: description => radiusHalfMassStellarDescription
     procedure :: unitsInSI   => radiusHalfMassStellarUnitsInSI
  end type nodePropertyExtractorRadiusHalfMassStellar

  interface nodePropertyExtractorRadiusHalfMassStellar
     !!{
     Constructors for the ``radiusHalfMassStellar'' output analysis class.
     !!}
     module procedure radiusHalfMassStellarConstructorParameters
     module procedure radiusHalfMassStellarConstructorInternal
  end interface nodePropertyExtractorRadiusHalfMassStellar

contains

  function radiusHalfMassStellarConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``radiusHalfMassStellar'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMassStellar)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(galacticStructureClass                    ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusHalfMassStellar(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function radiusHalfMassStellarConstructorParameters

  function radiusHalfMassStellarConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the ``radiusHalfMassStellar'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusHalfMassStellar)                        :: self
    class(galacticStructureClass                    ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function radiusHalfMassStellarConstructorInternal
  
  subroutine radiusHalfMassStellarDestructor(self)
    !!{
    Destructor for the ``radiusHalfMassStellar'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusHalfMassStellar), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine radiusHalfMassStellarDestructor

  double precision function radiusHalfMassStellarExtract(self,node,instance)
    !!{
    Implement a half-mass output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar
    implicit none
    class(nodePropertyExtractorRadiusHalfMassStellar), intent(inout)           :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    radiusHalfMassStellarExtract=self%galacticStructure_%radiusEnclosingMass(node,massFractional=0.5d0,massType=massTypeStellar)
    return
  end function radiusHalfMassStellarExtract


  function radiusHalfMassStellarName(self)
    !!{
    Return the name of the radiusHalfMassStellar property.
    !!}
    implicit none
    type (varying_string                            )                :: radiusHalfMassStellarName
    class(nodePropertyExtractorRadiusHalfMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassStellarName=var_str('radiusHalfMassStellar')
    return
  end function radiusHalfMassStellarName

  function radiusHalfMassStellarDescription(self)
    !!{
    Return a description of the radiusHalfMassStellar property.
    !!}
    implicit none
    type (varying_string                            )                :: radiusHalfMassStellarDescription
    class(nodePropertyExtractorRadiusHalfMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassStellarDescription=var_str('The stellar half-mass radius.')
    return
  end function radiusHalfMassStellarDescription

  double precision function radiusHalfMassStellarUnitsInSI(self)
    !!{
    Return the units of the radiusHalfMassStellar property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusHalfMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassStellarUnitsInSI=megaParsec
    return
  end function radiusHalfMassStellarUnitsInSI
