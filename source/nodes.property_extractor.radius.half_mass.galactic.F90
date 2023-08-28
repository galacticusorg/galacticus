!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a half-galactic mass radius output analysis property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusHalfMassGalactic">
   <description>A half-galactic mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusHalfMassGalactic
     !!{
     A half-galactic mass property extractor output analysis class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                radiusHalfMassGalacticDestructor
     procedure :: extract     => radiusHalfMassGalacticExtract
     procedure :: name        => radiusHalfMassGalacticName
     procedure :: description => radiusHalfMassGalacticDescription
     procedure :: unitsInSI   => radiusHalfMassGalacticUnitsInSI
  end type nodePropertyExtractorRadiusHalfMassGalactic

  interface nodePropertyExtractorRadiusHalfMassGalactic
     !!{
     Constructors for the ``radiusHalfMassGalactic'' output analysis class.
     !!}
     module procedure radiusHalfMassGalacticConstructorParameters
     module procedure radiusHalfMassGalacticConstructorInternal
  end interface nodePropertyExtractorRadiusHalfMassGalactic

contains

  function radiusHalfMassGalacticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``radiusHalfMassGalactic'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMassGalactic)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(galacticStructureClass                     ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusHalfMassGalactic(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function radiusHalfMassGalacticConstructorParameters

  function radiusHalfMassGalacticConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the ``radiusHalfMassGalactic'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusHalfMassGalactic)                        :: self
    class(galacticStructureClass                     ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function radiusHalfMassGalacticConstructorInternal
  
  subroutine radiusHalfMassGalacticDestructor(self)
    !!{
    Destructor for the ``radiusHalfMassGalactic'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine radiusHalfMassGalacticDestructor

  double precision function radiusHalfMassGalacticExtract(self,node,instance)
    !!{
    Implement a half-mass output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeGalactic
    implicit none
    class(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout), target   :: self
    type (treeNode                                   ), intent(inout), target   :: node
    type (multiCounter                               ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    radiusHalfMassGalacticExtract=self%galacticStructure_%radiusEnclosingMass(node,massFractional=0.5d0,massType=massTypeGalactic)
    return
  end function radiusHalfMassGalacticExtract


  function radiusHalfMassGalacticName(self)
    !!{
    Return the name of the radiusHalfMassGalactic property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusHalfMassGalacticName
    class(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassGalacticName=var_str('radiusHalfMassGalactic')
    return
  end function radiusHalfMassGalacticName

  function radiusHalfMassGalacticDescription(self)
    !!{
    Return a description of the radiusHalfMassGalactic property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusHalfMassGalacticDescription
    class(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassGalacticDescription=var_str('The galactic half-mass radius.')
    return
  end function radiusHalfMassGalacticDescription

  double precision function radiusHalfMassGalacticUnitsInSI(self)
    !!{
    Return the units of the radiusHalfMassGalactic property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassGalacticUnitsInSI=megaParsec
    return
  end function radiusHalfMassGalacticUnitsInSI
