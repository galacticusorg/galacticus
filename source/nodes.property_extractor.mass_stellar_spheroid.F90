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
Contains a module which implements a spheroid stellar mass output analysis property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassStellarSpheroid">
   <description>A spheroid stellar mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassStellarSpheroid
     !!{
     A stellar mass output analysis class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                massStellarSpheroidDestructor
     procedure :: extract     => massStellarSpheroidExtract
     procedure :: quantity    => massStellarSpheroidQuantity
     procedure :: name        => massStellarSpheroidName
     procedure :: description => massStellarSpheroidDescription
     procedure :: unitsInSI   => massStellarSpheroidUnitsInSI
  end type nodePropertyExtractorMassStellarSpheroid

  interface nodePropertyExtractorMassStellarSpheroid
     !!{
     Constructors for the ``massStellarSpheroid'' output analysis class.
     !!}
     module procedure massStellarSpheroidConstructorParameters
     module procedure massStellarSpheroidConstructorInternal
  end interface nodePropertyExtractorMassStellarSpheroid

contains

  function massStellarSpheroidConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massStellarSpheroid'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassStellarSpheroid)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(galacticStructureClass                  ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMassStellarSpheroid(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function massStellarSpheroidConstructorParameters

  function massStellarSpheroidConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the ``massStellarSpheroid'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorMassStellarSpheroid)                        :: self
    class(galacticStructureClass                  ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function massStellarSpheroidConstructorInternal
  
  subroutine massStellarSpheroidDestructor(self)
    !!{
    Destructor for the ``massStellarSpheroid'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassStellarSpheroid), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine massStellarSpheroidDestructor

  double precision function massStellarSpheroidExtract(self,node,instance)
    !!{
    Implement a stellar mass-weighted morphology output analysis.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, massTypeStellar, radiusLarge
    implicit none
    class           (nodePropertyExtractorMassStellarSpheroid), intent(inout)           :: self
    type            (treeNode                                ), intent(inout), target   :: node
    type            (multiCounter                            ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    massStellarSpheroidExtract=self%galacticStructure_%massEnclosed(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeSpheroid)
    return
  end function massStellarSpheroidExtract


  function massStellarSpheroidQuantity(self)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityMass
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: massStellarSpheroidQuantity
    class(nodePropertyExtractorMassStellarSpheroid     ), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarSpheroidQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massStellarSpheroidQuantity

  function massStellarSpheroidName(self)
    !!{
    Return the name of the massStellarSpheroid property.
    !!}
    implicit none
    type (varying_string                          )                :: massStellarSpheroidName
    class(nodePropertyExtractorMassStellarSpheroid), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarSpheroidName=var_str('massStellarSpheroid')
    return
  end function massStellarSpheroidName

  function massStellarSpheroidDescription(self)
    !!{
    Return a description of the massStellarSpheroid property.
    !!}
    implicit none
    type (varying_string                          )                :: massStellarSpheroidDescription
    class(nodePropertyExtractorMassStellarSpheroid), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarSpheroidDescription=var_str('The stellar mass of the galactic spheroid.')
    return
  end function massStellarSpheroidDescription

  double precision function massStellarSpheroidUnitsInSI(self)
    !!{
    Return the units of the massStellarSpheroid property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassStellarSpheroid), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarSpheroidUnitsInSI=massSolar
    return
  end function massStellarSpheroidUnitsInSI
