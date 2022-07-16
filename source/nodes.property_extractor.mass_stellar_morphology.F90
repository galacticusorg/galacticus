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
Contains a module which implements a stellar mass-weighted morphology output analysis property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassStellarMorphology">
   <description>A stellar mass-weighted morphology output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassStellarMorphology
     !!{
     A stellar mass output analysis class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                massStellarMorphologyDestructor
     procedure :: extract     => massStellarMorphologyExtract
     procedure :: name        => massStellarMorphologyName
     procedure :: description => massStellarMorphologyDescription
     procedure :: unitsInSI   => massStellarMorphologyUnitsInSI
  end type nodePropertyExtractorMassStellarMorphology

  interface nodePropertyExtractorMassStellarMorphology
     !!{
     Constructors for the ``massStellarMorphology'' output analysis class.
     !!}
     module procedure massStellarMorphologyConstructorParameters
     module procedure massStellarMorphologyConstructorInternal
  end interface nodePropertyExtractorMassStellarMorphology

contains

  function massStellarMorphologyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massStellarMorphology'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassStellarMorphology)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(galacticStructureClass                    ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMassStellarMorphology(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function massStellarMorphologyConstructorParameters

  function massStellarMorphologyConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the ``massStellarMorphology'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorMassStellarMorphology)                        :: self
    class(galacticStructureClass                    ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function massStellarMorphologyConstructorInternal
  
  subroutine massStellarMorphologyDestructor(self)
    !!{
    Destructor for the ``massStellarMorphology'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine massStellarMorphologyDestructor

  double precision function massStellarMorphologyExtract(self,node,instance)
    !!{
    Implement a stellar mass-weighted morphology output analysis.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, massTypeStellar, radiusLarge
    implicit none
    class           (nodePropertyExtractorMassStellarMorphology), intent(inout)           :: self
    type            (treeNode                                  ), intent(inout), target   :: node
    type            (multiCounter                              ), intent(inout), optional :: instance
    double precision                                                                      :: massStellarDisk, massStellarSpheroid
    !$GLC attributes unused :: self, instance

    massStellarDisk    =self%galacticStructure_%massEnclosed(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeDisk    )
    massStellarSpheroid=self%galacticStructure_%massEnclosed(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeSpheroid)
    if (massStellarDisk+massStellarSpheroid > 0.0d0) then
       massStellarMorphologyExtract=+  massStellarSpheroid &
            &                       /(                     &
            &                         +massStellarDisk     &
            &                         +massStellarSpheroid &
            &                       )
    else
       massStellarMorphologyExtract=+0.0d0
    end if
    return
  end function massStellarMorphologyExtract


  function massStellarMorphologyName(self)
    !!{
    Return the name of the massStellarMorphology property.
    !!}
    implicit none
    type (varying_string                            )                :: massStellarMorphologyName
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyName=var_str('morphologyMassStellarWeighted')
    return
  end function massStellarMorphologyName

  function massStellarMorphologyDescription(self)
    !!{
    Return a description of the massStellarMorphology property.
    !!}
    implicit none
    type (varying_string                            )                :: massStellarMorphologyDescription
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyDescription=var_str('The stellar mass-weighted mophology (bulge-to-total ratio) of the galaxy.')
    return
  end function massStellarMorphologyDescription

  double precision function massStellarMorphologyUnitsInSI(self)
    !!{
    Return the units of the massStellarMorphology property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyUnitsInSI=0.0d0
    return
  end function massStellarMorphologyUnitsInSI
