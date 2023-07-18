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
Contains a module which implements a radiusHalfMass property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusHalfMass">
   <description>
    A node property extractor which extracts the half-mass radius of the galaxy. The half-mass radius is output as {\normalfont
    \ttfamily [halfMassRadius]} (in Mpc).
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusHalfMass
     !!{
     A half-mass radius property extractor class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                radiusHalfMassDestructor
     procedure :: extract     => radiusHalfMassExtract
     procedure :: name        => radiusHalfMassName
     procedure :: description => radiusHalfMassDescription
     procedure :: unitsInSI   => radiusHalfMassUnitsInSI
  end type nodePropertyExtractorRadiusHalfMass

  interface nodePropertyExtractorRadiusHalfMass
     !!{
     Constructors for the ``radiusHalfMass'' output analysis class.
     !!}
     module procedure radiusHalfMassConstructorParameters
     module procedure radiusHalfMassConstructorInternal
  end interface nodePropertyExtractorRadiusHalfMass

contains

  function radiusHalfMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusHalfMass} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMass)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(galacticStructureClass             ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusHalfMass(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function radiusHalfMassConstructorParameters

  function radiusHalfMassConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiusHalfMass} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMass)                        :: self
    class(galacticStructureClass             ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function radiusHalfMassConstructorInternal
  
  subroutine radiusHalfMassDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiusHalfMass} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine radiusHalfMassDestructor

  double precision function radiusHalfMassExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar
    implicit none
    class(nodePropertyExtractorRadiusHalfMass), intent(inout), target   :: self
    type (treeNode                           ), intent(inout), target   :: node
    type (multiCounter                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    radiusHalfMassExtract=self%galacticStructure_%radiusEnclosingMass(node,massFractional=0.5d0,massType=massTypeStellar)
    return
  end function radiusHalfMassExtract

  function radiusHalfMassName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string                     )                :: radiusHalfMassName
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassName=var_str('radiusHalfMassStellar')
    return
  end function radiusHalfMassName

  function radiusHalfMassDescription(self)
    !!{
    Return a description of the radiusHalfMass property.
    !!}
    implicit none
    type (varying_string                     )                :: radiusHalfMassDescription
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassDescription=var_str('Radius enclosing half the galaxy stellar mass [Mpc].')
    return
  end function radiusHalfMassDescription

  double precision function radiusHalfMassUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassUnitsInSI=massSolar
    return
  end function radiusHalfMassUnitsInSI

