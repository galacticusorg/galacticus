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
Implements a half-galactic mass radius output analysis property extractor class.
!!}

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
   contains
     procedure :: extract     => radiusHalfMassGalacticExtract
     procedure :: name        => radiusHalfMassGalacticName
     procedure :: description => radiusHalfMassGalacticDescription
     procedure :: unitsInSI   => radiusHalfMassGalacticUnitsInSI
  end type nodePropertyExtractorRadiusHalfMassGalactic

  interface nodePropertyExtractorRadiusHalfMassGalactic
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusHalfMassGalactic} output analysis class.
     !!}
     module procedure radiusHalfMassGalacticConstructorParameters
  end interface nodePropertyExtractorRadiusHalfMassGalactic

contains

  function radiusHalfMassGalacticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusHalfMassGalactic} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMassGalactic)                :: self
    type (inputParameters                            ), intent(inout) :: parameters

    self=nodePropertyExtractorRadiusHalfMassGalactic()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusHalfMassGalacticConstructorParameters

  double precision function radiusHalfMassGalacticExtract(self,node,instance)
    !!{
    Implement a half-mass output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class(nodePropertyExtractorRadiusHalfMassGalactic), intent(inout), target   :: self
    type (treeNode                                   ), intent(inout), target   :: node
    type (multiCounter                               ), intent(inout), optional :: instance
    class(massDistributionClass                      )               , pointer  :: massDistribution_
    !$GLC attributes unused :: self, instance

    massDistribution_             => node             %massDistribution   (massType      =massTypeGalactic)
    radiusHalfMassGalacticExtract =  massDistribution_%radiusEnclosingMass(massFractional=0.5d0          )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
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
