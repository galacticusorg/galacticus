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
Implements a stellar mass effective radius node property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusEffectiveStellar">
   <description>A stellar mass effective radius node property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusEffectiveStellar
     !!{
     A stellar mass effective radius node property extractor class. That is, the projected radius (assuming a face-on galaxy)
     containing half the stellar mass is found.     
     !!}
     private
   contains
     procedure :: extract     => radiusEffectiveStellarExtract
     procedure :: name        => radiusEffectiveStellarName
     procedure :: description => radiusEffectiveStellarDescription
     procedure :: unitsInSI   => radiusEffectiveStellarUnitsInSI
  end type nodePropertyExtractorRadiusEffectiveStellar

  interface nodePropertyExtractorRadiusEffectiveStellar
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusEffectiveStellar} output analysis class.
     !!}
     module procedure radiusEffectiveStellarConstructorParameters
  end interface nodePropertyExtractorRadiusEffectiveStellar

contains

  function radiusEffectiveStellarConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusEffectiveStellar} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusEffectiveStellar)                :: self
    type (inputParameters                            ), intent(inout) :: parameters

    self=nodePropertyExtractorRadiusEffectiveStellar()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusEffectiveStellarConstructorParameters

  double precision function radiusEffectiveStellarExtract(self,node,instance)
    !!{
    Extract the stellar mass effective radius.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class(nodePropertyExtractorRadiusEffectiveStellar), intent(inout), target   :: self
    type (treeNode                                   ), intent(inout), target   :: node
    type (multiCounter                               ), intent(inout), optional :: instance
    class(massDistributionClass                      )               , pointer  :: massDistribution_
    !$GLC attributes unused :: self, instance

    massDistribution_             => node             %massDistribution              (massType      =massTypeStellar)
    radiusEffectiveStellarExtract =  massDistribution_%radiusCylindricalEnclosingMass(massFractional=0.5d0          )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function radiusEffectiveStellarExtract

  function radiusEffectiveStellarName(self)
    !!{
    Return the name of the radiusEffectiveStellar property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusEffectiveStellarName
    class(nodePropertyExtractorRadiusEffectiveStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusEffectiveStellarName=var_str('radiusEffectiveStellar')
    return
  end function radiusEffectiveStellarName

  function radiusEffectiveStellarDescription(self)
    !!{
    Return a description of the radiusEffectiveStellar property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusEffectiveStellarDescription
    class(nodePropertyExtractorRadiusEffectiveStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusEffectiveStellarDescription=var_str('The stellar mass effecive radius (i.e. the projected radius in the face-on galaxy containing half the stellar mass).')
    return
  end function radiusEffectiveStellarDescription

  double precision function radiusEffectiveStellarUnitsInSI(self)
    !!{
    Return the units of the radiusEffectiveStellar property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusEffectiveStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusEffectiveStellarUnitsInSI=megaParsec
    return
  end function radiusEffectiveStellarUnitsInSI
