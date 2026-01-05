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
Implements a half-stellar mass radius output analysis property extractor class.
!!}

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
   contains
     procedure :: extract     => radiusHalfMassStellarExtract
     procedure :: name        => radiusHalfMassStellarName
     procedure :: description => radiusHalfMassStellarDescription
     procedure :: unitsInSI   => radiusHalfMassStellarUnitsInSI
  end type nodePropertyExtractorRadiusHalfMassStellar

  interface nodePropertyExtractorRadiusHalfMassStellar
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusHalfMassStellar} output analysis class.
     !!}
     module procedure radiusHalfMassStellarConstructorParameters
  end interface nodePropertyExtractorRadiusHalfMassStellar

contains

  function radiusHalfMassStellarConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusHalfMassStellar} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMassStellar)                :: self
    type (inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorRadiusHalfMassStellar()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusHalfMassStellarConstructorParameters

  double precision function radiusHalfMassStellarExtract(self,node,instance)
    !!{
    Implement a half-mass output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class(nodePropertyExtractorRadiusHalfMassStellar), intent(inout), target   :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    class(massDistributionClass                     )               , pointer  :: massDistribution_
    !$GLC attributes unused :: self, instance

    massDistribution_            => node             %massDistribution   (massType      =massTypeStellar)
    radiusHalfMassStellarExtract =  massDistribution_%radiusEnclosingMass(massFractional=0.5d0          )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
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
