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
Implements a spheroid stellar mass output analysis property extractor class.
!!}

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
   contains
     procedure :: extract     => massStellarSpheroidExtract
     procedure :: quantity    => massStellarSpheroidQuantity
     procedure :: name        => massStellarSpheroidName
     procedure :: description => massStellarSpheroidDescription
     procedure :: unitsInSI   => massStellarSpheroidUnitsInSI
  end type nodePropertyExtractorMassStellarSpheroid

  interface nodePropertyExtractorMassStellarSpheroid
     !!{
     Constructors for the {\normalfont \ttfamily massStellarSpheroid} output analysis class.
     !!}
     module procedure massStellarSpheroidConstructorParameters
  end interface nodePropertyExtractorMassStellarSpheroid

contains

  function massStellarSpheroidConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massStellarSpheroid} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMassStellarSpheroid)                :: self
    type(inputParameters                         ), intent(inout) :: parameters

    self=nodePropertyExtractorMassStellarSpheroid()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massStellarSpheroidConstructorParameters

  double precision function massStellarSpheroidExtract(self,node,instance)
    !!{
    Implement a stellar mass-weighted morphology output analysis.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, massTypeStellar
    implicit none
    class(nodePropertyExtractorMassStellarSpheroid), intent(inout), target   :: self
    type (treeNode                                ), intent(inout), target   :: node
    type (multiCounter                            ), intent(inout), optional :: instance
    class(massDistributionClass                   )               , pointer  :: massDistribution_
    !$GLC attributes unused :: self, instance

    massDistribution_          => node             %massDistribution(massType=massTypeStellar,componentType=componentTypeSpheroid)
    massStellarSpheroidExtract =  massDistribution_%massTotal       (                                                            )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
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
