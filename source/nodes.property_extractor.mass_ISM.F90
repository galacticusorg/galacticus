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
Implements an ISM mass output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassISM">
   <description>An ISM mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassISM
     !!{
     A stellar mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => massISMExtract
     procedure :: quantity    => massISMQuantity
     procedure :: name        => massISMName
     procedure :: description => massISMDescription
     procedure :: unitsInSI   => massISMUnitsInSI
  end type nodePropertyExtractorMassISM

  interface nodePropertyExtractorMassISM
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassISM} output analysis class.
     !!}
     module procedure massISMConstructorParameters
  end interface nodePropertyExtractorMassISM

contains

  function massISMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassISM} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassISM)                :: self
    type (inputParameters             ), intent(inout) :: parameters

 
    self=nodePropertyExtractorMassISM()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massISMConstructorParameters

  double precision function massISMExtract(self,node,instance)
    !!{
    Implement a massISM output analysis.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , componentTypeSpheroid, massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class(nodePropertyExtractorMassISM), intent(inout), target   :: self
    type (treeNode                    ), intent(inout), target   :: node
    type (multiCounter                ), intent(inout), optional :: instance
    class(massDistributionClass       )               , pointer  :: massDistributionDisk, massDistributionSpheroid
    !$GLC attributes unused :: self, instance

    massDistributionDisk     =>  node                    %massDistribution(massType=massTypeGaseous,componentType=componentTypeDisk    )
    massDistributionSpheroid =>  node                    %massDistribution(massType=massTypeGaseous,componentType=componentTypeSpheroid)
    massISMExtract           =  +massDistributionDisk    %massTotal       (                                                            ) &
         &                      +massDistributionSpheroid%massTotal       (                                                            )
    !![
    <objectDestructor name="massDistributionDisk"    />
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
    return
  end function massISMExtract


  function massISMQuantity(self)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityMass
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: massISMQuantity
    class(nodePropertyExtractorMassISM                 ), intent(inout) :: self
    !$GLC attributes unused :: self

    massISMQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massISMQuantity

  function massISMName(self)
    !!{
    Return the name of the massISM property.
    !!}
    implicit none
    type (varying_string              )                :: massISMName
    class(nodePropertyExtractorMassISM), intent(inout) :: self
    !$GLC attributes unused :: self

    massISMName=var_str('massISM')
    return
  end function massISMName

  function massISMDescription(self)
    !!{
    Return a description of the massISM property.
    !!}
    implicit none
    type (varying_string              )                :: massISMDescription
    class(nodePropertyExtractorMassISM), intent(inout) :: self
    !$GLC attributes unused :: self

    massISMDescription=var_str('The mass of the interstellar medium in each galaxy.')
    return
  end function massISMDescription

  double precision function massISMUnitsInSI(self)
    !!{
    Return the units of the massISM property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassISM), intent(inout) :: self
    !$GLC attributes unused :: self

    massISMUnitsInSI=massSolar
    return
  end function massISMUnitsInSI
