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
Implements a cooling energy radiated property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorCoolingEnergyRadiated">
   <description>A cooling energy radiated property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorCoolingEnergyRadiated
     !!{
     A cooling energy radiated property extractor class.
     !!}
     private
     integer :: energyRadiatedID
   contains
     procedure :: extract     => coolingEnergyRadiatedExtract
     procedure :: name        => coolingEnergyRadiatedName
     procedure :: description => coolingEnergyRadiatedDescription
     procedure :: unitsInSI   => coolingEnergyRadiatedUnitsInSI
  end type nodePropertyExtractorCoolingEnergyRadiated

  interface nodePropertyExtractorCoolingEnergyRadiated
     !!{
     Constructors for the {\normalfont \ttfamily coolingEnergyRadiated} output analysis class.
     !!}
     module procedure coolingEnergyRadiatedConstructorParameters
     module procedure coolingEnergyRadiatedConstructorInternal
  end interface nodePropertyExtractorCoolingEnergyRadiated

contains

  function coolingEnergyRadiatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily coolingEnergyRadiated} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorCoolingEnergyRadiated)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
    
    self=nodePropertyExtractorCoolingEnergyRadiated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function coolingEnergyRadiatedConstructorParameters

  function coolingEnergyRadiatedConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily coolingEnergyRadiated} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorCoolingEnergyRadiated) :: self

    !![
    <addMetaProperty component="hotHalo" name="energyRadiatedBensonBower2010" id="self%energyRadiatedID" isEvolvable="yes"/>
    !!]
    return
  end function coolingEnergyRadiatedConstructorInternal

  double precision function coolingEnergyRadiatedExtract(self,node,instance)
    !!{
    Implement a cooling energy radiated property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodePropertyExtractorCoolingEnergyRadiated), intent(inout), target   :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    class(nodeComponentHotHalo                      ), pointer                 :: hotHalo
    !$GLC attributes unused :: instance

    hotHalo                      => node   %hotHalo                  (                     )
    coolingEnergyRadiatedExtract =  hotHalo%floatRank0MetaPropertyGet(self%energyRadiatedID)
    return
  end function coolingEnergyRadiatedExtract

  function coolingEnergyRadiatedName(self)
    !!{
    Return the name of the cooling energy radiated property.
    !!}
    implicit none
    type (varying_string                            )                :: coolingEnergyRadiatedName
    class(nodePropertyExtractorCoolingEnergyRadiated), intent(inout) :: self
    !$GLC attributes unused :: self

    coolingEnergyRadiatedName=var_str('hotHaloCoolingEnergyRadiated')
    return
  end function coolingEnergyRadiatedName

  function coolingEnergyRadiatedDescription(self)
    !!{
    Return a description of the cooling energy radiated property.
    !!}
    implicit none
    type (varying_string                            )                :: coolingEnergyRadiatedDescription
    class(nodePropertyExtractorCoolingEnergyRadiated), intent(inout) :: self
    !$GLC attributes unused :: self

    coolingEnergyRadiatedDescription=var_str('The energy radiated by the hot halo through cooling.')
    return
  end function coolingEnergyRadiatedDescription

  double precision function coolingEnergyRadiatedUnitsInSI(self)
    !!{
    Return the units of the cooling energy radiated property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    class(nodePropertyExtractorCoolingEnergyRadiated), intent(inout) :: self
    !$GLC attributes unused :: self

    coolingEnergyRadiatedUnitsInSI=ergs*gigaYear
    return
  end function coolingEnergyRadiatedUnitsInSI
