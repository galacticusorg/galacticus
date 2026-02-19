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
Implements a cold mode infall rate property extractor class.
!!}

  use :: Cooling_Cold_Mode_Infall_Rates, only : coldModeInfallRate, coldModeInfallRateClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRateInfallColdMode">
   <description>A cold mode infall rate property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRateInfallColdMode
     !!{
     A cold mode infall rate property extractor class.
     !!}
     private
     class(coldModeInfallRateClass), pointer :: coldModeInfallRate_ => null()
   contains
     final     ::                rateInfallColdModeDestructor
     procedure :: extract     => rateInfallColdModeExtract
     procedure :: name        => rateInfallColdModeName
     procedure :: description => rateInfallColdModeDescription
     procedure :: unitsInSI   => rateInfallColdModeUnitsInSI
  end type nodePropertyExtractorRateInfallColdMode

  interface nodePropertyExtractorRateInfallColdMode
     !!{
     Constructors for the \refClass{nodePropertyExtractorRateInfallColdMode} output analysis class.
     !!}
     module procedure rateInfallColdModeConstructorParameters
     module procedure rateInfallColdModeConstructorInternal
  end interface nodePropertyExtractorRateInfallColdMode

contains

  function rateInfallColdModeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRateInfallColdMode} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRateInfallColdMode)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(coldModeInfallRateClass                ), pointer       :: coldModeInfallRate_

    !![
    <objectBuilder class="coldModeInfallRate" name="coldModeInfallRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRateInfallColdMode(coldModeInfallRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coldModeInfallRate_"/>
    !!]
    return
  end function rateInfallColdModeConstructorParameters

  function rateInfallColdModeConstructorInternal(coldModeInfallRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorRateInfallColdMode} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRateInfallColdMode)                        :: self
    class(coldModeInfallRateClass                ), intent(in   ), target :: coldModeInfallRate_
    !![
    <constructorAssign variables="*coldModeInfallRate_"/>
    !!]

    return
  end function rateInfallColdModeConstructorInternal

  subroutine rateInfallColdModeDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorRateInfallColdMode} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self

    !![
    <objectDestructor name="self%coldModeInfallRate_"/>
    !!]
    return
  end subroutine rateInfallColdModeDestructor

  double precision function rateInfallColdModeExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily rateInfallColdMode} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRateInfallColdMode), intent(inout), target   :: self
    type (treeNode                               ), intent(inout), target   :: node
    type (multiCounter                           ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    rateInfallColdModeExtract=self%coldModeInfallRate_%infallRate(node)
    return
  end function rateInfallColdModeExtract

  function rateInfallColdModeName(self)
    !!{
    Return the name of the {\normalfont \ttfamily rateInfallColdMode} property.
    !!}
    implicit none
    type (varying_string                         )                :: rateInfallColdModeName
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeName=var_str('coldModeInfallRate')
    return
  end function rateInfallColdModeName

  function rateInfallColdModeDescription(self)
    !!{
    Return a description of the {\normalfont \ttfamily rateInfallColdMode} property.
    !!}
    implicit none
    type (varying_string                         )                :: rateInfallColdModeDescription
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeDescription=var_str('Rate of infall of cold mode material onto the galaxy [M☉/Gyr].')
    return
  end function rateInfallColdModeDescription

  double precision function rateInfallColdModeUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily rateInfallColdMode} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar
    implicit none
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeUnitsInSI=massSolar/gigaYear
    return
  end function rateInfallColdModeUnitsInSI


