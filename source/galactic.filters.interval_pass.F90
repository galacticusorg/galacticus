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
Implements an interval pass filter on any node property.
!!}

  use :: Node_Property_Extractors, only : nodePropertyExtractorScalar

  !![
  <galacticFilter name="galacticFilterIntervalPass">
   <description>an interval pass filter on any node property.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterIntervalPass
     !!{
     an interval pass galactic filter class on any node property.
     !!}
     private
     class           (nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
     double precision                                       :: thresholdLow                    , thresholdHigh
   contains
     final     ::           intervalPassDestructor
     procedure :: passes => intervalPassPasses
  end type galacticFilterIntervalPass

  interface galacticFilterIntervalPass
     !!{
     Constructors for the {\normalfont \ttfamily intervalPass} galactic filter class.
     !!}
     module procedure intervalPassConstructorParameters
     module procedure intervalPassConstructorInternal
  end interface galacticFilterIntervalPass

contains
  
  function intervalPassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily intervalPass} galactic filter class which takes a parameter set as input.
    !!}
    use :: Error                   , only : Error_Report
    use :: Input_Parameters        , only : inputParameter            , inputParameters
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass
    implicit none
    type            (galacticFilterIntervalPass)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (nodePropertyExtractorClass), pointer       :: nodePropertyExtractor_
    double precision                                            :: thresholdLow          , thresholdHigh
    
    !![
    <inputParameter>
      <name>thresholdLow</name>
      <source>parameters</source>
      <description>The low threshold value above which to pass.</description>
    </inputParameter>
    <inputParameter>
      <name>thresholdHigh</name>
      <source>parameters</source>
      <description>The high threshold value below which to pass.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=galacticFilterIntervalPass(thresholdLow,thresholdHigh,nodePropertyExtractor_)
       class default
       call Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function intervalPassConstructorParameters

  function intervalPassConstructorInternal(thresholdLow,thresholdHigh,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily intervalPass} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterIntervalPass )                        :: self
    class           (nodePropertyExtractorScalar), intent(in   ), target :: nodePropertyExtractor_
    double precision                             , intent(in   )         :: thresholdLow          , thresholdHigh
    !![
    <constructorAssign variables="thresholdLow, thresholdHigh, *nodePropertyExtractor_"/>
    !!]

    return
  end function intervalPassConstructorInternal

  subroutine intervalPassDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily intervalPass} galactic filter class.
    !!}
    implicit none
    type(galacticFilterIntervalPass), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine intervalPassDestructor

  logical function intervalPassPasses(self,node)
    !!{
    Implement an interval pass galactic filter.
    !!}
    implicit none
    class           (galacticFilterIntervalPass), intent(inout)         :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                                                    :: propertyValue

    propertyValue     =                  self%nodePropertyExtractor_%extract(node)
    intervalPassPasses= propertyValue >= self%thresholdLow                         &
         &             .and.                                                       &
         &              propertyValue <  self%thresholdHigh
    return
  end function intervalPassPasses
