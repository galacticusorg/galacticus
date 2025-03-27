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
  Implements a node property extractor class that combines star formation history masses and times.
  !!}

  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use :: Star_Formation_Histories  , only : starFormationHistoryClass
  use :: Output_Times              , only : outputTimesClass
                                                      
  !![
  <nodePropertyExtractor name="nodePropertyExtractorStarFormationHistory">
   <description>An output extractor property extractor class that combines star formation history masses and times.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorStarFormationHistory
     !!{
     An output extractor property extractor class that combines star formation history masses and times.
     !!}
     private
     class(starFormationHistoryClass   ), pointer :: starFormationHistory_ => null()
     class(outputTimesClass            ), pointer :: outputTimes_          => null()
     type (enumerationComponentTypeType)          :: component
   contains
     final :: starFormationHistoryDestructor
  end type nodePropertyExtractorStarFormationHistory

  interface nodePropertyExtractorStarFormationHistory
     !!{
     Constructors for the {\normalfont \ttfamily starFormationHistory} output extractor class.
     !!}
     module procedure starFormationHistoryConstructorParameters
     module procedure starFormationHistoryConstructorInternal
  end interface nodePropertyExtractorStarFormationHistory

contains

  function starFormationHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormationHistory} output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type(nodePropertyExtractorStarFormationHistory)                :: self
    type(inputParameters                          ), intent(inout) :: parameters
    class(starFormationHistoryClass               ), pointer       :: starFormationHistory_
    class(outputTimesClass                        ), pointer       :: outputTimes_
    type(varying_string                           )                :: component

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation history.</description>
    </inputParameter>
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    !!]
    self=nodePropertyExtractorStarFormationHistory(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),starFormationHistory_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"/>
    <objectDestructor name="outputTimes_"         />
    !!]
    return
  end function starFormationHistoryConstructorParameters

  function starFormationHistoryConstructorInternal(component,starFormationHistory_,outputTimes_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationHistory} output extractor property extractor class.
    !!}
    use :: Star_Formation_Histories, only : starFormationHistoryAgesFixedPerOutput
    use :: Error                   , only : Error_Report
    implicit none
    type (nodePropertyExtractorStarFormationHistory)                        :: self
    class(starFormationHistoryClass                ), intent(in   ), target :: starFormationHistory_
    class(outputTimesClass                         ), intent(in   ), target :: outputTimes_
    type (enumerationComponentTypeType             ), intent(in   )         :: component
    !![
    <constructorAssign variables="component, *starFormationHistory_, *outputTimes_"/>
    !!]
    
    allocate   (                                                  self%extractors                )
    allocate   (nodePropertyExtractorStarFormationHistoryMass  :: self%extractors     %extractor_)
    select type    (extractor_ => self%extractors     %extractor_)
    type is    (nodePropertyExtractorStarFormationHistoryMass )
       !![
       <referenceConstruct    isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorStarFormationHistoryMass (component,starFormationHistory_,outputTimes_)"/>
       !!]
    end select
    if (self%starFormationHistory_%ageDistribution() /= starFormationHistoryAgesFixedPerOutput) then
       allocate(                                                  self%extractors%next           )    
       allocate(nodePropertyExtractorStarFormationHistoryTimes :: self%extractors%next%extractor_)
       select type (extractor_ => self%extractors%next%extractor_)
       type is (nodePropertyExtractorStarFormationHistoryTimes)
          !![
	  <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorStarFormationHistoryTimes(component,starFormationHistory_             )"/>
          !!]
       end select
    end if
    return
  end function starFormationHistoryConstructorInternal

  subroutine starFormationHistoryDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationHistory} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorStarFormationHistory), intent(inout) :: self
    
    !![
    <objectDestructor name="self%starFormationHistory_"/>
    <objectDestructor name="self%outputTimes_"         />
    !!]
    return
  end subroutine starFormationHistoryDestructor
