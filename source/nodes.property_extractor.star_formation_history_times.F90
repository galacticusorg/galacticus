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
  Implements a property extractor class for the star formation history of a component.
  !!}
  
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use :: Star_Formation_Histories  , only : starFormationHistoryClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorStarFormationHistoryTimes">
    <description>A property extractor class for the star formation history tabulation times of a component.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorStarFormationHistoryTimes
     !!{
     A property extractor class for the star formation history tabulation times of a component.
     !!}
     private
     class(starFormationHistoryClass   ), pointer :: starFormationHistory_ => null()
     type (enumerationComponentTypeType)          :: component
   contains
     final     ::                 starFormationHistoryTimesDestructor
     procedure :: elementCount => starFormationHistoryTimesElementCount
     procedure :: extract      => starFormationHistoryTimesExtract
     procedure :: names        => starFormationHistoryTimesNames
     procedure :: descriptions => starFormationHistoryTimesDescriptions
     procedure :: unitsInSI    => starFormationHistoryTimesUnitsInSI
  end type nodePropertyExtractorStarFormationHistoryTimes
  
  interface nodePropertyExtractorStarFormationHistoryTimes
     !!{
     Constructors for the ``starFormationHistoryTimes'' output analysis class.
     !!}
     module procedure starFormationHistoryTimesConstructorParameters
     module procedure starFormationHistoryTimesConstructorInternal
  end interface nodePropertyExtractorStarFormationHistoryTimes
      
contains

  function starFormationHistoryTimesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormationHistoryTimes} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type (nodePropertyExtractorStarFormationHistoryTimes)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(starFormationHistoryClass                     ), pointer       :: starFormationHistory_
    type (varying_string                                )                :: component
    
    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation history.</description>
    </inputParameter>
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    !!]
    self=nodePropertyExtractorStarFormationHistoryTimes(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"/>
    !!]
    return
  end function starFormationHistoryTimesConstructorParameters

  function starFormationHistoryTimesConstructorInternal(component,starFormationHistory_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationHistoryTimes} property extractor class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNuclearStarCluster
    use :: Error                     , only : Error_Report
    implicit none
    type (nodePropertyExtractorStarFormationHistoryTimes)                        :: self
    type (enumerationComponentTypeType                  ), intent(in   )         :: component
    class(starFormationHistoryClass                     ), intent(in   ), target :: starFormationHistory_
    !![
    <constructorAssign variables="component, *starFormationHistory_"/>
    !!]
    
    if     (                                                                                                                          &
         &   component /= componentTypeDisk                                                                                           &
         &  .and.                                                                                                                     &
         &   component /= componentTypeSpheroid                                                                                       &
         &  .and.                                                                                                                     &
         &   component /= componentTypeNuclearStarCluster                                                                             &
         & ) call Error_Report("only 'disk', 'spheroid' and 'nuclearStarCluster' components are supported"//{introspection:location})    
    return
  end function starFormationHistoryTimesConstructorInternal

  subroutine starFormationHistoryTimesDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationHistoryTime} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorStarFormationHistoryTimes), intent(inout) :: self
    
    !![
    <objectDestructor name="self%starFormationHistory_"/>
    !!]
    return
  end subroutine starFormationHistoryTimesDestructor

  integer function starFormationHistoryTimesElementCount(self)
    !!{
    Return the number of elements in the {\normalfont \ttfamily starFormationHistoryTimes} property extractors.
    !!}
    implicit none
    class(nodePropertyExtractorStarFormationHistoryTimes), intent(inout) :: self

    starFormationHistoryTimesElementCount=1
    return
  end function starFormationHistoryTimesElementCount

  function starFormationHistoryTimesExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily starFormationHistoryTimes} property extractor.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNuclearStarCluster
    use :: Histories                 , only : history
    implicit none
    double precision                                                , dimension(:,:), allocatable :: starFormationHistoryTimesExtract
    class           (nodePropertyExtractorStarFormationHistoryTimes), intent(inout)               :: self
    type            (treeNode                                      ), intent(inout)               :: node
    type            (multiCounter                                  ), intent(inout) , optional    :: instance
    class           (nodeComponentDisk                             )                , pointer     :: disk
    class           (nodeComponentSpheroid                         )                , pointer     :: spheroid
    class           (nodeComponentNSC                              )                , pointer     :: nuclearStarCluster
    type            (history                                       )                              :: starFormationHistory
    double precision                                                , dimension(:  ), allocatable :: times 
    !$GLC attributes unused :: instance

    ! Get the relevant star formation history.
    select case (self%component%ID)
    case (componentTypeDisk              %ID)
       disk                 => node              %disk                ()
       starFormationHistory =  disk              %starFormationHistory()
    case (componentTypeSpheroid          %ID)
       spheroid             => node              %spheroid            ()
       starFormationHistory =  spheroid          %starFormationHistory()
    case (componentTypeNuclearStarCluster%ID)
       nuclearStarCluster   => node              %NSC                 ()
       starFormationHistory =  nuclearStarCluster%starFormationHistory()
    end select
    if (starFormationHistory%exists()) then
       times=self%starFormationHistory_%times(node=node,starFormationHistory=starFormationHistory,allowTruncation=.true.)
       allocate(starFormationHistoryTimesExtract(size(times),1))
       starFormationHistoryTimesExtract(:,1)=times
    else
       allocate(starFormationHistoryTimesExtract(0          ,1))
    end if
    return
  end function starFormationHistoryTimesExtract

  subroutine starFormationHistoryTimesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily starFormationHistoryTimes} properties.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class(nodePropertyExtractorStarFormationHistoryTimes), intent(inout)                             :: self
    type(varying_string                                 ), intent(inout), dimension(:) , allocatable :: names

    allocate(names(1))
    names(1)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)//"StarFormationHistoryTimes"
    return
  end subroutine starFormationHistoryTimesNames

  subroutine starFormationHistoryTimesDescriptions(self,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily starFormationHistoryTimes} property.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class(nodePropertyExtractorStarFormationHistoryTimes), intent(inout)                             :: self
    type (varying_string                                ), intent(inout), dimension(:) , allocatable :: descriptions

    allocate(descriptions(1))
    descriptions(1)="Star formation history tabulation times for the "//enumerationComponentTypeDecode(self%component,includePrefix=.false.)//" [Gyr]."
    return
  end subroutine starFormationHistoryTimesDescriptions

  function starFormationHistoryTimesUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily starFormationHistoryTimes} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                                , allocatable  , dimension(:) :: starFormationHistoryTimesUnitsInSI
    class           (nodePropertyExtractorStarFormationHistoryTimes), intent(inout)               :: self

    allocate(starFormationHistoryTimesUnitsInSI(1))
    starFormationHistoryTimesUnitsInSI(1)=gigaYear
    return
  end function starFormationHistoryTimesUnitsInSI
