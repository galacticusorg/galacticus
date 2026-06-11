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

  !!{RST
  Implements a node property extractor class that allows selection of output times at which to extract properties.
  !!}

  use :: Output_Times, only : outputTimesClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorOutputSelector" docformat="rst">
   <description>
   A wrapper property extractor that delegates extraction to one or more child ``nodePropertyExtractorClass`` objects but restricts output to a user-specified subset of output times. At each output time, the extractor checks whether that time matches one of the allowed output times (within a relative tolerance set by ``toleranceRelative``); non-matching times return zero-size datasets. This is useful when different properties need to be extracted at different output epochs without running separate simulations.
   </description>
   <linkedList type="multiExtractorList" variable="extractors" next="next" object="extractor_" objectType="nodePropertyExtractorClass"/>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorOutputSelector
     !!{RST
     A node output extractor class that allows selection of output times at which to extract properties.
     !!}
     private
     class           (outputTimesClass), pointer                   :: outputTimes_      => null()
     double precision                  , allocatable, dimension(:) :: times
     double precision                                              :: toleranceRelative
   contains
     !![
     <methods>
        <method method="initialize"  description="Initialize the object after construction."                                             />
        <method method="timeMatches" description="Return true if the current time matches a time for which we should extract properties."/>
     </methods>
     !!]
     final     ::                 outputSelectorDestructor
     procedure :: elementCount => outputSelectorElementCount
     procedure :: initialize   => outputSelectorInitialize
     procedure :: timeMatches  => outputSelectorTimeMatches
  end type nodePropertyExtractorOutputSelector

  interface nodePropertyExtractorOutputSelector
     !!{RST
     Constructors for the ``nodePropertyExtractorOutputSelector`` property extractor class.
     !!}
     module procedure outputSelectorConstructorParameters
     module procedure outputSelectorConstructorInternal
  end interface nodePropertyExtractorOutputSelector

contains

  function outputSelectorConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorOutputSelector`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorOutputSelector)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    type   (multiExtractorList                 ), pointer       :: extractor_
    integer                                                     :: i
    
    self      %extractors => null()
    extractor_            => null()
    do i=1,parameters%copiesCount('nodePropertyExtractor',zeroIfNotPresent=.true.)
       if (associated(extractor_)) then
          allocate(extractor_%next)
          extractor_ => extractor_%next
       else
          allocate(self%extractors)
          extractor_ => self%extractors
       end if
       !![
       <objectBuilder class="nodePropertyExtractor" name="extractor_%extractor_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <objectBuilder class="outputTimes" name="self%outputTimes_" source="parameters"/>
    <inputParameter docformat="rst">
      <name>toleranceRelative</name>
      <variable>self%toleranceRelative</variable>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The relative tolerance to accept when comparing times.
      </description>
    </inputParameter>
    !!]
    call self%initialize()
    !![
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    return
  end function outputSelectorConstructorParameters

  function outputSelectorConstructorInternal(extractors,outputTimes_,toleranceRelative) result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorOutputSelector`` property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorOutputSelector)                         :: self
    type            (multiExtractorList                 ), target , intent(in   ) :: extractors
    class           (outputTimesClass                   ), target , intent(in   ) :: outputTimes_
    double precision                                              , intent(in   ) :: toleranceRelative
    type            (multiExtractorList                 ), pointer                :: extractor_
    !![
    <constructorAssign variables="*outputTimes_, toleranceRelative"/>
    !!]
    
    self      %extractors => extractors
    extractor_            => extractors
    do while (associated(extractor_))
       !![
       <referenceCountIncrement owner="extractor_" object="extractor_"/>
       !!]
       extractor_ => extractor_%next
    end do
    call self%initialize()
    return
  end function outputSelectorConstructorInternal

  subroutine outputSelectorDestructor(self)
    !!{RST
    Destructor for the ``nodePropertyExtractorOutputSelector`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorOutputSelector), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine outputSelectorDestructor

  subroutine outputSelectorInitialize(self)
    !!{RST
    Initialize a ``outputSelector`` object with the list of times to select.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class  (nodePropertyExtractorOutputSelector), intent(inout) :: self
    integer(c_size_t                           )                :: i
    
    allocate(self%times(self%outputTimes_%count()))
    do i=1_c_size_t,self%outputTimes_%count()
       self%times(i)=self%outputTimes_%time(i)
    end do
    return
  end subroutine outputSelectorInitialize

  integer function outputSelectorElementCount(self,elementType,time)
    !!{RST
    Return the number of elements in the outputSelector property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorOutputSelector), intent(inout) :: self
    type            (enumerationElementTypeType         ), intent(in   ) :: elementType
    double precision                                     , intent(in   ) :: time

    if (self%timeMatches(time)) then
       outputSelectorElementCount=self%nodePropertyExtractorMulti%elementCount(elementType,time)
    else
       outputSelectorElementCount=0
    end if
    return
  end function outputSelectorElementCount

  logical function outputSelectorTimeMatches(self,time) result(matches)
    !!{RST
    Return true if the given ``time`` matches a time for which we should extract properties.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (nodePropertyExtractorOutputSelector), intent(inout) :: self
    double precision                                     , intent(in   ) :: time

    matches=any(Values_Agree(time,self%times,relTol=self%toleranceRelative))
    return
  end function outputSelectorTimeMatches
  
