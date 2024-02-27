!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a node property extractor class that allows selection of output times at which to extract properties.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorOutputSelector">
   <description>A node property extractor class that allows selection of output times at which to extract properties.</description>
   <linkedList type="multiExtractorList" variable="extractors" next="next" object="extractor_" objectType="nodePropertyExtractorClass"/>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorOutputSelector
     !!{
     A node output extractor class that allows selection of output times at which to extract properties.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_ => null()
     double precision                         , allocatable, dimension(:) :: redshifts                    , times
   contains
     !![
     <methods>
        <method method="timeMatches" description="Return true if the current time matches a time for which we should extract properties."/>
     </methods>
     !!]
     final     ::                 outputSelectorDestructor
     procedure :: elementCount => outputSelectorElementCount
     procedure :: timeMatches  => outputSelectorTimeMatches
  end type nodePropertyExtractorOutputSelector

  interface nodePropertyExtractorOutputSelector
     !!{
     Constructors for the ``outputSelector'' output extractor class.
     !!}
     module procedure outputSelectorConstructorParameters
     module procedure outputSelectorConstructorInternal
  end interface nodePropertyExtractorOutputSelector

contains

  function outputSelectorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``outputSelector'' output extractor property extractor class which takes a parameter set as input.
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
    <objectBuilder class="cosmologyFunctions" name="self%cosmologyFunctions_" source="parameters"/>
    !!]
    if (parameters%isPresent('redshifts')) then
       if (parameters%isPresent('times')) call Error_Report("only one of 'redshifts' and 'times' may be specified"//{introspection:location})
       allocate(self%redshifts(parameters%count('redshifts')))
       !![
       <inputParameter>
	 <name>redshifts</name>
	 <variable>self%redshifts</variable>
	 <description>A list of (space-separated) redshifts at which properties should be output.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
       allocate(self%times(size(self%redshifts)))
       do i=1,size(self%redshifts)
          self%times(i)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshifts(i)))
       end do
    else if (parameters%isPresent('times')) then
       if (parameters%isPresent('redshifts')) call Error_Report("only one of 'redshifts' and 'times' may be specified"//{introspection:location})
       allocate(self%times(parameters%count('times')))
       !![
       <inputParameter>
	 <name>times</name>
	 <variable>self%times</variable>
	 <description>A list of (space-separated) times at which properties should be output.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    else
       call Error_Report("either 'redshifts' or 'times' must be specified"//{introspection:location})
    end if
    !![
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    return
  end function outputSelectorConstructorParameters

  function outputSelectorConstructorInternal(extractors,cosmologyFunctions_,times) result(self)
    !!{
    Internal constructor for the ``outputSelector'' output extractor property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorOutputSelector)                              :: self
    type            (multiExtractorList                 ), target      , intent(in   ) :: extractors
    class           (cosmologyFunctionsClass            ), target      , intent(in   ) :: cosmologyFunctions_
    double precision                                     , dimension(:), intent(in   ) :: times
    type            (multiExtractorList                 ), pointer                     :: extractor_
    integer                                                                            :: i
    !![
    <constructorAssign variables="*cosmologyFunctions_, times"/>
    !!]
    
    self      %extractors => extractors
    extractor_            => extractors
    do while (associated(extractor_))
       !![
       <referenceCountIncrement owner="extractor_" object="extractor_"/>
       !!]
       extractor_ => extractor_%next
    end do
    allocate(self%redshifts(size(self%times)))
    do i=1,size(self%times)
       self%redshifts(i)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%times(i)))
    end do
   return
  end function outputSelectorConstructorInternal

  subroutine outputSelectorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily outputSelector} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorOutputSelector), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine outputSelectorDestructor

  integer function outputSelectorElementCount(self,elementType,time)
    !!{
    Return the number of elements in the outputSelectorple property extractors.
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
    !!{
    Return true if the given {\normalfont \ttfamily time} matches a time for which we should extract properties.
    !!}
    implicit none
    class           (nodePropertyExtractorOutputSelector), intent(inout) :: self
    double precision                                     , intent(in   ) :: time

    matches=any(time == self%times)
    return
  end function outputSelectorTimeMatches
  
