!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <outputTimes name="outputTimesList">
   <description>An output times class which simply reads a list of output times from a parameter.</description>
  </outputTimes>
  !!]
  type, extends(outputTimesClass) :: outputTimesList
     !!{
     Implementation of an output times class which reads a list of output times from a parameter.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_ => null()
     double precision                         , allocatable, dimension(:) :: times                        , redshifts
   contains
     final     ::                 listDestructor
     procedure :: count        => listCount
     procedure :: index        => listIndex
     procedure :: time         => listTime
     procedure :: redshift     => listRedshift
     procedure :: timeNext     => listTimeNext
     procedure :: timePrevious => listTimePrevious
  end type outputTimesList

  interface outputTimesList
     !!{
     Constructors for the {\normalfont \ttfamily list} output times class.
     !!}
     module procedure listConstructorParameters
     module procedure listConstructorInternal
  end interface outputTimesList

contains

  function listConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily list} output times class which takes a parameter set as input.
    !!}
    use :: Array_Utilities  , only : Array_Reverse
    use :: Input_Parameters , only : inputParameter, inputParameters
    use :: Memory_Management, only : allocateArray
    use :: Sorting          , only : sort
    implicit none
    type            (outputTimesList        )                            :: self
    type            (inputParameters        ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_
    double precision                         , allocatable, dimension(:) :: times
    integer         (c_size_t               )                            :: outputCount        , i

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    if      (parameters%isPresent('times'    )) then
       outputCount=parameters%count('times'    )
    else if (parameters%isPresent('redshifts')) then
       outputCount=parameters%count('redshifts')
    else
       outputCount=1_c_size_t
    end if
    call allocateArray(times,[outputCount])
    if (parameters%isPresent('times')) then
       !![
       <inputParameter>
         <name>times</name>
         <description>A list of (space-separated) times at which \glc\ results should be output. Times need not be in any particular order.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       call sort(times)
    else
       !![
       <inputParameter>
         <name>redshifts</name>
         <defaultValue>[0.0d0]</defaultValue>
         <variable>times</variable>
         <description>A list of (space-separated) redshifts at which \glc\ results should be output. Redshifts need not be in any particular order.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       call sort(times)
       times=Array_Reverse(times)
       do i=1,outputCount
          times(i)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(times(i)))
       end do
    end if
    self=outputTimesList(times,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function listConstructorParameters

  function listConstructorInternal(times,cosmologyFunctions_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily list} output times class which takes a parameter set as input.
    !!}
    use :: Memory_Management, only : allocateArray
    implicit none
    type            (outputTimesList        )                              :: self
    double precision                         , intent(in   ), dimension(:) :: times
    class           (cosmologyFunctionsClass), intent(in   ), target       :: cosmologyFunctions_
    integer         (c_size_t               )                              :: i
    !![
    <constructorAssign variables="times, *cosmologyFunctions_"/>
    !!]

    call allocateArray(self%redshifts,shape(times))
    do i=1,size(times)
       self%redshifts(i)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(times(i)))
    end do
    return
  end function listConstructorInternal

  subroutine listDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily list} output times class.
    !!}
    implicit none
    type(outputTimesList), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine listDestructor

  function listCount(self)
    !!{
    Return the number of outputs.
    !!}
    implicit none
    integer(c_size_t       )                :: listCount
    class  (outputTimesList), intent(inout) :: self

    listCount=size(self%times)
    return
  end function listCount

  double precision function listTime(self,indexOutput)
    !!{
    Returns the time of the output indexed by {\normalfont \ttfamily iOutput}.
    !!}
    implicit none
    class  (outputTimesList), intent(inout) :: self
    integer(c_size_t       ), intent(in   ) :: indexOutput

    if (indexOutput >=1 .and. indexOutput <= size(self%times)) then
       listTime=self%times(indexOutput)
    else
       listTime=-1.0d0
    end if
    return
  end function listTime

  double precision function listRedshift(self,indexOutput)
    !!{
    Returns the redshift of the output indexed by {\normalfont \ttfamily indexOutput}.
    !!}
    implicit none
    class  (outputTimesList), intent(inout) :: self
    integer(c_size_t       ), intent(in   ) :: indexOutput

    if (indexOutput >=1 .and. indexOutput <= size(self%times)) then
       listRedshift=self%redshifts(indexOutput)
    else
       listRedshift=-2.0d0
    end if
    return
  end function listRedshift

  function listIndex(self,time,findClosest)
    !!{
    Returns the index of the output given the corresponding time.
    !!}
    use :: Arrays_Search       , only : searchArray  , searchArrayClosest
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    integer         (c_size_t       )                          :: listIndex
    class           (outputTimesList), intent(inout)           :: self
    double precision                 , intent(in   )           :: time
    logical                          , intent(in   ), optional :: findClosest

    if (present(findClosest).and.findClosest) then
       listIndex=searchArrayClosest(self%times,time)
    else
       listIndex=searchArray            (self%times,time)
       if (Values_Differ(time,self%times(listIndex),relTol=1.0d-6)) &
            & call Error_Report('time does not correspond to an output'//{introspection:location})
    end if
    return
  end function listIndex

  double precision function listTimeNext(self,timeCurrent,indexOutput)
    !!{
    Returns the time of the next output after {\normalfont \ttfamily currentTime}.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (outputTimesList), intent(inout)           :: self
    double precision                 , intent(in   )           :: timeCurrent
    integer         (c_size_t       ), intent(  out), optional :: indexOutput
    integer         (c_size_t       )                          :: i

    ! If the current time exceeds the last output, return an unphysical value.
    if      (timeCurrent > self%times(size(self%times))) then
       listTimeNext=-1.0d0
       if (present(indexOutput)) indexOutput=-1
    else if (timeCurrent <  self%times(         1)) then
       listTimeNext=self%times(1)
       if (present(indexOutput)) indexOutput=+1
    else
       i=min(searchArray(self%times,timeCurrent)+1,size(self%times))
       listTimeNext=self%times(i)
       if (present(indexOutput)) indexOutput=i
    end if
    return
  end function listTimeNext

  double precision function listTimePrevious(self,timeCurrent)
    !!{
    Returns the time of the previous output prior to {\normalfont \ttfamily timeCurrent}.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (outputTimesList), intent(inout) :: self
    double precision                 , intent(in   ) :: timeCurrent

    if      (timeCurrent >  self%times(size(self%times))) then
       ! If the current time exceeds the last output, return the last output.
       listTimePrevious=self%times(size(self%times))
    else if (timeCurrent <= self%times(               1)) then
       ! If the current time preceeds the first output, return an unphysical value.
       listTimePrevious=-1.0d0
    else
       listTimePrevious=self%times(searchArray(self%times,timeCurrent))
    end if
    return
  end function listTimePrevious
