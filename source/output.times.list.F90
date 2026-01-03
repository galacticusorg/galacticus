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

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <outputTimes name="outputTimesList">
    <description>      
      An output times class which simply reads a list of output times specified via parameters. Times can be given as a
      (space-separated) list of actual cosmic times (in Gyr) via the {\normalfont \ttfamily [times]} parameter, or as a
      (space-separated) list of redshifts via the {\normalfont \ttfamily [redshifts]} parameter, or by a combination of the
      two. The {\normalfont \ttfamily [times]} parameter allows negative values which are interpreted as lookback times. For
      example, in a cosmological model where the universe is currently 13.8~Gyr old the following:
      \begin{verbatim}
      &lt;outputTimes value="list"&gt;
	&lt;redshifts value=" 0.0"/&gt;
	&lt;times     value="-1.0"/&gt;
      &lt;/outtputTimes&gt;
      \end{verbatim}
      will result in output at 12.8 and 13.8~Gyr.
    </description>
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
     Constructors for the \refClass{outputTimesList} output times class.
     !!}
     module procedure listConstructorParameters
     module procedure listConstructorInternal
  end interface outputTimesList

contains

  function listConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesList} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters , only : inputParameter, inputParameters
    use :: Sorting          , only : sort
    use :: Error            , only : Error_Report
    implicit none
    type            (outputTimesList        )                            :: self
    type            (inputParameters        ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_
    double precision                         , allocatable, dimension(:) :: times              , redshifts    , &
         &                                                                  timesFromRedshifts , timesCombined
    integer         (c_size_t               )                            :: i
    logical :: haveTimes, haveRedshifts

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    haveTimes    =parameters%isPresent('times'    )
    haveRedshifts=parameters%isPresent('redshifts')
    if (haveRedshifts.or..not.(haveRedshifts.or.haveTimes)) then
       allocate(redshifts         (max(1,parameters%count('redshifts',zeroIfNotPresent=.true.))))
       allocate(timesFromRedshifts(size(redshifts)                                             ))
       !![
       <inputParameter>
         <name>redshifts</name>
         <defaultValue>[0.0d0]</defaultValue>
         <description>A list of (space-separated) redshifts at which \glc\ results should be output. Redshifts need not be in any particular order.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       ! Convert redshifts to times.
       do i=1,size(redshifts)
          timesFromRedshifts(i)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshifts(i)))
       end do
    else
       allocate(timesFromRedshifts(0))
    end if
    if (haveTimes) then
       allocate(times(parameters%count('times')))
       !![
       <inputParameter>
         <name>times</name>
         <description>A list of (space-separated) times at which \glc\ results should be output. Times need not be in any particular order. Negative times are interpreted as look-back times.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       ! Convert any look-back times to actual times.
       do i=1,size(times)
          if (times(i) < 0.0d0) then
             times(i)=cosmologyFunctions_%cosmicTime(1.0d0)+times(i)
             if (times(i) < 0.0d0) call Error_Report('look-back time exceeds the age of the universe'//{introspection:location})
          end if
       end do
    else
       allocate(times(0))
    end if
    allocate(timesCombined(size(times)+size(timesFromRedshifts)))
    if (size(times             ) > 0) timesCombined(1            :size(times)                         )=times
    if (size(timesFromRedshifts) > 0) timesCombined(1+size(times):size(times)+size(timesFromRedshifts))=timesFromRedshifts
    call sort(timesCombined)
    self=outputTimesList(timesCombined,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function listConstructorParameters

  function listConstructorInternal(times,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \refClass{outputTimesList} output times class which takes a parameter set as input.
    !!}
    implicit none
    type            (outputTimesList        )                              :: self
    double precision                         , intent(in   ), dimension(:) :: times
    class           (cosmologyFunctionsClass), intent(in   ), target       :: cosmologyFunctions_
    integer         (c_size_t               )                              :: i
    !![
    <constructorAssign variables="times, *cosmologyFunctions_"/>
    !!]
    allocate(self%redshifts,mold=times)
    do i=1,size(times)
       self%redshifts(i)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(times(i)))
    end do
    return
  end function listConstructorInternal

  subroutine listDestructor(self)
    !!{
    Destructor for the \refClass{outputTimesList} output times class.
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
    use :: Arrays_Search       , only : searchArray   , searchArrayClosest
    use :: Display             , only : displayMessage, displayIndent     , displayUnindent
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ , Values_Agree
    implicit none
    integer         (c_size_t       )                          :: listIndex
    class           (outputTimesList), intent(inout)           :: self
    double precision                 , intent(in   )           :: time
    logical                          , intent(in   ), optional :: findClosest
    character       (len=16         )                          :: label

    if (present(findClosest).and.findClosest) then
       listIndex=searchArrayClosest(self%times,time)
    else
       listIndex=searchArray       (self%times,time)
       if (Values_Differ(time,self%times(listIndex),relTol=1.0d-6)) then
          ! Check for a match at the final time - because we use an array search function above the index of the final output is never returned.
          if (Values_Agree(time,self%times(size(self%times)),relTol=1.0d-6)) then
             listIndex=size(self%times)
          else
             call displayIndent('Output time matching:')
             write (label,'(f16.12)') time
             call displayMessage('Target time: '//label//' Gyr')
             call displayIndent('Available times:')
             do listIndex=1,size(self%times)
                write (label,'(f16.12)') self%times(listIndex)
                call displayMessage(label//' Gyr')
             end do
             call displayUnindent('')
             call displayUnindent('')
             call Error_Report('time does not correspond to an output'//{introspection:location})
          end if
       end if
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
    if      (timeCurrent >= self%times(size(self%times))) then
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

  double precision function listTimePrevious(self,timeCurrent,indexOutput)
    !!{
    Returns the time of the previous output prior to {\normalfont \ttfamily timeCurrent}.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (outputTimesList), intent(inout)           :: self
    double precision                 , intent(in   )           :: timeCurrent
    integer         (c_size_t       ), intent(  out), optional :: indexOutput
    integer         (c_size_t       )                          :: i

    if      (timeCurrent >  self%times(size(self%times))) then
       ! If the current time exceeds the last output, return the last output.
       listTimePrevious=self%times(size(self%times))
       if (present(indexOutput)) indexOutput=size(self%times)
    else if (timeCurrent <= self%times(               1)) then
       ! If the current time preceeds the first output, return an unphysical value.
       listTimePrevious=-1.0d0
       if (present(indexOutput)) indexOutput=-1
    else
       i               =searchArray(self%times,timeCurrent)
       listTimePrevious=self%times(i)
       if (present(indexOutput)) indexOutput=i
    end if
    return
  end function listTimePrevious
