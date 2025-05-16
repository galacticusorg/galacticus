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

  !![
  <outputTimes name="outputTimesUniformSpacingInTime">
   <description>An output times class which generates a set of output times spaced uniformly in time.</description>
  </outputTimes>
  !!]
  type, extends(outputTimesList) :: outputTimesUniformSpacingInTime
     !!{
     Implementation of an output times class which generates a set of output times spaced uniformly in time.
     !!}
     private
     double precision           :: timeMinimum    , timeMaximum    , &
          &                        redshiftMinimum, redshiftMaximum
     integer         (c_size_t) :: countTimes
  end type outputTimesUniformSpacingInTime

  interface outputTimesUniformSpacingInTime
     !!{
     Constructors for the \refClass{outputTimesUniformSpacingInTime} output times class.
     !!}
     module procedure uniformSpacingInTimeConstructorParameters
     module procedure uniformSpacingInTimeConstructorInternal
  end interface outputTimesUniformSpacingInTime

contains

  function uniformSpacingInTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesUniformSpacingInTime} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Error           , only : Error_Report
    implicit none
    type            (outputTimesUniformSpacingInTime)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    double precision                                                 :: timeMinimum        , timeMaximum    , &
         &                                                              redshiftMinimum    , redshiftMaximum
    integer         (c_size_t                       )                :: countTimes
    
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    if (parameters%isPresent('timeMinimum')) then
       if (parameters%isPresent('redshiftMaximum')) call Error_Report("can not specify both 'timeMinimum' and 'redshiftMaximum'"//{introspection:location})
       !![
       <inputParameter>
	 <name>timeMinimum</name>
	 <description>The minimum time at which to output. Negative times are interpreted as look-back times.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
       if (timeMinimum < 0.0d0) then
          timeMinimum=cosmologyFunctions_%cosmicTime(1.0d0)+timeMinimum
          if (timeMinimum < 0.0d0) call Error_Report('look-back time exceeds the age of the universe'//{introspection:location})
       end if
    else if (parameters%isPresent('redshiftMaximum')) then
       !![
       <inputParameter>
	 <name>redshiftMaximum</name>
	 <description>The maximum redshift at which to output.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
       timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    else
       call Error_Report("must specify either 'timeMinimum' or 'redshiftMaximum'"//{introspection:location})
    end if
    if (parameters%isPresent('timeMaximum')) then
       if (parameters%isPresent('redshiftMinimum')) call Error_Report("can not specify both 'timeMaximum' and 'redshiftMinimum'"//{introspection:location})
       !![
       <inputParameter>
	 <name>timeMaximum</name>
	 <description>The maximum time at which to output. Negative times are interpreted as look-back times.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
       if (timeMaximum < 0.0d0) then
          timeMaximum=cosmologyFunctions_%cosmicTime(1.0d0)+timeMaximum
          if (timeMaximum < 0.0d0) call Error_Report('look-back time exceeds the age of the universe'//{introspection:location})
       end if
    else if (parameters%isPresent('redshiftMinimum')) then
       !![
       <inputParameter>
	 <name>redshiftMinimum</name>
	 <description>The minimum redshift at which to output.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
       timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    else
       call Error_Report("must specify either 'timeMinimum' or 'redshiftMaximum'"//{introspection:location})
    end if
    !![
    <inputParameter>
      <name>countTimes</name>
      <description>The number of times at which to output.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=outputTimesUniformSpacingInTime(timeMinimum,timeMaximum,countTimes,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function uniformSpacingInTimeConstructorParameters

  function uniformSpacingInTimeConstructorInternal(timeMinimum,timeMaximum,countTimes,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{outputTimesUniformSpacingInTime} output times class.
    !!}
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLinear
    use :: Error           , only : Error_Report
    implicit none
    type            (outputTimesUniformSpacingInTime)                        :: self
    double precision                                 , intent(in   )         :: timeMinimum        , timeMaximum
    integer         (c_size_t                       ), intent(in   )         :: countTimes
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    integer         (c_size_t                       )                        :: i
    !![
    <constructorAssign variables="timeMinimum, timeMaximum, countTimes, *cosmologyFunctions_"/>
    !!]

    if (timeMinimum >= timeMaximum) call Error_Report('`timeMinimum` is after `timeMaximum`'//{introspection:location})
    allocate(self%times    (countTimes))
    allocate(self%redshifts(countTimes))
    self%times=Make_Range(timeMinimum,timeMaximum,int(countTimes),rangeTypeLinear)
    do i=1,countTimes
       self%redshifts(i)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%times(i)))
    end do
    self%redshiftMinimum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%timeMaximum))
    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%timeMinimum))
    return
  end function uniformSpacingInTimeConstructorInternal
