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
  Contains an output times class which uses a set of times spaced uniformly in $\log(1+z)$.
  !!}
  
  !![
  <outputTimes name="outputTimesLogarithmicSpacingInRedshift">
   <description>An output times class which uses a set of times spaced uniformly in $\log(1+z)$.</description>
  </outputTimes>
  !!]
  type, extends(outputTimesList) :: outputTimesLogarithmicSpacingInRedshift
     !!{
     Implementation of an output times class which uses a set of times spaced uniformly in $\log(1+z)$
     !!}
     private
     double precision           :: redshiftMinimum, redshiftMaximum
     integer         (c_size_t) :: countRedshifts
  end type outputTimesLogarithmicSpacingInRedshift

  interface outputTimesLogarithmicSpacingInRedshift
     !!{
     Constructors for the \refClass{outputTimesLogarithmicSpacingInRedshift} output times class.
     !!}
     module procedure logarithmicSpacingInRedshiftConstructorParameters
     module procedure logarithmicSpacingInRedshiftConstructorInternal
  end interface outputTimesLogarithmicSpacingInRedshift

contains

  function logarithmicSpacingInRedshiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesLogarithmicSpacingInRedshift} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputTimesLogarithmicSpacingInRedshift)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    double precision                                                         :: redshiftMinimum    , redshiftMaximum
    integer         (c_size_t                               )                :: countRedshifts

    !![
    <inputParameter>
      <name>redshiftMinimum</name>
      <description>The minimum redshift at which to output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <description>The maximum redshift at which to output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countRedshifts</name>
      <description>The number of redshifts at which to output.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=outputTimesLogarithmicSpacingInRedshift(redshiftMinimum,redshiftMaximum,countRedshifts,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function logarithmicSpacingInRedshiftConstructorParameters

  function logarithmicSpacingInRedshiftConstructorInternal(redshiftMinimum,redshiftMaximum,countRedshifts,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \refClass{outputTimesLogarithmicSpacingInRedshift} output times class which takes a parameter set as input.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (outputTimesLogarithmicSpacingInRedshift)                        :: self
    double precision                                         , intent(in   )         :: redshiftMinimum    , redshiftMaximum
    integer         (c_size_t                               ), intent(in   )         :: countRedshifts
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    integer         (c_size_t                               )                        :: i
    !![
    <constructorAssign variables="redshiftMinimum,redshiftMaximum,countRedshifts, *cosmologyFunctions_"/>
    !!]

    allocate(self%times    (countRedshifts))
    allocate(self%redshifts(countRedshifts))
    self%redshifts=Make_Range(1.0d0+redshiftMaximum,1.0d0+redshiftMinimum,int(countRedshifts),rangeTypeLogarithmic)-1.0d0
    do i=1,countRedshifts
       self%times(i)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshifts(i)))
    end do
    return
  end function logarithmicSpacingInRedshiftConstructorInternal
