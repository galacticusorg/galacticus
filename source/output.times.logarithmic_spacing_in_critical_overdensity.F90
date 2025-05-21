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
  Contains an output times class which uses a set of times spaced uniformly in $\log \delta_\mathrm{c}$.
  !!}

  use :: Cosmological_Density_Field, only : criticalOverdensityClass, cosmologicalMassVarianceClass
  
  !![
  <outputTimes name="outputTimesLogarithmicSpacingInCriticalOverdensity">
   <description>An output times class which uses a set of times spaced uniformly in $\log \delta_\mathrm{c}$.</description>
  </outputTimes>
  !!]
  type, extends(outputTimesList) :: outputTimesLogarithmicSpacingInCriticalOverdensity
     !!{
     Implementation of an output times class which uses a set of times spaced uniformly in $\log \delta_\mathrm{c}$
     !!}
     private
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     double precision                                         :: redshiftMinimum                    , redshiftMaximum
     integer         (c_size_t                     )          :: countTimes
   contains
     final :: logarithmicSpacingInCriticalOverdensityDestructor
  end type outputTimesLogarithmicSpacingInCriticalOverdensity

  interface outputTimesLogarithmicSpacingInCriticalOverdensity
     !!{
     Constructors for the \refClass{outputTimesLogarithmicSpacingInCriticalOverdensity} output times class.
     !!}
     module procedure logarithmicSpacingInCriticalOverdensityConstructorParameters
     module procedure logarithmicSpacingInCriticalOverdensityConstructorInternal
  end interface outputTimesLogarithmicSpacingInCriticalOverdensity

contains

  function logarithmicSpacingInCriticalOverdensityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesLogarithmicSpacingInCriticalOverdensity} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputTimesLogarithmicSpacingInCriticalOverdensity)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                          ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), pointer       :: cosmologicalMassVariance_
    double precision                                                                    :: redshiftMinimum          , redshiftMaximum
    integer         (c_size_t                                          )                :: countTimes

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
      <name>countTimes</name>
      <description>The number of times at which to output.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=outputTimesLogarithmicSpacingInCriticalOverdensity(redshiftMinimum,redshiftMaximum,countTimes,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function logarithmicSpacingInCriticalOverdensityConstructorParameters

  function logarithmicSpacingInCriticalOverdensityConstructorInternal(redshiftMinimum,redshiftMaximum,countTimes,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{outputTimesLogarithmicSpacingInCriticalOverdensity} output times class which takes a parameter set as input.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (outputTimesLogarithmicSpacingInCriticalOverdensity)                             :: self
    double precision                                                    , intent(in   )              :: redshiftMinimum                 , redshiftMaximum
    integer         (c_size_t                                          ), intent(in   )              :: countTimes
    class           (cosmologyFunctionsClass                           ), intent(in   ), target      :: cosmologyFunctions_
    class           (criticalOverdensityClass                          ), intent(in   ), target      :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target      :: cosmologicalMassVariance_
    double precision                                                    , dimension(:) , allocatable :: criticalOverdensities
    double precision                                                    , parameter                  :: massLarge                =1.0d15
    double precision                                                                                 :: timeMinimum                     , timeMaximum    , &
         &                                                                                              timeNow
    integer         (c_size_t                               )                                        :: i
    !![
    <constructorAssign variables="redshiftMinimum, redshiftMaximum, countTimes, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    allocate(self%times           (countTimes))
    allocate(self%redshifts       (countTimes))
    allocate(criticalOverdensities(countTimes))
    timeMinimum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeMaximum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    timeNow    =self%cosmologyFunctions_%cosmicTime(1.0d0                                                                )
    criticalOverdensities=Make_Range(                                                                               &
         &                           +self%criticalOverdensity_     %value       (time=timeMinimum               )  &
         &                           /self%cosmologicalMassVariance_%rootVariance(time=timeMinimum,mass=massLarge)  &
         &                           *self%cosmologicalMassVariance_%rootVariance(time=timeNow    ,mass=massLarge), &
         &                           +self%criticalOverdensity_     %value       (time=timeMaximum               )  &
         &                           /self%cosmologicalMassVariance_%rootVariance(time=timeMaximum,mass=massLarge)  &
         &                           *self%cosmologicalMassVariance_%rootVariance(time=timeNow    ,mass=massLarge), &
         &                           int(countTimes)                                                              , &
         &                           rangeTypeLogarithmic                                                           &
         &                          )
    do i=1,countTimes
       self%times    (i)=self%criticalOverdensity_%timeOfCollapse             (                                              criticalOverdensities(i) )
       self%redshifts(i)=self%cosmologyFunctions_ %redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%times                (i)))
    end do
    return
  end function logarithmicSpacingInCriticalOverdensityConstructorInternal

  subroutine logarithmicSpacingInCriticalOverdensityDestructor(self)
    !!{
    Destructor for the \refClass{outputTimesLogarithmicSpacingInCriticalOverdensity} output times class.
    !!}
    implicit none
    type(outputTimesLogarithmicSpacingInCriticalOverdensity), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine logarithmicSpacingInCriticalOverdensityDestructor
