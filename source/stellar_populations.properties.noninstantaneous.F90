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
  Implements a stellar population properties class based on the noninstantaneous recycling approximation.
  !!}

  use :: Output_Times                              , only : outputTimesClass
  use :: Stellar_Population_Selectors              , only : stellarPopulationSelectorClass
  use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass

  !![
  <stellarPopulationProperties name="stellarPopulationPropertiesNoninstantaneous">
   <description>
    A stellar population properties class based on the noninstantaneous recycling approximation---fully non-instantaneous
    recycling and metal enrichment are used. Recycling and metal production rates from simple stellar populations are computed,
    for any given \gls{imf}, from stellar evolution models. The rates of change are then:
    \begin{eqnarray}
     \dot{M}_\star &amp;=&amp; \phi - \int_0^t \phi(t^\prime) \dot{R}(t-t^\prime;Z_\mathrm{fuel}[t^\prime]) \d t^\prime, \\
     \dot{M}_\mathrm{fuel} &amp;=&amp; -\phi + \int_0^t \phi(t^\prime) \dot{R}(t-t^\prime;Z_\mathrm{fuel}[t]) \d t^\prime, \\
     \dot{M}_{\star,Z} &amp;=&amp; Z_\mathrm{fuel} \phi - \int_0^t \phi(t^\prime) Z_\mathrm{fuel}(t^\prime)
     \dot{R}(t-t^\prime;Z_\mathrm{fuel}[t^\prime]) \d t^\prime, \\
     \dot{M}_{\mathrm{fuel},Z} &amp;=&amp; -Z_\mathrm{fuel} \phi + \int_0^t \phi(t^\prime) \{ Z_\mathrm{fuel}(t^\prime)
     \dot{R}(t-t^\prime;Z_\mathrm{fuel}[t^\prime]) + \dot{p}(t-t^\prime;Z_\mathrm{fuel}[t^\prime]) \} \d t^\prime, \\
    \end{eqnarray}
    where $\dot{R}(t;Z)$ and $\dot{p}(t;Z)$ are the recycling and metal yield rates respectively from a stellar population of
    age $t$ and metallicity $Z$. The energy input rate is computed self-consistently from the star formation history.
   </description>
  </stellarPopulationProperties>
  !!]
  type, extends(stellarPopulationPropertiesClass) :: stellarPopulationPropertiesNoninstantaneous
     !!{
     A stellar population properties class based on the noninstantaneous recycling approximation.
     !!}
     private
     class (stellarPopulationSelectorClass             ), pointer :: stellarPopulationSelector_              => null()
     class (stellarPopulationBroadBandLuminositiesClass), pointer :: stellarPopulationBroadBandLuminosities_ => null()
     class (outputTimesClass                           ), pointer :: outputTimes_                            => null()
     ! Count of number of elements (plus total metals) that are to be tracked.
     integer                                                      :: elementsCount
     ! Count of the number of histories required by this implementation.
     integer                                                      :: historyCount_
     ! Indices for histories.
     integer                                                      :: recycledRateIndex                                , rateEnergyInputIndex     , &
          &                                                          returnedMetalRateBeginIndex                      , metalYieldRateBeginIndex , &
          &                                                          metalYieldRateEndIndex                           , returnedMetalRateEndIndex
     ! Number of times to store in histories.
     integer                                                      :: countHistoryTimes
   contains
     final     ::                  noninstantaneousDestructor
     procedure :: rates         => noninstantaneousRates
     procedure :: scales        => noninstantaneousScales
     procedure :: historyCount  => noninstantaneousHistoryCount
     procedure :: historyCreate => noninstantaneousHistoryCreate
  end type stellarPopulationPropertiesNoninstantaneous

  interface stellarPopulationPropertiesNoninstantaneous
     !!{
     Constructors for the \refClass{stellarPopulationPropertiesNoninstantaneous} stellar population class.
     !!}
     module procedure noninstantaneousConstructorParameters
     module procedure noninstantaneousConstructorInternal
  end interface stellarPopulationPropertiesNoninstantaneous

contains

  function noninstantaneousConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationPropertiesNoninstantaneous} stellar population properties class which takes a parameter list
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (stellarPopulationPropertiesNoninstantaneous)                :: self
    type   (inputParameters                            ), intent(inout) :: parameters
    class  (stellarPopulationSelectorClass             ), pointer       :: stellarPopulationSelector_
    class  (stellarPopulationBroadBandLuminositiesClass), pointer       :: stellarPopulationBroadBandLuminosities_
    class  (outputTimesClass                           ), pointer       :: outputTimes_
    integer                                                             :: countHistoryTimes

    !![
    <inputParameter>
      <name>countHistoryTimes</name>
      <defaultValue>10</defaultValue>
      <description>The number of times at which a galaxy's stellar properties history is stored.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="outputTimes"                            name="outputTimes_"                            source="parameters"/>
    <objectBuilder class="stellarPopulationSelector"              name="stellarPopulationSelector_"              source="parameters"/>
    <objectBuilder class="stellarPopulationBroadBandLuminosities" name="stellarPopulationBroadBandLuminosities_" source="parameters"/>
    !!]
    self=stellarPopulationPropertiesNoninstantaneous(countHistoryTimes,stellarPopulationSelector_,stellarPopulationBroadBandLuminosities_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                           />
    <objectDestructor name="stellarPopulationSelector_"             />
    <objectDestructor name="stellarPopulationBroadBandLuminosities_"/>
    !!]
   return
  end function noninstantaneousConstructorParameters

  function noninstantaneousConstructorInternal(countHistoryTimes,stellarPopulationSelector_,stellarPopulationBroadBandLuminosities_,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationPropertiesNoninstantaneous} stellar population properties class.
    !!}
    use :: Abundances_Structure, only : Abundances_Property_Count
    implicit none
    type   (stellarPopulationPropertiesNoninstantaneous)                        :: self
    integer                                             , intent(in   )         :: countHistoryTimes
    class  (stellarPopulationSelectorClass             ), intent(in   ), target :: stellarPopulationSelector_
    class  (stellarPopulationBroadBandLuminositiesClass), intent(in   ), target :: stellarPopulationBroadBandLuminosities_
    class  (outputTimesClass                           ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="countHistoryTimes, *stellarPopulationSelector_, *stellarPopulationBroadBandLuminosities_, *outputTimes_"/>
    !!]

    ! Get a count of the number of elements (plus total metals) that will be tracked.
    self%elementsCount=Abundances_Property_Count()
    ! Determine the number of histories that we must store. We require one for recycled mass, one for energy input and two
    ! (return rate and yield rate) for each element.
    self%historyCount_=2+2*self%elementsCount
    ! Establish indices in the array of histories where returned metal rates and yield rates will be stored.
    self%recycledRateIndex          =                                                   +1
    self%rateEnergyInputIndex       =self%recycledRateIndex                             +1
    self%returnedMetalRateBeginIndex=self%rateEnergyInputIndex                          +1
    self%returnedMetalRateEndIndex  =self%returnedMetalRateBeginIndex+self%elementsCount-1
    self%metalYieldRateBeginIndex   =self%returnedMetalRateEndIndex                     +1
    self%metalYieldRateEndIndex     =self%metalYieldRateBeginIndex   +self%elementsCount-1
    return
  end function noninstantaneousConstructorInternal

  subroutine noninstantaneousDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationPropertiesNoninstantaneous} stellar population properties class.
    !!}
    implicit none
    type(stellarPopulationPropertiesNoninstantaneous), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulationSelector_"             />
    <objectDestructor name="self%stellarPopulationBroadBandLuminosities_"/>
    <objectDestructor name="self%outputTimes_"                           />
    !!]
    return
  end subroutine noninstantaneousDestructor

  integer function noninstantaneousHistoryCount(self)
    !!{
    Returns the number of histories required by the noninstantaneous stellar populations properties class.
    !!}
    implicit none
    class(stellarPopulationPropertiesNoninstantaneous), intent(inout) :: self

    noninstantaneousHistoryCount=self%historyCount_
    return
  end function noninstantaneousHistoryCount

  subroutine noninstantaneousRates(self,rateStarFormation,abundancesFuel,component,node,history_,rateMassStellar,rateMassFuel,rateEnergyInput&
       &,rateAbundancesFuel,rateAbundancesStellar,rateLuminosityStellar,computeRateLuminosityStellar)
    !!{
    Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    !!}
    use            :: Abundances_Structure          , only : zeroAbundances
    use            :: Galacticus_Nodes              , only : nodeComponent         , nodeComponentBasic , treeNode
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Numerical_Interpolation       , only : interpolator
    use            :: Stellar_Luminosities_Structure, only : max                   , stellarLuminosities, zeroStellarLuminosities
    use            :: Stellar_Populations           , only : stellarPopulationClass
    implicit none
    class           (stellarPopulationPropertiesNoninstantaneous), intent(inout)                 :: self
    double precision                                             , intent(  out)                 :: rateEnergyInput              , rateMassFuel          , &
         &                                                                                          rateMassStellar
    type            (abundances                                 ), intent(inout)                 :: rateAbundancesFuel           , rateAbundancesStellar
    type            (stellarLuminosities                        ), intent(inout)                 :: rateLuminosityStellar
    double precision                                             , intent(in   )                 :: rateStarFormation
    type            (abundances                                 ), intent(in   )                 :: abundancesFuel
    class           (nodeComponent                              ), intent(in   )                 :: component
    type            (treeNode                                   ), intent(inout)                 :: node
    type            (history                                    ), intent(inout)                 :: history_
    logical                                                      , intent(in   )                 :: computeRateLuminosityStellar
    class           (nodeComponentBasic                         ), pointer                       :: basic
    class           (stellarPopulationClass                     ), pointer                       :: stellarPopulation_
    double precision                                             , dimension(self%elementsCount) :: fuelMetallicity              , fuelMetalsRateOfChange, &
         &                                                                                          metalReturnRate              , metalYieldRate        , &
         &                                                                                          stellarMetalsRateOfChange
    integer                                                                                      :: iElement
    integer         (c_size_t                                   )                                :: iHistory
    double precision                                                                             :: ageMaximum                   , ageMinimum            , &
         &                                                                                          currentTime                  , recyclingRate
    type            (interpolator                               )                                :: interpolator_

    ! If a history exists, compute rates.
    if (history_%exists()) then
       ! Get the current time.
       basic       => node %basic()
       currentTime =  basic%time ()
       ! Get interpolating factors in stellar population history.
       interpolator_=interpolator(history_%time)
       iHistory          =interpolator_%locate(currentTime)
       ! Get recycling, energy input, metal recycling and metal yield rates.
       recyclingRate  =history_%data(iHistory,self%          recycledRateIndex                               )
       rateEnergyInput=history_%data(iHistory,self%       rateEnergyInputIndex                               )
       metalReturnRate=history_%data(iHistory,self%returnedMetalRateBeginIndex:self%returnedMetalRateEndIndex)
       metalYieldRate =history_%data(iHistory,self%   metalYieldRateBeginIndex:self%   metalYieldRateEndIndex)
       ! Get the metallicity of the fuel supply.
       call abundancesFuel%serialize(fuelMetallicity)
       ! Set the stellar and fuel mass rates of change.
       rateMassStellar=+rateStarFormation-recyclingRate
       rateMassFuel   =-rateMassStellar
       ! Set the rates of change of the stellar and fuel metallicities.
       stellarMetalsRateOfChange=rateStarFormation*fuelMetallicity-metalReturnRate
       fuelMetalsRateOfChange   =-stellarMetalsRateOfChange+metalYieldRate
       call rateAbundancesStellar%deserialize(stellarMetalsRateOfChange)
       call rateAbundancesFuel   %deserialize(   fuelMetalsRateOfChange)
       ! Get the stellar population.
       stellarPopulation_ => self%stellarPopulationSelector_%select(rateStarFormation,abundancesFuel,component)
       ! Set luminosity rates of change.
       if (computeRateLuminosityStellar) call rateLuminosityStellar%setLuminosities(rateStarFormation,stellarPopulation_,self%stellarPopulationBroadBandLuminosities_,currentTime,abundancesFuel)
       ! Set rates of change in the stellar populations properties future history.
       do iHistory=1,size(history_%time)-1
          ! Find the age of the forming stellar population at the future time. We average over the time between successive timesteps
          ! to ensure that the rates will integrate to the correct values.
          ageMinimum=max(history_%time(iHistory  )-currentTime,0.0d0)
          ageMaximum=max(history_%time(iHistory+1)-currentTime,0.0d0)
          ! Check that it really is in the future and that the timestep over which the contribution to be made is non-zero.
          if (ageMaximum >= 0.0d0 .and. ageMaximum > ageMinimum) then
             ! Get the recycling rate.
             recyclingRate                                                                                =+stellarPopulation_%rateRecycling(abundancesFuel,ageMinimum,ageMaximum         ) &
                  &                                                                                        *rateStarFormation
             ! Accumulate the mass recycling rate from this population at the future time.
             history_      %data(iHistory,self%recycledRateIndex                                         )=+recyclingRate
             ! Get the (normalized) energy input rate.
             history_      %data(iHistory,self%rateEnergyInputIndex                                      )=+stellarPopulation_%rateEnergy   (abundancesFuel,ageMinimum,ageMaximum         ) &
                  &                                                                                        *rateStarFormation
             ! Accumulate the metal return rate from this population at the future time.
             history_      %data(iHistory,self%returnedMetalRateBeginIndex:self%returnedMetalRateEndIndex)=+recyclingRate   &
                  &                                                                                        *fuelMetallicity
             ! Loop over all elements (and total metallicity).
             do iElement=1,self%elementsCount
                ! Get the metal yield rate.
                history_   %data(iHistory,self%metalYieldRateBeginIndex+iElement-1                       )=+stellarPopulation_%rateYield    (abundancesFuel,ageMinimum,ageMaximum,iElement) &
                     &                                                                                     *rateStarFormation
             end do
          end if
       end do
    else
       ! No history exists - rates must all be zero.
       rateEnergyInput      =0.0d0
       rateMassFuel         =0.0d0
       rateMassStellar      =0.0d0
       rateAbundancesFuel   =zeroAbundances
       rateAbundancesStellar=zeroAbundances
       rateLuminosityStellar=zeroStellarLuminosities
    end if
    return
  end subroutine noninstantaneousRates

  subroutine noninstantaneousScales(self,massStellar,abundancesStellar,history_)
    !!{
    Set the scalings for error control on the absolute values of stellar population properties.
    !!}
    use :: Stellar_Feedback , only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (stellarPopulationPropertiesNoninstantaneous), intent(inout)                                :: self
    double precision                                             , intent(in   )                                :: massStellar
    type            (abundances                                 ), intent(in   )                                :: abundancesStellar
    type            (history                                    ), intent(inout)                                :: history_
    double precision                                             , parameter                                    :: massStellarMinimum       =1.0d0
    double precision                                             , parameter                                    :: abundancesStellarMinimum =1.0d0
    double precision                                                            , dimension(self%elementsCount) :: abundancesStellarUnpacked
    double precision                                             , allocatable  , dimension(:                 ) :: timeSteps
    integer                                                                                                     :: scaleIndex

    ! Get timesteps.
    call history_         %timeSteps(timeSteps                )
    ! Get abundances.
    call abundancesStellar%serialize(abundancesStellarUnpacked)
    ! Set scaling factors for recycled mass.
    history_   %data(:,self%recycledRateIndex                       )=max(massStellar                          ,massStellarMinimum      )                                       /timeSteps
    ! Set scaling factors for metal recycling rates.
    forall(scaleIndex=1:self%elementsCount)
       history_%data(:,self%returnedMetalRateBeginIndex-1+scaleIndex)=max(abundancesStellarUnpacked(scaleIndex),abundancesStellarMinimum)                                       /timeSteps
    end forall
    ! Set scaling factors for metal yield rates.
    forall(scaleIndex=1:self%elementsCount)
       history_%data(:,self%metalYieldRateBeginIndex   -1+scaleIndex)=max(abundancesStellarUnpacked(scaleIndex),abundancesStellarMinimum)                                       /timeSteps
    end forall
    ! Set scaling factors for energy input rates.
    history_   %data(:,self%rateEnergyInputIndex                    )=max(massStellar                          ,massStellarMinimum      )*feedbackEnergyInputAtInfinityCanonical/timeSteps
    ! Destroy temporary array.
    deallocate(timeSteps)
    return
  end subroutine noninstantaneousScales

  subroutine noninstantaneousHistoryCreate(self,node,history_)
    !!{
    Create any history required for storing stellar population properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic  , treeNode
    use :: Numerical_Ranges, only : rangeTypeLogarithmic
    implicit none
    class           (stellarPopulationPropertiesNoninstantaneous), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: node
    type            (history                                    ), intent(inout) :: history_
    class           (nodeComponentBasic                         ), pointer       :: basic
    double precision                                                             :: timeBegin, timeEnd

    ! Decide on start and end times for the history.
    basic     => node %basic()
    timeBegin =  basic%time ()
    timeEnd   =  self%outputTimes_%time(self%outputTimes_%count())
    ! Create the history.
    call history_%create(self%historyCount_,self%countHistoryTimes,timeBegin,timeEnd,rangeTypeLogarithmic)
    return
  end subroutine noninstantaneousHistoryCreate
