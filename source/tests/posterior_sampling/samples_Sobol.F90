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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests prior quantiles of model parameters and the Sobol posterior samples class.
!!}

program Test_Posterior_Samples_Sobol
  !!{RST
  Tests prior quantiles (``priorCumulative``) of model parameters, and the Sobol posterior samples class.
  !!}
  use :: Display                         , only : displayVerbositySet       , verbosityLevelStandard
  use :: Error                           , only : Error_Handler_Register
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Input_Parameters                , only : inputParameters
  use :: ISO_Varying_String              , only : var_str                   , char
  use :: Model_Parameters                , only : modelParameterClass       , modelParameterList
  use :: Models_Likelihoods_Constants    , only : logImpossible
  use :: Numerical_Quasi_Random_Sequences, only : gsl_qrng_sobol            , quasiRandomNumberGenerator
  use :: Posterior_Sampling_State        , only : posteriorSampleStateSimple
  use :: Numerical_Random_Numbers        , only : randomNumberGeneratorClass
  use :: Posterior_Sampling_State_Samples, only : posteriorSamplesClass     , posteriorSamplesSobol
  use :: Unit_Tests                      , only : Assert                    , Unit_Tests_Begin_Group, Unit_Tests_End_Group, &
       &                                          Unit_Tests_Finish
  implicit none
  integer                                     , parameter                  :: countSamples       =64, countTrials=9
  type            (inputParameters           )                             :: parameters         , parametersSimulation
  class           (randomNumberGeneratorClass), pointer                    :: randomNumberGenerator_
  class           (modelParameterClass       ), pointer                    :: modelParameter_
  class           (posteriorSamplesClass     ), pointer                    :: posteriorSamples_
  type            (modelParameterList        ), allocatable, dimension(:  ) :: modelParameters_
  type            (posteriorSampleStateSimple), allocatable, dimension(:  ) :: simulationStates
  type            (quasiRandomNumberGenerator)                             :: quasiRandomSequence
  double precision                            , allocatable, dimension(:,:) :: values             , quantiles
  double precision                            , allocatable, dimension(:  ) :: sequence
  integer                                     , allocatable, dimension(:  ) :: countPerBin
  double precision                                                          :: x                  , u          , &
       &                                                                       roundTripErrorMaximum
  integer                                                                   :: i                  , j          , &
       &                                                                       k                  , iDesign    , &
       &                                                                       countParameters
  logical                                                                   :: isBalanced         , isInPrior  , &
       &                                                                       isUnique

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call eventsHooksInitialize()
  call Unit_Tests_Begin_Group("Posterior samples: Sobol designs")
  parameters          =inputParameters(var_str('testSuite/parameters/posteriorSamplesSobol.xml'))
  parametersSimulation=parameters%subParameters('posteriorSampleSimulation')
  !![
  <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
  !!]
  ! Build the list of active model parameters.
  countParameters=parametersSimulation%copiesCount("modelParameter")
  allocate(modelParameters_(countParameters))
  do i=1,countParameters
     !![
     <objectBuilder class="modelParameter" name="modelParameter_" source="parametersSimulation" copy="i" />
     !!]
     modelParameters_(i)%modelParameter_ => modelParameter_
     !![
     <referenceCountIncrement owner="modelParameters_(i)" object="modelParameter_"/>
     <objectDestructor name="modelParameter_"/>
     !!]
  end do
  ! Prior quantiles: the cumulative probability must invert `priorInvert`, for each class of prior. Trial quantiles are
  ! spread across (0,1), including values close to each limit.
  call Unit_Tests_Begin_Group("Prior quantiles")
  do j=1,countParameters
     roundTripErrorMaximum=0.0d0
     do k=1,countTrials
        u                    =dble(k)/dble(countTrials+1)
        x                    =modelParameters_(j)%modelParameter_%priorInvert    (u)
        roundTripErrorMaximum=max(roundTripErrorMaximum,abs(modelParameters_(j)%modelParameter_%priorCumulative(x)-u))
     end do
     call Assert('priorCumulative(priorInvert(u)) = u: '//char(modelParameters_(j)%modelParameter_%name()),roundTripErrorMaximum,0.0d0,absTol=1.0d-10)
     call Assert('priorCumulative at prior minimum = 0: '//char(modelParameters_(j)%modelParameter_%name()),modelParameters_(j)%modelParameter_%priorCumulative(modelParameters_(j)%modelParameter_%priorMinimum()),0.0d0,absTol=1.0d-10)
  end do
  call Unit_Tests_End_Group()
  ! Designs: the first (randomly shifted) is built from the parameter file, the second (unshifted) directly.
  do iDesign=1,2
     select case (iDesign)
     case (1)
        !![
        <objectBuilder class="posteriorSamples" name="posteriorSamples_" source="parametersSimulation"/>
        !!]
     case (2)
        allocate(posteriorSamplesSobol :: posteriorSamples_)
        select type (posteriorSamples_)
        type is (posteriorSamplesSobol)
           !![
           <referenceConstruct object="posteriorSamples_" constructor="posteriorSamplesSobol(countSamples,.false.,randomNumberGenerator_)"/>
           !!]
        end select
     end select
     call posteriorSamples_%samples(simulationStates,modelParameters_)
     ! Recover the physical values and quantiles of each point.
     allocate(values   (countParameters,size(simulationStates)))
     allocate(quantiles(countParameters,size(simulationStates)))
     do i=1,size(simulationStates)
        values(:,i)=simulationStates(i)%get()
        do j=1,countParameters
           values   (j,i)=modelParameters_(j)%modelParameter_%unmap          (values(j,i))
           quantiles(j,i)=modelParameters_(j)%modelParameter_%priorCumulative(values(j,i))
        end do
     end do
     ! All points must lie within the support of each prior. (Support is tested through the prior density, since
     ! `priorMaximum` is undefined for a prior with no upper limit.)
     isInPrior=.true.
     do i=1,size(simulationStates)
        do j=1,countParameters
           if (modelParameters_(j)%modelParameter_%logPrior(values(j,i)) <= logImpossible) isInPrior=.false.
        end do
     end do
     ! Determine if the design is balanced - every one-dimensional projection has exactly one point in each interval
     ! [k/N,(k+1)/N).
     allocate(countPerBin(countSamples))
     isBalanced=.true.
     do j=1,countParameters
        countPerBin=0
        do i=1,countSamples
           k             =min(int(quantiles(j,i)*dble(countSamples))+1,countSamples)
           countPerBin(k)=countPerBin(k)+1
        end do
        if (any(countPerBin /= 1)) isBalanced=.false.
     end do
     deallocate(countPerBin)
     select case (iDesign)
     case (1)
        call Unit_Tests_Begin_Group("Randomly-shifted design")
        call Assert('number of points'                          ,size(simulationStates),countSamples)
        call Assert('points lie within the prior'               ,isInPrior             ,.true.      )
        call Assert('design is balanced'                        ,isBalanced            ,.true.      )
        ! Every point must be distinct from the unshifted sequence, confirming that a shift was applied.
        isUnique=.true.
        do j=1,countParameters
           if (any(abs(quantiles(j,:)-nint(quantiles(j,:)*dble(countSamples))/dble(countSamples)) < 1.0d-9)) isUnique=.false.
        end do
        call Assert('design is shifted off the unshifted grid'  ,isUnique              ,.true.      )
        call Unit_Tests_End_Group()
     case (2)
        call Unit_Tests_Begin_Group("Unshifted design")
        call Assert('number of points'                          ,size(simulationStates),countSamples)
        call Assert('points lie within the prior'               ,isInPrior             ,.true.      )
        ! The quantiles must be the GSL Sobol sequence itself (which omits the origin).
        allocate(sequence(countParameters))
        quasiRandomSequence=quasiRandomNumberGenerator(gsl_qrng_sobol,countDimensions=countParameters)
        roundTripErrorMaximum=0.0d0
        do i=1,countSamples
           call quasiRandomSequence%getVector(sequence)
           roundTripErrorMaximum=max(roundTripErrorMaximum,maxval(abs(quantiles(:,i)-sequence)))
        end do
        deallocate(sequence)
        call Assert('quantiles match the Sobol sequence'        ,roundTripErrorMaximum ,0.0d0       ,absTol=1.0d-10)
        call Unit_Tests_End_Group()
     end select
     deallocate(values          )
     deallocate(quantiles       )
     deallocate(simulationStates)
     !![
     <objectDestructor name="posteriorSamples_"/>
     !!]
  end do
  ! Clean up.
  !![
  <objectDestructor name="randomNumberGenerator_"/>
  !!]
  do i=1,countParameters
     !![
     <objectDestructor name="modelParameters_(i)%modelParameter_"/>
     !!]
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Posterior_Samples_Sobol
