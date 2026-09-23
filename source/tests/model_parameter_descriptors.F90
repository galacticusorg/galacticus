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
Contains a program which tests descriptors of model parameters and their priors.
!!}

program Test_Model_Parameter_Descriptors
  !!{RST
  Tests that the descriptor of a model parameter describes its prior and perturber under their own parameter names, that
  the descriptors of normal and log-normal distributions describe them correctly, and that a model parameter rebuilt from
  its descriptor has the same prior.
  !!}
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: Error             , only : Error_Handler_Register
  use :: Events_Hooks      , only : eventsHooksInitialize
  use :: Input_Parameters  , only : inputParameters
  use :: ISO_Varying_String, only : var_str            , varying_string        , char
  use :: Model_Parameters  , only : modelParameterClass
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                              , parameter                 :: countTrials   =9
  type            (inputParameters    )                            :: parameters          , parametersSimulation, &
       &                                                              descriptor          , descriptorParameter , &
       &                                                              descriptorPrior
  class           (modelParameterClass), pointer                   :: modelParameter_     , modelParameterRebuilt_
  type            (varying_string     )                            :: classPrior          , classPerturber
  character       (len=32             ), dimension(3  ), parameter :: classExpected=['logNormal','logNormal','normal   ']
  double precision                                                 :: x0                  , sigma               , &
       &                                                              limitLower          , differenceMaximum   , &
       &                                                              u                   , x
  integer                                                          :: i                   , k

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call eventsHooksInitialize()
  call Unit_Tests_Begin_Group("Model parameter descriptors")
  parameters          =inputParameters(var_str('testSuite/parameters/modelParameterDescriptors.xml'))
  parametersSimulation=parameters%subParameters('posteriorSampleSimulation')
  do i=1,parametersSimulation%copiesCount('modelParameter')
     !![
     <objectBuilder class="modelParameter" name="modelParameter_" source="parametersSimulation" copy="i"/>
     !!]
     call Unit_Tests_Begin_Group(char(modelParameter_%name()))
     descriptor=inputParameters()
     call modelParameter_%descriptor(descriptor,includeClass=.true.)
     descriptorParameter=descriptor%subParameters('modelParameter')
     ! The prior and perturber must be described under their own parameter names, not merged under their class name.
     call Assert('prior described as distributionFunction1DPrior'        ,descriptorParameter%isPresent('distributionFunction1DPrior'    ),.true. )
     call Assert('perturber described as distributionFunction1DPerturber',descriptorParameter%isPresent('distributionFunction1DPerturber'),.true. )
     call Assert('no entry under the class name'                         ,descriptorParameter%isPresent('distributionFunction1D'         ),.false.)
     call Assert('mapper described as operatorUnaryMapper'               ,descriptorParameter%isPresent('operatorUnaryMapper'            ),.true. )
     call descriptorParameter%value('distributionFunction1DPrior'    ,classPrior    )
     call descriptorParameter%value('distributionFunction1DPerturber',classPerturber)
     call Assert('prior class'    ,char(classPrior    ),trim(classExpected(i)))
     call Assert('perturber class',char(classPerturber),'cauchy'               )
     descriptorPrior=descriptorParameter%subParameters('distributionFunction1DPrior')
     select case (i)
     case (1)
        ! Log-normal specified by median and width: these, and the limits, must be described as given.
        call descriptorPrior%value('x0'        ,x0        )
        call descriptorPrior%value('sigma'     ,sigma     )
        call descriptorPrior%value('limitLower',limitLower)
        call Assert('x0'                        ,x0        ,250.0d0,relTol=1.0d-9)
        call Assert('sigma'                     ,sigma     ,  0.5d0,relTol=1.0d-9)
        call Assert('lower limit (not its log)' ,limitLower, 25.0d0,relTol=1.0d-9)
        call Assert('no mean or variance'       ,descriptorPrior%isPresent('mean').or.descriptorPrior%isPresent('variance'),.false.)
     case (2)
        ! Log-normal specified by mean and variance: described by the equivalent median and width, sigma^2 = ln(1+V/M^2) and
        ! x0 = M/sqrt(1+V/M^2).
        call descriptorPrior%value('x0'        ,x0        )
        call descriptorPrior%value('sigma'     ,sigma     )
        call descriptorPrior%value('limitUpper',limitLower)
        call Assert('x0 (from mean and variance)'   ,x0        ,3.0d0/sqrt(1.0d0+2.0d0/3.0d0**2),relTol=1.0d-9)
        call Assert('sigma (from mean and variance)',sigma     ,sqrt(log(1.0d0+2.0d0/3.0d0**2)),relTol=1.0d-9)
        call Assert('upper limit (not its log)'     ,limitLower,20.0d0                         ,relTol=1.0d-9)
     case (3)
        call Assert('no upper limit'                ,descriptorPrior%isPresent('limitUpper'),.false.)
     end select
     ! Rebuild the model parameter from its descriptor, and check that its prior is unchanged. A log-normal descriptor holding an
     ! invalid width can not be rebuilt (the constructor would raise a floating point exception, losing the output of this test),
     ! so is reported as a failure instead.
     if (i <= 2 .and. .not.(sigma > 1.0d-3 .and. sigma < 1.0d3)) then
        call Assert('prior rebuilt from descriptor is unchanged (descriptor width is invalid)',.false.,.true.)
     else
        !![
        <objectBuilder class="modelParameter" name="modelParameterRebuilt_" source="descriptor"/>
        !!]
        differenceMaximum=0.0d0
        do k=1,countTrials
           u                =dble(k)/dble(countTrials+1)
           x                =modelParameter_%priorInvert(u)
           differenceMaximum=max(differenceMaximum,abs(modelParameterRebuilt_%priorInvert(u)/x-1.0d0))
        end do
        call Assert('prior rebuilt from descriptor is unchanged',differenceMaximum,0.0d0,absTol=1.0d-8)
        !![
        <objectDestructor name="modelParameterRebuilt_"/>
        !!]
     end if
     !![
     <objectDestructor name="modelParameter_"/>
     !!]
     call Unit_Tests_End_Group()
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Model_Parameter_Descriptors
