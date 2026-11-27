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
Contains a program which tests parameter input.
!!}

program Test_Parameters
  !!{RST
  Test reading of input parameters.
  !!}
  use :: Cosmological_Density_Field, only : cosmologicalMassVariance, cosmologicalMassVarianceClass
  use :: Cosmology_Parameters      , only : cosmologyParameters     , cosmologyParametersClass
  use :: Cosmology_Functions       , only : cosmologyFunctions      , cosmologyFunctionsClass
  use :: Display                   , only : displayVerbositySet     , verbosityLevelStandard
  use :: IO_HDF5                   , only : hdf5File
  use :: ISO_Varying_String        , only : assignment(=)           , var_str                      , varying_string
  use :: Input_Parameters          , only : inputParameters         , inputParameter
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Kind_Numbers              , only : kind_int8
  use :: Unit_Tests                , only : Assert                  , Unit_Tests_Begin_Group       , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (hdf5File                     )              :: outputFile
  type            (varying_string               )              :: parameterFile            , parameterValue
  class           (cosmologyParametersClass     ), pointer     :: cosmologyParameters_
  class           (cosmologyFunctionsClass      ), pointer     :: cosmologyFunctions_
  class           (cosmologicalMassVarianceClass), pointer     :: cosmologicalMassVariance_
  type            (inputParameters              ), target      :: testParameters
  type            (inputParameters              ), allocatable :: wrapper1                 ,  wrapper2
  type            (inputParameter               ), pointer     :: testParameter
  double precision                                             :: valueNumerical
  integer                                                      :: valueInteger
  integer         (kind_int8                    )              :: valueLong
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call Functions_Global_Set()
  ! Open an output file.
  outputFile=hdf5File("testSuite/outputs/testParameters.hdf5",overWrite=.true.)
  parameterFile  ='testSuite/parameters/testsParameters.xml'
  testParameters=inputParameters(parameterFile,outputParametersGroup=outputFile)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Parameter input")
  ! Test retrieval of cosmology parameters (simple).
  !![
  <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="testParameters"/>
  !!]
  !![
  <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="testParameters"/>
  !!]
  call Unit_Tests_Begin_Group("Retrieve cosmological parameters (simple)")
  call Assert('Ωₘ  ',cosmologyParameters_%OmegaMatter    (), 0.27250d0,relTol=1.0d-6)
  call Assert('Ωb  ',cosmologyParameters_%OmegaBaryon    (), 0.04550d0,relTol=1.0d-6)
  call Assert('ΩΛ  ',cosmologyParameters_%OmegaDarkEnergy(), 0.72750d0,relTol=1.0d-6)
  call Assert('H₀  ',cosmologyParameters_%HubbleConstant (),70.20000d0,relTol=1.0d-6)
  call Assert('TCMB',cosmologyParameters_%temperatureCMB (), 2.72548d0,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Test retrieval of cosmological mass variance through a reference.
  !![
  <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="testParameters"/>
  !!]
  call Unit_Tests_Begin_Group("Parameter referencing")
  call Assert('σ₈ via reference'          ,cosmologicalMassVariance_%sigma8     (                          ),0.912d0,relTol=1.0d-6)
  call Assert('Test presence of reference',testParameters           %isPresent  ('cosmologicalMassVariance'),.true.               )
  call Assert('Test count of references'  ,testParameters           %copiesCount('cosmologicalMassVariance'),1                    )
  call Unit_Tests_End_Group()
  ! Test conditionals.
  call Unit_Tests_Begin_Group("Parameter evaluation")
  call testParameters%value('active1',valueNumerical)
  call Assert('conditional value',valueNumerical,0.2d0,absTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Test adding, retrieving, resetting, reading a parameter.
  call Unit_Tests_Begin_Group("Parameter adding")
  call testParameters%addParameter('addedParameter','qwertyuiop'  )
  call testParameters%value       ('addedParameter',parameterValue)
  call Assert('added parameter exists'                      ,testParameters%isPresent('addedParameter'),.true.               )
  call Assert('added parameter value is correct'            ,parameterValue                            ,var_str('qwertyuiop'))
  call testParameters%reset()
  call Assert('added parameter no longer exists after reset',testParameters%isPresent('addedParameter'),.false.              )
  call testParameters%addParameter('addedParameter','asdfghjkl'   )
  call testParameters%value       ('addedParameter',parameterValue)
  call Assert('re-added parameter exists'                   ,testParameters%isPresent('addedParameter'),.true.               )
  call Assert('re-added parameter value is correct'         ,parameterValue                            ,var_str('asdfghjkl' ))
  call Unit_Tests_End_Group()
  ! Test evaluation.
  call Unit_Tests_Begin_Group("Parameter evaluation")
  call testParameters%value('fixedValue'      ,valueNumerical)
  call Assert('fixed value'                        ,valueNumerical,+1.234000000d0                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue1'   ,valueNumerical)
  call Assert('derived value'                      ,valueNumerical,+1.234000000d1                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue2'   ,valueNumerical)
  call Assert('derived value [recursive]'          ,valueNumerical,-8.264825587d0                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue6'   ,valueNumerical)
  call Assert('derived value [min function]'       ,valueNumerical,+1.000000000d0                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue7'   ,valueNumerical)
  call Assert('derived value [max function]'       ,valueNumerical,+2.000000000d0                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue8'   ,valueNumerical)
  call Assert('derived value [default]'            ,valueNumerical,+1.234500000d3                                       ,absTol=1.0d-6)
  call testParameters%value('derivedValue11'  ,valueNumerical)
  call Assert('derived value [deep]'               ,valueNumerical,-5.264825587d0                                       ,absTol=1.0d-6)
  call testParameters%value('derivedInteger1' ,valueInteger)
  call Assert('derived value [integer]'            ,valueInteger  ,13                                                                 )
  call testParameters%value('derivedInteger2' ,valueInteger)
  call Assert('derived value [integer; default]'   ,valueInteger  ,21                                                                 )
  call testParameters%value('derivedLong1'    ,valueLong   )
  call Assert('derived value [long integer]'       ,valueLong     ,6000000000000_kind_int8                                            )
  call testParameters%value('derivedInteger3' ,valueInteger)
  call Assert('derived value [integer; rounded]'   ,valueInteger  ,29                                                                 )
  call testParameters%value('derivedInteger4',valueInteger)
  call Assert('derived value [integer; truncated]',valueInteger   ,1                                                                  )
  call testParameters%value('derivedText1'   ,parameterValue)
  call Assert('derived value [text]'              ,parameterValue ,var_str('mouse_filteredPower_1.00_0.200000E+01_0006')              )
  allocate(wrapper1)
  allocate(wrapper2)
  wrapper1=testParameters%subParameters('wrapper1')
  wrapper2=wrapper1      %subParameters('wrapper2')
  call wrapper1%value('derivedValue4',valueNumerical)
  call Assert('derived value [relative]'             ,valueNumerical, 33.42d0,absTol=1.0d-6)
  call wrapper2%value('derivedValue5',valueNumerical)
  call Assert('derived value [relative]'             ,valueNumerical,118.08d0,absTol=1.0d-6)
  deallocate(wrapper2)
  deallocate(wrapper1)
  !! Reset the value of the fixed parameter.
  testParameter => testParameters%node('fixedValue')
  call testParameter %set  (2.468d0)
  call testParameters%reset(       )
  call testParameters%value('fixedValue'   ,valueNumerical)
  call Assert('fixed value [post-reset]'             ,valueNumerical,+2.468000000d0,absTol=1.0d-6)
  call testParameters%value('derivedValue1',valueNumerical)
  call Assert('derived value [post-reset]'           ,valueNumerical,+2.468000000d1,absTol=1.0d-6)
  call testParameters%value('derivedValue2',valueNumerical)
  call Assert('derived value [recursive; post-reset]',valueNumerical,-1.344442852d2,absTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Test the model fact registry.
  call Unit_Tests_Begin_Group("Model facts")
  call Assert('fact not provided initially'                        ,testParameters%modelFactIsProvided('test:fact:a'),.false.)
  ! Facts recorded through sub-parameters are recorded for the whole model (the registry is held by the shared document).
  allocate(wrapper1)
  wrapper1=testParameters%subParameters('wrapper1')
  call wrapper1      %modelFactProvide('test:fact:a')
  call wrapper1      %modelFactProvide('test:fact:a')
  call Assert('fact provided through sub-parameters is seen by root',testParameters%modelFactIsProvided('test:fact:a'),.true. )
  call Assert('fact provided through sub-parameters is seen by self',wrapper1      %modelFactIsProvided('test:fact:a'),.true. )
  call Assert('other fact not provided'                            ,testParameters%modelFactIsProvided('test:fact:b'),.false.)
  ! A requirement which is met passes validation (an unmet one would be a fatal error).
  call testParameters%modelFactRequire('test:fact:a','a test requires fact a','provide fact a','test location')
  call testParameters%modelFactRequire('test:fact:a','a test requires fact a','provide fact a','test location')
  call testParameters%modelFactsValidate()
  call Assert('met requirement passes validation'                  ,.true.                                            ,.true. )
  ! Exclusive claims repeated from one location (as by each OpenMP thread), and single claims on distinct facts, pass validation
  ! (two claims on one fact from different locations would be a fatal error).
  call wrapper1      %modelFactExclusive('test:exclusive:a','claimant one','location one')
  call testParameters%modelFactExclusive('test:exclusive:a','claimant one','location one')
  call testParameters%modelFactExclusive('test:exclusive:b','claimant two','location two')
  call testParameters%modelFactsValidate()
  call Assert('repeated and distinct exclusive claims pass validation',.true.                                         ,.true. )
  ! Resetting the parameters (as is done at the end of each task) must not clear the facts, since a later task (e.g. in a
  ! taskMulti) belongs to the same model.
  deallocate(wrapper1)
  call testParameters%reset()
  call Assert('fact survives parameter reset'                      ,testParameters%modelFactIsProvided('test:fact:a'),.true. )
  ! Resetting the facts (as is done between models built from the same document, e.g. by likelihoods) clears facts and
  ! requirements, so that validation of the next model is not affected by the previous one.
  call testParameters%modelFactsReset()
  call Assert('fact cleared by fact reset'                         ,testParameters%modelFactIsProvided('test:fact:a'),.false.)
  call testParameters%modelFactsValidate()
  call Assert('requirement cleared by fact reset'                  ,.true.                                            ,.true. )
  ! Exclusive claims are also cleared, so a claim on the same fact from a different location (as by a different object in the
  ! next model) does not conflict.
  call testParameters%modelFactExclusive('test:exclusive:a','claimant three','location three')
  call testParameters%modelFactsValidate()
  call Assert('exclusive claims cleared by fact reset'             ,.true.                                            ,.true. )
  ! Facts can be recorded again after a reset.
  call testParameters%modelFactProvide('test:fact:a')
  call Assert('fact can be provided again after reset'             ,testParameters%modelFactIsProvided('test:fact:a'),.true. )
  call testParameters%modelFactsReset()
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Close down.
  call testParameters%reset  ()
  call testParameters%destroy()
  !![
  <objectDestructor name="cosmologyParameters_"     />
  <objectDestructor name="cosmologyFunctions_"     />
  <objectDestructor name="cosmologicalMassVariance_"/>
  !!]
end program Test_Parameters
