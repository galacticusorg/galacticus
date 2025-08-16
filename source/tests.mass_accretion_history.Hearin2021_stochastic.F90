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
Contains a program which tests the Hearin (2021) stochastic halo mass formation history.
!!}

program Test_Hearin2021_Stochastic_MAH
  !!{
  Tests the Hearin (2021) halo mass formation history algorithm.
  !!}
  use, intrinsic :: ISO_C_Binding                            , only : c_long
  use            :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryHearin2021Stochastic
  use            :: Display                                  , only : displayVerbositySet                                   , verbosityLevelStandard
  use            :: Events_Hooks                             , only : eventsHooksInitialize
  use            :: File_Utilities                           , only : Count_Lines_in_File
  use            :: Functions_Global_Utilities               , only : Functions_Global_Set
  use            :: Input_Paths                              , only : inputPath                                             , pathTypeExec
  use            :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize                          , nodeComponentBasic               , treeNode , mergerTree
  use            :: Input_Parameters                         , only : inputParameters
  use            :: ISO_Varying_String                       , only : char
  use            :: Node_Components                          , only : Node_Components_Initialize                            , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Numerical_Random_Numbers                 , only : randomNumberGeneratorGSL
  use            :: Unit_Tests                               , only : Assert                                                , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (mergerTree                                            ), target                          :: treeEarly                                      , treeLate
  class           (nodeComponentBasic                                    ), pointer                         :: basicEarly                                     , basicLate
  type            (darkMatterHaloMassAccretionHistoryHearin2021Stochastic)                                  :: darkMatterHaloMassAccretionHistoryEarlyForming_, darkMatterHaloMassAccretionHistoryLateForming_
  double precision                                                        , dimension(:      ), allocatable :: massLateTarget                                 , massLateRecovered                             , &
       &                                                                                                       massEarlyTarget                                , massEarlyRecovered
  double precision                                                                                          :: time                                           , timeLogarithmic
  integer                                                                                                   :: i                                              , referenceFile                                 , &
       &                                                                                                       countTimes
  type            (inputParameters                                       )                                  :: parameters
  character       (len=1024                                              )                                  :: line
  double precision                                                        , dimension(2,2,3  ), parameter   :: means  =reshape(                                    &
       &                                                                                                                       [                                   &
       !                                                                                                                       Each row is the mean of the stated parameter for:
       !                                                                                                                         early forming/low mass : early forming/high mass : late forming/low mass : late forming/high mass
       !                                                                                                                       u_early
       &                                                                                                                       +0.70d0, +3.50d0, +0.50d0, +2.81d0, &
       !                                                                                                                       u_late
       &                                                                                                                       -0.40d0, -0.40d0, -3.05d0, -1.65d0, &
       !                                                                                                                       log₁₀(t₀)
       &                                                                                                                       -0.39d0, +0.91d0, +0.16d0, +1.90d0  &
       &                                                                                                                       ]                                 , &
       &                                                                                                                       [2,2,3]                             &
       &                                                                                                                      )
  double precision                                                        , dimension(2,2,3,3), parameter  :: cholesky=reshape(                                                         &
       &                                                                                                                       [                                                        &
       &                                                                                                                        -huge(0.0d0), -huge(0.0d0), -huge(0.0d0), -huge(0.0d0), &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                        -huge(0.0d0), -huge(0.0d0), -huge(0.0d0), -huge(0.0d0), &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                              0.0d0 ,       0.0d0 ,       0.0d0 ,       0.0d0 , &
       &                                                                                                                        -huge(0.0d0), -huge(0.0d0), -huge(0.0d0), -huge(0.0d0)  &
       &                                                                                                                       ]                                                      , &
       &                                                                                                                       [2,2,3,3]                                                &
       &                                                                                                                      )
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Hearin 2021 mass accretion history algorithm")
  ! Test Correa et al. 2015 algorithm.
  parameters=inputParameters()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Construct early-forming tree.
  treeEarly %index                     =  1
  treeEarly %initializedUntil          =  0.0d0
  treeEarly %isTreeInitialized         =  .false.
  treeEarly %event                     => null()
  treeEarly %firstTree                 => treeEarly
  treeEarly %nodeBase                  => treeNode                (                 )
  treeEarly %nodeBase        %hostTree => treeEarly
  basicEarly                           => treeEarly%nodeBase%basic(autoCreate=.true.)
  allocate(randomNumberGeneratorGSL :: treeEarly%randomNumberGenerator_)
  select type (randomNumberGenerator_ => treeEarly%randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=8322_c_long)
  end select
  call treeEarly%randomNumberGenerator_%seedSet   (seed=treeEarly%index,offset=.true.)
  call treeEarly%properties            %initialize(                                  )
  ! Construct late-forming tree.
  treeLate  %index=2
  treeLate  %initializedUntil          =  0.0d0
  treeLate  %isTreeInitialized         =  .false.
  treeLate  %event                     => null()
  treeLate  %firstTree                 => treeLate
  treeLate  %nodeBase                  => treeNode                (                 )
  treeLate  %nodeBase        %hostTree => treeLate
  basicLate                            => treeLate %nodeBase%basic(autoCreate=.true.)
  allocate(randomNumberGeneratorGSL :: treeLate %randomNumberGenerator_)
  select type (randomNumberGenerator_ => treeLate %randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=8322_c_long)
  end select
  call treeLate %randomNumberGenerator_%seedSet   (seed=treeLate %index,offset=.true.)
  call treeLate %properties            %initialize(                                  )
  ! Build mass accretion history objects.
  !![
  <referenceConstruct object="darkMatterHaloMassAccretionHistoryLateForming_" >
   <constructor>
    darkMatterHaloMassAccretionHistoryHearin2021Stochastic(                           &amp;
     &amp;                                                 rateRollOver    =3.5d0   , &amp;
     &amp;                                                 fractionLateLow =1.0d0   , &amp;
     &amp;                                                 fractionLateHigh=1.0d0   , &amp;
     &amp;                                                 means           =means   , &amp;
     &amp;                                                 cholesky        =cholesky  &amp;
     &amp;                                                )
   </constructor>
  </referenceConstruct> 
  <referenceConstruct object="darkMatterHaloMassAccretionHistoryEarlyForming_">
   <constructor>
    darkMatterHaloMassAccretionHistoryHearin2021Stochastic(                           &amp;
     &amp;                                                 rateRollOver    =3.5d0   , &amp;
     &amp;                                                 fractionLateLow =0.0d0   , &amp;
     &amp;                                                 fractionLateHigh=0.0d0   , &amp;
     &amp;                                                 means           =means   , &amp;
     &amp;                                                 cholesky        =cholesky  &amp;
     &amp;                                                )
   </constructor>
  </referenceConstruct> 
  !!]
  ! Set node properties.
  call basicEarly%massSet(10.0d0**12.5d0)
  call basicLate %massSet(10.0d0**12.5d0)
  call basicEarly%timeSet(13.8d0        )
  call basicLate %timeSet(13.8d0        )
  ! Open the reference file produced by Andrew Hearin's "diffmah" code and read the parameter values used.
  countTimes=Count_Lines_In_File(char(inputPath(pathTypeExec))//'testSuite/data/hearin2021MAHMean.txt',comment_char='#')
  open(newUnit=referenceFile,file=char(inputPath(pathTypeExec))//'testSuite/data/hearin2021MAHMean.txt',status='old',form='formatted')
  do i=1,7
     read (referenceFile,'(a)') line
  end do
  allocate(massLateTarget    (countTimes))
  allocate(massEarlyTarget   (countTimes))
  allocate(massLateRecovered (countTimes))
  allocate(massEarlyRecovered(countTimes))
  do i=1,countTimes
     read (referenceFile,*) timeLogarithmic,massEarlyTarget(i),massLateTarget(i)
     time                 =10.0d0**timeLogarithmic
     massEarlyRecovered(i)=log10(darkMatterHaloMassAccretionHistoryEarlyForming_%mass(treeEarly%nodeBase,time))
     massLateRecovered (i)=log10(darkMatterHaloMassAccretionHistoryLateForming_ %mass(treeLate %nodeBase,time))
  end do
  call Assert('mass accretion history (early-forming)',massEarlyTarget,massEarlyRecovered,relTol=1.0d-6)
  call Assert('mass accretion history (late-forming)' ,massLateTarget ,massLateRecovered ,relTol=1.0d-6)
  ! End unit tests.
  close(referenceFile)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Hearin2021_Stochastic_MAH
