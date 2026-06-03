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

!!{
Contains a program which tests the Hearin (2021) halo mass formation history.
!!}

program Test_Hearin2021_MAH
  !!{
  Tests the Hearin (2021) halo mass formation history algorithm.
  !!}
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryHearin2021
  use :: Display                                  , only : displayVerbositySet                         , verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: File_Utilities                           , only : Count_Lines_in_File
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Input_Paths                              , only : inputPath                                   , pathTypeExec
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize                , nodeComponentBasic               , treeNode
  use :: Input_Parameters                         , only : inputParameters
  use :: ISO_Varying_String                       , only : char                                        , varying_string                   , operator(//)
  use :: Node_Components                          , only : Node_Components_Initialize                  , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                               , only : Assert                                      , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none  
  type            (treeNode                                    ), pointer                   :: node
  class           (nodeComponentBasic                          ), pointer                   :: basic
  type            (darkMatterHaloMassAccretionHistoryHearin2021)                            :: darkMatterHaloMassAccretionHistory_
  double precision                                              , dimension(:), allocatable :: massTarget                         , massRecovered
  double precision                                                                          :: time                               , powerLawIndexEarly, &
       &                                                                                       powerLawIndexLate                  , rateRollover      , &
       &                                                                                       logMassMaximum                     , logTimeMaximum    , &
       &                                                                                       massMaximum                        , timeMaximum
  integer                                                                                   :: i                                  , referenceFile     , &
       &                                                                                       countTimes
  type            (inputParameters                             )                            :: parameters
  character       (len=1024                                    )                            :: line
  type            (varying_string                              )                            :: referenceFileName
  
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
  ! Create a node.
  node  => treeNode      (                 )
  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Open the reference file produced by Andrew Hearin's "diffmah" code and read the parameter values used.
  referenceFilename=inputPath(pathTypeExec)//'testSuite/data/diffmah_massAccretionHistory.txt'
  countTimes=Count_Lines_In_File(referenceFilename,comment_char='#')
  open(newUnit=referenceFile,file=char(referenceFilename),status='old',form='formatted')
  read (referenceFile,'(a)') line
  read (line(index(line,"=")+1:len(line)),*) powerLawIndexEarly
  read (referenceFile,'(a)') line
  read (line(index(line,"=")+1:len(line)),*) powerLawIndexLate
  read (referenceFile,'(a)') line
  read (line(index(line,"=")+1:len(line)),*) rateRollover
  read (referenceFile,'(a)') line
  read (line(index(line,"=")+1:len(line)),*) logMassMaximum
  read (referenceFile,'(a)') line
  read (line(index(line,"=")+1:len(line)),*) logTimeMaximum
  read (referenceFile,'(a)') line
  read (referenceFile,'(a)') line
  timeMaximum=10.0d0**logTimeMaximum
  massMaximum=10.0d0**logMassMaximum
  ! Get required objects.
  !![
  <referenceConstruct object="darkMatterHaloMassAccretionHistory_">
   <constructor>
    darkMatterHaloMassAccretionHistoryHearin2021(                                        &amp;
     &amp;                                        powerLawIndexEarly=powerLawIndexEarly, &amp;
     &amp;                                        powerLawIndexLate=powerLawIndexLate  , &amp;
     &amp;                                        rateRollOver     =rateRollover       , &amp;
     &amp;                                        timeMaximum      =timeMaximum          &amp;
     &amp;                                       )
   </constructor>
  </referenceConstruct> 
  !!]
  ! Set node properties.
  call basic%massSet(massMaximum)
  call basic%timeSet(timeMaximum)
  allocate(massTarget   (countTimes))
  allocate(massRecovered(countTimes))
  do i=1,countTimes
     read (referenceFile,*) time,massTarget(i)
     massRecovered(i)=darkMatterHaloMassAccretionHistory_%mass(node,time)
  end do
  call Assert('mass accretion history',massTarget,massRecovered,relTol=1.0d-3)
  ! End unit tests.
  close(referenceFile)
  call node%destroy()
  deallocate(node,massTarget,massRecovered)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Hearin2021_MAH
