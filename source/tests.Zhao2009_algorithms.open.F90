!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a program which tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in an
!% open Universe.

program Test_Zhao2009_Open
  !% Tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in an open
  !% Universe. Comparisons are made to the \href{http://202.127.29.4/dhzhao/mandc_calculator.htm}{``{\normalfont \ttfamily mandc}''} Note that
  !% comparison tolerances are relatively large since we have not attempted to match details (such as critical density
  !% calculation) with ``{\normalfont \ttfamily mandc}''.
  use :: Cosmology_Functions                      , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistory , darkMatterHaloMassAccretionHistoryClass
  use :: Dark_Matter_Profiles_Concentration       , only : darkMatterProfileConcentration     , darkMatterProfileConcentrationClass
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: File_Utilities                           , only : Count_Lines_in_File
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Display                       , only : Galacticus_Verbosity_Level_Set     , verbosityStandard
  use :: Galacticus_Function_Classes_Destroys     , only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize       , nodeComponentBasic                     , treeNode
  use :: Galacticus_Paths                         , only : galacticusPath                     , pathTypeExec
  use :: ISO_Varying_String                       , only : varying_string                     , assignment(=)                          , operator(//)                       , char
  use :: Input_Parameters                         , only : inputParameters
  use :: Node_Components                          , only : Node_Components_Initialize         , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: String_Handling                          , only : operator(//)
  use :: Unit_Tests                               , only : Assert                             , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                               )                         , pointer :: node
  class           (nodeComponentBasic                     )                         , pointer :: basic
  integer                                                  , dimension(1), parameter          :: logarithmicHaloMasses           =[12]
  double precision                                         , dimension(1), parameter          :: concentrationDifferenceTolerance=[3.6d-2], timeDifferenceTolerance=[2.3d-2]
  class           (cosmologyFunctionsClass                )                         , pointer :: cosmologyFunctions_
  class           (darkMatterProfileConcentrationClass    )                         , pointer :: darkMatterProfileConcentration_
  class           (darkMatterHaloMassAccretionHistoryClass)                         , pointer :: darkMatterHaloMassAccretionHistory_
  type            (varying_string                         )                                   :: fileName                                 , message                         , &
       &                                                                                         parameterFile
  integer                                                                                     :: dataLinesInFile                          , fUnit                           , &
       &                                                                                         iLine                                    , iMass                           , &
       &                                                                                         totalLinesInFile
  double precision                                                                            :: concentrationDifferenceMaximum           , haloMass                        , &
       &                                                                                         ourConcentration                         , ourTime                         , &
       &                                                                                         redshift                                 , theirConcentration              , &
       &                                                                                         theirTime                                , timeDifferenceMaximum
  type            (inputParameters                        )                                   :: parameters

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: open cosmology")

  ! Test Zhao et al. 2009 algorithms in an open universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Zhao2009Algorithms/open.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)

  ! Create a node.
  node  => treeNode      (                 )

  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Get the default cosmology functions object.
  cosmologyFunctions_                 => cosmologyFunctions                ()
  ! Get the default concentrations object.
  darkMatterProfileConcentration_     => darkMatterProfileConcentration    ()
  ! Get the default mass accretion history object.
  darkMatterHaloMassAccretionHistory_ => darkMatterHaloMassAccretionHistory()

  ! Loop over halo masses to test.
  do iMass=1,size(logarithmicHaloMasses)

     ! Count lines in the "mandc" comparison file.
     fileName=char(galacticusPath(pathTypeExec))//'testSuite/data/zhao2009MassAccretionHistories/mandcoutputOpenlgM'
     fileName=fileName//logarithmicHaloMasses(iMass)//'.data'
     totalLinesInFile=Count_Lines_in_File(fileName    )
     dataLinesInFile =Count_Lines_in_File(fileName,'#')-1

     ! Discard file header.
     open(newunit=fUnit,file=char(fileName),status='old',form='formatted')
     do iLine=1,(totalLinesInFile-dataLinesInFile)
        read (fUnit,*)
     end do

     ! Initialize maximum differences to zero.
     timeDifferenceMaximum         =0.0d0
     concentrationDifferenceMaximum=0.0d0

     ! Read all data lines from the comparison file.
     do iLine=1,dataLinesInFile
        read (fUnit,*) redshift,haloMass,theirConcentration

        ! Compute the corresponding cosmological time.
        theirTime=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))

        ! Set the mass and time of the original node.
        call basic%massSet(10.0d0**logarithmicHaloMasses(iMass))
        call basic%timeSet(cosmologyFunctions_%cosmicTime(1.0d0))

        ! Get the time corresponding to the current halo mass.
        ourTime=darkMatterHaloMassAccretionHistory_%time(node,haloMass)

        ! Set the node mass and time to the current values.
        call basic%massSet(haloMass )
        call basic%timeSet(theirTime)

        ! Get the corresponding halo concentration.
        ourConcentration=darkMatterProfileConcentration_%concentration(node)

        ! Compute the difference between our values and the comparison values.
        timeDifferenceMaximum         =max(                                                             &
             &                              timeDifferenceMaximum                                       &
             &                             ,abs(ourTime         -theirTime         )/theirTime          &
             &                             )
        concentrationDifferenceMaximum=max(                                                             &
             &                              concentrationDifferenceMaximum                              &
             &                             ,abs(ourConcentration-theirConcentration)/theirConcentration &
             &                            )

     end do
     close(fUnit)

     ! Perform the tests.
     message='10^'
     message=message//logarithmicHaloMasses(iMass)//' M⊙ halo mass accretion history'
     call Assert(char(message),timeDifferenceMaximum         ,0.0d0,absTol=timeDifferenceTolerance         (iMass))
     message='10^'
     message=message//logarithmicHaloMasses(iMass)//' M⊙ halo concentration history'
     call Assert(char(message),concentrationDifferenceMaximum,0.0d0,absTol=concentrationDifferenceTolerance(iMass))

  end do

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Galacticus_Function_Classes_Destroy()

end program Test_Zhao2009_Open
