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

!% Contains a program which tests the \cite{diemer_universal_2014} halo concentration algorithm.

program Test_DiemerKravtsov2014_Concentration
  !% Tests the \cite{diemer_universal_2014} halo concentration algorithm. Values of concentration were taken from their website\footnote{File no longer available---was downloaded from {\normalfont \ttfamily http://www.benediktdiemer.com/wp-content/uploads/2014/07/Concentration\_WMAP7\_median.txt}}.
  use :: Cosmological_Density_Field          , only : cosmologicalMassVariance                        , cosmologicalMassVarianceClass    , criticalOverdensity                , criticalOverdensityClass
  use :: Cosmology_Functions                 , only : cosmologyFunctions                              , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParameters                             , cosmologyParametersClass         , hubbleUnitsLittleH
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationDiemerKravtsov2014
  use :: Display                             , only : displayVerbositySet                             , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : File_Exists
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Error                    , only : Galacticus_Error_Report
  use :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                    , nodeComponentBasic               , treeNode
  use :: Galacticus_Paths                    , only : galacticusPath                                  , pathTypeExec
  use :: ISO_Varying_String                  , only : assignment(=)                                   , char                             , operator(//)                       , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize                      , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Power_Spectra                       , only : powerSpectrum                                   , powerSpectrumClass
  use :: System_Command                      , only : System_Command_Do
  use :: Unit_Tests                          , only : Assert                                          , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                        ), pointer :: node
  class           (nodeComponentBasic                              ), pointer :: basic
  class           (cosmologyFunctionsClass                         ), pointer :: cosmologyFunctions_
  class           (cosmologyParametersClass                        ), pointer :: cosmologyParameters_
  class           (criticalOverdensityClass                        ), pointer :: criticalOverdensity_
  class           (cosmologicalMassVarianceClass                   ), pointer :: cosmologicalMassVariance_
  class           (powerSpectrumClass                              ), pointer :: powerSpectrum_
  type            (darkMatterProfileConcentrationDiemerKravtsov2014)          :: darkMatterProfileConcentration_
  type            (varying_string                                  )          :: parameterFile
  type            (inputParameters                                 )          :: parameters
  double precision                                                            :: ourConcentration               , differenceFractional, &
       &                                                                         concentration                  , mass                , &
       &                                                                         redshift                       , nu                  , &
       &                                                                         differenceFractionalMaximum
  integer                                                                     :: referenceUnit                  , ioStatus            , &
       &                                                                         i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("DiemerKravtsov2014 halo concentration algorithm")
  ! Test DiemerKravtsov2014 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/DiemerKravtsov2014HaloConcentration/testParameters.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)

  ! Get the data file if we don't have it.
  if (.not.File_Exists(galacticusPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt")) then
     call System_Command_Do(                                                                           &
          &                 "wget http://www.benediktdiemer.com/wp-content/uploads/cM_WMAP7.txt -O "// &
          &                 galacticusPath(pathTypeExec)                                            // &
          &                 "testSuite/data/diemerKravtsov2014Concentration.txt"                       &
          &                )
     if (.not.File_Exists(galacticusPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt")) &
          & call Galacticus_Error_Report('unable to retrieve reference dataset'//{introspection:location})
  end if
  ! Create a node.
  node                            => treeNode                                        (                 )
  ! Get the basic component.
  basic                           => node%basic                                      (autoCreate=.true.)
  ! Get required objects.
  cosmologyFunctions_             => cosmologyFunctions                              (                 )
  cosmologyParameters_            => cosmologyParameters                             (                 )
  criticalOverdensity_            => criticalOverdensity                             (                 )
  cosmologicalMassVariance_       => cosmologicalMassVariance                        (                 )
  powerSpectrum_                  => powerSpectrum                                   (                 )
  darkMatterProfileConcentration_ =  darkMatterProfileConcentrationDiemerKravtsov2014(                           &
       &                                                                              0.69d0                   , &
       &                                                                              6.58d0                   , &
       &                                                                              1.37d0                   , &
       &                                                                              6.82d0                   , &
       &                                                                              1.42d0                   , &
       &                                                                              1.12d0                   , &
       &                                                                              1.69d0                   , &
       &                                                                              0.00d0                   , &
       &                                                                              cosmologyFunctions_      , &
       &                                                                              cosmologyParameters_     , &
       &                                                                              criticalOverdensity_     , &
       &                                                                              cosmologicalMassVariance_, &
       &                                                                              powerSpectrum_             &
       &                                                                             )
  ! Read the reference file.
  differenceFractionalMaximum=0.0d0
  open(newUnit=referenceUnit,file=char(galacticusPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt"),status='old',form='formatted',iostat=ioStatus)
  do i=1,7
     read (referenceUnit,*,ioStat=ioStatus) ! Skip header.
  end do
  do while (ioStatus == 0)
     read (referenceUnit,*,ioStat=ioStatus) redshift,nu,mass,concentration
     if (ioStatus /= 0) exit
     ! Set the time for the node.
     call basic%timeSet            (cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     ! Set the mass of the original node (Diemer & Kravtsov masses are in units of M☉/h, so
     ! we convert from that system).
     call basic%massSet(mass/cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH))
     ! Compute and compare concentration at z=0.
     ourConcentration           =darkMatterProfileConcentration_%concentration(node)
     differenceFractional       =abs(ourConcentration-concentration)/concentration
     differenceFractionalMaximum=max(differenceFractionalMaximum,differenceFractional)
  end do
  ! Assert that the maximum fractional difference is not too large. The ~3% differences are
  ! presumably because we don't use precisely the same transfer function as do Diemer &
  ! Kravtsov.
  call Assert("Halo concentration in WMAP7 reference model",differenceFractionalMaximum,0.0d0,absTol=2.9d-2)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Galacticus_Function_Classes_Destroy()
end program Test_DiemerKravtsov2014_Concentration
