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
Contains a program which tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in a
dark energy Universe.
!!}

program Test_Zhao2009_Dark_Energy
  !!{
  Tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in a dark energy
  Universe. Comparisons are made to the \href{http://202.127.29.4/dhzhao/mandc_calculator.htm}{``{\normalfont \ttfamily mandc}''} Note that
  comparison tolerances are relatively large since we have not attempted to match details (such as critical density
  calculation) with ``{\normalfont \ttfamily mandc}''.
  !!}
  use :: Cosmological_Density_Field               , only : criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                      , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                     , only : cosmologyParametersSimple                                   , hubbleUnitsLittleH
  use :: Dark_Matter_Particles                    , only : darkMatterParticleCDM
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryZhao2009
  use :: Dark_Matter_Profiles_Concentration       , only : darkMatterProfileConcentrationZhao2009
  use :: Display                                  , only : displayVerbositySet                                         , verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: File_Utilities                           , only : Count_Lines_in_File
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize                                , nodeComponentBasic                     , treeNode
  use :: Input_Paths                              , only : inputPath                                                   , pathTypeExec
  use :: ISO_Varying_String                       , only : assignment(=)                                               , char                                   , operator(//)                       , varying_string              , &
       &                                                   var_str
  use :: Input_Parameters                         , only : inputParameters
  use :: Linear_Growth                            , only : linearGrowthCollisionlessMatter
  use :: Node_Components                          , only : Node_Components_Initialize                                  , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: String_Handling                          , only : operator(//)                                                , String_Superscript
  use :: Power_Spectra_Primordial                 , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred     , only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions          , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                       , only : transferFunctionBBKS
  use :: Unit_Tests                               , only : Assert                                                      , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                    ), pointer                 :: node
  class           (nodeComponentBasic                                          ), pointer                 :: basic
  integer                                                                       , dimension(2), parameter :: logarithmicHaloMasses              =[12,15]
  double precision                                                              , dimension(2), parameter :: concentrationDifferenceTolerance   =[3.1d-2,7.9d-4], timeDifferenceTolerance=[2.5d-2,2.0d-2]
  type            (cosmologyParametersSimple                                   )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )                          :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower                       )                          :: cosmologicalMassVariance_
  type            (linearGrowthCollisionlessMatter                             )                          :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat                           )                          :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             )                          :: powerSpectrumPrimordial_
  type            (transferFunctionBBKS                                        )                          :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    )                          :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                       )                          :: darkMatterParticle_
  type            (darkMatterProfileConcentrationZhao2009                      )                          :: darkMatterProfileConcentration_
  type            (darkMatterHaloMassAccretionHistoryZhao2009                  )                          :: darkMatterHaloMassAccretionHistory_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                          :: criticalOverdensity_
  type            (varying_string                                              )                          :: fileName                                           , message                         , &
       &                                                                                                     parameterFile
  integer                                                                                                 :: dataLinesInFile                                    , fUnit                           , &
       &                                                                                                     iLine                                              , iMass                           , &
       &                                                                                                     totalLinesInFile
  double precision                                                                                        :: concentrationDifferenceMaximum                     , haloMass                        , &
       &                                                                                                     ourConcentration                                   , ourTime                         , &
       &                                                                                                     redshift                                           , theirConcentration              , &
       &                                                                                                     theirTime                                          , timeDifferenceMaximum
  type            (inputParameters                                             )                          :: parameters
  character       (len=2                                                       )                          :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: dark energy cosmology")
  ! Test Zhao et al. 2009 algorithms in a dark energy universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Zhao2009Algorithms.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node  => treeNode      (                 )
  basic => node    %basic(autoCreate=.true.)
  ! Construct required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                         &amp;
     &amp;                                                        OmegaMatter                        =  0.25d0                           , &amp;
     &amp;                                                        OmegaBaryon                        =  0.00d0                           , &amp;
     &amp;                                                        OmegaDarkEnergy                    =  0.75d0                           , &amp;
     &amp;                                                        temperatureCMB                     =  2.70d0                           , &amp;
     &amp;                                                        HubbleConstant                     =100.00d0                             &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                               (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                              (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                              (                                                                         &amp;
     &amp;                                                        index_                             =+1.0d0                             , &amp;
     &amp;                                                        running                            =+0.0d0                             , &amp;
     &amp;                                                        runningRunning                     =+0.0d0                             , &amp;
     &amp;                                                        wavenumberReference                =+1.0d0                             , &amp;
     &amp;                                                        runningSmallScalesOnly             =.false.                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionBBKS                                        (                                                                          &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                     (                                                                         &amp;
     &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordial_           , &amp;
     &amp;                                                        transferFunction_                  =transferFunction_                  , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                            (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                        (                                                                         &amp;
     &amp;                                                        sigma8                             =0.8d0                              , &amp;
     &amp;                                                        tolerance                          =1.0d-4                             , &amp;
     &amp;                                                        toleranceTopHat                    =1.0d-4                             , &amp;
     &amp;                                                        nonMonotonicIsFatal                =.true.                             , &amp;
     &amp;                                                        monotonicInterpolation             =.false.                            , &amp;
     &amp;                                                        truncateAtParticleHorizon          =.false.                            , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_, &amp;
     &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                          &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        tableStore                         =.true.                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloMassAccretionHistory_">
   <constructor>
    darkMatterHaloMassAccretionHistoryZhao2009                        (                                                                    &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        criticalOverdensity_               =criticalOverdensity_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileConcentration_"    >
   <constructor>
    darkMatterProfileConcentrationZhao2009                       (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        darkMatterHaloMassAccretionHistory_=darkMatterHaloMassAccretionHistory_  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  ! Loop over halo masses to test.
  do iMass=1,size(logarithmicHaloMasses)
     ! Count lines in the "mandc" comparison file.
     fileName=inputPath(pathTypeExec)//'testSuite/data/zhao2009MassAccretionHistories/mandcoutputDarkEnergylgM'
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
     write (label,'(i2)') logarithmicHaloMasses(iMass)
     message=var_str("10")//String_Superscript(trim(adjustl(label)))//" M⊙ halo mass accretion history"
     call Assert(char(message),timeDifferenceMaximum         ,0.0d0,absTol=timeDifferenceTolerance         (iMass))
     message=var_str("10")//String_Superscript(trim(adjustl(label)))//" M⊙ halo concentration history"
     call Assert(char(message),concentrationDifferenceMaximum,0.0d0,absTol=concentrationDifferenceTolerance(iMass))
  end do
  ! End unit tests.
  call node%destroy()
  deallocate(node)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Zhao2009_Dark_Energy
