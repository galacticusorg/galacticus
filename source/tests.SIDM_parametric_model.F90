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

!+    Contributions to this file made by: Niusha Ahvazi

!!{RST
Tests of the inputs to the SIDM parametric model of :cite:t:`yang_parametric_2024` (implemented by the :galacticus-class:`nodeOperatorSIDMParametric` class). A cold dark matter (CDM) halo mass accretion history, extracted from an N-body simulation, is read from file; for each step the maximum circular velocity and virial radius of the corresponding NFW halo are computed and compared against the values tabulated in the simulation. This validates the CDM-side quantities from which the parametric SIDM solution is built.

The reference data files ``testSuite/data/SIDM/data_799_cdm_NFW.txt`` and ``testSuite/data/SIDM/data_799_cdm.txt`` are tabulated from an N-body simulation and committed alongside this test.
!!}

program Tests_SIDM_Parametric_Model
  use :: Cosmological_Density_Field                , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                      , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Mass_Accretion_Histories , only : darkMatterHaloMassAccretionHistoryCorrea2015
  use :: Dark_Matter_Halo_Scales                   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Particles                     , only : darkMatterParticleCDM                             , darkMatterParticleSIDMVelocityDependent
  use :: Dark_Matter_Profiles_Concentration        , only : darkMatterProfileConcentrationFixed
  use :: Dark_Matter_Profiles_DMO                  , only : darkMatterProfileDMONFW
  use :: Display                                   , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Error                                     , only : Error_Report
  use :: Events_Hooks                              , only : eventsHooksInitialize
  use :: File_Utilities                            , only : Count_Lines_In_File
  use :: Functions_Global_Utilities                , only : Functions_Global_Set
  use :: Galacticus_Nodes                          , only : mergerTree                                        , nodeClassHierarchyInitialize           , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                                    treeNode                                          , treeNodeList
  use :: Input_Parameters                          , only : inputParameters
  use :: Linear_Growth                             , only : linearGrowthCollisionlessMatter
  use :: Mass_Distributions                        , only : massDistributionClass
  use :: Node_Components                           , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Nodes_Operators                           , only : nodeOperatorSIDMParametric
  use :: Numerical_Constants_Prefixes              , only : kilo
  use :: Power_Spectra_Primordial                  , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred      , only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions           , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                        , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                                , only : Assert                                            , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Virial_Density_Contrast                   , only : virialDensityContrastBryanNorman1998
  implicit none
  ! Hubble parameter (h=H₀/100 km/s/Mpc) of the simulation from which the reference data were extracted.
  double precision                                                    , parameter                 :: hubbleParameter      =0.7d0
  ! Reference data files (CDM halo mass accretion history from an N-body simulation).
  character       (len=*                                             ), parameter                 :: fileNameNFW          ='testSuite/data/SIDM/data_799_cdm_NFW.txt', &
       &                                                                                             fileNameCDM          ='testSuite/data/SIDM/data_799_cdm.txt'
  type            (inputParameters                                   )                            :: parameters
  ! Cosmology, power spectrum, and structure-formation objects.
  type            (cosmologyParametersSimple                         ), pointer                   :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                    ), pointer                   :: cosmologyFunctions_
  type            (virialDensityContrastBryanNorman1998              ), pointer                   :: virialDensityContrast_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer                   :: darkMatterHaloScale_
  type            (darkMatterProfileDMONFW                           ), pointer                   :: darkMatterProfileDMONFW_
  type            (nodeOperatorSIDMParametric                        ), pointer                   :: nodeOperator_
  type            (darkMatterParticleCDM                             )                            :: darkMatterParticleCDM_
  type            (darkMatterParticleSIDMVelocityDependent           )                            :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter                   )                            :: linearGrowth_
  type            (powerSpectrumPrimordialPowerLaw                   )                            :: powerSpectrumPrimordial_
  type            (powerSpectrumPrimordialTransferredSimple          )                            :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumWindowFunctionTopHat                 )                            :: powerSpectrumWindowFunction_
  type            (transferFunctionEisensteinHu1999                  )                            :: transferFunction_
  type            (cosmologicalMassVarianceFilteredPower             )                            :: cosmologicalMassVariance_
  type            (darkMatterHaloMassAccretionHistoryCorrea2015      )                            :: darkMatterHaloMassAccretionHistory_
  type            (darkMatterProfileConcentrationFixed               )                            :: darkMatterProfileConcentration_
  ! Tree, nodes, and components.
  type            (mergerTree                                        )                            :: tree
  type            (treeNodeList                                      ), allocatable, dimension(:) :: nodes
  class           (nodeComponentBasic                                ), pointer                   :: basic
  class           (nodeComponentDarkMatterProfile                    ), pointer                   :: darkMatterProfile
  class           (massDistributionClass                             ), pointer                   :: massDistribution_
  double precision                                                    , parameter                 :: radiusMaximumFactorNFW=2.1626d0 ! The ratio of R_max to R_s in an NFW profile.
  ! Reference data and computed quantities.
  double precision                                                    , allocatable, dimension(:) :: expansionFactor          , massVirialData    , radiusVirialData             , &
       &                                                                                             radiusVelocityMaximumData, timeData          , radiusScale                  , &
       &                                                                                             velocityMaximumData      , velocityMaximum   , radiusVirial
  integer                                                                                         :: countData                , unitFile          , ioStatus            , i
  double precision                                                                                :: column1                  , column2           , column3             , column4, &
       &                                                                                             column5                  , column6           , column7             , column8

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM parametric model: CDM accretion-history inputs")
  ! Read the parameter file and initialize the node-component hierarchy.
  parameters=inputParameters('testSuite/parameters/nodes/nodes_SIDM_parametric.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the cosmology, power spectrum, and structure-formation object stack. These are the dependencies required to construct
  ! the SIDM parametric node operator, and to compute NFW halo properties for comparison with the simulation data.
  allocate(cosmologyParameters_    )
  allocate(cosmologyFunctions_     )
  allocate(virialDensityContrast_  )
  allocate(darkMatterHaloScale_    )
  allocate(darkMatterProfileDMONFW_)
  allocate(nodeOperator_           )
  cosmologyParameters_                =cosmologyParametersSimple                         (                                                                              &
       &                                                                                  OmegaMatter                            =0.2815d0                            , &
       &                                                                                  OmegaBaryon                            =0.0465d0                            , &
       &                                                                                  OmegaDarkEnergy                        =0.7185d0                            , &
       &                                                                                  temperatureCMB                         =2.78d0                              , &
       &                                                                                  HubbleConstant                         =100.0d0*hubbleParameter               &
       &                                                                                 )
  cosmologyFunctions_                 =cosmologyFunctionsMatterLambda                    (                                                                              &
       &                                                                                  cosmologyParameters_                   =cosmologyParameters_                  &
       &                                                                                 )
  virialDensityContrast_              =virialDensityContrastBryanNorman1998              (                                                                              &
       &                                                                                  allowUnsupportedCosmology              =.false.                             , &
       &                                                                                  cosmologyParameters_                   =cosmologyParameters_                , &
       &                                                                                  cosmologyFunctions_                    =cosmologyFunctions_                   &
       &                                                                                 )
  darkMatterHaloScale_                =darkMatterHaloScaleVirialDensityContrastDefinition(                                                                              &
       &                                                                                  cosmologyParameters_                   =cosmologyParameters_                , &
       &                                                                                  cosmologyFunctions_                    =cosmologyFunctions_                 , &
       &                                                                                  virialDensityContrast_                 =virialDensityContrast_                &
       &                                                                                 )
  darkMatterProfileDMONFW_            =darkMatterProfileDMONFW                           (                                                                              &
       &                                                                                  velocityDispersionUseSeriesExpansion   =.false.                             , &
       &                                                                                  darkMatterHaloScale_                   =darkMatterHaloScale_                  &
       &                                                                                 )
  darkMatterParticleCDM_              =darkMatterParticleCDM                             (                                                                              &
       &                                                                                 )
  darkMatterParticle_                 =darkMatterParticleSIDMVelocityDependent           (                                                                              &
       &                                                                                  velocityCharacteristic                 = 24.3289794155754d0                 , &
       &                                                                                  sigma0                                 =147.1008800000000d0                 , &
       &                                                                                  darkMatterParticle_                    =darkMatterParticleCDM_                &
       &                                                                                 )
  linearGrowth_                       =linearGrowthCollisionlessMatter                   (                                                                              &
       &                                                                                  cosmologyParameters_                    =cosmologyParameters_               , &
       &                                                                                  cosmologyFunctions_                     =cosmologyFunctions_                  &
       &                                                                                 )
  powerSpectrumWindowFunction_        =powerSpectrumWindowFunctionTopHat                 (                                                                              &
       &                                                                                  cosmologyParameters_                    =cosmologyParameters_                 &
       &                                                                                 )
  powerSpectrumPrimordial_            =powerSpectrumPrimordialPowerLaw                   (                                                                              &
       &                                                                                  index_                                  =0.971d0                            , &
       &                                                                                  running                                 =0.000d0                            , &
       &                                                                                  runningRunning                          =0.000d0                            , &
       &                                                                                  wavenumberReference                     =1.000d0                            , &
       &                                                                                  runningSmallScalesOnly                  =.false.                              &
       &                                                                                 )
  transferFunction_                   =transferFunctionEisensteinHu1999                  (                                                                              &
       &                                                                                  neutrinoNumberEffective                 =3.046d0                            , &
       &                                                                                  neutrinoMassSummed                      =0.0d0                              , &
       &                                                                                  darkMatterParticle_                     =darkMatterParticleCDM_             , &
       &                                                                                  cosmologyParameters_                    =cosmologyParameters_               , &
       &                                                                                  cosmologyFunctions_                     =cosmologyFunctions_                  &
       &                                                                                 )
  powerSpectrumPrimordialTransferred_ =powerSpectrumPrimordialTransferredSimple          (                                                                              &
       &                                                                                  powerSpectrumPrimordial_                =powerSpectrumPrimordial_           , &
       &                                                                                  transferFunction_                       =transferFunction_                  , &
       &                                                                                  linearGrowth_                           =linearGrowth_                        &
       &                                                                                 )
  cosmologicalMassVariance_           =cosmologicalMassVarianceFilteredPower             (                                                                              &
       &                                                                                  sigma8                                  =0.82d+00                           , &
       &                                                                                  tolerance                               =4.00d-06                           , &
       &                                                                                  toleranceTopHat                         =1.00d-06                           , &
       &                                                                                  rootVarianceLogarithmicGradientTolerance=1.00d-12                           , &
       &                                                                                  nonMonotonicIsFatal                     =.true.                             , &
       &                                                                                  monotonicInterpolation                  =.false.                            , &
       &                                                                                  truncateAtParticleHorizon               =.false.                            , &
       &                                                                                  storeTabulations                        =.true.                             , &
       &                                                                                  integrationFailureIsFatal               =.true.                             , &
       &                                                                                  cosmologyParameters_                    =cosmologyParameters_               , &
       &                                                                                  cosmologyFunctions_                     =cosmologyFunctions_                , &
       &                                                                                  linearGrowth_                           =linearGrowth_                      , &
       &                                                                                  powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_, &
       &                                                                                  powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_         &
       &                                                                                 )
  darkMatterHaloMassAccretionHistory_ =darkMatterHaloMassAccretionHistoryCorrea2015      (                                                                              &
       &                                                                                  cosmologyFunctions_                     =cosmologyFunctions_                , &
       &                                                                                  linearGrowth_                           =linearGrowth_                      , &
       &                                                                                  cosmologicalMassVariance_               =cosmologicalMassVariance_            &
       &                                                                                 )
  darkMatterProfileConcentration_     =darkMatterProfileConcentrationFixed               (                                                                              &
       &                                                                                  concentration_                          =5.0d0                              , &
       &                                                                                  virialDensityContrast_                  =virialDensityContrast_             , &
       &                                                                                  darkMatterProfileDMO_                   =darkMatterProfileDMONFW_             &
       &                                                                                 )
  nodeOperator_                       =nodeOperatorSIDMParametric                        (                                                                              &
       &                                                                                  alpha                                   =2.0d0                              , &
       &                                                                                  C                                       =0.75d0                             , &
       &                                                                                  darkMatterParticle_                     =darkMatterParticle_                , &
       &                                                                                  darkMatterHaloMassAccretionHistory_     =darkMatterHaloMassAccretionHistory_, &
       &                                                                                  darkMatterProfileDMO_                   =darkMatterProfileDMONFW_           , &
       &                                                                                  cosmologyFunctions_                     =cosmologyFunctions_                , &
       &                                                                                  cosmologyParameters_                    =cosmologyParameters_               , &
       &                                                                                  darkMatterHaloScale_                    =darkMatterHaloScale_               , &
       &                                                                                  virialDensityContrast_                  =virialDensityContrast_             , &
       &                                                                                  darkMatterProfileConcentration_         =darkMatterProfileConcentration_      &
       &                                                                                 )
  ! Read the reference CDM data. The number of data lines is determined from the file (excluding the single header line). Both
  ! files are tabulated from earliest to latest time; we fill our arrays in reverse so that index 1 corresponds to the final
  ! (base) halo and increasing index corresponds to its progenitors.
  countData=Count_Lines_In_File(fileNameNFW)-1
  allocate(expansionFactor          (countData))
  allocate(massVirialData           (countData))
  allocate(radiusVirialData         (countData))
  allocate(radiusVelocityMaximumData(countData))
  allocate(timeData                 (countData))
  allocate(radiusScale              (countData))
  allocate(velocityMaximumData      (countData))
  allocate(velocityMaximum          (countData))
  allocate(radiusVirial             (countData))
  allocate(nodes                    (countData))
  ! Read the NFW halo properties: expansion factor, virial mass, virial radius, and radius of maximum circular velocity.
  open(newUnit=unitFile,file=fileNameNFW,status='old',action='read',iostat=ioStatus)
  if (ioStatus /= 0) call Error_Report('unable to open reference data file "'//fileNameNFW//'"'//{introspection:location})
  read(unitFile,*) ! Skip the header line.
  do i=countData,1,-1
     read(unitFile,*) column1,column2,column3,column4
     expansionFactor          (i)=column1
     massVirialData           (i)=column2                   /hubbleParameter          ! [M☉] (remove factor of h).
     radiusVirialData         (i)=column3*expansionFactor(i)/hubbleParameter          ! [kpc] (comoving, /h -> physical).
     radiusVelocityMaximumData(i)=column4*expansionFactor(i)/hubbleParameter          ! [kpc]
     timeData                 (i)=cosmologyFunctions_%cosmicTime(expansionFactor(i))
     radiusScale              (i)=radiusVelocityMaximumData(i)/radiusMaximumFactorNFW ! NFW scale radius [kpc].
  end do
  close(unitFile)
  ! Read the maximum circular velocity (column 7) from the CDM data file.
  open(newUnit=unitFile,file=fileNameCDM,status='old',action='read',iostat=ioStatus)
  if (ioStatus /= 0) call Error_Report('unable to open reference data file "'//fileNameCDM//'"'//{introspection:location})
  read(unitFile,*) ! Skip the header line.
  do i=countData,1,-1
     read(unitFile,*) column1,column2,column3,column4,column5,column6,column7,column8
     velocityMaximumData(i)=column7 ! [km/s].
  end do
  close(unitFile)
  ! Create the nodes and link them into a single (smooth-accretion) branch.
  do i=1,countData
     nodes(i)%node => treeNode()
  end do
  do i=1,countData-1
     nodes(i  )%node%firstChild => nodes(i+1)%node
     nodes(i+1)%node%parent     => nodes(i  )%node
  end do
  tree%nodeBase => nodes(1)%node
  ! For each node, set the basic and dark matter profile properties from the data and compute the NFW maximum circular velocity
  ! and virial radius.
  do i=1,countData
     basic             => nodes(i)%node%basic            (autoCreate=.true.)
     darkMatterProfile => nodes(i)%node%darkMatterProfile(autoCreate=.true.)
     call basic            %massSet(massVirialData(i))
     call basic            %timeSet(timeData      (i))
     call darkMatterProfile%scaleSet(radiusScale  (i)/kilo)                                       ! kpc -> Mpc.
     massDistribution_    => darkMatterProfileDMONFW_%get                         (nodes(i)%node)
     velocityMaximum  (i) =  massDistribution_       %velocityRotationCurveMaximum(             ) ! [km/s].
     radiusVirial     (i) =  darkMatterHaloScale_    %radiusVirial                (nodes(i)%node) ! [Mpc].
  end do
  ! Compare the computed NFW properties against the simulation reference values, requiring agreement to better than 1%.
  call Assert("NFW maximum circular velocity vs simulation",velocityMaximum     ,velocityMaximumData,relTol=1.0d-2)
  call Assert("NFW virial radius vs simulation"            ,radiusVirial   *kilo,radiusVirialData   ,relTol=1.0d-2)
  ! Destroy the nodes.
  do i=1,countData
     call nodes(i)%node%destroy()
  end do
  ! Clean up.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_Parametric_Model
