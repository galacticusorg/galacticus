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
Contains a program which tests concentration models.
!!}

program Test_Concentrations
  !!{
  Tests concentration models.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower             , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationClass               , darkMatterProfileConcentrationDiemerKravtsov2014            , darkMatterProfileConcentrationDuttonMaccio2014, darkMatterProfileConcentrationLudlow2016Fit  , &
       &                                              darkMatterProfileConcentrationPrada2011           , duttonMaccio2014FitTypeNFWVirial                            , duttonMaccio2014FitTypeNFWCritical200         , darkMatterProfileConcentrationDiemerJoyce2019
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profile_Scales          , only : darkMatterProfileScaleRadiusConcentration         , darkMatterProfileScaleRadiusLudlow2016Analytic
  use :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMONFW
  use :: Display                             , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : Count_Lines_in_File
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Calculations_Resets                 , only : Calculations_Reset
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                      , nodeComponentBasic                                          , treeNode                                      , nodeClassHierarchyFinalize                   , &
       &                                              nodeComponentDarkMatterProfile
  use :: Input_Paths                         , only : inputPath                                         , pathTypeExec
  use :: ISO_Varying_String                  , only : assignment(=)                                     , char                                                        , operator(//)                                  , operator(==)                                 , &
          &                                           varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize                           , Node_Components_Thread_Uninitialize           , Node_Components_Uninitialize
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: String_Handling                     , only : String_Split_Words
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1998
  use :: Virial_Density_Contrast             , only : virialDensityContrastClass                        , virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastFixed, fixedDensityTypeCritical
  use :: Unit_Tests                          , only : Assert                                            , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group                          , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                    ), pointer                             :: node
  class           (nodeComponentBasic                                          ), pointer                             :: basic
  class           (nodeComponentDarkMatterProfile                              ), pointer                             :: darkMatterProfile
  class           (darkMatterProfileConcentrationClass                         ), pointer                             :: darkMatterProfileConcentration_
  class           (virialDensityContrastClass                                  ), pointer                             :: virialDensityContrast_
  type            (cosmologyParametersSimple                                   ), pointer                             :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              ), pointer                             :: cosmologyFunctionsMatterLambda_
  type            (linearGrowthCollisionlessMatter                             ), pointer                             :: linearGrowthCollisionlessMatter_
  type            (cosmologicalMassVarianceFilteredPower                       ), pointer                             :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionTopHat                           ), pointer                             :: powerSpectrumWindowFunctionTopHat_
  type            (powerSpectrumPrimordialPowerLaw                             ), pointer                             :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionEisensteinHu1998                            ), pointer                             :: transferFunctionEisensteinHu1998_
  type            (powerSpectrumPrimordialTransferredSimple                    ), pointer                             :: powerSpectrumPrimordialTransferredSimple_
  type            (powerSpectrumStandard                                       ), pointer                             :: powerSpectrumStandard_
  type            (darkMatterParticleCDM                                       ), pointer                             :: darkMatterParticleCDM_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer                             :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (darkMatterProfileScaleRadiusConcentration                   ), pointer                             :: darkMatterProfileScaleRadiusConcentration_
  type            (darkMatterProfileScaleRadiusLudlow2016Analytic              ), pointer                             :: darkMatterProfileScaleRadiusLudlow2016Analytic_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition          ), pointer                             :: darkMatterHaloScaleVirialDensityContrastDefinition_
  type            (darkMatterProfileDMONFW                                     ), pointer                             :: darkMatterProfileDMONFW_
  type            (varying_string                                              )                                      :: parameterFile
  integer                                                                       , parameter                           :: countModels                                                   =7
  type            (varying_string                                              ), dimension(countModels)              :: modelName                                                       , modelLabel         , &
       &                                                                                                                 modelDensityContrast
  double precision                                                              , dimension(countModels)              :: modelTolerance
  double precision                                                              , dimension(2          )              :: redshifts
  type            (inputParameters                                             ), pointer                             :: parameters
  character       (len= 128                                                    ), dimension(4          )              :: columns
  character       (len=   3                                                    )                                      :: label
  double precision                                                              , dimension(:          ), allocatable :: mass                                                            , concentrationTarget, &
       &                                                                                                                 concentration
  character       (len=1024                                                    )                                      :: line
  double precision                                                                                                    :: OmegaMatter                                                     , OmegaBaryon        , &
       &                                                                                                                 OmegaDarkEnergy                                                 , temperatureCMB     , &
       &                                                                                                                 powerSpectrumIndex                                              , HubbleConstant     , &
       &                                                                                                                 neutrinoNumberEffective                                         , sigma8             , &
       &                                                                                                                 redshift
  integer                                                                                                             :: colossusFile                                                    , status             , &
       &                                                                                                                 countMasses                                                     , i                  , &
       &                                                                                                                 iModel                                                          , iRedshift

  ! Initialize.
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Specify all models to run.
  modelName           (1)='Diemer & Kravtsov (2015)'
  modelLabel          (1)='diemer15_orig_200c'
  modelDensityContrast(1)='200c'
  modelTolerance      (1)=1.2d-1 ! Use a relatively large tolerance here as this model depends on the slope of the power spectrum,
                                 ! which becomes oscillatory due to BAOs, leading to small but significant differences relative to
                                 ! Colossus.
  modelName           (2)='Dutton et al. (2014; 200c)'
  modelLabel          (2)='dutton14_200c'
  modelDensityContrast(2)='200c'
  modelTolerance      (2)=1.0d-3
  modelName           (3)='Dutton et al. (2014; vir)'
  modelLabel          (3)='dutton14_vir'
  modelDensityContrast(3)='vir'
  modelTolerance      (3)=1.0d-3
  modelName           (4)='Ludlow et al. (2016)'
  modelLabel          (4)='ludlow16_200c'
  modelDensityContrast(4)='200c'
  modelTolerance      (4)=3.5d-2
  modelName           (5)='Prada et al. (2012)'
  modelLabel          (5)='prada12_200c'
  modelDensityContrast(5)='200c'
  modelTolerance      (5)=1.0d-2
  modelName           (6)='Diemer & Joyce (2019)'
  modelLabel          (6)='diemer19_200c'
  modelDensityContrast(6)='200c'
  modelTolerance      (6)=1.0d-2
  modelName           (7)='Diemer & Joyce (2019; vir [converted])'
  modelLabel          (7)='diemer19_vir'
  modelDensityContrast(7)='vir'
  modelTolerance      (7)=1.0d-2
  ! Set redshifts.
  redshifts=[0.0d0,1.0d0]
  ! Iterate over models.
  call Unit_Tests_Begin_Group("Concentration algorithms")
  do iRedshift=1,2
     redshift=redshifts(iRedshift)
     write (label,'(f3.1)') redshift
        call Unit_Tests_Begin_Group("Redshift: z = "//trim(label))
        do iModel=1,countModels
           allocate(parameters)
           parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/concentrations_'//modelDensityContrast(iModel)//'.xml'
           parameters   =inputParameters(parameterFile)
           call nodeClassHierarchyInitialize     (parameters)
           call Node_Components_Initialize       (parameters)
           call Node_Components_Thread_Initialize(parameters)
           ! Read the Colossus target data (and parameters used) from file.
           countMasses=Count_Lines_in_File(inputPath(pathTypeExec)//'testSuite/data/concentrationsColossus/'//char(modelLabel(iModel))//'_z'//trim(label)//'.txt','#')
           i=0
           allocate(mass(               countMasses))
           allocate(concentration      (countMasses))
           allocate(concentrationTarget(countMasses))
           open(newUnit=colossusFile,file=char(inputPath(pathTypeExec))//'testSuite/data/concentrationsColossus/'//char(modelLabel(iModel))//'_z'//trim(label)//'.txt',status='old',form='formatted')
           do while (.true.)
              read (colossusFile,'(a)',ioStat=status) line
              if (status /= 0) exit
              if (line(1:1) == "#") then
                 call String_Split_Words(columns,line," ")
                 if (columns(3) == "=") then
                    select case (trim(columns(2)))
                    case ('OmegaMatter'             )
                       read (columns(4),*) OmegaMatter
                    case ('OmegaBaryon'             )
                       read (columns(4),*) OmegaBaryon
                    case ('OmegaDarkEnergy'         )
                       read (columns(4),*) OmegaDarkEnergy
                    case ('temperatureCMB'          )
                       read (columns(4),*) temperatureCMB
                    case ('index'                   )
                       read (columns(4),*) powerSpectrumIndex
                    case ('HubbleConstant'          )
                       read (columns(4),*) HubbleConstant
                    case ('effectiveNumberNeutrinos')
                       read (columns(4),*) neutrinoNumberEffective
                    case ('sigma_8'                 )
                       read (columns(4),*) sigma8
                    end select
                 end if
              else
                 i=i+1
                 read (line,*) mass(i),concentrationTarget(i)
              end if
           end do
           close(colossusFile)
           ! Convert masses from h⁻¹M☉ to M☉.
           mass=mass/(HubbleConstant/100.0d0)
           ! Construct all required objects.
           allocate(cosmologyParametersSimple_                                   )
           allocate(cosmologyFunctionsMatterLambda_                              )
           allocate(linearGrowthCollisionlessMatter_                             )
           allocate(cosmologicalMassVarianceFilteredPower_                       )
           allocate(powerSpectrumWindowFunctionTopHat_                           )
           allocate(powerSpectrumPrimordialPowerLaw_                             )
           allocate(transferFunctionEisensteinHu1998_                            )
           allocate(powerSpectrumPrimordialTransferredSimple_                    )
           allocate(powerSpectrumStandard_                                       )
           allocate(darkMatterParticleCDM_                                       )
           allocate(criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
           allocate(darkMatterHaloScaleVirialDensityContrastDefinition_          )
           allocate(darkMatterProfileScaleRadiusConcentration_                   )
           allocate(darkMatterProfileScaleRadiusLudlow2016Analytic_              )
           allocate(darkMatterProfileDMONFW_                                     )
           !![
           <referenceConstruct object="darkMatterParticleCDM_"                                        >
	     <constructor>
	       darkMatterParticleCDM                                        (                                                                                    &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="cosmologyParametersSimple_"                                    >
	     <constructor>
	       cosmologyParametersSimple                                    (                                                                                    &amp;
                &amp;                                                        OmegaMatter                             =OmegaMatter                              , &amp;
                &amp;                                                        OmegaBaryon                             =OmegaBaryon                              , &amp;
                &amp;                                                        OmegaDarkEnergy                         =OmegaDarkEnergy                          , &amp;
                &amp;                                                        temperatureCMB                          =temperatureCMB                           , &amp;
                &amp;                                                        HubbleConstant                          =HubbleConstant                             &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="cosmologyFunctionsMatterLambda_"                               >
	     <constructor>
	       cosmologyFunctionsMatterLambda                               (                                                                                    &amp;
                &amp;                                                        cosmologyParameters_                    =cosmologyParametersSimple_                 &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="linearGrowthCollisionlessMatter_"                              >
	     <constructor>
	       linearGrowthCollisionlessMatter                              (                                                                                    &amp;
                &amp;                                                        cosmologyParameters_                    =cosmologyParametersSimple_               , &amp;
                &amp;                                                        cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_            &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="powerSpectrumPrimordialPowerLaw_"                              >
	     <constructor>
	       powerSpectrumPrimordialPowerLaw                              (                                                                                    &amp;
                &amp;                                                        index_                                  =powerSpectrumIndex                       , &amp;
                &amp;                                                        running                                 =+0.0d0                                   , &amp;
                &amp;                                                        runningRunning                          =+0.0d0                                   , &amp;
                &amp;                                                        wavenumberReference                     =+1.0d0                                   , &amp;
                &amp;                                                        runningSmallScalesOnly                  =.false.                                    &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="transferFunctionEisensteinHu1998_"                             >
	     <constructor>
	       transferFunctionEisensteinHu1998                             (                                                                                    &amp;
                &amp;                                                        darkMatterParticle_                     =darkMatterParticleCDM_                   , &amp;
                &amp;                                                        cosmologyParameters_                    =cosmologyParametersSimple_               , &amp;
                &amp;                                                        cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_            &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="powerSpectrumPrimordialTransferredSimple_"                     >
	     <constructor>
	       powerSpectrumPrimordialTransferredSimple                     (                                                                                    &amp;
                &amp;                                                        powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_         , &amp;
                &amp;                                                        transferFunction_                       =transferFunctionEisensteinHu1998_        , &amp;
                &amp;                                                        linearGrowth_                           =linearGrowthCollisionlessMatter_           &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="powerSpectrumWindowFunctionTopHat_"                            >
	     <constructor>
	       powerSpectrumWindowFunctionTopHat                            (                                                                                    &amp;
                &amp;                                                        cosmologyParameters_                    =cosmologyParametersSimple_                 &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="cosmologicalMassVarianceFilteredPower_"                        >
	     <constructor>
	       cosmologicalMassVarianceFilteredPower                        (                                                                                    &amp;
                &amp;                                                        sigma8                                  =sigma8                                   , &amp;
                &amp;                                                        tolerance                               =1.0d-4                                   , &amp;
                &amp;                                                        toleranceTopHat                         =1.0d-4                                   , &amp;
	        &amp;                                                        rootVarianceLogarithmicGradientTolerance=1.0d-9                                   , &amp;
	        &amp;                                                        integrationFailureIsFatal               =.true.                                   , &amp;
	        &amp;                                                        storeTabulations                        =.true.                                   , &amp;
                &amp;                                                        nonMonotonicIsFatal                     =.true.                                   , &amp;
                &amp;                                                        monotonicInterpolation                  =.false.                                  , &amp;
                &amp;                                                        truncateAtParticleHorizon               =.false.                                  , &amp;
                &amp;                                                        cosmologyParameters_                    =cosmologyParametersSimple_               , &amp;
                &amp;                                                        cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_          , &amp;
                &amp;                                                        linearGrowth_                           =linearGrowthCollisionlessMatter_         , &amp;
                &amp;                                                        powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_, &amp;
                &amp;                                                        powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionTopHat_         &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="powerSpectrumStandard_"                                        >
	     <constructor>
	       powerSpectrumStandard                                        (                                                                                    &amp;
                &amp;                                                        cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_   , &amp;
                &amp;                                                        powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_  &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_">
	     <constructor>
	       criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                                     &amp;
                &amp;                                                        linearGrowth_                           =linearGrowthCollisionlessMatter_         , &amp;
                &amp;                                                        cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_          , &amp;
                &amp;                                                        cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_   , &amp;
                &amp;                                                        darkMatterParticle_                     =darkMatterParticleCDM_                   , &amp;
                &amp;                                                        tableStore                              =.true.                                     &amp;
                &amp;                                                       )
	     </constructor>
	   </referenceConstruct>
	   !!]
           ! Construct the concentration object for this test.
           select case (char(modelName(iModel)))
           case ('Diemer & Kravtsov (2015)'              )
              allocate(darkMatterProfileConcentrationDiemerKravtsov2014 :: darkMatterProfileConcentration_)
           case ('Dutton et al. (2014; 200c)'            )
              allocate(darkMatterProfileConcentrationDuttonMaccio2014   :: darkMatterProfileConcentration_)
           case ('Dutton et al. (2014; vir)'             )
              allocate(darkMatterProfileConcentrationDuttonMaccio2014   :: darkMatterProfileConcentration_)
           case ('Ludlow et al. (2016)'                  )
              allocate(darkMatterProfileConcentrationLudlow2016Fit      :: darkMatterProfileConcentration_)
           case ('Prada et al. (2012)'                   )
              allocate(darkMatterProfileConcentrationPrada2011          :: darkMatterProfileConcentration_)
           case ('Diemer & Joyce (2019)'                 )
              allocate(darkMatterProfileConcentrationDiemerJoyce2019    :: darkMatterProfileConcentration_)
           case ('Diemer & Joyce (2019; vir [converted])')
              allocate(darkMatterProfileConcentrationDiemerJoyce2019    :: darkMatterProfileConcentration_)
           end select
           select type (darkMatterProfileConcentration_)
           type is (darkMatterProfileConcentrationDiemerKravtsov2014)
              !![
              <referenceConstruct object="darkMatterProfileConcentration_">
		<constructor>
		  darkMatterProfileConcentrationDiemerKravtsov2014(                                                                                          &amp;
		   &amp;                                           kappa                    =0.69d0                                                        , &amp;
		   &amp;                                           phi0                     =6.58d0                                                        , &amp;
		   &amp;                                           phi1                     =1.37d0                                                        , &amp;
		   &amp;                                           eta0                     =6.82d0                                                        , &amp;
		   &amp;                                           eta1                     =1.42d0                                                        , &amp;
		   &amp;                                           alpha                    =1.12d0                                                        , &amp;
		   &amp;                                           beta                     =1.69d0                                                        , &amp;
		   &amp;                                           scatter                  =0.00d0                                                        , &amp;
		   &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
		   &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		   &amp;                                           criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_ , &amp;
		   &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                        , &amp;
		   &amp;                                           powerSpectrum_           =powerSpectrumStandard_                                          &amp;
		   &amp;                                          )
		</constructor>
              </referenceConstruct>
              !!]
           type is (darkMatterProfileConcentrationDuttonMaccio2014)
              if (modelName(iModel) == 'Dutton et al. (2014; vir)') then
                 !![
		 <referenceConstruct object="darkMatterProfileConcentration_">
		   <constructor>
		     darkMatterProfileConcentrationDuttonMaccio2014  (                                                                                          &amp;
		      &amp;                                           fitType                  =duttonMaccio2014FitTypeNFWVirial                              , &amp;
		      &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		      &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                                 &amp;
		      &amp;                                          )
		   </constructor>
		 </referenceConstruct>
                 !!]
              else if (modelName(iModel) == 'Dutton et al. (2014; 200c)') then
                 !![
		 <referenceConstruct object="darkMatterProfileConcentration_">
		   <constructor>
		     darkMatterProfileConcentrationDuttonMaccio2014  (                                                                                          &amp;
		      &amp;                                           fitType                  =duttonMaccio2014FitTypeNFWCritical200                         , &amp;
		      &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		      &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                                 &amp;
		      &amp;                                          )
		   </constructor>
		 </referenceConstruct>
                 !!]
              end if
           type is (darkMatterProfileConcentrationLudlow2016Fit)
              !![
              <referenceConstruct object="darkMatterProfileConcentration_">
		<constructor>
		  darkMatterProfileConcentrationLudlow2016Fit     (                                                                                          &amp;
		   &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
		   &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		   &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                          &amp;
		   &amp;                                          )
		</constructor>
              </referenceConstruct>
              !!]
           type is (darkMatterProfileConcentrationPrada2011)
              !![
              <referenceConstruct object="darkMatterProfileConcentration_">
		<constructor>
		  darkMatterProfileConcentrationPrada2011         (                                                                                          &amp;
		   &amp;                                           A                        =2.881d0                                                       , &amp;
		   &amp;                                           B                        =1.257d0                                                       , &amp;
		   &amp;                                           C                        =1.022d0                                                       , &amp;
		   &amp;                                           D                        =0.060d0                                                       , &amp;
		   &amp;                                           C0                       =3.681d0                                                       , &amp;
		   &amp;                                           C1                       =5.033d0                                                       , &amp;
		   &amp;                                           X0                       =0.424d0                                                       , &amp;
		   &amp;                                           X1                       =0.526d0                                                       , &amp;
		   &amp;                                           inverseSigma0            =1.047d0                                                       , &amp;
		   &amp;                                           inverseSigma1            =1.646d0                                                       , &amp;
		   &amp;                                           alpha                    =6.948d0                                                       , &amp;
		   &amp;                                           beta                     =7.386d0                                                       , &amp;
		   &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
		   &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		   &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                          &amp;
		   &amp;                                          )
		</constructor>
              </referenceConstruct>
              !!]
           type is (darkMatterProfileConcentrationDiemerJoyce2019)
              !![
              <referenceConstruct object="darkMatterProfileConcentration_">
		<constructor>
		  darkMatterProfileConcentrationDiemerJoyce2019   (                                                                                          &amp;
		   &amp;                                           kappa                    =0.41d0                                                        , &amp;
		   &amp;                                           a0                       =2.45d0                                                        , &amp;
		   &amp;                                           a1                       =1.82d0                                                        , &amp;
		   &amp;                                           b0                       =3.20d0                                                        , &amp;
		   &amp;                                           b1                       =2.30d0                                                        , &amp;
		   &amp;                                           cAlpha                   =0.21d0                                                        , &amp;
		   &amp;                                           scatter                  =0.00d0                                                        , &amp;
		   &amp;                                           truncateConcentration    =.false.                                                       , &amp;
		   &amp;                                           includeUpturn            =.true.                                                        , &amp;
		   &amp;                                           truncateUpturn           =.false.                                                       , &amp;
		   &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
		   &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
		   &amp;                                           criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_ , &amp;
		   &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                        , &amp;
		   &amp;                                           linearGrowth_            =linearGrowthCollisionlessMatter_                                &amp;
		   &amp;                                          )
		</constructor>
              </referenceConstruct>
              !!]
           end select
           ! Construct the virial density contrast object.
           select case (char(modelDensityContrast(iModel)))
           case ("200c")
              allocate(virialDensityContrastFixed                                     :: virialDensityContrast_)
              select type (virialDensityContrast_)
              type is (virialDensityContrastFixed                                    )
                 !![
		 <referenceConstruct object="virialDensityContrast_">
		   <constructor>
		     virialDensityContrastFixed                                    (                                                            &amp;
		      &amp;                                                         densityContrastValue      =200.0d0                        , &amp;
		      &amp;                                                         densityType               =fixedDensityTypeCritical       , &amp;
		      &amp;                                                         turnAroundOverVirialRadius=  2.0d0                        , &amp;
		      &amp;                                                         cosmologyParameters_      =cosmologyParametersSimple_     , &amp;
		      &amp;                                                         cosmologyFunctions_       =cosmologyFunctionsMatterLambda_  &amp;
		      &amp;                                                        )
		   </constructor>
		 </referenceConstruct>
                 !!]
              end select
           case ("vir" )
              allocate(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt :: virialDensityContrast_)
              select type (virialDensityContrast_)
              type is (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)
                 !![
		 <referenceConstruct object="virialDensityContrast_">
		   <constructor>
		     virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                            &amp;
		      &amp;                                                         tableStore                =.true.                         , &amp;
		      &amp;                                                         cosmologyFunctions_       =cosmologyFunctionsMatterLambda_  &amp;
		      &amp;                                                        )
		   </constructor>
		 </referenceConstruct>
                 !!]
              end select
           end select
           ! Build the scale calculator.
           !![
	   <referenceConstruct object="darkMatterHaloScaleVirialDensityContrastDefinition_">
	     <constructor>
	       darkMatterHaloScaleVirialDensityContrastDefinition(                                                                                          &amp;
                &amp;                                             cosmologyParameters_                =cosmologyParametersSimple_                         , &amp;
                &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsMatterLambda_                    , &amp;
                &amp;                                             virialDensityContrast_              =virialDensityContrast_                               &amp;
                &amp;                                            )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="darkMatterProfileDMONFW_">
	     <constructor>
	       darkMatterProfileDMONFW                           (                                                                                          &amp;
                &amp;                                             velocityDispersionUseSeriesExpansion=.true.                                             , &amp;
                &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleVirialDensityContrastDefinition_  &amp;
                &amp;                                            )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="darkMatterProfileScaleRadiusConcentration_">
	     <constructor>
	       darkMatterProfileScaleRadiusConcentration         (                                                                                          &amp;
	        &amp;                                             correctForConcentrationDefinition   =.true.                                             , &amp;
	        &amp;                                             useMeanConcentration                =.true.                                             , &amp;
                &amp;                                             cosmologyParameters_                =cosmologyParametersSimple_                         , &amp;
                &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsMatterLambda_                    , &amp;
                &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleVirialDensityContrastDefinition_, &amp;
                &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFW_                           , &amp;
                &amp;                                             virialDensityContrast_              =virialDensityContrast_                             , &amp;
                &amp;                                             darkMatterProfileConcentration_     =darkMatterProfileConcentration_                      &amp;
                &amp;                                            )
	     </constructor>
	   </referenceConstruct>
	   <referenceConstruct object="darkMatterProfileScaleRadiusLudlow2016Analytic_">
	     <constructor>
	       darkMatterProfileScaleRadiusLudlow2016Analytic    (                                                                                                    &amp;
	        &amp;                                             C                                   =650.00d0                                                     , &amp;
	        &amp;                                             f                                   =  0.02d0                                                     , &amp;
                &amp;                                             cosmologyParameters_                =cosmologyParametersSimple_                                   , &amp;
                &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsMatterLambda_                              , &amp;
                &amp;                                             darkMatterProfileScaleRadius_       =darkMatterProfileScaleRadiusConcentration_                   , &amp;
                &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleVirialDensityContrastDefinition_          , &amp;
                &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFW_                                     , &amp;
                &amp;                                             virialDensityContrast_              =virialDensityContrast_                                       , &amp;
	        &amp;                                             criticalOverdensity_                =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &amp;
	        &amp;                                             cosmologicalMassVariance_           =cosmologicalMassVarianceFilteredPower_                       , &amp;
                &amp;                                             linearGrowth_                       =linearGrowthCollisionlessMatter_                               &amp;
                &amp;                                            )
	     </constructor>
	   </referenceConstruct>
           !!]
           ! Create a node and assign a time corresponding to z=0.
           node             => treeNode                   (                 )
           basic             => node    %basic            (autoCreate=.true.)
           darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
           call basic%timeSet            (                                                                      &
                &                         cosmologyFunctionsMatterLambda_%cosmicTime                 (          &
                &                         cosmologyFunctionsMatterLambda_%expansionFactorFromRedshift (         &
                &                                                                                      redshift &
                &                                                                                     )         &
                &                                                                                    )          &
                &                        )
           call basic%timeLastIsolatedSet(                                                                      &
                &                         basic                          %time                       (          &
                &                                                                                    )          &
                &                        )
           ! Iterate over masses evaluating concentration.
           do i=1,countMasses
              call basic%massSet(mass(i))
              call Calculations_Reset(node)
              concentration(i)=+darkMatterHaloScaleVirialDensityContrastDefinition_%radiusVirial(node) &
                   &           /darkMatterProfileScaleRadiusConcentration_         %radius      (node)
           end do
           ! Assert the result.
           call Assert(char(modelName(iModel)),concentration,concentrationTarget,relTol=modelTolerance(iModel))
           ! For the Ludlow et al. (2016) model, also check the direct (non-fitting function) calculation.
           if (modelLabel(iModel) == "ludlow16_200c") then
              ! Iterate over masses evaluating concentration.
              do i=1,countMasses
                 call basic%massSet(mass(i))
                 call Calculations_Reset(node)
                 concentration(i)=+darkMatterHaloScaleVirialDensityContrastDefinition_%radiusVirial(node) &
                      &           /darkMatterProfileScaleRadiusLudlow2016Analytic_    %radius      (node)
              end do
              ! Assert the result.
              call Assert(char(modelName(iModel))//" [direct calculation]",concentration,concentrationTarget,relTol=modelTolerance(iModel))
          end if
           ! Clean up.
           call parameters%reset  ()
           call parameters%destroy()
           deallocate(parameters         )
           deallocate(mass               )
           deallocate(concentration      )
           deallocate(concentrationTarget)
           !![
	   <objectDestructor name="darkMatterProfileConcentration_"                              />
	   <objectDestructor name="cosmologyParametersSimple_"                                   />
	   <objectDestructor name="cosmologyFunctionsMatterLambda_"                              />
	   <objectDestructor name="linearGrowthCollisionlessMatter_"                             />
	   <objectDestructor name="cosmologicalMassVarianceFilteredPower_"                       />
	   <objectDestructor name="powerSpectrumWindowFunctionTopHat_"                           />
	   <objectDestructor name="powerSpectrumPrimordialPowerLaw_"                             />
	   <objectDestructor name="transferFunctionEisensteinHu1998_"                            />
	   <objectDestructor name="powerSpectrumPrimordialTransferredSimple_"                    />
	   <objectDestructor name="powerSpectrumStandard_"                                       />
	   <objectDestructor name="darkMatterParticleCDM_"                                       />
	   <objectDestructor name="criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_"/>
	   <objectDestructor name="darkMatterProfileScaleRadiusConcentration_"                   />
	   <objectDestructor name="darkMatterHaloScaleVirialDensityContrastDefinition_"          />
	   <objectDestructor name="darkMatterProfileDMONFW_"                                     />
	   <objectDestructor name="virialDensityContrast_"                                       />
           !!]
           call Node_Components_Thread_Uninitialize()
           call Node_Components_Uninitialize       ()
           call nodeClassHierarchyFinalize         ()
        end do
        call Unit_Tests_End_Group()
     end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Concentrations
