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

!% Contains a program which tests concentration models.

program Test_Concentrations
  !% Tests concentration models.
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower   , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationClass     , darkMatterProfileConcentrationDiemerKravtsov2014             , darkMatterProfileConcentrationDuttonMaccio2014, darkMatterProfileConcentrationLudlow2016Fit, &
          &                                           darkMatterProfileConcentrationPrada2011
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : Count_Lines_in_File
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Calculations_Resets      , only : Galacticus_Calculations_Reset
  use :: Galacticus_Display                  , only : Galacticus_Verbosity_Level_Set          , verbosityStandard
  use :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize            , nodeComponentBasic                                           , treeNode
  use :: Galacticus_Paths                    , only : galacticusPath                          , pathTypeExec
  use :: ISO_Varying_String                  , only : assignment(=)                           , char                                                         , operator(//)                                  , operator(==)                               , &
          &                                           varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize              , Node_Components_Thread_Initialize                            , Node_Components_Thread_Uninitialize           , Node_Components_Uninitialize
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: String_Handling                     , only : String_Split_Words
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group                                       , Unit_Tests_End_Group                          , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                    ), pointer                             :: node
  class           (nodeComponentBasic                                          ), pointer                             :: basic
  class           (darkMatterProfileConcentrationClass                         ), pointer                             :: darkMatterProfileConcentration_
  type            (cosmologyParametersSimple                                   )                                      :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )                                      :: cosmologyFunctionsMatterLambda_
  type            (linearGrowthCollisionlessMatter                             )                                      :: linearGrowthCollisionlessMatter_
  type            (cosmologicalMassVarianceFilteredPower                       )                                      :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionTopHat                           )                                      :: powerSpectrumWindowFunctionTopHat_
  type            (powerSpectrumPrimordialPowerLaw                             )                                      :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionEisensteinHu1999                            )                                      :: transferFunctionEisensteinHu1999_
  type            (powerSpectrumPrimordialTransferredSimple                    )                                      :: powerSpectrumPrimordialTransferredSimple_
  type            (powerSpectrumStandard                                       )                                      :: powerSpectrumStandard_
  type            (darkMatterParticleCDM                                       )                                      :: darkMatterParticleCDM_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                                      :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (varying_string                                              )                                      :: parameterFile
  integer                                                                       , parameter                           :: countModels                                                   =5
  type            (varying_string                                              ), dimension(countModels)              :: modelName                                                       , modelLabel         , &
       &                                                                                                                 modelDensityContrast
  double precision                                                              , dimension(countModels)              :: modelTolerance
  type            (inputParameters                                             ), pointer                             :: parameters
  character       (len= 128                                                    ), dimension(4          )              :: columns
  double precision                                                              , dimension(:          ), allocatable :: mass                                                            , concentrationTarget, &
       &                                                                                                                 concentration
  character       (len=1024                                                    )                                      :: line
  double precision                                                                                                    :: OmegaMatter                                                     , OmegaBaryon        , &
       &                                                                                                                 OmegaDarkEnergy                                                 , temperatureCMB     , &
       &                                                                                                                 powerSpectrumIndex                                              , HubbleConstant     , &
       &                                                                                                                 neutrinoNumberEffective                                         , sigma8
  integer                                                                                                             :: colossusFile                                                    , status             , &
       &                                                                                                                 countMasses                                                     , i                  , &
       &                                                                                                                 iModel

  ! Initialize.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Specify all models to run.
  modelName           (1)='Diemer & Kravtsov (2015)'
  modelLabel          (1)='diemer15_orig_200c'
  modelDensityContrast(1)='200c'
  modelTolerance      (1)=1.5d-2
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
  ! Iterate over models.
  call Unit_Tests_Begin_Group("Concentration algorithms")
  do iModel=1,countModels
     allocate(parameters)
     parameterFile=galacticusPath(pathTypeExec)//'testSuite/parameters/concentrations_'//modelDensityContrast(iModel)//'.xml'
     parameters   =inputParameters(parameterFile)
     call parameters%markGlobal()
     call nodeClassHierarchyInitialize     (parameters)
     call Node_Components_Initialize       (parameters)
     call Node_Components_Thread_Initialize(parameters)
     ! Read the Colossus target data (and parameters used) from file.
     countMasses=Count_Lines_in_File(galacticusPath(pathTypeExec)//'testSuite/data/concentrationsColossus/'//char(modelLabel(iModel))//'.txt','#')
     i=0
     allocate(mass(               countMasses))
     allocate(concentration      (countMasses))
     allocate(concentrationTarget(countMasses))
     open(newUnit=colossusFile,file='testSuite/data/concentrationsColossus/'//char(modelLabel(iModel))//'.txt',status='old',form='formatted')
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
     !# <referenceConstruct object="darkMatterParticleCDM_"                                        >
     !#  <constructor>
     !#   darkMatterParticleCDM                                        (                                                                               &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="cosmologyParametersSimple_"                                    >
     !#  <constructor>
     !#   cosmologyParametersSimple                                    (                                                                               &amp;
     !#    &amp;                                                        OmegaMatter                        =OmegaMatter                              , &amp;
     !#    &amp;                                                        OmegaBaryon                        =OmegaBaryon                              , &amp;
     !#    &amp;                                                        OmegaDarkEnergy                    =OmegaDarkEnergy                          , &amp;
     !#    &amp;                                                        temperatureCMB                     =temperatureCMB                           , &amp;
     !#    &amp;                                                        HubbleConstant                     =HubbleConstant                             &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="cosmologyFunctionsMatterLambda_"                               >
     !#  <constructor>
     !#   cosmologyFunctionsMatterLambda                               (                                                                               &amp;
     !#    &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_                 &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="linearGrowthCollisionlessMatter_"                              >
     !#  <constructor>
     !#   linearGrowthCollisionlessMatter                              (                                                                               &amp;
     !#    &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
     !#    &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_            &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="powerSpectrumPrimordialPowerLaw_"                              >
     !#  <constructor>
     !#   powerSpectrumPrimordialPowerLaw                              (                                                                               &amp;
     !#    &amp;                                                        index                              =powerSpectrumIndex                       , &amp;
     !#    &amp;                                                        running                            =+0.0d0                                   , &amp;
     !#    &amp;                                                        wavenumberReference                =+1.0d0                                     &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="transferFunctionEisensteinHu1999_"                             >
     !#  <constructor>
     !#   transferFunctionEisensteinHu1999                             (                                                                               &amp;
     !#    &amp;                                                        neutrinoNumberEffective            =neutrinoNumberEffective                  , &amp;
     !#    &amp;                                                        neutrinoMassSummed                 =0.0d0                                    , &amp;
     !#    &amp;                                                        darkMatterParticle_                =darkMatterParticleCDM_                   , &amp;
     !#    &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
     !#    &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_            &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="powerSpectrumPrimordialTransferredSimple_"                     >
     !#  <constructor>
     !#   powerSpectrumPrimordialTransferredSimple                     (                                                                               &amp;
     !#    &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw_         , &amp;
     !#    &amp;                                                        transferFunction_                  =transferFunctionEisensteinHu1999_        , &amp;
     !#    &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_           &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="powerSpectrumWindowFunctionTopHat_"                            >
     !#  <constructor>
     !#   powerSpectrumWindowFunctionTopHat                            (                                                                               &amp;
     !#    &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_                 &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="cosmologicalMassVarianceFilteredPower_"                        >
     !#  <constructor>
     !#   cosmologicalMassVarianceFilteredPower                        (                                                                               &amp;
     !#    &amp;                                                        sigma8                             =sigma8                                   , &amp;
     !#    &amp;                                                        tolerance                          =1.0d-4                                   , &amp;
     !#    &amp;                                                        toleranceTopHat                    =1.0d-4                                   , &amp;
     !#    &amp;                                                        nonMonotonicIsFatal                =.true.                                   , &amp;
     !#    &amp;                                                        monotonicInterpolation             =.false.                                  , &amp;
     !#    &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
     !#    &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_          , &amp;
     !#    &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_         , &amp;
     !#    &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple_, &amp;
     !#    &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionTopHat_         &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="powerSpectrumStandard_"                                        >
     !#  <constructor>
     !#   powerSpectrumStandard                                        (                                                                               &amp;
     !#    &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower_   , &amp;
     !#    &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple_  &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     !# <referenceConstruct object="criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_">
     !#  <constructor>
     !#   criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                               &amp;
     !#    &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_         , &amp;
     !#    &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_          , &amp;
     !#    &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower_   , &amp;
     !#    &amp;                                                        darkMatterParticle_                =darkMatterParticleCDM_                   , &amp;
     !#    &amp;                                                        tableStore                         =.true.                                     &amp;
     !#    &amp;                                                       )
     !#  </constructor>
     !# </referenceConstruct>
     ! Construct the concentration object for this test.
     select case (char(modelName(iModel)))
     case ('Diemer & Kravtsov (2015)'  )
        allocate(darkMatterProfileConcentrationDiemerKravtsov2014 :: darkMatterProfileConcentration_)
     case ('Dutton et al. (2014; 200c)')
        allocate(darkMatterProfileConcentrationDuttonMaccio2014   :: darkMatterProfileConcentration_)
     case ('Dutton et al. (2014; vir)' )
        allocate(darkMatterProfileConcentrationDuttonMaccio2014   :: darkMatterProfileConcentration_)
     case ('Ludlow et al. (2016)'      )
        allocate(darkMatterProfileConcentrationLudlow2016Fit      :: darkMatterProfileConcentration_)
     case ('Prada et al. (2012)'       )
        allocate(darkMatterProfileConcentrationPrada2011          :: darkMatterProfileConcentration_)
     end select
     select type (darkMatterProfileConcentration_)
     type is (darkMatterProfileConcentrationDiemerKravtsov2014)
        !#    <referenceConstruct object="darkMatterProfileConcentration_">
        !#     <constructor>
        !#      darkMatterProfileConcentrationDiemerKravtsov2014(                                                                                          &amp;
        !#       &amp;                                           kappa                    =0.69d0                                                        , &amp;
        !#       &amp;                                           phi0                     =6.58d0                                                        , &amp;
        !#       &amp;                                           phi1                     =1.37d0                                                        , &amp;
        !#       &amp;                                           eta0                     =6.82d0                                                        , &amp;
        !#       &amp;                                           eta1                     =1.42d0                                                        , &amp;
        !#       &amp;                                           alpha                    =1.12d0                                                        , &amp;
        !#       &amp;                                           beta                     =1.69d0                                                        , &amp;
        !#       &amp;                                           scatter                  =0.00d0                                                        , &amp;
        !#       &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
        !#       &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
        !#       &amp;                                           criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &amp;
        !#       &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                        , &amp;
        !#       &amp;                                           powerSpectrum_           =powerSpectrumStandard_                                          &amp;
        !#       &amp;                                          )
        !#     </constructor>
        !#    </referenceConstruct>
     type is (darkMatterProfileConcentrationDuttonMaccio2014)
        if (modelName(iModel) == 'Dutton et al. (2014; vir)') then
           !# <referenceConstruct object="darkMatterProfileConcentration_">
           !#  <constructor>
           !#   darkMatterProfileConcentrationDuttonMaccio2014  (                                                                                          &amp;
           !#    &amp;                                           fitType                  ='nfwVirial'                                                   , &amp;
           !#    &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
           !#    &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                                 &amp;
           !#    &amp;                                          )
           !#  </constructor>
           !# </referenceConstruct>
        else if (modelName(iModel) == 'Dutton et al. (2014; 200c)') then
           !# <referenceConstruct object="darkMatterProfileConcentration_">
           !#  <constructor>
           !#   darkMatterProfileConcentrationDuttonMaccio2014  (                                                                                          &amp;
           !#    &amp;                                           fitType                  ='nfwCritical200'                                              , &amp;
           !#    &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
           !#    &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                                 &amp;
           !#    &amp;                                          )
           !#  </constructor>
           !# </referenceConstruct>
        end if
     type is (darkMatterProfileConcentrationLudlow2016Fit)
        !#    <referenceConstruct object="darkMatterProfileConcentration_">
        !#     <constructor>
        !#      darkMatterProfileConcentrationLudlow2016Fit     (                                                                                          &amp;
        !#       &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
        !#       &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
        !#       &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                          &amp;
        !#       &amp;                                          )
        !#     </constructor>
        !#    </referenceConstruct>
     type is (darkMatterProfileConcentrationPrada2011)
        !#    <referenceConstruct object="darkMatterProfileConcentration_">
        !#     <constructor>
        !#      darkMatterProfileConcentrationPrada2011         (                                                                                          &amp;
        !#       &amp;                                           A                        =2.881d0                                                       , &amp;
        !#       &amp;                                           B                        =1.257d0                                                       , &amp;
        !#       &amp;                                           C                        =1.022d0                                                       , &amp;
        !#       &amp;                                           D                        =0.060d0                                                       , &amp;
        !#       &amp;                                           C0                       =3.681d0                                                       , &amp;
        !#       &amp;                                           C1                       =5.033d0                                                       , &amp;
        !#       &amp;                                           X0                       =0.424d0                                                       , &amp;
        !#       &amp;                                           X1                       =0.526d0                                                       , &amp;
        !#       &amp;                                           inverseSigma0            =1.047d0                                                       , &amp;
        !#       &amp;                                           inverseSigma1            =1.646d0                                                       , &amp;
        !#       &amp;                                           alpha                    =6.948d0                                                       , &amp;
        !#       &amp;                                           beta                     =7.386d0                                                       , &amp;
        !#       &amp;                                           cosmologyFunctions_      =cosmologyFunctionsMatterLambda_                               , &amp;
        !#       &amp;                                           cosmologyParameters_     =cosmologyParametersSimple_                                    , &amp;
        !#       &amp;                                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                          &amp;
        !#       &amp;                                          )
        !#     </constructor>
        !#    </referenceConstruct>
     end select
     ! Create a node and assign a time corresponding to z=0.
     node  => treeNode      (                 )
     basic => node    %basic(autoCreate=.true.)
     call basic%timeSet            (                                                                   &
          &                         cosmologyFunctionsMatterLambda_%cosmicTime                 (       &
          &                         cosmologyFunctionsMatterLambda_%expansionFactorFromRedshift (      &
          &                                                                                      0.0d0 &
          &                                                                                     )      &
          &                                                                                    )       &
          &                        )
     call basic%timeLastIsolatedSet(                                                                   &
          &                         basic                          %time                       (       &
          &                                                                                    )       &
          &                        )
     ! Iterate over masses evaluating concentration.
     do i=1,countMasses
        call basic%massSet(mass(i))
        call Galacticus_Calculations_Reset(node)
        concentration(i)=darkMatterProfileConcentration_%concentration(node)
     end do
     ! Assert the result.
     call Assert(char(modelName(iModel)),concentration,concentrationTarget,relTol=modelTolerance(iModel))
     ! Clean up.
     call parameters%reset  ()
     call parameters%destroy()
     deallocate(parameters)
     deallocate(mass               )
     deallocate(concentration      )
     deallocate(concentrationTarget)
     !# <objectDestructor name="darkMatterProfileConcentration_"/>
     call Node_Components_Thread_Uninitialize()
     call Node_Components_Uninitialize       ()
     call Galacticus_Function_Classes_Destroy()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Concentrations
