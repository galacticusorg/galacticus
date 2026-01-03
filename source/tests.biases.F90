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
Contains a program which tests halo bias models.
!!}

program Test_Biases
  !!{
  Tests concentration models.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower   , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Biases             , only : darkMatterHaloBiasClass                 , darkMatterHaloBiasPressSchechter                            , darkMatterHaloBiasSheth2001        , darkMatterHaloBiasTinker2010
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : Count_Lines_in_File
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Calculations_Resets                 , only : Calculations_Reset
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize            , nodeComponentBasic                                          , treeNode                           , nodeClassHierarchyFinalize
  use :: Input_Paths                         , only : inputPath                               , pathTypeExec
  use :: ISO_Varying_String                  , only : assignment(=)                           , char                                                        , operator(//)                       , operator(==)                , &
          &                                           varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize              , Node_Components_Thread_Initialize                           , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: String_Handling                     , only : String_Split_Words
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Virial_Density_Contrast             , only : virialDensityContrastFixed              , fixedDensityTypeCritical
  implicit none
  type            (treeNode                                                    ), pointer                             :: node
  class           (nodeComponentBasic                                          ), pointer                             :: basic
  class           (darkMatterHaloBiasClass                                     ), pointer                             :: darkMatterHaloBias_
  type            (virialDensityContrastFixed                                  ), pointer                             :: virialDensityContrast_
  type            (cosmologyParametersSimple                                   ), pointer                             :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              ), pointer                             :: cosmologyFunctionsMatterLambda_
  type            (linearGrowthCollisionlessMatter                             ), pointer                             :: linearGrowthCollisionlessMatter_
  type            (cosmologicalMassVarianceFilteredPower                       ), pointer                             :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionTopHat                           ), pointer                             :: powerSpectrumWindowFunctionTopHat_
  type            (powerSpectrumPrimordialPowerLaw                             ), pointer                             :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionEisensteinHu1999                            ), pointer                             :: transferFunctionEisensteinHu1999_
  type            (powerSpectrumPrimordialTransferredSimple                    ), pointer                             :: powerSpectrumPrimordialTransferredSimple_
  type            (darkMatterParticleCDM                                       ), pointer                             :: darkMatterParticleCDM_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer                             :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (varying_string                                              )                                      :: parameterFile
  integer                                                                       , parameter                           :: countModels                                                  =3
  type            (varying_string                                              ), dimension(countModels)              :: modelName                                                      , modelLabel    , &
       &                                                                                                                 modelDensityContrast
  double precision                                                              , dimension(countModels)              :: modelTolerance
  type            (inputParameters                                             ), pointer                             :: parameters
  character       (len= 128                                                    ), dimension(4          )              :: columns
  double precision                                                              , dimension(:          ), allocatable :: mass                                                           , biasTarget    , &
       &                                                                                                                 bias                                                           , redshift
  character       (len=1024                                                    )                                      :: line
  double precision                                                                                                    :: OmegaMatter                                                    , OmegaBaryon   , &
       &                                                                                                                 OmegaDarkEnergy                                                , temperatureCMB, &
       &                                                                                                                 powerSpectrumIndex                                             , HubbleConstant, &
       &                                                                                                                 neutrinoNumberEffective                                        , sigma8
  integer                                                                                                             :: colossusFile                                                   , status        , &
       &                                                                                                                 countMasses                                                    , i             , &
       &                                                                                                                 iModel

  ! Initialize.
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Specify all models to run.
  modelName           (1)='Press-Schechter'
  modelLabel          (1)='cole89_NA'
  modelDensityContrast(1)='NA'
  modelTolerance      (1)=2.0d-2
  modelName           (2)='Sheth et al. (2001)'
  modelLabel          (2)='sheth01_NA'
  modelDensityContrast(2)='NA'
  modelTolerance      (2)=1.0d-2
  modelName           (3)='Tinker et al. (2010; 200c)'
  modelLabel          (3)='tinker10_200c'
  modelDensityContrast(3)='200c'
  modelTolerance      (3)=3.0d-2
  ! Iterate over models.
  call Unit_Tests_Begin_Group("Bias algorithms")
  do iModel=1,countModels
     allocate(parameters)
     parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/haloBias_'//modelDensityContrast(iModel)//'.xml'
     parameters   =inputParameters(parameterFile)
     call nodeClassHierarchyInitialize     (parameters)
     call Node_Components_Initialize       (parameters)
     call Node_Components_Thread_Initialize(parameters)
     ! Read the Colossus target data (and parameters used) from file.
     countMasses=Count_Lines_in_File(inputPath(pathTypeExec)//'testSuite/data/haloBiasesColossus/'//char(modelLabel(iModel))//'.txt','#')
     i=0
     allocate(mass      (countMasses))
     allocate(redshift  (countMasses))
     allocate(bias      (countMasses))
     allocate(biasTarget(countMasses))
     open(newUnit=colossusFile,file='testSuite/data/haloBiasesColossus/'//char(modelLabel(iModel))//'.txt',status='old',form='formatted')
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
           read (line,*) mass(i),redshift(i),biasTarget(i)
        end if
     end do
     close(colossusFile)
     ! Convert masses from h⁻¹M☉ to M☉.
     mass=mass/(HubbleConstant/100.0d0)
     ! Construct all required objects.
     allocate(virialDensityContrast_                                       )
     allocate(cosmologyParametersSimple_                                   )
     allocate(cosmologyFunctionsMatterLambda_                              )
     allocate(linearGrowthCollisionlessMatter_                             )
     allocate(cosmologicalMassVarianceFilteredPower_                       )
     allocate(powerSpectrumWindowFunctionTopHat_                           )
     allocate(powerSpectrumPrimordialPowerLaw_                             )
     allocate(transferFunctionEisensteinHu1999_                            )
     allocate(powerSpectrumPrimordialTransferredSimple_                    )
     allocate(darkMatterParticleCDM_                                       )
     allocate(criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
     !![
     <referenceConstruct object="darkMatterParticleCDM_"                                        >
      <constructor>
       darkMatterParticleCDM                                        (                                                                               &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="cosmologyParametersSimple_"                                    >
      <constructor>
       cosmologyParametersSimple                                    (                                                                               &amp;
        &amp;                                                        OmegaMatter                        =OmegaMatter                              , &amp;
        &amp;                                                        OmegaBaryon                        =OmegaBaryon                              , &amp;
        &amp;                                                        OmegaDarkEnergy                    =OmegaDarkEnergy                          , &amp;
        &amp;                                                        temperatureCMB                     =temperatureCMB                           , &amp;
        &amp;                                                        HubbleConstant                     =HubbleConstant                             &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="cosmologyFunctionsMatterLambda_"                               >
      <constructor>
       cosmologyFunctionsMatterLambda                               (                                                                               &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_                 &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="linearGrowthCollisionlessMatter_"                              >
      <constructor>
       linearGrowthCollisionlessMatter                              (                                                                               &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
        &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_            &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="powerSpectrumPrimordialPowerLaw_"                              >
      <constructor>
       powerSpectrumPrimordialPowerLaw                              (                                                                               &amp;
        &amp;                                                        index_                             =powerSpectrumIndex                       , &amp;
        &amp;                                                        running                            =+0.0d0                                   , &amp;
        &amp;                                                        runningRunning                     =+0.0d0                                   , &amp;
        &amp;                                                        wavenumberReference                =+1.0d0                                   , &amp;
        &amp;                                                        runningSmallScalesOnly             =.false.                                    &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="transferFunctionEisensteinHu1999_"                             >
      <constructor>
       transferFunctionEisensteinHu1999                             (                                                                               &amp;
        &amp;                                                        neutrinoNumberEffective            =neutrinoNumberEffective                  , &amp;
        &amp;                                                        neutrinoMassSummed                 =0.0d0                                    , &amp;
        &amp;                                                        darkMatterParticle_                =darkMatterParticleCDM_                   , &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
        &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_            &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="powerSpectrumPrimordialTransferredSimple_"                     >
      <constructor>
       powerSpectrumPrimordialTransferredSimple                     (                                                                               &amp;
        &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw_         , &amp;
        &amp;                                                        transferFunction_                  =transferFunctionEisensteinHu1999_        , &amp;
        &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_           &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="powerSpectrumWindowFunctionTopHat_"                            >
      <constructor>
       powerSpectrumWindowFunctionTopHat                            (                                                                               &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_                 &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="cosmologicalMassVarianceFilteredPower_"                        >
      <constructor>
       cosmologicalMassVarianceFilteredPower                        (                                                                               &amp;
        &amp;                                                        sigma8                             =sigma8                                   , &amp;
        &amp;                                                        tolerance                          =1.0d-4                                   , &amp;
        &amp;                                                        toleranceTopHat                    =1.0d-4                                   , &amp;
        &amp;                                                        nonMonotonicIsFatal                =.true.                                   , &amp;
        &amp;                                                        monotonicInterpolation             =.false.                                  , &amp;
        &amp;                                                        truncateAtParticleHorizon          =.false.                                  , &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
        &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_          , &amp;
        &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_         , &amp;
        &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple_, &amp;
        &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionTopHat_         &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_">
      <constructor>
       criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                               &amp;
        &amp;                                                        linearGrowth_                      =linearGrowthCollisionlessMatter_         , &amp;
        &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_          , &amp;
        &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower_   , &amp;
        &amp;                                                        darkMatterParticle_                =darkMatterParticleCDM_                   , &amp;
        &amp;                                                        tableStore                         =.true.                                     &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     <referenceConstruct object="virialDensityContrast_">
      <constructor>
       virialDensityContrastFixed                                   (                                                                               &amp;
        &amp;                                                        densityContrastValue               =200.0d0                                  , &amp;
        &amp;                                                        densityType                        =fixedDensityTypeCritical                 , &amp;
        &amp;                                                        turnAroundOverVirialRadius         =  2.0d0                                  , &amp;
        &amp;                                                        cosmologyParameters_               =cosmologyParametersSimple_               , &amp;
        &amp;                                                        cosmologyFunctions_                =cosmologyFunctionsMatterLambda_            &amp;
        &amp;                                                       )
      </constructor>
     </referenceConstruct>
     !!]
     ! Construct the concentration object for this test.
     select case (char(modelName(iModel)))
     case ('Press-Schechter')
        allocate(darkMatterHaloBiasPressSchechter :: darkMatterHaloBias_)
     case ('Sheth et al. (2001)')
        allocate(darkMatterHaloBiasSheth2001      :: darkMatterHaloBias_)
     case ('Tinker et al. (2010; 200c)')
        allocate(darkMatterHaloBiasTinker2010     :: darkMatterHaloBias_)
     end select
     select type (darkMatterHaloBias_)
     type is (darkMatterHaloBiasPressSchechter)
        !![
        <referenceConstruct object="darkMatterHaloBias_">
         <constructor>
          darkMatterHaloBiasPressSchechter(                                                                                         &amp;
           &amp;                           criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &amp;
           &amp;                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                         &amp;
           &amp;                          )
         </constructor>
        </referenceConstruct>
        !!]
     type is (darkMatterHaloBiasSheth2001)
        !![
        <referenceConstruct object="darkMatterHaloBias_">
         <constructor>
          darkMatterHaloBiasSheth2001     (                                                                                         &amp;
              &amp;                        criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &amp;
              &amp;                        cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                         &amp;
           &amp;                          )
         </constructor>
        </referenceConstruct>
        !!]
     type is (darkMatterHaloBiasTinker2010)
        !![
        <referenceConstruct object="darkMatterHaloBias_">
         <constructor>
          darkMatterHaloBiasTinker2010    (                                                                                         &amp;
           &amp;                           criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &amp;
           &amp;                           cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                       , &amp;
           &amp;                           virialDensityContrast_   =virialDensityContrast_                                         &amp;
           &amp;                          )
         </constructor>
        </referenceConstruct>
        !!]
     end select
     ! Create a node.
     node  => treeNode      (                 )
     basic => node    %basic(autoCreate=.true.)
     ! Iterate over masses evaluating bias.
     do i=1,countMasses
        call basic%timeSet            (                                                                         &
             &                         cosmologyFunctionsMatterLambda_%cosmicTime                 (             &
             &                         cosmologyFunctionsMatterLambda_%expansionFactorFromRedshift (            &
             &                                                                                      redshift(i) &
             &                                                                                     )            &
             &                                                                                    )             &
             &                        )
        call basic%timeLastIsolatedSet(                                                                         &
             &                         basic                          %time                       (             &
             &                                                                                    )             &
             &                        )
        call basic%massSet            (                                                                         &
             &                                                                                      mass    (i) &
             &                        )
        call Calculations_Reset(node)
        bias(i)=darkMatterHaloBias_%bias(node)
     end do
     ! Assert the result.
     call Assert(char(modelName(iModel)),bias,biasTarget,relTol=modelTolerance(iModel))
     ! Clean up.
     call parameters%reset  ()
     call parameters%destroy()
     deallocate(parameters)
     deallocate(mass      )
     deallocate(redshift  )
     deallocate(bias      )
     deallocate(biasTarget)
     !![
     <objectDestructor name="darkMatterHaloBias_"                                          />
     <objectDestructor name="virialDensityContrast_"                                       />
     <objectDestructor name="cosmologyParametersSimple_"                                   />
     <objectDestructor name="cosmologyFunctionsMatterLambda_"                              />
     <objectDestructor name="linearGrowthCollisionlessMatter_"                             />
     <objectDestructor name="cosmologicalMassVarianceFilteredPower_"                       />
     <objectDestructor name="powerSpectrumWindowFunctionTopHat_"                           />
     <objectDestructor name="powerSpectrumPrimordialPowerLaw_"                             />
     <objectDestructor name="transferFunctionEisensteinHu1999_"                            />
     <objectDestructor name="powerSpectrumPrimordialTransferredSimple_"                    />
     <objectDestructor name="darkMatterParticleCDM_"                                       />
     <objectDestructor name="criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_"/>
     !!]
     call Node_Components_Thread_Uninitialize()
     call Node_Components_Uninitialize       ()
     call nodeClassHierarchyFinalize         ()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Biases
