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
Contains a program which tests that environmental-averaging of the peak-background split Press-Schechter mass function behaves as
expected.
!!}

program Tests_Halo_Mass_Function_Environmental_Average
  !!{
  Tests that environmental-averaging of the peak-background split Press-Schechter mass function behaves as expected.
  !!}
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, cosmologicalMassVarianceFilteredPower  , haloEnvironmentNormal                  , criticalOverdensityPeakBackgroundSplit, &
       &                                              cosmologicalMassVariancePeakBackgroundSplit
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple                                    , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                                          , verbosityLevelStandard
  use :: Excursion_Sets_Barriers             , only : excursionSetBarrierCriticalOverdensity
  use :: Excursion_Sets_First_Crossings      , only : excursionSetFirstCrossingFarahiMidpoint                      , excursionSetFirstCrossingLinearBarrier
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : treeNode                                                     , mergerTree                            , nodeComponentBasic                     , nodeClassHierarchyInitialize
  use :: Halo_Mass_Functions                 , only : haloMassFunctionPressSchechter                               , haloMassFunctionEnvironmentAveraged
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize                                   , Node_Components_Thread_Initialize     , Node_Components_Thread_Uninitialize    , Node_Components_Uninitialize
  use :: Numerical_Ranges                    , only : Make_Range, rangeTypeLogarithmic
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionCAMB
  use :: Unit_Tests                          , only : Assert                                                      , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                                                                              :: i
  double precision                                                                                     :: time
  integer                                                                       , parameter            :: massCount                           =13
  double precision                                                              , parameter            :: massMinimum                         =1.0d6, massMaximum                    =1.0d12
  double precision                                                              , dimension(massCount) :: mass                                      , massFunction                          , &
       &                                                                                                  massFunctionEnvironmentAveraged
  type            (cosmologyParametersSimple                                   ), pointer              :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              ), pointer              :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower                       ), pointer              :: cosmologicalMassVariance_
  type            (cosmologicalMassVariancePeakBackgroundSplit                 ), pointer              :: cosmologicalMassVarianceConditioned_
  type            (linearGrowthCollisionlessMatter                             ), pointer              :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat                           ), pointer              :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             ), pointer              :: powerSpectrumPrimordial_
  type            (transferFunctionCAMB                                        ), pointer              :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    ), pointer              :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                       ), pointer              :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer              :: criticalOverdensity_
  type            (criticalOverdensityPeakBackgroundSplit                      ), pointer              :: criticalOverdensityConditioned_
  type            (haloMassFunctionPressSchechter                              ), pointer              :: haloMassFunction_                   , haloMassFunctionConditioned_
  type            (haloMassFunctionEnvironmentAveraged                         ), pointer              :: haloMassFunctionEnvironmentAveraged_
  type            (haloEnvironmentNormal                                       ), pointer              :: haloEnvironment_
  type            (excursionSetBarrierCriticalOverdensity                      ), pointer              :: excursionSetBarrier_                , excursionSetBarrierConditioned_
  type            (excursionSetFirstCrossingLinearBarrier                      ), pointer              :: excursionSetFirstCrossing_          , excursionSetFirstCrossingConditioned_
  type            (treeNode                                                    ), pointer              :: node
  type            (mergerTree                                                  )                       :: tree
  class           (nodeComponentBasic                                          ), pointer              :: basic
  type            (inputParameters                                             )                       :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize.
  parameters=inputParameters()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Halo mass function: environment-averaged")
  ! Construct required objects.
  allocate(cosmologyParameters_                 )
  allocate(cosmologyFunctions_                  )
  allocate(cosmologicalMassVariance_            )
  allocate(cosmologicalMassVarianceConditioned_ )
  allocate(linearGrowth_                        )
  allocate(powerSpectrumWindowFunction_         )
  allocate(powerSpectrumPrimordial_             )
  allocate(transferFunction_                    )
  allocate(powerSpectrumPrimordialTransferred_  )
  allocate(darkMatterParticle_                  )
  allocate(criticalOverdensity_                 )
  allocate(criticalOverdensityConditioned_      )
  allocate(haloMassFunction_                    )
  allocate(haloMassFunctionConditioned_         )
  allocate(haloMassFunctionEnvironmentAveraged_ )
  allocate(haloEnvironment_                     )
  allocate(excursionSetBarrier_                 )
  allocate(excursionSetBarrierConditioned_      )
  allocate(excursionSetFirstCrossing_           )
  allocate(excursionSetFirstCrossingConditioned_)
  !![
  <referenceConstruct object="darkMatterParticle_"                 >
   <constructor>
    darkMatterParticleCDM                                       (                                                                           &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"                >
   <constructor>
    cosmologyParametersSimple                                   (                                                                           &amp;
     &amp;                                                       OmegaMatter                        = 0.238d0                             , &amp;
     &amp;                                                       OmegaBaryon                        = 0.045d0                             , &amp;
     &amp;                                                       OmegaDarkEnergy                    = 0.762d0                             , &amp;
     &amp;                                                       temperatureCMB                     = 2.700d0                             , &amp;
     &amp;                                                       HubbleConstant                     =70.000d0                               &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                 >
   <constructor>
    cosmologyFunctionsMatterLambda                              (                                                                           &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  !!]
  time=cosmologyFunctions_%cosmicTime(1.0d0)
  !![  
  <referenceConstruct object="linearGrowth_"                       >
   <constructor>
    linearGrowthCollisionlessMatter                             (                                                                           &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                    &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"            >
   <constructor>
    powerSpectrumPrimordialPowerLaw                             (                                                                           &amp;
     &amp;                                                       index_                             =+0.951d0                             , &amp;
     &amp;                                                       running                            =+0.000d0                             , &amp;
     &amp;                                                       runningRunning                     =+0.000d0                             , &amp;
     &amp;                                                       wavenumberReference                =+1.000d0                             , &amp;
     &amp;                                                       runningSmallScalesOnly             =.false.                                &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                   >
   <constructor>
    transferFunctionCAMB                                        (                                                                           &amp;
     &amp;                                                       redshift                           =0.000d0                              , &amp;
     &amp;                                                       cambCountPerDecade                 =0                                    , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                  , &amp;
     &amp;                                                       darkMatterParticle_                =darkMatterParticle_                    &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_" >
   <constructor>
    powerSpectrumPrimordialTransferredSimple                    (                                                                           &amp;
     &amp;                                                       powerSpectrumPrimordial_           =powerSpectrumPrimordial_             , &amp;
     &amp;                                                       transferFunction_                  =transferFunction_                    , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                          &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"        >
   <constructor>
    powerSpectrumWindowFunctionTopHat                           (                                                                           &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"           >
   <constructor>
    cosmologicalMassVarianceFilteredPower                       (                                                                           &amp;
     &amp;                                                       sigma8                             =0.74d+0                              , &amp;
     &amp;                                                       tolerance                          =1.00d-4                              , &amp;
     &amp;                                                       toleranceTopHat                    =1.00d-4                              , &amp;
     &amp;                                                       nonMonotonicIsFatal                =.true.                               , &amp;
     &amp;                                                       monotonicInterpolation             =.false.                              , &amp;
     &amp;                                                       truncateAtParticleHorizon          =.false.                              , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                  , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                        , &amp;
     &amp;                                                       powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_  , &amp;
     &amp;                                                       powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_           &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"                 >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                           &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                        , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                  , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_            , &amp;
     &amp;                                                       darkMatterParticle_                =darkMatterParticle_                  , &amp;
     &amp;                                                       tableStore                         =.true.                                 &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloEnvironment_"                     >
   <constructor>
    haloEnvironmentNormal                                       (                                                                           &amp;
     &amp;                                                       time                               =time                                 , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                  , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_            , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                        , &amp;
     &amp;                                                       criticalOverdensity_               =criticalOverdensity_                 , &amp;
     &amp;                                                       massEnvironment                    =1.0d13                                 &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVarianceConditioned_" >
   <constructor>
    cosmologicalMassVariancePeakBackgroundSplit                 (                                                                           &amp;
     &amp;                                                       haloEnvironment_                   =haloEnvironment_                     , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_            , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                    &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensityConditioned_"      >
   <constructor>
    criticalOverdensityPeakBackgroundSplit                      (                                                                           &amp;
     &amp;                                                       criticalOverdensity_               =criticalOverdensity_                 , &amp;
     &amp;                                                       haloEnvironment_                   =haloEnvironment_                     , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                  , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_            , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                          &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="excursionSetBarrier_"                 >
   <constructor>
     excursionSetBarrierCriticalOverdensity                     (                                                                           &amp;
      &amp;                                                      criticalOverdensity_               =criticalOverdensity_                 , &amp;
      &amp;                                                      cosmologicalMassVariance_          =cosmologicalMassVariance_              &amp;
      &amp;                                                     )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="excursionSetBarrierConditioned_"      >
   <constructor>
     excursionSetBarrierCriticalOverdensity                     (                                                                           &amp;
      &amp;                                                      criticalOverdensity_               =criticalOverdensityConditioned_      , &amp;
      &amp;                                                      cosmologicalMassVariance_          =cosmologicalMassVarianceConditioned_   &amp;
      &amp;                                                     )
   </constructor>
  </referenceConstruct>  
  <referenceConstruct object="excursionSetFirstCrossing_"           >
   <constructor>
     excursionSetFirstCrossingLinearBarrier                     (                                                                           &amp;
      &amp;                                                      fractionalTimeStep                 =0.01d0                               , &amp;
      &amp;                                                      excursionSetBarrier_               =excursionSetBarrier_                 , &amp;
      &amp;                                                      cosmologicalMassVariance_          =cosmologicalMassVariance_              &amp;
      &amp;                                                     )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="excursionSetFirstCrossingConditioned_">
   <constructor>
     excursionSetFirstCrossingLinearBarrier                     (                                                                           &amp;
      &amp;                                                      fractionalTimeStep                 =0.01d0                               , &amp;
      &amp;                                                      excursionSetBarrier_               =excursionSetBarrierConditioned_      , &amp;
      &amp;                                                      cosmologicalMassVariance_          =cosmologicalMassVarianceConditioned_   &amp;
      &amp;                                                     )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloMassFunction_"                    >
   <constructor>
    haloMassFunctionPressSchechter                              (                                                                           &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_            , &amp;
     &amp;                                                       excursionSetFirstCrossing_         =excursionSetFirstCrossing_             &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloMassFunctionConditioned_"         >
   <constructor>
    haloMassFunctionPressSchechter                              (                                                                           &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                 , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVarianceConditioned_ , &amp;
     &amp;                                                       excursionSetFirstCrossing_         =excursionSetFirstCrossingConditioned_  &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloMassFunctionEnvironmentAveraged_" >
   <constructor>
    haloMassFunctionEnvironmentAveraged                         (                                                                           &amp;
     &amp;                                                       includeUnoccupiedVolume            =.true.                               , &amp;
     &amp;                                                       haloMassFunctionConditioned_       =haloMassFunctionConditioned_         , &amp;
     &amp;                                                       haloMassFunctionUnconditioned_     =haloMassFunction_                    , &amp;
     &amp;                                                       haloEnvironment_                   =haloEnvironment_                     , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  !!]
  mass          =  Make_Range(massMinimum,massMaximum,massCount,rangeTypeLogarithmic)
  node          => treeNode  (hostTree  =tree   )
  tree%nodeBase => node
  basic         => node%basic(autoCreate=.true.)
  call tree %properties%initialize         (    )
  call basic           %timeSet            (time)
  call basic           %timeLastIsolatedSet(time)
  do i=1,massCount
     call basic%massSet(mass(i))
     massFunction                   (i)=haloMassFunction_                   %differential(time,mass(i),node)
     massFunctionEnvironmentAveraged(i)=haloMassFunctionEnvironmentAveraged_%differential(time,mass(i),node)
  end do
  call Assert('⟨n(M|δ)⟩ₑₙᵥ = n(M) for Press-Schechter mass function',massFunction,massFunctionEnvironmentAveraged,relTol=1.0d-4)
  ! Clean up.
  call node%destroy()
  deallocate(node)
  !![
  <objectDestructor name="cosmologyParameters_"                 />
  <objectDestructor name="cosmologyFunctions_"                  />
  <objectDestructor name="cosmologicalMassVariance_"            />
  <objectDestructor name="cosmologicalMassVarianceConditioned_" />
  <objectDestructor name="linearGrowth_"                        />
  <objectDestructor name="powerSpectrumWindowFunction_"         />
  <objectDestructor name="powerSpectrumPrimordial_"             />
  <objectDestructor name="transferFunction_"                    />
  <objectDestructor name="powerSpectrumPrimordialTransferred_"  />
  <objectDestructor name="darkMatterParticle_"                  />
  <objectDestructor name="criticalOverdensity_"                 />
  <objectDestructor name="criticalOverdensityConditioned_"      />
  <objectDestructor name="haloMassFunction_"                    />
  <objectDestructor name="haloMassFunctionConditioned_"         />
  <objectDestructor name="haloMassFunctionEnvironmentAveraged_" />
  <objectDestructor name="haloEnvironment_"                     />
  <objectDestructor name="excursionSetBarrier_"                 />
  <objectDestructor name="excursionSetBarrierConditioned_"      />
  <objectDestructor name="excursionSetFirstCrossing_"           />
  <objectDestructor name="excursionSetFirstCrossingConditioned_"/>
  !!]
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Tests_Halo_Mass_Function_Environmental_Average
