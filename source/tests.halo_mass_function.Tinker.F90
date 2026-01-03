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
Contains a program which tests the \cite{tinker_towardhalo_2008} mass function by comparing to Jeremy Tinker's
\href{http://web.archive.org/web/20220313133942/https://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.
!!}

program Tests_Halo_Mass_Function_Tinker
  !!{
  Tests the \cite{tinker_towardhalo_2008} mass function by comparing to Jeremy Tinker's
  \href{http://web.archive.org/web/20220313133942/https://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.
  !!}
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple                                   , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                                         , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : Count_Lines_In_File
  use :: Halo_Mass_Functions                 , only : haloMassFunctionTinker2008
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                                      , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast             , only : virialDensityContrastFixed                                  , fixedDensityTypeMean
  implicit none
  integer                                                                                                   :: fUnit                              , i           , &
       &                                                                                                       massCount
  double precision                                                                                          :: time
  double precision                                                              , allocatable, dimension(:) :: mass                               , massFunction, &
       &                                                                                                       massFunctionTinker
  logical                                                                       , allocatable, dimension(:) :: success
  type            (cosmologyParametersSimple                                   )                            :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )                            :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower                       )                            :: cosmologicalMassVariance_
  type            (linearGrowthCollisionlessMatter                             )                            :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat                           )                            :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             )                            :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                            )                            :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    )                            :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                       )                            :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                            :: criticalOverdensity_
  type            (virialDensityContrastFixed                                  )                            :: virialDensityContrast_
  type            (haloMassFunctionTinker2008                                  )                            :: haloMassFunction_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Halo mass function: Tinker et al. (2008)")
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
     &amp;                                                        OmegaMatter                        = 0.30d0                            , &amp;
     &amp;                                                        OmegaBaryon                        = 0.04d0                            , &amp;
     &amp;                                                        OmegaDarkEnergy                    = 0.70d0                            , &amp;
     &amp;                                                        temperatureCMB                     = 2.70d0                            , &amp;
     &amp;                                                        HubbleConstant                     =71.00d0                              &amp;
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
    transferFunctionEisensteinHu1999                             (                                                                         &amp;
     &amp;                                                        neutrinoNumberEffective            =3.046d0                            , &amp;
     &amp;                                                        neutrinoMassSummed                 =0.000d0                            , &amp;
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
     &amp;                                                        sigma8                             =0.9d0                              , &amp;
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
  <referenceConstruct object="virialDensityContrast_"             >
   <constructor>
    virialDensityContrastFixed                                   (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        densityType                        =fixedDensityTypeMean               , &amp;
     &amp;                                                        densityContrastValue               =200.0d0                            , &amp;
     &amp;                                                        turnAroundOverVirialRadius         =  2.0d0                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloMassFunction_"                  >
   <constructor>
    haloMassFunctionTinker2008                                   (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        virialDensityContrast_             =virialDensityContrast_               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  time=cosmologyFunctions_%cosmicTime(1.0d0)
  ! Determine number of masses in reference data file and allocate arrays.
  massCount=Count_Lines_In_File('testSuite/data/haloMassFunction/tinker.txt')
  allocate(mass              (massCount))
  allocate(massFunction      (massCount))
  allocate(massFunctionTinker(massCount))
  allocate(success           (massCount))
  ! Ensure that critical density and critical overdensity for collapse are consistent with values used in our input file to
  ! Tinker's code.
  call Assert('critical density consistency'                 ,cosmologyParameters_%densityCritical(    )/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2,2.7751950000000000d11,relTol=1.0d-6)
  call Assert('critical overdensity for collapse consistency',criticalOverdensity_%value          (time)                                                           ,1.6755779626281502d00,relTol=1.0d-6)
  ! Compute mass function for each reference mass.
  open(newUnit=fUnit,file='testSuite/data/haloMassFunction/tinker.txt',status='old',form='formatted')
  do i=1,massCount
     read (fUnit,*) mass(i),massFunctionTinker(i)
     mass              (i)=mass              (i)/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
     massFunctionTinker(i)=massFunctionTinker(i)*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**4
     massFunction      (i)=haloMassFunction_%differential(time,mass(i))
  end do
  ! Assert that our mass function agrees with the reference data.
  success=                                                       &
       &   abs(1.0d0-massFunction/massFunctionTinker) < 1.0d-02  &
       &  .or.                                                   &
       &   abs(      massFunction-massFunctionTinker) < 1.0d-20
  call Assert(                           &
       &      'halo mass mass function', &
       &      all(success),              &
       &      .true.                     &
       &     )
  ! Clean up memory.
  deallocate(mass              )
  deallocate(massFunction      )
  deallocate(massFunctionTinker)
  deallocate(success           )
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Halo_Mass_Function_Tinker
