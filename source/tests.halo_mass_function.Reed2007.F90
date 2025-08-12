!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a program which tests the \cite{reed_halo_2007} mass function by comparing to Darren Reed's
\href{http://icc.dur.ac.uk/Research/PublicDownloads/genmf_readme.html}{\textsc{genmf}} code.
!!}

program Tests_Halo_Mass_Function_Reed2007
  !!{
  Tests the \cite{reed_halo_2007} mass function by comparing to Darren Reed's
  \href{http://icc.dur.ac.uk/Research/PublicDownloads/genmf_readme.html}{\textsc{genmf}} code.
  !!}
  use :: Cosmological_Density_Field          , only : criticalOverdensityFixed                , cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple               , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : Count_Lines_In_File
  use :: Halo_Mass_Functions                 , only : haloMassFunctionReed2007
  use :: IO_HDF5                             , only : ioHDF5AccessInitialize
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionFile
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: ISO_Varying_String                  , only : var_str                                 , operator(//)                         , char                , varying_string
  use :: String_Handling                     , only : operator(//)
  implicit none
  integer                                                                               :: fUnit                              , i           , &
       &                                                                                   massCount
  double precision                                                                      :: time                               , redshift
  integer                                                                               :: iRedshift
  character       (len=64                                  )                            :: label
  double precision                                          , allocatable, dimension(:) :: mass                               , massFunction, &
       &                                                                                   massFunctionReed
  logical                                                   , allocatable, dimension(:) :: success
  type            (cosmologyParametersSimple               )                            :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )                            :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower   )                            :: cosmologicalMassVariance_
  type            (linearGrowthCollisionlessMatter         )                            :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat       )                            :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw         )                            :: powerSpectrumPrimordial_
  type            (transferFunctionFile                    )                            :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple)                            :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                   )                            :: darkMatterParticle_
  type            (criticalOverdensityFixed                )                            :: criticalOverdensity_
  type            (haloMassFunctionReed2007                )                            :: haloMassFunction_
  type            (varying_string                          )                            :: fileName

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Halo mass function: Reed et al. (2007)")
  ! Construct required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                                                     &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                                                     &amp;
     &amp;                                                        OmegaMatter                        = 0.23800d0                                                     , &amp;
     &amp;                                                        OmegaBaryon                        = 0.00000d0                                                     , &amp;
     &amp;                                                        OmegaDarkEnergy                    = 0.76200d0                                                     , &amp;
     &amp;                                                        temperatureCMB                     = 2.72548d0                                                     , &amp;
     &amp;                                                        HubbleConstant                     =70.00000d0                                                       &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                               (                                                                                                     &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                             &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                              (                                                                                                     &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                           , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                              (                                                                                                     &amp;
     &amp;                                                        index_                             =+0.951d0                                                       , &amp;
     &amp;                                                        running                            =+0.000d0                                                       , &amp;
     &amp;                                                        runningRunning                     =+0.000d0                                                       , &amp;
     &amp;                                                        wavenumberReference                =+1.000d0                                                       , &amp;
     &amp;                                                        runningSmallScalesOnly             =.false.                                                          &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionFile                                        (                                                                                                      &amp;
     &amp;                                                        fileName                           ='testSuite/data/haloMassFunction/reed2007TransferFunction.hdf5', &amp;
     &amp;                                                        redshift                           =0.000d0                                                        , &amp;
     &amp;                                                        acceptNegativeValues               =.false.                                                        , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                           , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                     (                                                                                                     &amp;
     &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordial_                                       , &amp;
     &amp;                                                        transferFunction_                  =transferFunction_                                              , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                                                    &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                            (                                                                                                     &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                             &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                        (                                                                                                     &amp;
     &amp;                                                        sigma8                             =0.74d+0                                                        , &amp;
     &amp;                                                        tolerance                          =1.00d-4                                                        , &amp;
     &amp;                                                        toleranceTopHat                    =1.00d-4                                                        , &amp;
     &amp;                                                        nonMonotonicIsFatal                =.true.                                                         , &amp;
     &amp;                                                        monotonicInterpolation             =.false.                                                        , &amp;
     &amp;                                                        truncateAtParticleHorizon          =.false.                                                        , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                           , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                                            , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                                                  , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_                            , &amp;
     &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_                                     &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensityFixed(                                                                                                                                          &amp;
     &amp;                                                        criticalOverdensity_               =1.686d0                                                        , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                                                  , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                                            , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_                                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="haloMassFunction_"                  >
   <constructor>
    haloMassFunctionReed2007                                    (                                                                                                     &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                                          , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_                                     , &amp;
     &amp;                                                        criticalOverdensity_               =criticalOverdensity_                                            &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  ! Iterate over redshifts.
  do iRedshift=1,2
     select case (iRedshift)
     case (1)
        redshift= 0.0d0
     case (2)
        redshift=10.0d0
     end select
     time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
     ! Determine number of masses in reference data file and allocate arrays.
     massCount=Count_Lines_In_File(char(var_str('testSuite/data/haloMassFunction/reed2007MassFunction_z')//int(redshift)//'.txt'))-1
     allocate(mass            (massCount))
     allocate(massFunction    (massCount))
     allocate(massFunctionReed(massCount))
     allocate(success         (massCount))
     ! Compute mass function for each reference mass.
     fileName=var_str('testSuite/data/haloMassFunction/reed2007MassFunction_z')//int(redshift)//'.txt'
     open(newUnit=fUnit,file=char(fileName),status='old',form='formatted')
     read (fUnit,*)
     do i=1,massCount
        read (fUnit,*) mass(i),massFunctionReed(i)
        mass            (i)=(10.0d0**mass            (i))/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
        massFunctionReed(i)=(10.0d0**massFunctionReed(i))*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**3/log(10.0d0)/mass(i)
        massFunction    (i)=haloMassFunction_%differential(time,mass(i))
     end do
     ! Assert that our mass function agrees with the reference data.
     success=                                                     &
          &   abs(1.0d0-massFunction/massFunctionReed) < 1.0d-02  &
          &  .or.                                                 &
          &   abs(      massFunction-massFunctionReed) < 1.0d-20
     write (label,'(a,f4.1)')  'halo mass mass function at z=',redshift
     call Assert(               &
          &      trim(label  ), &
          &      all (success), &
          &      .true.         &
          &     )
     ! Clean up memory.
     deallocate(mass            )
     deallocate(massFunction    )
     deallocate(massFunctionReed)
     deallocate(success         )
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Halo_Mass_Function_Reed2007
