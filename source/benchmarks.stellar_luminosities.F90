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
Contains a program to benchmark stellar population luminosity calculations.
!!}

program Benchmark_Stellar_Populations_Luminosities
  !!{
  Benchmarking of stellar population luminosity calculations.
  !!}
  use, intrinsic :: ISO_Fortran_Env                           , only : output_unit
  use            :: Abundances_Structure                      , only : abundances                                    , metallicityTypeLinearByMassSolar
  use            :: Cosmology_Functions                       , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters                      , only : cosmologyParametersSimple
  use            :: Display                                   , only : displayVerbositySet                           , verbosityLevelWorking
  use            :: Events_Hooks                              , only : eventsHooksInitialize
  use            :: Input_Paths                               , only : inputPath                                     , pathTypeDataDynamic                      , pathTypeDataStatic
  use            :: ISO_Varying_String                        , only : char                                          , operator(//)                             , var_str
  use            :: Input_Parameters                          , only : inputParameters
  use            :: Instruments_Filters                       , only : Filter_Get_Index
  use            :: Kind_Numbers                              , only : kind_int8
  use            :: Stellar_Astrophysics                      , only : stellarAstrophysics                           , stellarAstrophysicsFile
  use            :: Stellar_Astrophysics_Tracks               , only : stellarTracksFile
  use            :: Stellar_Astrophysics_Winds                , only : stellarWindsLeitherer1992
  use            :: Stellar_Feedback                          , only : stellarFeedbackStandard
  use            :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesStandard
  use            :: Stellar_Population_Spectra                , only : stellarPopulationSpectraFile
  use            :: Stellar_Population_Spectra_Postprocess    , only : stellarPopulationSpectraPostprocessorIdentity , stellarPopulationSpectraPostprocessorList
  use            :: Stellar_Populations                       , only : stellarPopulationStandard
  use            :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001
  use            :: Supernovae_Population_III                 , only : supernovaePopulationIIIHegerWoosley2002
  use            :: Supernovae_Type_Ia                        , only : supernovaeTypeIaNagashima2005
  implicit none
  type            (inputParameters                               ), target                                 :: parameters
  integer                                                         , parameter                              :: filterCount                               =  137  , populationCount      =20, &
       &                                                                                                      trialCount                                =10000
  double precision                                                , parameter                              :: redshiftStep                              =0.1d0
  type            (stellarPopulationSpectraPostprocessorList     ), dimension(filterCount*populationCount) :: stellarPopulationSpectraPostprocessorList_
  type            (stellarPopulationSpectraPostprocessorIdentity ), target                                 :: stellarPopulationSpectraPostprocessor_
  double precision                                                , dimension(filterCount*populationCount) :: redshift                                          , age                     , &
       &                                                                                                      luminosity
  integer                                                         , dimension(filterCount*populationCount) :: luminosityIndex                                   , filterIndex
  type            (abundances                                    )                                         :: abundances_
  type            (initialMassFunctionChabrier2001               )                                         :: initialMassFunction_
  type            (stellarAstrophysicsFile                       )                                         :: stellarAstrophysics_
  type            (stellarPopulationStandard                     )                                         :: stellarPopulation_
  type            (stellarFeedbackStandard                       )                                         :: stellarFeedback_
  type            (stellarTracksFile                             )                                         :: stellarTracks_
  type            (stellarWindsLeitherer1992                     )                                         :: stellarWinds_
  type            (supernovaeTypeIaNagashima2005                 )                                         :: supernovaeTypeIa_
  type            (supernovaePopulationIIIHegerWoosley2002       )                                         :: supernovaePopulationIII_
  type            (stellarPopulationSpectraFile                  )                                         :: stellarPopulationSpectra_
  type            (cosmologyParametersSimple                     )                                         :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                )                                         :: cosmologyFunctions_
  type            (stellarPopulationBroadBandLuminositiesStandard)                                         :: stellarPopulationBroadBandLuminosities_
  integer                                                                                                  :: trial                                             , i
  integer         (kind=kind_int8                                )                                         :: countStart                                        , countEnd                , &
       &                                                                                                      countRate
  integer         (kind=kind_int8                                ), dimension(trialCount                 ) :: trialTime
  character       (len =   3                                     )                                         :: units
  double precision                                                                                         :: timeMean                                          , timeStandardDeviation   , &
       &                                                                                                      timeMeanError

  parameters=inputParameters()
  call displayVerbositySet(verbosityLevelWorking)
  call eventsHooksInitialize()
  ! Construct cosmology and stellar populations.
  cosmologyParameters_                  =cosmologyParametersSimple                      (                                                                                                                                                   &
       &                                                                                 OmegaMatter                          = 0.3d0                                                                                                     , &
       &                                                                                 OmegaBaryon                          = 0.0d0                                                                                                     , &
       &                                                                                 OmegaDarkEnergy                      = 0.7d0                                                                                                     , &
       &                                                                                 temperatureCMB                       = 2.7d0                                                                                                     , &
       &                                                                                 HubbleConstant                       =70.0d0                                                                                                       &
       &                                                                                )
  cosmologyFunctions_                   =cosmologyFunctionsMatterLambda                 (                                                                                                                                                   &
       &                                                                                 cosmologyParameters_                 =cosmologyParameters_                                                                                         &
       &                                                                                 )
  initialMassFunction_                  =initialMassFunctionChabrier2001                (                                                                                                                                                   &
       &                                                                                 massLower                            =+0.10d0                                                                                                    , &
       &                                                                                 massTransition                       =+1.00d0                                                                                                    , &
       &                                                                                 massUpper                            =+1.25d2                                                                                                    , &
       &                                                                                 exponent                             =-2.30d0                                                                                                    , &
       &                                                                                 massCharacteristic                   =+0.08d0                                                                                                    , &
       &                                                                                 sigma                                =+0.69d0                                                                                                      &
       &                                                                                )
  stellarAstrophysics_                  =stellarAstrophysicsFile                        (                                                                                                                                                   &
       &                                                                                 fileName                            =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/stellarPropertiesCompilationStandard.xml'           &
       &                                                                                )
  stellarTracks_                        =stellarTracksFile                              (                                                                                                                                                   &
       &                                                                                 fileName                            =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Stellar_Tracks_Padova.hdf5'                         &
       &                                                                                )
  stellarWinds_                         =stellarWindsLeitherer1992                      (                                                                                                                                                   &
       &                                                                                 stellarTracks_                      =stellarTracks_                                                                                                &
       &                                                                                )
  supernovaeTypeIa_                     =supernovaeTypeIaNagashima2005                  (                                                                                                                                                   &
       &                                                                                 stellarAstrophysics_                =stellarAstrophysics_                                                                                          &
       &                                                                                )
  supernovaePopulationIII_              =supernovaePopulationIIIHegerWoosley2002        (                                                                                                                                                   &
       &                                                                                 stellarAstrophysics_                =stellarAstrophysics_                                                                                          &
       &                                                                                )
  stellarFeedback_                      =stellarFeedbackStandard                        (                                                                                                                                                   &
       &                                                                                 initialMassForSupernovaeTypeII       =8.0d00                                                                                                     , &
       &                                                                                 supernovaEnergy                      =1.0d51                                                                                                     , &
       &                                                                                 supernovaeTypeIa_                    =supernovaeTypeIa_                                                                                          , &
       &                                                                                 supernovaePopulationIII_             =supernovaePopulationIII_                                                                                   , &
       &                                                                                 stellarWinds_                        =stellarWinds_                                                                                              , &
       &                                                                                 stellarAstrophysics_                 =stellarAstrophysics_                                                                                         &
       &                                                                                )
  stellarPopulationSpectra_             =stellarPopulationSpectraFile                   (                                                                                                                                                   &
       &                                                                                 forceZeroMetallicity                 =.false.                                                                                                    , &
       &                                                                                 fileName                             =char(inputPath(pathTypeDataStatic))//'stellarPopulations/SSP_Spectra_BC2003_lowResolution_imfSalpeter.hdf5'  &
       &                                                                                )
  stellarPopulation_                    =stellarPopulationStandard                      (                                                                                                                                                   &
       &                                                                                 instantaneousRecyclingApproximation  =.false.                                                                                                    , &
       &                                                                                 instantaneousYieldApproximation      =.false.                                                                                                    , &
       &                                                                                 instantaneousEnergyInputApproximation=.false.                                                                                                    , &
       &                                                                                 massLongLived                        =1.0d0                                                                                                      , &
       &                                                                                 ageEffective                         =1.0d1                                                                                                      , &
       &                                                                                 initialMassFunction_                 =initialMassFunction_                                                                                       , &
       &                                                                                 stellarAstrophysics_                 =stellarAstrophysics_                                                                                       , &
       &                                                                                 stellarFeedback_                     =stellarFeedback_                                                                                           , &
       &                                                                                 supernovaeTypeIa_                    =supernovaeTypeIa_                                                                                          , &
       &                                                                                 stellarPopulationSpectra_            =stellarPopulationSpectra_                                                                                    &
       &                                                                                )
  stellarPopulationSpectraPostprocessor_=stellarPopulationSpectraPostprocessorIdentity  (                                                                                                                                                   &
       &                                                                                )
  stellarPopulationBroadBandLuminosities_=stellarPopulationBroadBandLuminositiesStandard(                                                                                                                                                   &
       &                                                                                 integrationToleranceRelative         =4.0d-3                                                                                                     , &
       &                                                                                 integrationToleranceDegrade          =.false.                                                                                                    , &
       &                                                                                 maximumAgeExceededIsFatal            =.true.                                                                                                     , &
       &                                                                                 storeToFile                          =.true.                                                                                                     , &
       &                                                                                 storeDirectory                       =inputPath(pathTypeDataDynamic)//'stellarPopulations'                                                         &
       &                                                                                )
  ! Initialize filters and metallicities.
  call abundances_%metallicitySet(1.0d0,metallicityTypeLinearByMassSolar)
  do i=1,filterCount*populationCount
     luminosityIndex                           (i)                                        =  i
     stellarPopulationSpectraPostprocessorList_(i)%stellarPopulationSpectraPostprocessor_ => stellarPopulationSpectraPostprocessor_
  end do
  do i=1,populationCount
     filterIndex((i-1)*filterCount+1:i*filterCount)=                                                                                                                                                                           &
          &      [                                                                                                                                                                                                             &
          &       Filter_Get_Index(var_str("2MASS_H"              )),Filter_Get_Index(var_str("2MASS_J"              )),Filter_Get_Index(var_str("2MASS_Ks"             )),Filter_Get_Index(var_str("ACSWFC_f435w"         )), &
          &       Filter_Get_Index(var_str("ACSWFC_f606w"         )),Filter_Get_Index(var_str("ACSWFC_f814w"         )),Filter_Get_Index(var_str("ACSWFC_f850lp"        )),Filter_Get_Index(var_str("Buser_B"              )), &
          &       Filter_Get_Index(var_str("Buser_U"              )),Filter_Get_Index(var_str("Buser_V"              )),Filter_Get_Index(var_str("DES_Y"                )),Filter_Get_Index(var_str("DES_g"                )), &
          &       Filter_Get_Index(var_str("DES_i"                )),Filter_Get_Index(var_str("DES_r"                )),Filter_Get_Index(var_str("DES_u"                )),Filter_Get_Index(var_str("DES_z"                )), &
          &       Filter_Get_Index(var_str("Euclid_H"             )),Filter_Get_Index(var_str("Euclid_J"             )),Filter_Get_Index(var_str("Euclid_VIS"           )),Filter_Get_Index(var_str("Euclid_Y"             )), &
          &       Filter_Get_Index(var_str("F814W"                )),Filter_Get_Index(var_str("FourStar_FS118"       )),Filter_Get_Index(var_str("FourStar_FS209"       )),Filter_Get_Index(var_str("FourStar_H"           )), &
          &       Filter_Get_Index(var_str("FourStar_Hlong"       )),Filter_Get_Index(var_str("FourStar_Hshort"      )),Filter_Get_Index(var_str("FourStar_J"           )),Filter_Get_Index(var_str("FourStar_J1"          )), &
          &       Filter_Get_Index(var_str("FourStar_J2"          )),Filter_Get_Index(var_str("FourStar_J3"          )),Filter_Get_Index(var_str("FourStar_Ks"          )),Filter_Get_Index(var_str("Galex_FUV"            )), &
          &       Filter_Get_Index(var_str("Galex_NUV"            )),Filter_Get_Index(var_str("HeliumContinuum"      )),Filter_Get_Index(var_str("Herschel_PACS_100um"  )),Filter_Get_Index(var_str("Herschel_PACS_160um"  )), &
          &       Filter_Get_Index(var_str("Herschel_PACS_70um"   )),Filter_Get_Index(var_str("Herschel_SPIRE_PSL"   )),Filter_Get_Index(var_str("Herschel_SPIRE_PSM"   )),Filter_Get_Index(var_str("Herschel_SPIRE_PSW"   )), &
          &       Filter_Get_Index(var_str("IR_K"                 )),Filter_Get_Index(var_str("JWST_MIRI_f1130w"     )),Filter_Get_Index(var_str("JWST_MIRI_f1280w"     )),Filter_Get_Index(var_str("JWST_MIRI_f1500w"     )), &
          &       Filter_Get_Index(var_str("JWST_MIRI_f1800w"     )),Filter_Get_Index(var_str("JWST_MIRI_f2100w"     )),Filter_Get_Index(var_str("JWST_MIRI_f2550w"     )),Filter_Get_Index(var_str("JWST_MIRI_f560w"      )), &
          &       Filter_Get_Index(var_str("JWST_MIRI_f770w"      )),Filter_Get_Index(var_str("JWST_NIRCAM_f070w"    )),Filter_Get_Index(var_str("JWST_NIRCAM_f090w"    )),Filter_Get_Index(var_str("JWST_NIRCAM_f115w"    )), &
          &       Filter_Get_Index(var_str("JWST_NIRCAM_f150w"    )),Filter_Get_Index(var_str("JWST_NIRCAM_f200w"    )),Filter_Get_Index(var_str("LSST_g"               )),Filter_Get_Index(var_str("LSST_i"               )), &
          &       Filter_Get_Index(var_str("LSST_r"               )),Filter_Get_Index(var_str("LSST_u"               )),Filter_Get_Index(var_str("LSST_y"               )),Filter_Get_Index(var_str("LSST_y4"              )), &
          &       Filter_Get_Index(var_str("LSST_z"               )),Filter_Get_Index(var_str("Lyc"                  )),Filter_Get_Index(var_str("OxygenContinuum"      )),Filter_Get_Index(var_str("RGO_B"                )), &
          &       Filter_Get_Index(var_str("RGO_I"                )),Filter_Get_Index(var_str("RGO_V"                )),Filter_Get_Index(var_str("SDSS_g"               )),Filter_Get_Index(var_str("SDSS_g_1.3airmass"    )), &
          &       Filter_Get_Index(var_str("SDSS_i"               )),Filter_Get_Index(var_str("SDSS_i_1.3airmass"    )),Filter_Get_Index(var_str("SDSS_r"               )),Filter_Get_Index(var_str("SDSS_r_1.3airmass"    )), &
          &       Filter_Get_Index(var_str("SDSS_u"               )),Filter_Get_Index(var_str("SDSS_z"               )),Filter_Get_Index(var_str("SDSS_z_1.3airmass"    )),Filter_Get_Index(var_str("Spitzer_IRAC_Channel1")), &
          &       Filter_Get_Index(var_str("Spitzer_IRAC_Channel2")),Filter_Get_Index(var_str("Spitzer_IRAC_Channel3")),Filter_Get_Index(var_str("Spitzer_IRAC_Channel4")),Filter_Get_Index(var_str("Spitzer_IRS_16um"     )), &
          &       Filter_Get_Index(var_str("Spitzer_IRS_22um"     )),Filter_Get_Index(var_str("Spitzer_MIPS_160um"   )),Filter_Get_Index(var_str("Spitzer_MIPS_24um"    )),Filter_Get_Index(var_str("Spitzer_MIPS_70um"    )), &
          &       Filter_Get_Index(var_str("SuprimeCam_B"         )),Filter_Get_Index(var_str("SuprimeCam_I"         )),Filter_Get_Index(var_str("SuprimeCam_IA427"     )),Filter_Get_Index(var_str("SuprimeCam_IA464"     )), &
          &       Filter_Get_Index(var_str("SuprimeCam_IA484"     )),Filter_Get_Index(var_str("SuprimeCam_IA505"     )),Filter_Get_Index(var_str("SuprimeCam_IA527"     )),Filter_Get_Index(var_str("SuprimeCam_IA574"     )), &
          &       Filter_Get_Index(var_str("SuprimeCam_IA624"     )),Filter_Get_Index(var_str("SuprimeCam_IA679"     )),Filter_Get_Index(var_str("SuprimeCam_IA709"     )),Filter_Get_Index(var_str("SuprimeCam_IA738"     )), &
          &       Filter_Get_Index(var_str("SuprimeCam_IA767"     )),Filter_Get_Index(var_str("SuprimeCam_IA827"     )),Filter_Get_Index(var_str("SuprimeCam_R"         )),Filter_Get_Index(var_str("SuprimeCam_V"         )), &
          &       Filter_Get_Index(var_str("SuprimeCam_Y"         )),Filter_Get_Index(var_str("SuprimeCam_gPrime"    )),Filter_Get_Index(var_str("SuprimeCam_iPrime"    )),Filter_Get_Index(var_str("SuprimeCam_rPrime"    )), &
          &       Filter_Get_Index(var_str("SuprimeCam_zPrime"    )),Filter_Get_Index(var_str("UKIRT_H"              )),Filter_Get_Index(var_str("UKIRT_J"              )),Filter_Get_Index(var_str("UKIRT_K"              )), &
          &       Filter_Get_Index(var_str("VIRCAM_H"             )),Filter_Get_Index(var_str("VIRCAM_J"             )),Filter_Get_Index(var_str("VIRCAM_Ks"            )),Filter_Get_Index(var_str("VIRCAM_NB118"         )), &
          &       Filter_Get_Index(var_str("VIRCAM_NB980"         )),Filter_Get_Index(var_str("VIRCAM_NB990"         )),Filter_Get_Index(var_str("VIRCAM_Y"             )),Filter_Get_Index(var_str("VIRCAM_Z"             )), &
          &       Filter_Get_Index(var_str("WFC3IR_f105w"         )),Filter_Get_Index(var_str("WFC3IR_f125w"         )),Filter_Get_Index(var_str("WFC3IR_f160w"         )),Filter_Get_Index(var_str("WFCAM_H"              )), &
          &       Filter_Get_Index(var_str("WFCAM_J"              )),Filter_Get_Index(var_str("WFCAM_K"              )),Filter_Get_Index(var_str("WFCAM_Y"              )),Filter_Get_Index(var_str("WFCAM_Z"              )), &
          &       Filter_Get_Index(var_str("Roman_F062"           )),Filter_Get_Index(var_str("Roman_F087"           )),Filter_Get_Index(var_str("Roman_F106"           )),Filter_Get_Index(var_str("Roman_F129"           )), &
          &       Filter_Get_Index(var_str("Roman_F146"           )),Filter_Get_Index(var_str("Roman_F158"           )),Filter_Get_Index(var_str("Roman_F184"           )),Filter_Get_Index(var_str("Roman_F213"           )), &
          &       Filter_Get_Index(var_str("WIRCAM_K"             )),Filter_Get_Index(var_str("bJ"                   )),Filter_Get_Index(var_str("xRayFull"             )),Filter_Get_Index(var_str("xRayHard"             )), &
          &       Filter_Get_Index(var_str("xRaySoft"             ))                                                                                                                                                           &
          &      ]
     redshift((i-1)*filterCount+1:i*filterCount)=                                                 dble(i-1)*redshiftStep
     age     ((i-1)*filterCount+1:i*filterCount)=+cosmologyFunctions_%cosmicTime                 (                        &
          &                                       cosmologyFunctions_%expansionFactorFromRedshift (                       &
          &                                                                                        0.0d0                  &
          &                                                                                       )                       &
          &                                                                                      )                        &
          &                                      -cosmologyFunctions_%cosmicTime                 (                        &
          &                                       cosmologyFunctions_%expansionFactorFromRedshift (                       &
          &                                                                                        dble(i-1)*redshiftStep &
          &                                                                                       )                       &
          &                                                                                      )
  end do
  ! Get clock units.
  call System_Clock(countStart,countRate)
  select case (countRate)
  case (      1000)
     units="ms"
  case (   1000000)
     units="μs"
  case (1000000000)
     units="ns"
  end select
  ! Begin trials.
  do trial=1,trialCount
     call System_Clock(countStart)
     luminosity=stellarPopulationBroadBandLuminosities_%luminosities(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessorList_,stellarPopulation_,abundances_,age,redshift)
     call System_Clock(countEnd  )
     trialTime(trial)=+countEnd   &
          &           -countStart
  end do
  ! Compute mean and standard deviation of benchmark time. Skip the first two evaluations as they include start-up overheads.
  timeMean             =+       dble(sum(trialTime(3:trialCount)   ))/dble(trialCount-2)
  timeStandardDeviation=+sqrt(+(dble(sum(trialTime(3:trialCount)**2))/dble(trialCount-2)-timeMean**2) &
       &                      *                                       dble(trialCount-2)              &
       &                      /                                       dble(trialCount-3)              &
       &                     )
  timeMeanError        =+timeStandardDeviation                                                        &
       &                /sqrt(                                        dble(trialCount-2)              &
       &                     )
  ! Report benchmark information.
  write (output_unit,'(a,1x,a,1x,a,f12.1,1x,f12.1,1x,a1,a2,a1)') 'BENCHMARK','stellarPopulationInterpolation','"Stellar population interpolation"',timeMean,timeMeanError,'"',trim(adjustl(units)),'"'
  ! Clean up.
  call parameters%destroy()
end program Benchmark_Stellar_Populations_Luminosities


