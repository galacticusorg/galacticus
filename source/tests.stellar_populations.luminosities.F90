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
Contains a program to test stellar population luminosities.
!!}

program Test_Stellar_Populations_Luminosities
  !!{
  Tests of stellar population luminosities.
  !!}
  use :: Abundances_Structure                      , only : abs                                           , abundances                               , max
  use :: Cosmology_Functions                       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                      , only : cosmologyParametersSimple
  use :: Display                                   , only : displayVerbositySet                           , verbosityLevelWorking
  use :: Events_Hooks                              , only : eventsHooksInitialize
  use :: File_Utilities                            , only : File_Exists                                   , Directory_Make
  use :: Input_Paths                               , only : inputPath                                     , pathTypeDataDynamic                      , pathTypeDataStatic
  use :: ISO_Varying_String                        , only : char                                          , operator(//)                             , var_str
  use :: Input_Parameters                          , only : inputParameters
  use :: Instruments_Filters                       , only : Filter_Get_Index
  use :: Stellar_Astrophysics                      , only : stellarAstrophysics                           , stellarAstrophysicsFile
  use :: Stellar_Astrophysics_Tracks               , only : stellarTracksFile
  use :: Stellar_Astrophysics_Winds                , only : stellarWindsLeitherer1992
  use :: Stellar_Feedback                          , only : stellarFeedbackStandard
  use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesStandard
  use :: Stellar_Population_Spectra                , only : stellarPopulationSpectraFile
  use :: Stellar_Population_Spectra_Postprocess    , only : stellarPopulationSpectraPostprocessorIdentity , stellarPopulationSpectraPostprocessorList
  use :: Stellar_Populations                       , only : stellarPopulationStandard
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001
  use :: Supernovae_Population_III                 , only : supernovaePopulationIIIHegerWoosley2002
  use :: Supernovae_Type_Ia                        , only : supernovaeTypeIaNagashima2005
  use :: System_Download                           , only : download
  use :: Unit_Tests                                , only : Assert                                        , Unit_Tests_Begin_Group                   , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (inputParameters                               ), target        :: parameters
  type            (stellarPopulationSpectraPostprocessorList     ), dimension( 2) :: stellarPopulationSpectraPostprocessorList_
  type            (stellarPopulationSpectraPostprocessorIdentity ), target        :: stellarPopulationSpectraPostprocessor_
  double precision                                                , dimension( 2) :: redshift                                       , age                                        , &
       &                                                                             luminosity                                     , magnitude
  double precision                                                , dimension(15) :: columns
  integer                                                         , dimension( 2) :: luminosityIndex                                , filterIndex
  type            (abundances                                    )                :: abundances_
  type            (initialMassFunctionChabrier2001               )                :: initialMassFunction_
  type            (stellarAstrophysicsFile                       )                :: stellarAstrophysics_
  type            (stellarPopulationStandard                     )                :: stellarPopulation_
  type            (stellarFeedbackStandard                       )                :: stellarFeedback_
  type            (stellarTracksFile                             )                :: stellarTracks_
  type            (stellarWindsLeitherer1992                     )                :: stellarWinds_
  type            (supernovaeTypeIaNagashima2005                 )                :: supernovaeTypeIa_
  type            (supernovaePopulationIIIHegerWoosley2002       )                :: supernovaePopulationIII_
  type            (stellarPopulationSpectraFile                  )                :: stellarPopulationSpectra_
  type            (cosmologyParametersSimple                     )                :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                )                :: cosmologyFunctions_
  type            (stellarPopulationBroadBandLuminositiesStandard)                :: stellarPopulationBroadBandLuminosities_
  character       (len=1024                                      )                :: line
  integer                                                                         :: bc2003                                         , i                                          , &
       &                                                                             status
  double precision                                                                :: metallicity                                    , OmegaMatter                                , &
       &                                                                             OmegaDarkEnergy                                , HubbleConstant                             , &
       &                                                                             distanceModulus                                , differenceMaximumMagnitudeAbsoluteRestFrame, &
       &                                                                             differenceMaximumMagnitudeAbsoluteObservedFrame, differenceMaximumMagnitudeApparent

  parameters=inputParameters()
  call eventsHooksInitialize()
  call displayVerbositySet  (verbosityLevelWorking)
  ! Ensure that we have the required stellar population spectra file.
  call Directory_Make(char(inputPath(pathTypeDataDynamic))//'stellarPopulations/SSP_Spectra_BC2003_lowResolution_imfSalpeter.hdf5')
  if (.not.File_Exists(char(inputPath(pathTypeDataDynamic))//'stellarPopulations/SSP_Spectra_BC2003_lowResolution_imfSalpeter.hdf5'))  &
       & call download(                                                                                                                &
       &               "https://drive.google.com/uc?export=download&id=1DI52tMO4PEN-eGk79-0w2BHu9yaEcFMp"                            , &
       &               char(inputPath(pathTypeDataDynamic))//'stellarPopulations/SSP_Spectra_BC2003_lowResolution_imfSalpeter.hdf5'    &
       &              )
  ! Begin parsing the Bruzual & Charlot model file.
  open(newunit=bc2003,file=char(inputPath(pathTypeDataStatic))//"stellarPopulations/bc2003_lr_m72_salp_ssp.magnitude_F121.txt",status='old',form='formatted',iostat=status)
  do i=1,6
     read (bc2003,'(a)',iostat=status) line
  end do
  i=index(line,"Z=")
  read (line(i+2:),*)  metallicity
  do i=1,22
     read (bc2003,'(a)',iostat=status) line
  end do
  i=index(line,"Ho =  ")
  read (line(i+6:),*)  HubbleConstant
  i=index(line,"Omega = ")
  read (line(i+8:),*)  OmegaMatter
  i=index(line,"Omega_lambda = ")
  read (line(i+15:),*)  OmegaDarkEnergy
  do while (line(1:6) /= "#0.000")
     read (bc2003,'(a)',iostat=status) line
  end do
  ! Construct cosmology and stellar populations.
  cosmologyParameters_                   =cosmologyParametersSimple                     (                                                                                                                                                   &
       &                                                                                 OmegaMatter                          =OmegaMatter                                                                                                , &
       &                                                                                 OmegaBaryon                          =0.0d0                                                                                                      , &
       &                                                                                 OmegaDarkEnergy                      =OmegaDarkEnergy                                                                                            , &
       &                                                                                 temperatureCMB                       =2.7d0                                                                                                      , &
       &                                                                                 HubbleConstant                       =HubbleConstant                                                                                               &
       &                                                                                )
  cosmologyFunctions_                    =cosmologyFunctionsMatterLambda                (                                                                                                                                                   &
       &                                                                                 cosmologyParameters_                 =cosmologyParameters_                                                                                         &
       &                                                                                 )
  initialMassFunction_                   =initialMassFunctionChabrier2001               (                                                                                                                                                   &
       &                                                                                 massLower                            =+0.10d0                                                                                                    , &
       &                                                                                 massTransition                       =+1.00d0                                                                                                    , &
       &                                                                                 massUpper                            =+1.25d2                                                                                                    , &
       &                                                                                 exponent                             =-2.30d0                                                                                                    , &
       &                                                                                 massCharacteristic                   =+0.08d0                                                                                                    , &
       &                                                                                 sigma                                =+0.69d0                                                                                                      &
       &                                                                                )
  stellarAstrophysics_                   =stellarAstrophysicsFile                       (                                                                                                                                                   &
       &                                                                                 fileName                             =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/stellarPropertiesCompilationStandard.xml'          &
       &                                                                                )
  stellarTracks_                         =stellarTracksFile                             (                                                                                                                                                   &
       &                                                                                 fileName                             =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Stellar_Tracks_Padova.hdf5'                        &
       &                                                                                )
  stellarWinds_                          =stellarWindsLeitherer1992                     (                                                                                                                                                   &
       &                                                                                 stellarTracks_                       =stellarTracks_                                                                                               &
       &                                                                                )
  supernovaeTypeIa_                      =supernovaeTypeIaNagashima2005                 (                                                                                                                                                   &
       &                                                                                 stellarAstrophysics_                 =stellarAstrophysics_                                                                                         &
       &                                                                                )
  supernovaePopulationIII_               =supernovaePopulationIIIHegerWoosley2002       (                                                                                                                                                   &
       &                                                                                 stellarAstrophysics_                 =stellarAstrophysics_                                                                                         &
       &                                                                                )
  stellarFeedback_                       =stellarFeedbackStandard                       (                                                                                                                                                   &
       &                                                                                 initialMassForSupernovaeTypeII       =8.0d00                                                                                                     , &
       &                                                                                 supernovaEnergy                      =1.0d51                                                                                                     , &
       &                                                                                 supernovaeTypeIa_                    =supernovaeTypeIa_                                                                                          , &
       &                                                                                 supernovaePopulationIII_             =supernovaePopulationIII_                                                                                   , &
       &                                                                                 stellarWinds_                        =stellarWinds_                                                                                              , &
       &                                                                                 stellarAstrophysics_                 =stellarAstrophysics_                                                                                         &
       &                                                                                )
  stellarPopulationSpectra_              =stellarPopulationSpectraFile                  (                                                                                                                                                   &
       &                                                                                 forceZeroMetallicity                 =.false.                                                                                                    , &
       &                                                                                 fileName                             =char(inputPath(pathTypeDataStatic))//'stellarPopulations/SSP_Spectra_BC2003_lowResolution_imfSalpeter.hdf5'  &
       &                                                                                )
  stellarPopulation_                     =stellarPopulationStandard                     (                                                                                                                                                   &
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
  stellarPopulationSpectraPostprocessor_ =stellarPopulationSpectraPostprocessorIdentity (                                                                                                                                                   &
       &                                                                                )
  stellarPopulationBroadBandLuminosities_=stellarPopulationBroadBandLuminositiesStandard(                                                                                                                                                   &
       &                                                                                 integrationToleranceRelative         =4.0d-3                                                                                                     , &
       &                                                                                 integrationToleranceDegrade          =.false.                                                                                                    , &
       &                                                                                 maximumAgeExceededIsFatal            =.true.                                                                                                     , &
       &                                                                                 storeToFile                          =.true.                                                                                                     , &
       &                                                                                 storeDirectory                       =inputPath(pathTypeDataDynamic)//'stellarPopulations'                                                         &
       &                                                                                )
  ! Create a list of stellar population spectra postprocessors.
  stellarPopulationSpectraPostprocessorList_(1)%stellarPopulationSpectraPostprocessor_ => stellarPopulationSpectraPostprocessor_
  stellarPopulationSpectraPostprocessorList_(2)%stellarPopulationSpectraPostprocessor_ => stellarPopulationSpectraPostprocessor_
  ! Initialize filters and metallicities.
  call abundances_%metallicitySet(metallicity)
  luminosityIndex=[                          1                   ,                          1                   ]
  filterIndex    =[Filter_Get_Index(var_str('SDSS_g_1.3airmass')),Filter_Get_Index(var_str('SDSS_g_1.3airmass'))]
  ! Initialize variables to track maximum differences from Bruzual & Charlot results.
  differenceMaximumMagnitudeAbsoluteRestFrame    =0.0d0
  differenceMaximumMagnitudeAbsoluteObservedFrame=0.0d0
  differenceMaximumMagnitudeApparent             =0.0d0
  ! Compute magnitudes at each age tabulated in the Bruzual & Charlot model file.
  do while (status == 0)
     read (bc2003,'(a)',iostat=status) line
     if (status     /= 0    ) exit ! Finish if we reach the end of the file.
     read (line,*) columns
     if (columns(3) == 0.0d0) exit ! Finish if we reach zero age.
     ! Construct band properties and evaluate luminosities. Increment the observed frame luminosity index by 1 to ensure it is treated as a new band.
     luminosityIndex(2)=luminosityIndex(2)+1
     redshift          =[0.0d0,columns(1)]
     age               =       columns(3)
     luminosity        =stellarPopulationBroadBandLuminosities_%luminosities(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessorList_,stellarPopulation_,abundances_,age,redshift)
     magnitude         =-2.5d0*log10(      luminosity) &
          &             -2.5d0*log10(1.0d0+redshift  )   ! Bruzual & Charlot observed frame absolute magnitudes include the factor for compression of photon frequencies.
     distanceModulus   =+25.0d0                                                              &
          &             + 5.0d0                                                              &
          &             *log10(                                                              &
          &                    cosmologyFunctions_%distanceLuminosity         (              &
          &                    cosmologyFunctions_%cosmicTime                  (             &
          &                    cosmologyFunctions_%expansionFactorFromRedshift  (            &
          &                                                                      redshift(2) &
          &                                                                     )            &
          &                                                                    )             &
          &                                                                   )              &
          &                   )
     ! Compute maximum difference from Bruzual & Charlot.
     differenceMaximumMagnitudeAbsoluteRestFrame    =max(differenceMaximumMagnitudeAbsoluteRestFrame    ,abs(magnitude(1)                -columns( 9)))
     differenceMaximumMagnitudeAbsoluteObservedFrame=max(differenceMaximumMagnitudeAbsoluteObservedFrame,abs(magnitude(2)                -columns(11)))
     differenceMaximumMagnitudeApparent             =max(differenceMaximumMagnitudeApparent             ,abs(magnitude(2)+distanceModulus-columns(12)))
  end do
  close(bc2003)
  call Unit_Tests_Begin_Group("Stellar population luminosities"              )
  call Unit_Tests_Begin_Group("Comparison with Bruzual & Charlot (2003) code")
  call Assert('rest-frame absolute magnitude'    ,differenceMaximumMagnitudeAbsoluteRestFrame    ,0.0d0,absTol=0.03d0)
  call Assert('observed-frame absolute magnitude',differenceMaximumMagnitudeAbsoluteObservedFrame,0.0d0,absTol=0.03d0)
  call Assert('apparent magnitude'               ,differenceMaximumMagnitudeApparent             ,0.0d0,absTol=0.05d0)
  call Unit_Tests_End_Group  ()
  call Unit_Tests_End_Group  ()
  call Unit_Tests_Finish     ()
  call parameters%destroy()

end program Test_Stellar_Populations_Luminosities


