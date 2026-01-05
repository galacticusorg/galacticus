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
Contains a program to test stellar populations.
!!}

program Test_Stellar_Populations
  !!{
  Tests of stellar populations.
  !!}
  use :: Abundances_Structure                      , only : abundances
  use :: Display                                   , only : displayVerbositySet                    , verbosityLevelWorking
  use :: Events_Hooks                              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities                , only : Functions_Global_Set
  use :: Input_Paths                               , only : inputPath                              , pathTypeDataStatic
  use :: ISO_Varying_String                        , only : char
  use :: Input_Parameters                          , only : inputParameters
  use :: Galacticus_Nodes                          , only : nodeClassHierarchyInitialize
  use :: Node_Components                           , only : Node_Components_Initialize
  use :: Numerical_Constants_Astronomical          , only : metallicitySolar
  use :: Stellar_Astrophysics                      , only : stellarAstrophysics                    , stellarAstrophysicsFile
  use :: Stellar_Astrophysics_Tracks               , only : stellarTracksFile
  use :: Stellar_Astrophysics_Winds                , only : stellarWindsLeitherer1992
  use :: Stellar_Feedback                          , only : stellarFeedbackStandard
  use :: Stellar_Population_Spectra                , only : stellarPopulationSpectraFSPS
  use :: Stellar_Populations                       , only : stellarPopulationStandard
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001
  use :: Supernovae_Population_III                 , only : supernovaePopulationIIIHegerWoosley2002
  use :: Supernovae_Type_Ia                        , only : supernovaeTypeIaNagashima2005
  use :: Unit_Tests                                , only : Assert                                 , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                         , parameter :: ageMinimum               =0.0d0, ageMaximum=10.0d0
  type            (inputParameters                        ), target    :: parameters
  type            (abundances                             )            :: abundances_
  type            (initialMassFunctionChabrier2001        )            :: initialMassFunction_
  type            (stellarAstrophysicsFile                )            :: stellarAstrophysics_
  type            (stellarPopulationStandard              )            :: stellarPopulation_
  type            (stellarFeedbackStandard                )            :: stellarFeedback_
  type            (stellarTracksFile                      )            :: stellarTracks_
  type            (stellarWindsLeitherer1992              )            :: stellarWinds_
  type            (supernovaeTypeIaNagashima2005          )            :: supernovaeTypeIa_
  type            (supernovaePopulationIIIHegerWoosley2002)            :: supernovaePopulationIII_
  type            (stellarPopulationSpectraFSPS           )            :: stellarPopulationSpectra_
  double precision                                                     :: recycledMass                   , yieldMetals

  call displayVerbositySet(verbosityLevelWorking)
  parameters=inputParameters()
  call Functions_Global_Set        (          )
  call eventsHooksInitialize       (          )
  call nodeClassHierarchyInitialize(parameters)
  call Node_Components_Initialize  (parameters)
  call abundances_%metallicitySet  (metallicitySolar)
  initialMassFunction_     =initialMassFunctionChabrier2001        (                                                                                                                                           &
       &                                                            massLower                            =+0.10d0                                                                                            , &
       &                                                            massTransition                       =+1.00d0                                                                                            , &
       &                                                            massUpper                            =+1.25d2                                                                                            , &
       &                                                            exponent                             =-2.30d0                                                                                            , &
       &                                                            massCharacteristic                   =+0.08d0                                                                                            , &
       &                                                            sigma                                =+0.69d0                                                                                              &
       &                                                           )
  stellarAstrophysics_     =stellarAstrophysicsFile                (                                                                                                                                           &
       &                                                            fileName                             =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/stellarPropertiesCompilationStandard.xml'  &
       &                                                           )
  stellarTracks_           =stellarTracksFile                      (                                                                                                                                           &
       &                                                            fileName                             =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Stellar_Tracks_Padova.hdf5'                &
       &                                                           )
  stellarWinds_            =stellarWindsLeitherer1992              (                                                                                                                                           &
       &                                                            stellarTracks_                       =stellarTracks_                                                                                       &
       &                                                           )
  supernovaeTypeIa_        =supernovaeTypeIaNagashima2005          (                                                                                                                                           &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_                                                                                 &
       &                                                           )
  supernovaePopulationIII_ =supernovaePopulationIIIHegerWoosley2002(                                                                                                                                           &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_                                                                                 &
       &                                                           )
  stellarFeedback_         =stellarFeedbackStandard                (                                                                                                                                           &
       &                                                            initialMassForSupernovaeTypeII       =8.0d00                                                                                             , &
       &                                                            supernovaEnergy                      =1.0d51                                                                                             , &
       &                                                            supernovaeTypeIa_                    =supernovaeTypeIa_                                                                                  , &
       &                                                            supernovaePopulationIII_             =supernovaePopulationIII_                                                                           , &
       &                                                            stellarWinds_                        =stellarWinds_                                                                                      , &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_                                                                                 &
       &                                                           )
  stellarPopulationSpectra_=stellarPopulationSpectraFSPS           (                                                                                                                                           &
       &                                                            forceZeroMetallicity                 =.false.                                                                                            , &
       &                                                            initialMassFunction_                 =initialMassFunction_                                                                                 &
       &                                                           )
  stellarPopulation_       =stellarPopulationStandard              (                                                                                                                                           &
       &                                                            instantaneousRecyclingApproximation  =.false.                                                                                            , &
       &                                                            instantaneousYieldApproximation      =.false.                                                                                            , &
       &                                                            instantaneousEnergyInputApproximation=.false.                                                                                            , &
       &                                                            massLongLived                        =1.0d0                                                                                              , &
       &                                                            ageEffective                         =1.0d1                                                                                              , &
       &                                                            initialMassFunction_                 =initialMassFunction_                                                                               , &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_                                                                               , &
       &                                                            stellarFeedback_                     =stellarFeedback_                                                                                   , &
       &                                                            supernovaeTypeIa_                    =supernovaeTypeIa_                                                                                  , &
       &                                                            stellarPopulationSpectra_            =stellarPopulationSpectra_                                                                            &
       &                                                           )
  call Unit_Tests_Begin_Group("Stellar population functions")
  ! Compute the recycled mass, yield etc. by getting the mean production rate in our time interval and multiplying by the size of
  ! the interval. This is compared to a previously computed value for the Chabrier (2001) IMF.
  recycledMass=stellarPopulation_%recycledFractionInstantaneous()
  call Assert('recycled fraction',recycledMass,4.6d-1,relTol=1.0d-2)
  yieldMetals =stellarPopulation_%yieldInstantaneous    ()
  call Assert('metal yield'      ,yieldMetals ,3.7d-2,relTol=1.0d-2)
  call Unit_Tests_End_Group  ()
  call Unit_Tests_Finish     ()
end program Test_Stellar_Populations


