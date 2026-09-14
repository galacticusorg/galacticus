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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test the HII region absorption stellar population spectra postprocessor.
!!}

program Test_HII_Region_Absorption
  !!{RST
  Tests the :galacticus-class:`stellarPopulationSpectraPostprocessorHIIRegionAbsorption` postprocessor: that it applies
  the escape fraction shortward of the Lyman limit and nothing longward of it, that it does not depend on redshift, and
  that it composes with an age window in a sequence without spoiling the sharpness of that window.
  !!}
  use :: Display                               , only : displayVerbositySet                      , verbosityLevelStandard
  use :: HII_Region_Escape_Fraction            , only : hiiRegionEscapeFractionFixed
  use :: Numerical_Constants_Atomic            , only : lymanSeriesLimitWavelengthHydrogen_atomic
  use :: Stellar_Population_Spectra_Postprocess, only : postprocessorList                        , stellarPopulationSpectraPostprocessorHIIRegionAbsorption, stellarPopulationSpectraPostprocessorRecent, stellarPopulationSpectraPostprocessorSequence
  use :: Unit_Tests                            , only : Assert                                   , Unit_Tests_Begin_Group                                  , Unit_Tests_End_Group                       , Unit_Tests_Finish
  implicit none
  ! An escape fraction of 0.1 for HII regions younger than 10 Myr, and unity thereafter.
  double precision                                                          , parameter :: escapeFraction          =0.1d0 , ageLimit  =1.0d-2
  double precision                                                          , parameter :: ageYoung                =5.0d-3, ageOld    =2.0d-2
  type            (hiiRegionEscapeFractionFixed                            ), pointer   :: hiiRegionEscapeFraction_
  type            (stellarPopulationSpectraPostprocessorHIIRegionAbsorption), pointer   :: absorption_
  type            (stellarPopulationSpectraPostprocessorRecent             ), pointer   :: recent_
  type            (stellarPopulationSpectraPostprocessorSequence           ), pointer   :: sequence_
  type            (postprocessorList                                       ), pointer   :: members
  double precision                                                                      :: ageMinimum                     , ageMaximum

  call displayVerbositySet(verbosityLevelStandard)
  allocate(hiiRegionEscapeFraction_)
  allocate(absorption_             )
  allocate(recent_                 )
  allocate(sequence_               )
  !![
  <referenceConstruct object="hiiRegionEscapeFraction_" constructor="hiiRegionEscapeFractionFixed                            (escapeFraction,ageLimit )"/>
  <referenceConstruct object="absorption_"              constructor="stellarPopulationSpectraPostprocessorHIIRegionAbsorption(hiiRegionEscapeFraction_)"/>
  <referenceConstruct object="recent_"                  constructor="stellarPopulationSpectraPostprocessorRecent             (ageYoung                )"/>
  !!]
  call Unit_Tests_Begin_Group("HII region absorption postprocessor")

  call Unit_Tests_Begin_Group("Multiplier")
  call Assert("ionizing light from a young population is reduced to the escape fraction",absorption_%multiplier(        5.0d2                                    ,ageYoung,0.0d0),escapeFraction,relTol=1.0d-12)
  call Assert("ionizing light just shortward of the Lyman limit is reduced"             ,absorption_%multiplier(0.999d0*lymanSeriesLimitWavelengthHydrogen_atomic,ageYoung,0.0d0),escapeFraction,relTol=1.0d-12)
  call Assert("light at the Lyman limit is untouched"                                   ,absorption_%multiplier(        lymanSeriesLimitWavelengthHydrogen_atomic,ageYoung,0.0d0),1.0d0         ,relTol=1.0d-12)
  call Assert("non-ionizing light from a young population is untouched"                 ,absorption_%multiplier(        1.5d3                                    ,ageYoung,0.0d0),1.0d0         ,relTol=1.0d-12)
  call Assert("ionizing light from an old population escapes"                           ,absorption_%multiplier(        5.0d2                                    ,ageOld  ,0.0d0),1.0d0         ,relTol=1.0d-12)
  call Assert("the multiplier does not depend on redshift"                              ,absorption_%multiplier(        5.0d2                                    ,ageYoung,3.0d0),escapeFraction,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Flags")
  call absorption_%ageRange(ageMinimum,ageMaximum)
  call Assert("is not redshift dependent"      ,absorption_%isRedshiftDependent(),.false.)
  call Assert("reports a sharp age window"     ,absorption_%ageWindowIsSharp   (),.true. )
  call Assert("reports an unbounded age range" ,ageMinimum == 0.0d0 .and. ageMaximum == huge(0.0d0),.true.)
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("In a sequence with an age window")
  allocate(members     )
  allocate(members%next)
  members     %postprocessor_ => absorption_
  members%next%postprocessor_ => recent_
  !![
  <referenceConstruct object="sequence_" constructor="stellarPopulationSpectraPostprocessorSequence(members)"/>
  !!]
  call sequence_%ageRange(ageMinimum,ageMaximum)
  call Assert("the window remains sharp"                                ,sequence_%ageWindowIsSharp(),.true.)
  call Assert("the window is that of the age window"                   ,[ageMinimum,ageMaximum],[0.0d0,ageYoung],absTol=1.0d-15)
  call Assert("young ionizing light is reduced to the escape fraction" ,sequence_%multiplier(5.0d2,0.5d0*ageYoung,0.0d0),escapeFraction,relTol=1.0d-12)
  call Assert("ionizing light older than the window is suppressed"     ,sequence_%multiplier(5.0d2,2.0d0*ageYoung,0.0d0),0.0d0         ,absTol=1.0d-15)
  call Unit_Tests_End_Group()

  call Unit_Tests_End_Group()
  !![
  <objectDestructor name="sequence_"               />
  <objectDestructor name="recent_"                 />
  <objectDestructor name="absorption_"             />
  <objectDestructor name="hiiRegionEscapeFraction_"/>
  !!]
  call Unit_Tests_Finish()

end program Test_HII_Region_Absorption
