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

!% Contains a program which tests cooling function functionality.

program Test_Cooling_Functions
  !% Tests cooling function functionality.
  use :: Abundances_Structure            , only : abundances                             , metallicityTypeLinearByMassSolar
  use :: Chemical_Abundances_Structure   , only : chemicalAbundances                     , zeroChemicalAbundances
  use :: Chemical_States                 , only : chemicalStateAtomicCIECloudy
  use :: Cooling_Functions               , only : coolingFunctionSummation                        , coolingFunctionAtomicCIECloudy  , coolingFunctionCMBCompton, coolantList
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Display                         , only : displayVerbositySet                    , verbosityLevelStandard
  use :: Numerical_Constants_Astronomical, only : gigaYear
  use :: Numerical_Constants_Physical    , only : boltzmannsConstant
  use :: Numerical_Constants_Units       , only : ergs
  use :: Radiation_Fields                , only : radiationFieldCosmicMicrowaveBackground
  use :: Unit_Tests                      , only : Assert                                 , Unit_Tests_Begin_Group          , Unit_Tests_End_Group     , Unit_Tests_Finish
  implicit none
  type            (coolantList                            ), pointer :: coolants
  type            (coolingFunctionCMBCompton              ), target  :: coolingFunctionCMBCompton_
  type            (coolingFunctionSummation               )          :: coolingFunctionSummation_
  type            (coolingFunctionAtomicCIECloudy)                   :: coolingFunctionAtomicCIECloudy_
  type            (cosmologyParametersSimple              )          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda         )          :: cosmologyFunctions_
  type            (chemicalStateAtomicCIECloudy           )          :: chemicalState_
  type            (abundances                             )          :: gasAbundances
  type            (chemicalAbundances                     )          :: chemicalDensities
  type            (radiationFieldCosmicMicrowaveBackground)          :: radiation
  double precision                                                   :: numberDensityHydrogen          , temperature           , &
       &                                                                coolantSummation               , coolantCMBCompton     , &
       &                                                                timescaleCooling               , coolantAtomicCIECloudy

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cooling functions")
  ! Construct cooling functions.
  !# <referenceConstruct object="cosmologyParameters_"           >
  !#  <constructor>
  !#   cosmologyParametersSimple     (                                            &amp;
  !#    &amp;                         OmegaMatter         = 0.28120d0           , &amp;
  !#    &amp;                         OmegaBaryon         = 0.04611d0           , &amp;
  !#    &amp;                         OmegaDarkEnergy     = 0.71880d0           , &amp;
  !#    &amp;                         temperatureCMB      = 2.72548d0           , &amp;
  !#    &amp;                         HubbleConstant      =69.70000d0             &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="cosmologyFunctions_"            >
  !#  <constructor>
  !#   cosmologyFunctionsMatterLambda(                                            &amp;
  !#    &amp;                          cosmologyParameters_=cosmologyParameters_  &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="chemicalState_" >
  !#  <constructor>
  !#   chemicalStateAtomicCIECloudy  (                                            &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="coolingFunctionCMBCompton_"     >
  !#  <constructor>
  !#   coolingFunctionCMBCompton     (                                           &amp;
  !#    &amp;                          chemicalState_      =chemicalState_       &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="coolingFunctionAtomicCIECloudy_">
  !#  <constructor>
  !#   coolingFunctionAtomicCIECloudy(                                           &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  allocate(coolants     )
  allocate(coolants%next)
  coolants     %coolingFunction => coolingFunctionCMBCompton_
  coolants%next%coolingFunction => coolingFunctionCMBCompton_
  !# <referenceConstruct object="coolingFunctionSummation_"      >
  !#  <constructor>
  !#   coolingFunctionSummation      (                                           &amp;
  !#    &amp;                         coolants             =coolants             &amp;
  !#    &amp;                        )
  !#  </constructor>
  !# </referenceConstruct>
  ! Define plasma conditions.
  numberDensityHydrogen=1.0d-4 ! cm⁻³
  temperature          =1.0d+6 ! K
  chemicalDensities    =zeroChemicalAbundances
  radiation            =radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
  call radiation    %       timeSet(1.0d0                                                  )
  call gasAbundances%metallicitySet(1.0d0,metallicityType= metallicityTypeLinearByMassSolar)
  ! Summed cooling function should be twice the CMB Compton cooling function.
  coolantSummation       =  coolingFunctionSummation_ %coolingFunction                   (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunction                   (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('summation'        ,coolantSummation,2.0d0*coolantCMBCompton,relTol=1.0d-6)
  ! Logarithmic derivatives of summed cooling function should equal those of the CMB Compton cooling function.
  coolantSummation       =  coolingFunctionSummation_ %coolingFunctionDensityLogSlope    (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunctionDensityLogSlope    (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('density slope'    ,coolantSummation,      coolantCMBCompton,relTol=1.0d-6)
  ! Logarithmic derivatives of summed cooling function should equal those of the CMB Compton cooling function.
  coolantSummation       =  coolingFunctionSummation_ %coolingFunctionTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunctionTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('temperature slope',coolantSummation,      coolantCMBCompton,relTol=1.0d-6)
  ! Begin repeatibility tests.
  call Unit_Tests_Begin_Group("Repeatibility")
  ! Compute timescale for Compton cooling off of CMB at the present day.
  call radiation%timeSet(13.8d0)
  coolantCMBCompton=+coolingFunctionCMBCompton_%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  timescaleCooling =+1.5d0                 &
       &            *numberDensityHydrogen &
       &            *boltzmannsConstant    &
       &            *temperature           &
       &            /ergs                  &
       &            /coolantCMBCompton     &
       &            /gigaYear
  call Assert('CMB Compton cooling timescale at z=0',timescaleCooling,976.65342064763729d0,relTol=1.0d-6)
  ! Compute Cloudy cooling time.
  coolantAtomicCIECloudy=+coolingFunctionAtomicCIECloudy_%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('Cloudy CIE cooling function',coolantAtomicCIECloudy,1.4155d-30,relTol=1.0d-6)
  call Unit_Tests_End_Group       ()
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()
end program Test_Cooling_Functions
