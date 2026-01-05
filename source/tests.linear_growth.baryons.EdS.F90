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

program Tests_Linear_Growth_EdS_Baryons
  !!{
  Tests linear growth calculations.
  !!}
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Display                   , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Intergalactic_Medium_State, only : intergalacticMediumStateSimple
  use :: Linear_Growth             , only : componentDarkMatter           , linearGrowthBaryonsDarkMatter
  use :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group       , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                , dimension(5), parameter :: redshift                 =[0.0d0,1.0d0,3.0d0,9.0d0,30.0d0]
  type            (cosmologyParametersSimple     )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda)                          :: cosmologyFunctions_
  type            (linearGrowthBaryonsDarkMatter )                          :: linearGrowth_
  type            (intergalacticMediumStateSimple)                          :: intergalacticMediumState_
  character       (len=1024                      )                          :: message
  integer                                                                   :: i
  double precision                                                          :: expansionFactor                                           , linearGrowthFactor        , &
       &                                                                       exponent                                                  , linearGrowthFactorExpected

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Construct model.
  cosmologyParameters_=cosmologyParametersSimple          (                                                                 &
       &                                                   OmegaMatter                          = 1.00d0                  , &
       &                                                   OmegaBaryon                          = 0.05d0                  , &
       &                                                   OmegaDarkEnergy                      = 0.00d0                  , &
       &                                                   temperatureCMB                       = 2.78d0                  , &
       &                                                   HubbleConstant                       =70.00d0                    &
       &                                                  )
  cosmologyFunctions_ =cosmologyFunctionsMatterLambda     (                                                                 &
       &                                                                                         cosmologyParameters_       &
       &                                                  )
  intergalacticMediumState_=intergalacticMediumStateSimple(                                                                 &
       &                                                   reionizationRedshift                 = 8.00d0                  , &
       &                                                   reionizationTemperature              = 1.00d4                  , &
       &                                                   preReionizationTemperature           = 1.00d4                  , &
       &                                                   cosmologyFunctions_                  =cosmologyFunctions_      , &
       &                                                   cosmologyParameters_                 =cosmologyParameters_       &
       &                                                  )
  linearGrowth_       =linearGrowthBaryonsDarkMatter      (                                                                 &
       &                                                   redshiftInitial                      =100.0d0                  , &
       &                                                   redshiftInitialDelta                 =  1.0d0                  , &
       &                                                   cambCountPerDecade                   =  0                      , &
       &                                                   darkMatterOnlyInitialConditions      =.false.                  , &
       &                                                   cosmologyParameters_                 =cosmologyParameters_     , &
       &                                                   cosmologyParametersInitialConditions_=cosmologyParameters_     , &
       &                                                   cosmologyFunctions_                  =cosmologyFunctions_      , &
       &                                                   intergalacticMediumState_            =intergalacticMediumState_  &
       &                                                  )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: Einstein-de Sitter with baryons")
  ! Iterate over redshifts.
  do i=1,size(redshift)
     expansionFactor   =cosmologyFunctions_%expansionFactorFromRedshift(                                  redshift       (i)                              )
     ! For an Einstein-de Sitter universe with dark matter fraction, fₓ, in the large wavenumber limit (such that baryons do not
     ! cluster) the growing mode in the dark matter grows as δ ∝ a^{3p/2} with p=[√(1+24fₓ)-1]/6
     linearGrowthFactor=linearGrowth_      %value                      (wavenumber=1.0d+3,expansionFactor=expansionFactor   ,component=componentDarkMatter)
     exponent          =+(                                            &
          &               +sqrt(                                      &
          &                     + 1.0d0                               &
          &                     +24.0d0                               &
          &                     *(                                    &
          &                       +cosmologyParameters_%OmegaMatter() &
          &                       -cosmologyParameters_%OmegaBaryon() &
          &                      )                                    &
          &                     /  cosmologyParameters_%OmegaMatter() &
          &                    )                                      &
          &               -1.0d0                                      &
          &              )                                            &
          &             /6.0d0
     linearGrowthFactorExpected=+expansionFactor**(1.5d0*exponent)
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [k ⟶ ∞ ; z = ",redshift(i),"]"
     call Assert(trim(message),linearGrowthFactor,linearGrowthFactorExpected,relTol=1.0d-3)
     ! For an Einstein-de Sitter universe in the small wavenumber limit (such that baryons do cluster) the growing mode in the
     ! dark matter grows as δ ∝ a.
     linearGrowthFactor=linearGrowth_      %value                      (wavenumber=1.0d-3,expansionFactor=expansionFactor   ,component=componentDarkMatter)
     exponent          =+(                                            &
          &               +sqrt(                                      &
          &                     + 1.0d0                               &
          &                     +24.0d0                               &
          &                     *(                                    &
          &                       +cosmologyParameters_%OmegaMatter() &
          &                       -cosmologyParameters_%OmegaBaryon() &
          &                      )                                    &
          &                     /  cosmologyParameters_%OmegaMatter() &
          &                    )                                      &
          &               -1.0d0                                      &
          &              )                                            &
          &             /6.0d0
     linearGrowthFactorExpected=+expansionFactor
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [k ⟶ 0 ; z = ",redshift(i),"]"
     call Assert(trim(message),linearGrowthFactor,linearGrowthFactorExpected,relTol=1.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Linear_Growth_EdS_Baryons
