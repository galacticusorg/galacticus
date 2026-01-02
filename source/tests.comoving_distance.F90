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

program Tests_Comoving_Distance
  !!{
  Tests comoving distance calculations for various universes. Distances calculated using Python
  \href{http://www.astro.ucla.edu/~wright/CC.python}{implementation} of Ned Wright's cosmology
  calculator.
  !!}
  use :: Cosmology_Functions             , only : cosmologyFunctions       , cosmologyFunctionsClass, cosmologyFunctionsMatterLambda
  use :: Cosmology_Functions_Options     , only : distanceTypeComoving
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple, hubbleUnitsTime
  use :: Display                         , only : displayVerbositySet      , verbosityLevelStandard
  use :: Unit_Tests                      , only : Assert                   , Unit_Tests_Begin_Group , Unit_Tests_End_Group          , Unit_Tests_Finish
  use :: Numerical_Constants_Physical    , only : speedLight
  use :: Numerical_Constants_Astronomical, only : megaParsec               , gigaYear
  implicit none
  double precision                                , dimension(8), parameter :: redshift                               =[0.1000000d0, 1.0000000d0, 3.0000000d0, 9.0000000d0,30.0000000d0,100.0000000d0,300.0000000d0,1000.0000000d0]
  double precision                                , dimension(8), target    :: distanceEdS                            =[2.7903130d0,17.5614328d0,29.9792345d0,40.9979133d0,49.1894753d0, 53.9916697d0, 56.4991337d0,  58.0449076d0]
  double precision                                , dimension(8), target    :: distanceOpen                           =[2.8365460d0,19.5660920d0,36.2487160d0,53.3136690d0,67.2862090d0, 75.8567320d0, 80.4014770d0,  83.2171590d0]
  double precision                                , dimension(8), target    :: distanceCosmologicalConstant           =[2.9291813d0,23.1267935d0,44.4897529d0,64.4722404d0,79.4221370d0, 88.1893579d0, 92.7669523d0,  95.5884181d0]
  double precision                                , dimension(:), pointer   :: distance_
  class           (cosmologyFunctionsClass       )              , pointer   :: cosmologyFunctions_
  type            (cosmologyParametersSimple     )                          :: cosmologyParametersCosmologicalConstant                                                                                                             , cosmologyParametersOpen, &
       &                                                                       cosmologyParametersEdS
  type            (cosmologyFunctionsMatterLambda)              , target    :: cosmologyFunctionsCosmologicalConstant                                                                                                              , cosmologyFunctionsOpen , &
       &                                                                       cosmologyFunctionsEdS
  character       (len=1024                      )                          :: message
  integer                                                                   :: i                                                                                                                                                   , iExpansion
  double precision                                                          :: distance                                                                                                                                            , distanceModulus        , &
       &                                                                       time                                                                                                                                                , timeLookup

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Comoving distance")
  ! Cosmology functions for an Einstein-de Sitter Universe.
  !![
  <referenceConstruct object="cosmologyParametersEdS"                 >
   <constructor>
    cosmologyParametersSimple     (                                                              &amp;
     &amp;                         OmegaMatter         =1.00d0                                 , &amp;
     &amp;                         OmegaBaryon         =0.00d0                                 , &amp;
     &amp;                         OmegaDarkEnergy     =0.00d0                                 , &amp;
     &amp;                         temperatureCMB      =2.78d0                                 , &amp;
     &amp;                         HubbleConstant      =1.00d4                                   &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsEdS"                  >
   <constructor>
    cosmologyFunctionsMatterLambda(                                                              &amp;
     &amp;                         cosmologyParameters_=cosmologyParametersEdS                   &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Cosmology functions for an open Universe.
  !![
  <referenceConstruct object="cosmologyParametersOpen"                >
   <constructor>
    cosmologyParametersSimple     (                                                              &amp;
     &amp;                         OmegaMatter         =0.30d0                                 , &amp;
     &amp;                         OmegaBaryon         =0.00d0                                 , &amp;
     &amp;                         OmegaDarkEnergy     =0.00d0                                 , &amp;
     &amp;                         temperatureCMB      =2.78d0                                 , &amp;
     &amp;                         HubbleConstant      =1.00d4                                   &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsOpen"                 >
   <constructor>
    cosmologyFunctionsMatterLambda(                                                              &amp;
     &amp;                         cosmologyParameters_=cosmologyParametersOpen                  &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Cosmology functions for a cosmological constant Universe.
  !![
  <referenceConstruct object="cosmologyParametersCosmologicalConstant">
   <constructor>
    cosmologyParametersSimple     (                                                              &amp;
     &amp;                         OmegaMatter         =0.30d0                                 , &amp;
     &amp;                         OmegaBaryon         =0.00d0                                 , &amp;
     &amp;                         OmegaDarkEnergy     =0.70d0                                 , &amp;
     &amp;                         temperatureCMB      =2.78d0                                 , &amp;
     &amp;                         HubbleConstant      =1.00d4                                   &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsCosmologicalConstant" >
   <constructor>
    cosmologyFunctionsMatterLambda(                                                              &amp;
     &amp;                         cosmologyParameters_=cosmologyParametersCosmologicalConstant  &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Iterate over cosmologies.
  do i=1,3
     select case (i)
     case (1)
        call Unit_Tests_Begin_Group("Einstein-de Sitter"   )
        cosmologyFunctions_ => cosmologyFunctionsEdS
        distance_           => distanceEdS
        ! For Einstein-de Sitter test the calculation of the comoving distance to the particle horizon. This is
        ! ∫_0^t c dt' / a(t') = 2c/H(t)
        call Assert(                                                                                                       &
             &             "Comoving distance to particle horizon"                                                       , &
             &             cosmologyFunctions_    %distanceParticleHorizonComoving(cosmologyFunctions_%cosmicTime(1.0d0)), &
             &             +2.0d0                                                                                          &
             &             *speedLight                                                                                     &
             &             *gigaYear                                                                                       &
             &             /megaParsec                                                                                     &
             &             /cosmologyParametersEdS%HubbleConstant                 (hubbleUnitsTime                      ), &
             &      relTol=+1.0d-5                                                                                         &
             &     )
     case (2)
        call Unit_Tests_Begin_Group("Open universe"        )
        cosmologyFunctions_ => cosmologyFunctionsOpen
        distance_           => distanceOpen
     case (3)
        call Unit_Tests_Begin_Group("Cosmological constant")
        cosmologyFunctions_ => cosmologyFunctionsCosmologicalConstant
        distance_           => distanceCosmologicalConstant        
     end select
     do iExpansion=1,size(redshift)
        time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion)))
        distance=cosmologyFunctions_%distanceComoving(time)
        write (message,'(a,f6.1,a)') "comoving distance [        z=",redshift(iExpansion),"    ]"
        call Assert(trim(message),distance,distance_(iExpansion),relTol=1.0d-3)
        timeLookup=cosmologyFunctions_%timeAtDistanceComoving(distance)
        write (message,'(a,f6.1,a)') "cosmic time       [distance =",distance            ," Mpc]"
        call Assert(trim(message),timeLookup,time,relTol=1.0d-3)
        distance=cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshift(iExpansion))
        write (message,'(a,f6.1,a)') "comoving distance [direct; z=",redshift(iExpansion),"    ]"
        call Assert(trim(message),distance,distance_(iExpansion),relTol=1.0d-3)
        distanceModulus=25.0d0+5.0d0*log10(distance*(1.0d0+redshift(iExpansion)))
        distance=cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,distanceModulus=distanceModulus)
        write (message,'(a,f6.1,a)') "comoving distance [        D=",distanceModulus     ,"    ]"
        call Assert(trim(message),distance,distance_(iExpansion),relTol=1.0d-3)
     end do
     call Unit_Tests_End_Group()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Comoving_Distance
