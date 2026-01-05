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

program Tests_Linear_Growth_Cosmological_Constant
  !!{
  Tests linear growth calculations for a cosmological constant Universe. Growth rates are compared to calculations taken from
  Andrew Hamilton's \href{https://jila.colorado.edu/~ajsh/growl/index.html}{growl}.
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Cosmology_Functions , only : cosmologyFunctionsMatterLambda
  use :: Display             , only : displayVerbositySet           , verbosityLevelStandard
  use :: Linear_Growth       , only : componentDarkMatter           , linearGrowthCollisionlessMatter, normalizeMatterDominated
  use :: Unit_Tests          , only : Assert                        , Unit_Tests_Begin_Group         , Unit_Tests_End_Group    , Unit_Tests_Finish
  implicit none
  double precision                                 , dimension(8), parameter :: redshift              =[0.000d0,1.0000d0,3.0000d0,9.0d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision                                 , dimension(8), parameter :: growthFactorDarkEnergy=[0.7789810167707876d0,0.9531701355446482d0,0.9934824792317063d0,0.9995762227500181d0,0.9999857599010360d0,0.9999995882349219d0,0.9999999844434028d0,0.9999999995770291d0]
  type            (cosmologyParametersSimple      )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda )                          :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter)                          :: linearGrowth_
  character       (len=1024                       )                          :: message
  integer                                                                    :: iExpansion
  double precision                                                           :: expansionFactor                                                                                                                                                                                 , linearGrowthFactor

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: cosmological constant cosmology")
  ! Test growth factor in a dark energy universe.
  !![
  <referenceConstruct object="cosmologyParameters_">
   <constructor>
    cosmologyParametersSimple      (                                           &amp;
     &amp;                          OmegaMatter         = 0.30d0             , &amp;
     &amp;                          OmegaBaryon         = 0.00d0             , &amp;
     &amp;                          OmegaDarkEnergy     = 0.70d0             , &amp;
     &amp;                          temperatureCMB      = 2.78d0             , &amp;
     &amp;                          HubbleConstant      =73.00d0               &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_" >
   <constructor>
    cosmologyFunctionsMatterLambda (                                           &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_  &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"       >
   <constructor>
    linearGrowthCollisionlessMatter(                                           &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_, &amp;
     &amp;                          cosmologyFunctions_ =cosmologyFunctions_   &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  !!]
  do iExpansion=1,size(redshift)
     expansionFactor   =cosmologyFunctions_%expansionFactorFromRedshift(                                                      &
          &                                                                              redshift               (iExpansion)  &
          &                                                            )
     linearGrowthFactor=linearGrowth_      %value                      (                                                      &
          &                                                             expansionFactor=expansionFactor                     , &
          &                                                             component      =componentDarkMatter                 , &
          &                                                             normalize      =normalizeMatterDominated              &
          &                                                            )                                                      &
          &                                                            /expansionFactor
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
     call Assert(trim(message),linearGrowthFactor,growthFactorDarkEnergy(iExpansion),relTol=1.0d-3)
  end do

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Linear_Growth_Cosmological_Constant
