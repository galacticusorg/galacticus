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

program Tests_Linear_Growth_Cosmological_Constant
  !!{RST
  Tests linear growth calculations for a cosmological constant Universe. Growth rates are compared to calculations taken from Andrew Hamilton's `growl <https://jila.colorado.edu/~ajsh/growl/index.html>`_.

  Also tests that the tabulated growth factor is independent of the order in which epochs are requested. Two objects are presented
  with the same epochs in opposite orders, each value being recorded as soon as it is requested, and the results are required to
  agree *exactly*. The epochs reach beyond the initial tabulation at both ends, so that values recorded before the tabulation is
  extended (to earlier or later epochs) are compared with values recorded after. (When the initial condition was imposed at the
  earliest tabulated epoch, extending the tabulation to earlier epochs rebuilt it, changing every value slightly - see issue
  #1780.)
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Cosmology_Functions , only : cosmologyFunctionsMatterLambda
  use :: Display             , only : displayVerbositySet           , verbosityLevelStandard
  use :: Linear_Growth       , only : componentDarkMatter           , linearGrowthCollisionlessMatter, normalizeMatterDominated, linearGrowthNonClusteringBaryonsDarkMatter
  use :: Unit_Tests          , only : Assert                        , Unit_Tests_Begin_Group         , Unit_Tests_End_Group    , Unit_Tests_Finish
  implicit none
  double precision                                 , dimension(8), parameter :: redshift              =[0.000d0,1.0000d0,3.0000d0,9.0d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision                                 , dimension(8), parameter :: growthFactorDarkEnergy=[0.7789810167707876d0,0.9531701355446482d0,0.9934824792317063d0,0.9995762227500181d0,0.9999857599010360d0,0.9999995882349219d0,0.9999999844434028d0,0.9999999995770291d0]
  ! Epochs (in Gyr) at which the determinism of the tabulation is tested.
  double precision                                 , dimension(6), parameter :: timesDeterminism      =[3.1d-9,3.7d-3,0.83d0,2.9d0,9.1d0,47.0d0]
  type            (cosmologyParametersSimple      )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda )                          :: cosmologyFunctions_                                                                                                                                                                             , cosmologyFunctionsAscending_, &
       &                                                                        cosmologyFunctionsDescending_
  type            (linearGrowthCollisionlessMatter)                          :: linearGrowth_                                                                                                                                                                                   , linearGrowthAscending_      , &
       &                                                                        linearGrowthDescending_
  double precision                                 , dimension(6)            :: resultAscending                                                                                                                                                                                 , resultDescending
  character       (len=1024                       )                          :: message
  integer                                                                    :: iExpansion                                                                                                                                                                                      , i
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

  ! Test that the growth factor is independent of the order in which epochs are requested. Each growth object is given its own
  ! cosmology functions object, so that the tabulation of the expansion factor is also built up in opposite orders.
  !![
  <referenceConstruct object="cosmologyFunctionsAscending_" >
   <constructor>
    cosmologyFunctionsMatterLambda (                                                    &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_           &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsDescending_">
   <constructor>
    cosmologyFunctionsMatterLambda (                                                    &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_           &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowthAscending_"       >
   <constructor>
    linearGrowthCollisionlessMatter(                                                    &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_         , &amp;
     &amp;                          cosmologyFunctions_ =cosmologyFunctionsAscending_   &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowthDescending_"      >
   <constructor>
    linearGrowthCollisionlessMatter(                                                    &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_         , &amp;
     &amp;                          cosmologyFunctions_ =cosmologyFunctionsDescending_  &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  !!]
  do i=1,size(timesDeterminism)
     resultAscending (i)=linearGrowthAscending_ %value(time=timesDeterminism(                         i))
     resultDescending(i)=linearGrowthDescending_%value(time=timesDeterminism(size(timesDeterminism)+1-i))
  end do
  call Assert("growth factor is independent of the order in which epochs are requested"                         ,resultAscending,resultDescending(size(resultDescending):1:-1),absTol=0.0d0)
  do i=1,size(timesDeterminism)
     resultAscending (i)=linearGrowthAscending_ %value(time=timesDeterminism(i),normalize=normalizeMatterDominated)
     resultDescending(i)=linearGrowthDescending_%value(time=timesDeterminism(i),normalize=normalizeMatterDominated)
  end do
  call Assert("growth factor normalized to the matter dominated solution is independent of the order of requests",resultAscending,resultDescending                                ,absTol=0.0d0)

  ! Repeat the order test for the growth factor of dark matter with non-clustering baryons, which is tabulated in the same way.
  block
    type   (cosmologyFunctionsMatterLambda            ) :: cosmologyFunctionsBaryonsAscending_, cosmologyFunctionsBaryonsDescending_
    type   (linearGrowthNonClusteringBaryonsDarkMatter) :: linearGrowthBaryonsAscending_      , linearGrowthBaryonsDescending_
    !![
    <referenceConstruct object="cosmologyFunctionsBaryonsAscending_" >
     <constructor>
      cosmologyFunctionsMatterLambda            (                                                           &amp;
       &amp;                                     cosmologyParameters_=cosmologyParameters_                  &amp;
       &amp;                                    )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsBaryonsDescending_">
     <constructor>
      cosmologyFunctionsMatterLambda            (                                                           &amp;
       &amp;                                     cosmologyParameters_=cosmologyParameters_                  &amp;
       &amp;                                    )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="linearGrowthBaryonsAscending_"       >
     <constructor>
      linearGrowthNonClusteringBaryonsDarkMatter(                                                           &amp;
       &amp;                                     cosmologyParameters_=cosmologyParameters_                , &amp;
       &amp;                                     cosmologyFunctions_ =cosmologyFunctionsBaryonsAscending_   &amp;
       &amp;                                    )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="linearGrowthBaryonsDescending_"      >
     <constructor>
      linearGrowthNonClusteringBaryonsDarkMatter(                                                           &amp;
       &amp;                                     cosmologyParameters_=cosmologyParameters_                , &amp;
       &amp;                                     cosmologyFunctions_ =cosmologyFunctionsBaryonsDescending_  &amp;
       &amp;                                    )
     </constructor>
    </referenceConstruct>
    !!]
    do i=1,size(timesDeterminism)
       resultAscending (i)=linearGrowthBaryonsAscending_ %value(time=timesDeterminism(                         i))
       resultDescending(i)=linearGrowthBaryonsDescending_%value(time=timesDeterminism(size(timesDeterminism)+1-i))
    end do
    call Assert("growth factor with non-clustering baryons is independent of the order in which epochs are requested",resultAscending,resultDescending(size(resultDescending):1:-1),absTol=0.0d0)
  end block

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Linear_Growth_Cosmological_Constant
