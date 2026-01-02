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

program Tests_Linear_Growth_Dark_Energy
  !!{
  Tests linear growth calculations for a dark energy Universe. Growth rates are compared to Figure 1 of Linder and Jenkins
  (2003; MNRAS; 346; 573; http://adsabs.harvard.edu/abs/2003MNRAS.346..573L).
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Cosmology_Functions , only : cosmologyFunctionsMatterDarkEnergy
  use :: Display             , only : displayVerbositySet               , verbosityLevelStandard
  use :: Linear_Growth       , only : componentDarkMatter               , linearGrowthCollisionlessMatter, normalizeMatterDominated
  use :: Unit_Tests          , only : Assert                            , Unit_Tests_Begin_Group         , Unit_Tests_End_Group    , Unit_Tests_Finish
  implicit none
  double precision                                    , dimension(13), parameter :: redshift              =[0.000000d0,0.052632d0,0.149425d0,0.265823d0,0.449275d0,0.666667d0,1.000000d0,1.325580d0,1.857140d0,2.846150d0,4.555560d0,8.090910d0,17.867900d0]
  double precision                                    , dimension(13), parameter :: growthFactorDarkEnergy=[0.73d0,0.75d0,0.78d0,0.81d0,0.85d0,0.88d0,0.92d0,0.94d0,0.96d0,0.98d0,0.99d0,1.00d0,1.00d0]
  type            (cosmologyParametersSimple         )                           :: cosmologyParameters_
  type            (cosmologyFunctionsMatterDarkEnergy)                           :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter   )                           :: linearGrowth_
  character       (len=1024                          )                           :: message
  integer                                                                        :: iExpansion
  double precision                                                               :: expansionFactor                                                                                                                                                        , linearGrowthFactor

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: dark energy cosmology")
  ! Test growth factor in a dark energy universe.
  !![
  <referenceConstruct object="cosmologyParameters_">
   <constructor>
    cosmologyParametersSimple      (                                                     &amp;
     &amp;                          OmegaMatter                   = 0.30d0             , &amp;
     &amp;                          OmegaBaryon                   = 0.00d0             , &amp;
     &amp;                          OmegaDarkEnergy               = 0.70d0             , &amp;
     &amp;                          temperatureCMB                = 2.78d0             , &amp;
     &amp;                          HubbleConstant                =73.00d0               &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_" >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                  &amp;
     &amp;                             cosmologyParameters_       =cosmologyParameters_, &amp;
     &amp;                             darkEnergyEquationOfStateW0=-0.8d0              , &amp;
     &amp;                             darkEnergyEquationOfStateW1=+0.0d0                &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"       >
   <constructor>
    linearGrowthCollisionlessMatter(                                                     &amp;
     &amp;                          cosmologyParameters_          =cosmologyParameters_, &amp;
     &amp;                          cosmologyFunctions_           =cosmologyFunctions_   &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  !!]
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowthFactor=linearGrowth_%value(expansionFactor=expansionFactor,component=componentDarkMatter,normalize=normalizeMatterDominated)/expansionFactor
     write (message,'(a,f7.2,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
    call Assert(trim(message),linearGrowthFactor,growthFactorDarkEnergy(iExpansion),relTol=5.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Linear_Growth_Dark_Energy
