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

program Tests_Linear_Growth_Dark_Energy
  !% Tests linear growth calculations for a dark energy Universe. Growth rates are compared to Figure 1 of Linder and Jenkins
  !% (2003; MNRAS; 346; 573; http://adsabs.harvard.edu/abs/2003MNRAS.346..573L).
  use :: Cosmology_Functions, only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: ISO_Varying_String , only : assignment(=)      , varying_string
  use :: Input_Parameters   , only : inputParameters
  use :: Linear_Growth      , only : componentDarkMatter, linearGrowth           , linearGrowthClass   , normalizeMatterDominated
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                         , dimension(13), parameter :: redshift                 =[0.000000d0,0.052632d0,0.149425d0,0.265823d0,0.449275d0,0.666667d0,1.000000d0,1.325580d0,1.857140d0,2.846150d0,4.555560d0,8.090910d0,17.867900d0]
  double precision                         , dimension(13), parameter :: growthFactorDarkEnergy   =[0.73d0,0.75d0,0.78d0,0.81d0,0.85d0,0.88d0,0.92d0,0.94d0,0.96d0,0.98d0,0.99d0,1.00d0,1.00d0]
  class           (cosmologyFunctionsClass), pointer                  :: cosmologyFunctions_
  class           (linearGrowthClass      ), pointer                  :: linearGrowth_
  type            (varying_string         )                           :: parameterFile
  character       (len=1024               )                           :: message
  integer                                                             :: iExpansion
  double precision                                                    :: expansionFactor                                                                                                                                                            , linearGrowthFactor
  type            (inputParameters        )                           :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: dark energy cosmology")

  ! Test growth factor in a dark energy universe.
  parameterFile='testSuite/parameters/linearGrowth/darkEnergy.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get the default cosmology functions object.
  cosmologyFunctions_ => cosmologyFunctions()
  linearGrowth_       => linearGrowth      ()
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
