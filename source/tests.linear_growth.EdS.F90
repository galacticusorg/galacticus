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

program Tests_Linear_Growth_EdS
  !% Tests linear growth calculations.
  use :: Cosmology_Functions, only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: ISO_Varying_String , only : assignment(=)      , varying_string
  use :: Input_Parameters   , only : inputParameters
  use :: Linear_Growth      , only : componentDarkMatter, linearGrowth           , linearGrowthClass
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                         , dimension(8), parameter :: redshift                 =[0.0d0,1.0d0,3.0d0,9.0d0,30.0d0,100.0d0,300.0d0,1000.0d0]
  class           (cosmologyFunctionsClass), pointer                 :: cosmologyFunctions_
  class           (linearGrowthClass      ), pointer                 :: linearGrowth_
  type            (varying_string         )                          :: parameterFile
  character       (len=1024               )                          :: message
  integer                                                            :: iExpansion
  double precision                                                   :: expansionFactor                                                                    , linearGrowthFactor
  type            (inputParameters        )                          :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: Einstein-de Sitter")

  ! Test growth factor in an Einstein-de Sitter universe. Growth factor should equal the expansion factor.
  parameterFile='testSuite/parameters/linearGrowth/EdS.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get the default cosmology functions object.
  cosmologyFunctions_ => cosmologyFunctions()
  linearGrowth_       => linearGrowth      ()
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowthFactor=linearGrowth_%value(expansionFactor=expansionFactor,component=componentDarkMatter)
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
     call Assert(trim(message),linearGrowthFactor,expansionFactor,relTol=1.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Linear_Growth_EdS
