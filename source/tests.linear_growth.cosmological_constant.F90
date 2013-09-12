!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  !% Tests linear growth calculations for a cosmological constant Universe. Growth rates are compared to calculations taken from
  !% Andrew Hamilton's "growl" code available at: http://casa.colorado.edu/~ajsh/growl/
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Linear_Growth
  use Cosmology_Functions
  use Memory_Management
  implicit none
  double precision                         , dimension(8), parameter :: redshift                 =[0.000d0,1.0000d0,3.0000d0,9.0d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision                         , dimension(8), parameter :: growthFactorDarkEnergy   =[0.7789810167707876d0,0.9531701355446482d0,0.9934824792317063d0,0.9995762227500181d0,0.9999857599010360d0,0.9999995882349219d0,0.9999999844434028d0,0.9999999995770291d0]
  class           (cosmologyFunctionsClass), pointer                 :: cosmologyFunctionsDefault
  type            (varying_string         )                          :: parameterFile
  character       (len=1024               )                          :: message
  integer                                                            :: iExpansion
  double precision                                                   :: expansionFactor                                                                                                                                                                                    , linearGrowth

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.linear_growth.cosmological_constant.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: cosmological constant cosmology")

  ! Test growth factor in a dark energy universe.
  parameterFile='testSuite/parameters/linearGrowth/cosmologicalConstant.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions()
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctionsDefault%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowth=Linear_Growth_Factor(aExpansion=expansionFactor,component=linearGrowthComponentDarkMatter,normalize=normalizeMatterDominated)/expansionFactor
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
     call Assert(trim(message),linearGrowth,growthFactorDarkEnergy(iExpansion),relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Linear_Growth_Cosmological_Constant
