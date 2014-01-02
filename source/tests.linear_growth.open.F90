!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

program Tests_Linear_Growth_Open
  !% Tests linear growth calculations for an open Universe. Growth rates are compared to calculations taken from:
  !% http://www.icosmos.co.uk/index.html
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Linear_Growth
  use Cosmology_Functions
  implicit none
  double precision                         , dimension(8), parameter :: redshift                 =[0.0000d0,1.0000d0,3.0000d0,9.0000d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision                         , dimension(8), parameter :: growthFactorOpen         =[0.4568354614082405d0,0.6176697062100953d0,0.7581803450845095d0,0.8844205217773703d0,0.9590358011003045d0,0.9869986440083428d0,0.9955930852515837d0,0.9986700650973155d0]
  class           (cosmologyFunctionsClass), pointer                 :: cosmologyFunctionsDefault
  type            (varying_string         )                          :: parameterFile
  character       (len=1024               )                          :: message
  integer                                                            :: iExpansion
  double precision                                                   :: expansionFactor                                                                                                                                                                                    , linearGrowth

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: open cosmology")

  ! Test growth factor in an open universe.
  parameterFile='testSuite/parameters/linearGrowth/open.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions()
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctionsDefault%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowth=Linear_Growth_Factor(aExpansion=expansionFactor,component=linearGrowthComponentDarkMatter,normalize=normalizeMatterDominated)/expansionFactor
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
     call Assert(trim(message),linearGrowth,growthFactorOpen(iExpansion),relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Linear_Growth_Open
