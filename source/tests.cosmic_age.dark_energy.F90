!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

program Tests_Cosmic_Age_Dark_Energy
  !% Tests cosmic age calculations for a dark energy Universe. Ages calculated using Python implementation of Ned Wright's
  !% cosmology calculator available from: http://www.astro.ucla.edu/~wright/CC.python
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Memory_Management
  implicit none
  double precision, parameter, dimension(8) :: redshift     =[0.0000d0,1.0000d0,3.0000d0,9.0000d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision, parameter, dimension(8) :: ageDarkEnergy=[0.0942699818d0,0.0402619685d0,0.0147878493d0,0.0037621019d0&
       &,0.0006895264d0,0.0001172509d0,0.0000227901d0,0.0000037578d0]
  type(varying_string)                      :: parameterFile
  character(len=1024)                       :: message
  integer                                   :: iExpansion
  double precision                          :: age,expansionFactor

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.cosmic_age.dark_energy.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cosmic age: dark energy cosmology")

  ! Test cosmic age in a dark energy universe.
  parameterFile='testSuite/parameters/cosmicAge/darkEnergy.xml'
  call Input_Parameters_File_Open(parameterFile)
  do iExpansion=1,size(redshift)
     expansionFactor=Expansion_Factor_From_Redshift(redshift(iExpansion))
     age=Cosmology_Age(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageDarkEnergy(iExpansion),relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Cosmic_Age_Dark_Energy
