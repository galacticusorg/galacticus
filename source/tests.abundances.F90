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

!% Contains a program to test abundances objects functions.

program Test_Abundances
  !% Test abundances objects.
  use Unit_Tests
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Abundances_Structure
  implicit none
  type            (varying_string)               :: parameterFile
  type            (abundances    )               :: abundances1    , abundances2, abundances3
  double precision                , dimension(5) :: abundancesArray

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.abundances.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Abundances objects functions")

  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/abundances/testAbundances.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Initialize abundances.
  call abundances1%deserialize([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0])
  call abundances2%deserialize([5.0d0,4.0d0,3.0d0,2.0d0,1.0d0])
  abundances3=max(abundances1,abundances2)
  call abundances3%serialize(abundancesArray)
  call Assert('max() function',abundancesArray,[5.0d0,4.0d0,3.0d0,4.0d0,5.0d0])

  ! Close the input parameter file.
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Abundances
