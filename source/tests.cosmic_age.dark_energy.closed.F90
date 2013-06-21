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

program Tests_Cosmic_Age_Dark_Energy_Closed
  !% Tests cosmic age calculations for a dark energy Universe with no dark energy but more than the critical density of matter
  !% such that it recollapses. Ages calculated using Python implementation of Ned Wright's cosmology calculator available from:
  !% http://www.astro.ucla.edu/~wright/CC.python
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Cosmological_Parameters
  use Numerical_Constants_Math
  use Memory_Management
  implicit none
  double precision                , dimension(8), parameter :: redshift                =[0.0000d0,1.0000d0,3.0000d0,9.0000d0,30.000000d0,100.0000d0,300.000000d0,1000.000d0]
  double precision                , dimension(8), parameter :: ageClosed               =[3.436956d-02,8.612654d-03,2.775283d-03,6.703704d-04,1.204878d-04,2.036305d-05,3.950938d-06,6.510672d-07]
  type            (varying_string)                          :: parameterFile
  character       (len=1024      )                          :: message
  integer                                                   :: iExpansion
  double precision                                          :: age                                                                                                                               , expansionFactor, &
       &                                                       expansionFactorSymmetric                                                                                                          , timeTurnaround

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.cosmic_age.dark_energy.closed.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cosmic age: dark energy solver, closed cosmology")

  ! Open the parameter file.
  parameterFile='testSuite/parameters/cosmicAge/darkEnergyClosed.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Compute the time of maximum expansion for the Universe. In this simple, OmegaM=10, OmegaDE=0 Universe this is analytically
  ! calculable and equals:
  timeTurnaround=(5.0d0/27.0d0)*Pi/H_0_invGyr()

  ! Test cosmic age in a dark energy universe.
  do iExpansion=1,size(redshift)
     expansionFactor=Expansion_Factor_From_Redshift(redshift(iExpansion))
     age=Cosmology_Age(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age        [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageClosed(iExpansion),relTol=1.0d-3)
     expansionFactorSymmetric=Expansion_Factor(2.0d0*timeTurnaround-age)
     write (message,'(a,f6.1,a)') "collapse symmetry [z=",redshift(iExpansion),"]"
     call Assert(trim(message),expansionFactor,expansionFactorSymmetric,absTol=0.03d0)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Cosmic_Age_Dark_Energy_Closed
