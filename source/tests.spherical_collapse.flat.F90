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

!% Contains a program which tests spherical collapse calculations for a flat Universe.

program Tests_Spherical_Collapse_Flat
  !% Tests spherical collapse calculations for a flat Universe. Compares results to the fitting formula of
  !% \cite{bryan_statistical_1998}.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Memory_Management
  use Virial_Density_Contrast
  use Numerical_Constants_Math
  implicit none
  double precision         , dimension(7) :: redshift=[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  type     (varying_string)               :: parameterFile
  character(len=1024      )               :: message
  integer                                 :: iExpansion
  double precision                        :: age,expansionFactor,densityContrast,x,bryanNormanFit

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.spherical_collapse.flat.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: flat cosmology")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/flat.xml'
  call Input_Parameters_File_Open(parameterFile)
  do iExpansion=1,size(redshift)
     expansionFactor=Expansion_Factor_From_Redshift(redshift(iExpansion))
     age            =Cosmology_Age(expansionFactor)
     densityContrast=Halo_Virial_Density_Contrast(age)
     x              =Omega_Matter_Total(age)-1.0d0
     bryanNormanFit =(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/Omega_Matter_Total(age)
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast [z=",redshift(iExpansion),";Ωₘ=",Omega_Matter_Total(age),"]"
     call Assert(trim(message),densityContrast,bryanNormanFit,relTol=1.1d-2)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Spherical_Collapse_Flat
