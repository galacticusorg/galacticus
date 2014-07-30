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

!% Contains a program which tests spherical collapse calculations for an open Universe.

program Tests_Spherical_Collapse_Open
  !% Tests spherical collapse calculations for an open Universe. Compares results to the fitting formula of
  !% \cite{bryan_statistical_1998}.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Memory_Management
  use Virial_Density_Contrast
  use Numerical_Constants_Math
  implicit none
  double precision                         , dimension(7) :: redshift                 =[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  class           (cosmologyFunctionsClass), pointer      :: cosmologyFunctionsDefault
  class           (virialDensityContrastClass), pointer      :: virialDensityContrast_
  type            (varying_string         )               :: parameterFile
  character       (len=1024               )               :: message
  integer                                                 :: iExpansion
  double precision                                        :: age                                                                     , bryanNormanFit , &
       &                                                     densityContrast                                                         , expansionFactor, &
       &                                                     x

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.spherical_collapse.open.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: open cosmology")

  ! Test spherical collapse in an open universe.
  parameterFile='testSuite/parameters/sphericalCollapse/open.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions   ()
  virialDensityContrast_    => virialDensityContrast()
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctionsDefault%expansionFactorFromRedshift(redshift(iExpansion))
     age            =cosmologyFunctionsDefault%cosmicTime(expansionFactor)
     densityContrast=virialDensityContrast_%densityContrast(age)
     x              =cosmologyFunctionsDefault%omegaMatterEpochal(age)-1.0d0
     bryanNormanFit =(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)/cosmologyFunctionsDefault%omegaMatterEpochal(age)
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctionsDefault%omegaMatterEpochal(age),"]"
     call Assert(trim(message),densityContrast,bryanNormanFit,relTol=1.0d-2)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Spherical_Collapse_Open
