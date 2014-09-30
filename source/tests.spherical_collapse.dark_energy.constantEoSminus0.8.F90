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

!% Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using a flat,
!% $\omega=-0.8$ cosmology.

program Tests_Spherical_Collapse_Dark_Energy_Omega_Zero_Point_Eight
  !% Tests spherical collapse calculations for a dark energy Universe, specifically using a flat, $\omega=-0.8$
  !% cosmology. Compares results to points read from Figure~6 of \cite{horellou_dark_2005} using \href{http://datathief.org/}{\sc
  !% DataThief}.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Cosmology_Functions
  use Virial_Density_Contrast
  implicit none
  double precision                            , dimension(3) :: redshift                     =[0.00d0,1.00d0,2.00d0]
  double precision                            , dimension(3) :: virialDensityContrastExpected=[367.81d0,217.63d0,192.72d0]
  class           (cosmologyFunctionsClass   ), pointer      :: cosmologyFunctionsDefault
  class           (virialDensityContrastClass), pointer      :: virialDensityContrast_
  double precision                            , parameter    :: massDummy                    =1.0d0
  type            (varying_string            )               :: parameterFile
  character       (len=1024                  )               :: message
  integer                                                    :: iExpansion
  double precision                                           :: age                                                       , expansionFactor, &
       &                                                        virialDensityContrastActual

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.spherical_collapse.dark_energy.constantEoSminus0.8.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (ω=-0.8 cosmology)")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/darkEnergy.constantEoSminus0.8.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions   ()
  virialDensityContrast_    => virialDensityContrast()
  do iExpansion=1,size(redshift)
     expansionFactor            =cosmologyFunctionsDefault%expansionFactorFromRedshift(redshift(iExpansion))
     age                        =cosmologyFunctionsDefault%cosmicTime(expansionFactor)
     virialDensityContrastActual=virialDensityContrast_   %densityContrast(massDummy,age)
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctionsDefault%omegaMatterEpochal(age),"]"
     call Assert(trim(message),virialDensityContrastActual,virialDensityContrastExpected(iExpansion),relTol=5.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Spherical_Collapse_Dark_Energy_Omega_Zero_Point_Eight
