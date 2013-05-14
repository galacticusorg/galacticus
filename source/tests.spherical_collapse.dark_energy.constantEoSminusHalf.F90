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

!% Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using a flat,
!% $\omega=-1/2$ cosmology.

program Tests_Spherical_Collapse_Dark_Energy_Omega_Half
  !% Tests spherical collapse calculations for a dark energy Universe, specifically using a flat, $\omega=-1/2$
  !% cosmology. Compares results to the fitting function of
  !% \citeauthor{weinberg_constraining_2003}~(\citeyear{weinberg_constraining_2003}; eqn.~18).
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Memory_Management
  use Numerical_Constants_Math
  use Critical_Overdensity
  use Linear_Growth
  implicit none
  double precision         , dimension(7) :: redshift=[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  type     (varying_string)               :: parameterFile
  character(len=1024      )               :: message
  integer                                 :: iExpansion
  double precision                        :: age,expansionFactor,criticalOverdensity,criticalOverdensityExpected,omega,alpha
  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.spherical_collapse.dark_energy.constantEoSminusHalf.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (ω=-½ cosmology)")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/darkEnergy.constantEoSminusHalf.xml'
  call Input_Parameters_File_Open(parameterFile)
  do iExpansion=1,size(redshift)
     expansionFactor            =Expansion_Factor_From_Redshift(redshift(iExpansion))
     age                        =Cosmology_Age(expansionFactor)
     criticalOverdensity        =Critical_Overdensity_for_Collapse(age)
     omega                      =Cosmology_Dark_Energy_Equation_Of_State(age)
     alpha                      =0.353d0*omega**4+1.044d0*omega**3+1.128d0*omega**2+0.555d0*omega+0.131d0
     criticalOverdensityExpected=(3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)*(1.0d0+alpha*log10(Omega_Matter_Total(age)))/Linear_Growth_Factor(age)
     write (message,'(a,f6.1,a,f6.4,a)') "critical density for collapse [z=",redshift(iExpansion),";Ωₘ=",Omega_Matter_Total(age),"]"
     call Assert(trim(message),criticalOverdensity,criticalOverdensityExpected,relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Spherical_Collapse_Dark_Energy_Omega_Half
