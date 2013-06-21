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

!% Contains a program that tests power spectrum calculations.

program Tests_Power_Spectrum
  !% Tests power spectrum calculations.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Power_Spectra
  use Numerical_Constants_Math
  use Cosmological_Parameters
  implicit none
  double precision                , parameter :: radiusNormalization=8.0d0        !  Radius for sigma(M) normalization in Mpc/h.
  type            (varying_string)            :: parameterFile
  double precision                            :: mass                     , ratio                                                , &
       &                                         sigma

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.power_spectrum.size')
  ! Read parameters.
  parameterFile='testSuite/parameters/powerSpectrum.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Power spectra")
  ! Test that sigma_8 is correctly recovered.
  mass=(4.0d0*Pi/3.0d0)*Omega_Matter()*Critical_Density()*(radiusNormalization/Little_H_0())**3
  call Assert('σ₈ consistency',Cosmological_Mass_Root_Variance(mass),sigma_8(),relTol=1.0d-6)
  ! Test that sigma(M) scales as expected.
  ratio=Cosmological_Mass_Root_Variance(1.0d10)/Cosmological_Mass_Root_Variance(1.0d12)
  call Assert('σ(M) scaling',ratio,100.0d0**((-1.0d0+3.0d0)/6.0d0),relTol=1.0d-6)
  ! Test power spectrum normalization. For a power-law n=-1 power spectrum, the integral over k^2 P(k) W(k)^2/2 Pi^2 can be
  ! computed analytically and is equal to 9/[8 Pi^2]. We can therefore express this integral as 2 Pi P(k_8) 9/[8 Pi^2] / R_8^3
  ! where R_8=(8/h)Mpc and k_8=2 Pi/R_8.
  sigma=sqrt(Power_Spectrum(2.0d0*Pi*Little_H_0()/radiusNormalization)*9.0d0/4.0d0/Pi/(radiusNormalization/Little_H_0())**3)
  call Assert('P(k) normalization',sigma,sigma_8(),relTol=1.0d-6)
  ! Test that power spectrum scales as expected.
  ratio=Power_Spectrum(1.0d0)/Power_Spectrum(0.1d0)
  call Assert('P(k) scaling',ratio,0.1d0,relTol=1.0d-6)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
  ! Close the parameter file.
  call Input_Parameters_File_Close
end program Tests_Power_Spectrum
