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

program Tests_Sigma
  !% Tests
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Power_Spectra
  use Numerical_Ranges
  use Cosmology_Parameters
  use Numerical_Constants_Math
  implicit none
  type            (varying_string          )                       :: parameterFile
  integer                                   , parameter            :: massCount              =10
  double precision                          , parameter            :: massMaximum            =1.0d15, massMinimum  =1.0d6
  double precision                          , dimension(massCount) :: mass                          , massFromSigma      , &
       &                                                              sigma
  class           (cosmologyParametersClass), pointer              :: thisCosmologyParameters
  integer                                                          :: iMass
  double precision                                                 :: mass8                         , radius8            , &
       &                                                              sigma8

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.sigma.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Power spectrum: σ(M)")

  ! Open the parameter file.
  parameterFile='parameters.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Create an array of masses.
  mass=Make_Range(massMinimum,massMaximum,massCount,rangeType=rangeTypeLogarithmic)

  ! Check that converting from mass to sigma and back to mass gives consistent answers.
  do iMass=1,massCount
     sigma        (iMass)=Cosmological_Mass_Root_Variance    (mass (iMass))
     massFromSigma(iMass)=Mass_from_Cosmolgical_Root_Variance(sigma(iMass))
  end do
  call Assert('M -> σ(M) -> M conversion loop',mass,massFromSigma,relTol=1.0d-3)

  ! Get the default cosmology.
  thisCosmologyParameters => cosmologyParameters()

  ! Compute the mass corresponding to 8Mpc/h.
  radius8=8.0d0/thisCosmologyParameters%HubbleConstant(unitsLittleH)
  mass8=4.0d0*Pi*thisCosmologyParameters%densityCritical()*thisCosmologyParameters%OmegaMatter()*radius8**3/3.0d0
  sigma8=Cosmological_Mass_Root_Variance(mass8)
  call Assert('σ₈ equals specified value',sigma8,sigma_8(),relTol=1.0d-6)

  ! Close the input parameter file.
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Sigma
