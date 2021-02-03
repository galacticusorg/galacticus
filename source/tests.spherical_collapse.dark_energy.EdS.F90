!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using an Einstein-de
!% Sitter cosmology (e.g. \citealt{kitayama_semianalytic_1996}; eqn.~A2).

program Tests_Spherical_Collapse_Dark_Energy_EdS
  !% Tests spherical collapse calculations for a dark energy Universe, specifically using an Einstein-de
  !% Sitter cosmology. Compares results to the analytic solution.
  use :: Cosmological_Density_Field, only : criticalOverdensity  , criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctions   , cosmologyFunctionsClass
  use :: Display                   , only : displayVerbositySet  , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: ISO_Varying_String        , only : assignment(=)        , varying_string
  use :: Input_Parameters          , only : inputParameters
  use :: Numerical_Constants_Math  , only : Pi
  use :: Unit_Tests                , only : Assert               , Unit_Tests_Begin_Group    , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast   , only : virialDensityContrast, virialDensityContrastClass
  implicit none
  double precision                            , dimension(7) :: redshift                   =[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  class           (cosmologyFunctionsClass   ), pointer      :: cosmologyFunctions_
  class           (virialDensityContrastClass), pointer      :: virialDensityContrast_
  class           (criticalOverdensityClass  ), pointer      :: criticalOverdensity_
  double precision                            , parameter    :: massDummy                  =1.0d0
  type            (varying_string            )               :: parameterFile
  character       (len=1024                  )               :: message
  integer                                                    :: iExpansion
  double precision                                           :: age                                                                       , criticalOverdensityValue     , &
       &                                                        criticalOverdensityExpected                                               , expansionFactor              , &
       &                                                        virialDensityContrastActual                                               , virialDensityContrastExpected
  type            (inputParameters           )               :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (Einstein-de Sitter cosmology)")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/darkEnergy.EdS.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get the default cosmology functions object.
  cosmologyFunctions_    => cosmologyFunctions   ()
  virialDensityContrast_ => virialDensityContrast()
  criticalOverdensity_   => criticalOverdensity  ()
  do iExpansion=1,size(redshift)
     expansionFactor            =cosmologyFunctions_ %expansionFactorFromRedshift(redshift       (iExpansion))
     age                        =cosmologyFunctions_ %cosmicTime                 (expansionFactor            )
     criticalOverdensityValue   =criticalOverdensity_%value                      (age                        )
     criticalOverdensityExpected=3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0
     write (message,'(a,f6.1,a,f6.4,a)') "critical density for collapse [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),criticalOverdensityValue,criticalOverdensityExpected,relTol=3.0d-5)
     virialDensityContrastActual  =virialDensityContrast_%densityContrast(massDummy,age)
     virialDensityContrastExpected=18.0d0*Pi**2
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast       [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),virialDensityContrastActual,virialDensityContrastExpected,relTol=3.0d-5)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Dark_Energy_EdS
