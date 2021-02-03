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

!% Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using an open
!% cosmology.

program Tests_Spherical_Collapse_Dark_Energy_Open
  !% Tests spherical collapse calculations for a dark energy Universe, specifically using an open cosmology. Compares results to
  !% the analytic solution (e.g. \citealt{kitayama_semianalytic_1996}; eqn.~A4).
  use :: Cosmological_Density_Field, only : criticalOverdensity  , criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctions   , cosmologyFunctionsClass
  use :: Display                   , only : displayVerbositySet  , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: ISO_Varying_String        , only : assignment(=)        , varying_string
  use :: Input_Parameters          , only : inputParameters
  use :: Linear_Growth             , only : linearGrowth         , linearGrowthClass
  use :: Numerical_Constants_Math  , only : Pi
  use :: Unit_Tests                , only : Assert               , Unit_Tests_Begin_Group  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                          , dimension(7) :: redshift                   =[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  class           (cosmologyFunctionsClass ), pointer      :: cosmologyFunctions_
  class           (criticalOverdensityClass), pointer      :: criticalOverdensity_
  type            (varying_string          )               :: parameterFile
  character       (len=1024                )               :: message
  integer                                                  :: iExpansion
  double precision                                         :: age                                                                       , criticalOverdensityValue, &
       &                                                      criticalOverdensityExpected                                               , etaf                    , &
       &                                                      expansionFactor
  type            (inputParameters         )               :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (open cosmology)")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/darkEnergy.open.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get the default cosmology functions object.
  cosmologyFunctions_  => cosmologyFunctions ()
  criticalOverdensity_ => criticalOverdensity()
  do iExpansion=1,size(redshift)
     expansionFactor            =cosmologyFunctions_ %expansionFactorFromRedshift(redshift       (iExpansion))
     age                        =cosmologyFunctions_ %cosmicTime                 (expansionFactor            )
     criticalOverdensityValue   =criticalOverdensity_%value                      (age                        )
     etaf                       =acosh(2.0d0/cosmologyFunctions_%omegaMatterEpochal(age)-1.0d0)
     criticalOverdensityExpected=1.5d0*(3.0d0*sinh(etaf)*(sinh(etaf)-etaf)/(cosh(etaf)-1.0d0)**2-2.0d0)*(1.0d0+(2.0d0*Pi/(sinh(etaf)-etaf))**(2.0d0/3.0d0))
     write (message,'(a,f6.1,a,f6.4,a)') "critical density for collapse [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),criticalOverdensityValue,criticalOverdensityExpected,relTol=2.1d-4)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Dark_Energy_Open
