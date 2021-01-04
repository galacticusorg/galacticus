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

!% Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using a flat,
!% $\omega=-0.6$ cosmology.

program Tests_Spherical_Collapse_Dark_Energy_Omega_Zero_Point_Six
  !% Tests spherical collapse calculations for a dark energy Universe, specifically using a flat, $\omega=-0.6$
  !% cosmology. Compares results to points read from Figure~6 of \cite{horellou_dark_2005} using
  !% \href{http://datathief.org/}{\normalfont \scshape DataThief}.
  use :: Cosmology_Functions    , only : cosmologyFunctions            , cosmologyFunctionsClass
  use :: Events_Hooks           , only : eventsHooksInitialize
  use :: Galacticus_Display     , only : Galacticus_Verbosity_Level_Set, verbosityStandard
  use :: ISO_Varying_String     , only : varying_string                , assignment(=)
  use :: Input_Parameters       , only : inputParameters
  use :: Unit_Tests             , only : Assert                        , Unit_Tests_Begin_Group    , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast, only : virialDensityContrast         , virialDensityContrastClass
  implicit none
  double precision                            , dimension(3) :: redshift                     =[0.00d0,1.00d0,2.00d0]
  double precision                            , dimension(3) :: virialDensityContrastExpected=[390.44d0,241.35d0,208.17d0]
  double precision                            , parameter    :: massDummy                    =1.0d0
  class           (cosmologyFunctionsClass   ), pointer      :: cosmologyFunctions_
  class           (virialDensityContrastClass), pointer      :: virialDensityContrast_
  type            (varying_string            )               :: parameterFile
  character       (len=1024                  )               :: message
  integer                                                    :: iExpansion
  double precision                                           :: age                                                       , expansionFactor, &
       &                                                        virialDensityContrastActual
  type            (inputParameters           )               :: parameters

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (ω=-0.6 cosmology)")

  ! Test spherical collapse in a flat universe.
  parameterFile='testSuite/parameters/sphericalCollapse/darkEnergy.constantEoSminus0.6.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get the default cosmology functions object.
  cosmologyFunctions_    => cosmologyFunctions   ()
  virialDensityContrast_ => virialDensityContrast()
  do iExpansion=1,size(redshift)
     expansionFactor            =cosmologyFunctions_%expansionFactorFromRedshift(redshift       (iExpansion))
     age                        =cosmologyFunctions_%cosmicTime                 (expansionFactor            )
     virialDensityContrastActual=virialDensityContrast_%densityContrast               (massDummy,age              )
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),virialDensityContrastActual,virialDensityContrastExpected(iExpansion),relTol=1.0d-2)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Dark_Energy_Omega_Zero_Point_Six
