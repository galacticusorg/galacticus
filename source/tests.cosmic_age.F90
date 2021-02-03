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

program Tests_Cosmic_Age
  !% Tests cosmic age calculations for various Universers. Ages calculated using Python
  !% \href{http://www.astro.ucla.edu/~wright/CC.python}{implementation} of Ned Wright's cosmology calculator.
  use :: Cosmology_Functions     , only : cosmologyFunctions       , cosmologyFunctionsClass, cosmologyFunctionsMatterDarkEnergy, cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters    , only : cosmologyParametersSimple, hubbleUnitsTime
  use :: Display                 , only : displayVerbositySet      , verbosityLevelStandard
  use :: ISO_Varying_String      , only : assignment(=)            , varying_string
  use :: Input_Parameters        , only : inputParameters
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                   , Unit_Tests_Begin_Group , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  double precision                                    , dimension(8), parameter :: redshift                                      =[0.0000000000d+0,1.0000000000d+0,3.0000000000d+0,9.0000000000d+0,3.0000000000d+1,1.0000000000d+2,3.0000000000d+2,1.0000000000d+3]
  double precision                                    , dimension(8), parameter :: ageEdS                                        =[0.6518682071d+0,0.2304700578d+0,0.0814833670d+0,0.0206137656d+0,0.0037766710d+0,0.0006421713d+0,0.0001248044d+0,0.0000205705d+0]
  double precision                                    , dimension(8), parameter :: ageOpen                                       =[0.0790841462d+0,0.0327062977d+0,0.0128686687d+0,0.0035287732d+0,0.0006745641d+0,0.0001164482d+0,0.0000227374d+0,0.0000037552d+0]
  double precision                                    , dimension(8), parameter :: ageCosmologicalConstant                       =[0.0942699818d+0,0.0402619685d+0,0.0147878493d+0,0.0037621019d+0,0.0006895264d+0,0.0001172509d+0,0.0000227901d+0,0.0000037578d+0]
  double precision                                    , dimension(8), parameter :: ageClosed                                     =[3.4369560000d-2,8.6126540000d-3,2.7752830000d-3,6.7037040000d-4,1.2048780000d-4,2.0363050000d-5,3.9509380000d-6,6.5106720000d-7]
  class           (cosmologyFunctionsClass           ), pointer                 :: cosmologyFunctions_
  type            (cosmologyParametersSimple         )                          :: cosmologyParametersClosed                                                                                                                                                       , cosmologyParametersCosmologicalConstant         , &
       &                                                                           cosmologyParametersOpen
  type            (cosmologyFunctionsMatterLambda    )                          :: cosmologyFunctionsCosmologicalConstant                                                                                                                                          , cosmologyFunctionsOpen
  type            (cosmologyFunctionsMatterDarkEnergy)                          :: cosmologyFunctionsDarkEnergyClosed                                                                                                                                              , cosmologyFunctionsDarkEnergyCosmologicalConstant, &
       &                                                                           cosmologyFunctionsDarkEnergyOmegaMinusOneThird
  type            (varying_string                    )                          :: parameterFile
  character       (len=1024                          )                          :: message
  type            (inputParameters                   )                          :: parameters
  integer                                                                       :: iExpansion
  double precision                                                              :: age                                                                                                                                                                             , expansionFactor                                 , &
       &                                                                           expansionFactorSymmetric                                                                                                                                                        , timeTurnaround

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cosmic age")
  ! Cosmology functions for in an Einstein-de Sitter universe. For this case, we use the default settings.
  parameterFile='testSuite/parameters/cosmicAge/EdS.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  cosmologyFunctions_ => cosmologyFunctions()
  ! Define other cosmological parameters.
  cosmologyParametersOpen                =cosmologyParametersSimple( 0.3d0,0.0d0,0.0d0,2.78d0,10000.0d0)
  cosmologyParametersCosmologicalConstant=cosmologyParametersSimple( 0.3d0,0.0d0,0.7d0,2.78d0,10000.0d0)
  cosmologyParametersClosed              =cosmologyParametersSimple(10.0d0,0.0d0,0.0d0,2.78d0,10000.0d0)
  ! Cosmology functions for an open Universe.
  cosmologyFunctionsOpen                          =cosmologyFunctionsMatterLambda    (cosmologyParametersOpen                                   )
  ! Cosmology functions for a cosmological constant Universe.
  cosmologyFunctionsCosmologicalConstant          =cosmologyFunctionsMatterLambda    (cosmologyParametersCosmologicalConstant                   )
  ! Cosmology functions for a cosmological constant Universe using dark energy method.
  cosmologyFunctionsDarkEnergyCosmologicalConstant=cosmologyFunctionsMatterDarkEnergy(cosmologyParametersCosmologicalConstant,-1.0d0      ,0.0d0)
  ! Cosmology functions for a closed Universe using dark energy method.
  cosmologyFunctionsDarkEnergyClosed              =cosmologyFunctionsMatterDarkEnergy(cosmologyParametersClosed              , 0.0d0      ,0.0d0)
  ! Cosmology functions for a closed Universe using dark energy method.
  cosmologyFunctionsDarkEnergyOmegaMinusOneThird  =cosmologyFunctionsMatterDarkEnergy(cosmologyParametersCosmologicalConstant,-1.0d0/3.0d0,0.0d0)
  ! Compute the time of maximum expansion for the Universe. In this simple, OmegaM=10, OmegaDE=0 Universe this is analytically
  ! calculable and equals:
  timeTurnaround=(5.0d0/27.0d0)*Pi/cosmologyParametersClosed%HubbleConstant(hubbleUnitsTime)
  ! Evaluate ages for matter + cosmological constant universes.
  call Unit_Tests_Begin_Group("Matter + Cosmological Constant")
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift       (iExpansion))
     age            =cosmologyFunctions_%cosmicTime                 (expansionFactor            )
     write (message,'(a,f6.1,a)') "cosmic age: EdS                       [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageEdS                 (iExpansion),relTol=1.0d-3)
     age=cosmologyFunctionsOpen                %cosmicTime(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age: open                      [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageOpen                (iExpansion),relTol=1.0d-3)
     age=cosmologyFunctionsCosmologicalConstant%cosmicTime(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age: cosmological constant     [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageCosmologicalConstant(iExpansion),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! Evaluate ages for a dark energy Universe. We consider a few simple cases for which external solutions are available:
  !  * A cosmological constant;
  !  * A closed universe with no dark energy (to test the solving of re-collapsing universes);
  !  * A flat universe with equation of state -1/3 which behaves like the equiavlent open universe.
  call Unit_Tests_Begin_Group("Matter + Dark Energy")
  do iExpansion=1,size(redshift)
     expansionFactor         =cosmologyFunctionsDarkEnergyCosmologicalConstant%expansionFactorFromRedshift(redshift(iExpansion))
     age                     =cosmologyFunctionsDarkEnergyCosmologicalConstant%cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: cosmological constant     [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageCosmologicalConstant(iExpansion),relTol=1.0d-3)
     age                     =cosmologyFunctionsDarkEnergyClosed              %cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: closed                    [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageClosed              (iExpansion),relTol=1.0d-3)
     expansionFactorSymmetric=cosmologyFunctionsDarkEnergyClosed              %expansionFactor(2.0d0*timeTurnaround-age)
     write (message,'(a,f6.1,a)') "cosmic age: closed: collapse symmetry [z=",redshift(iExpansion),"]"
     call Assert(trim(message),expansionFactor,expansionFactorSymmetric,absTol=0.03d0)
     age                     =cosmologyFunctionsDarkEnergyOmegaMinusOneThird  %cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: w^{-1/3}                  [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageOpen                (iExpansion),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Cosmic_Age
