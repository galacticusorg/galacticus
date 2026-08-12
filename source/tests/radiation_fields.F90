!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

program Tests_Radiation_Fields
  !!{RST
  Tests the tabulation of rate coefficients---the integral of a radiation field over a cross section---which a radiation field
  depending on time alone builds and extends as it is asked for further times.

  The tests are of the *determinism* of that tabulation: the rate coefficient returned for a given time must not depend on which
  other times have been asked for. Three fields are used, differing only in what they are asked and in what order: one is asked
  for its rate coefficients from the earliest time to the latest, one from the latest to the earliest, and a third is built
  afresh for each time and asked once, so that its tabulation spans only the neighbourhood of that time. All three must agree
  exactly, which requires that the times tabulated lie on an absolute lattice and that the values already computed be carried
  over exactly when the tabulation is extended.

  The cosmic microwave background is used as the radiation field, so the rate coefficient must also decrease monotonically with
  time as the universe expands and the background cools---a check that the values are placed at the times they were computed
  for.
  !!}
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Display                   , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyInitialize  , treeNode
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize    , Node_Components_Thread_Initialize
  use :: Radiation_Fields          , only : crossSectionFunctionTemplate  , radiationFieldCosmicMicrowaveBackground
  use :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Times at which the rate coefficient is evaluated, in Gyr. These span a decade, so that the tabulation must be extended
  ! several times, and include exact points of the lattice on which it is built - a value carried over to the wrong place shows
  ! itself most clearly at a time which needs no interpolation.
  double precision                                         , dimension(5), parameter :: time                =[1.0d0,2.0d0,3.16227766016838d0,6.0d0,10.0d0]
  ! Wavelength range over which to integrate, in Angstroms. This spans the peak of the cosmic microwave background over the
  ! whole range of times used here.
  double precision                                         , dimension(2), parameter :: wavelengthRange     =[1.0d6,1.0d9]
  double precision                                         , dimension(5)            :: rateAscending                                                     , rateDescending     , &
       &                                                                                rateSingle
  type            (cosmologyParametersSimple              )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda         )                          :: cosmologyFunctions_
  type            (radiationFieldCosmicMicrowaveBackground)                          :: radiationAscending                                                , radiationDescending
  ! One field per time, each constructed once and asked a single question - these must not be re-assigned in a loop, as the
  ! assignment would finalize the previous field and so release the cosmology functions object which they all share.
  type            (radiationFieldCosmicMicrowaveBackground), dimension(5)            :: radiationSingle
  procedure       (crossSectionFunctionTemplate           ), pointer                 :: crossSection_
  type            (treeNode                               ), pointer                 :: node
  type            (inputParameters                        )                          :: parameters
  character       (len=1024                               )                          :: message
  integer                                                                            :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Build a node - one must be passed to the radiation field, although a field which depends on time alone can not make any use
  ! of it.
  parameters=inputParameters()
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node          => treeNode()
  crossSection_ => crossSection
  ! Construct the cosmology.
  cosmologyParameters_=cosmologyParametersSimple     (                                           &
       &                                              OmegaMatter         = 0.2815d0           , &
       &                                              OmegaBaryon         = 0.0465d0           , &
       &                                              OmegaDarkEnergy     = 0.7185d0           , &
       &                                              temperatureCMB      = 2.7255d0           , &
       &                                              HubbleConstant      =69.3000d0             &
       &                                             )
  cosmologyFunctions_ =cosmologyFunctionsMatterLambda(                                           &
       &                                              cosmologyParameters_=cosmologyParameters_  &
       &                                             )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Radiation fields")
  ! Evaluate the rate coefficient in each direction of time, and for a field built afresh for each time.
  radiationAscending =radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
  radiationDescending=radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
  do i=1,size(time)
     call radiationAscending %timeSet(time(i))
     rateAscending (i)=radiationAscending %integrateOverCrossSection(wavelengthRange,crossSection_,node)
  end do
  do i=size(time),1,-1
     call radiationDescending%timeSet(time(i))
     rateDescending(i)=radiationDescending%integrateOverCrossSection(wavelengthRange,crossSection_,node)
  end do
  do i=1,size(time)
     radiationSingle(i)=radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
     call radiationSingle(i)%timeSet(time(i))
     rateSingle     (i)=radiationSingle(i)%integrateOverCrossSection(wavelengthRange,crossSection_,node)
  end do
  call Assert("rate coefficient is independent of the order in which times are requested",rateAscending,rateDescending,absTol=0.0d0)
  call Assert("rate coefficient is independent of the range of times tabulated"          ,rateAscending,rateSingle    ,absTol=0.0d0)
  do i=2,size(time)
     write (message,'(a,f6.2,a)') "rate coefficient decreases as the background cools [t = ",time(i)," Gyr]"
     call Assert(trim(message),rateAscending(i) < rateAscending(i-1),.true.)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  call node%destroy()
  deallocate(node)

contains

  double precision function crossSection(wavelength)
    !!{RST
    A cross section, in cm\ :math:`^2`, used to test the tabulation of rate coefficients. The precise form is unimportant - what
    matters is that it varies over the range of wavelengths integrated - so a simple power law is used, normalized to a typical
    atomic cross section at a wavelength of :math:`10^7` Angstroms.
    !!}
    implicit none
    double precision, intent(in   ) :: wavelength

    crossSection=+1.0d-18      &
         &       *(            &
         &         +wavelength &
         &         /1.0d+7     &
         &        )**(-3.0d0)
    return
  end function crossSection

end program Tests_Radiation_Fields
