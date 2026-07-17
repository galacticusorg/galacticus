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

!!{RST
Contains a module which accumulates and reports statistics on the outcome of the
:cite:t:`johnson_random_2021` dark matter profile scale radius energy model.
!!}

module Johnson2021_Statistics
  !!{RST
  Accumulates and reports statistics on the outcome of the :cite:t:`johnson_random_2021` dark matter profile scale radius
  energy model.

  The energy model fails whenever it computes a positive total energy for a halo, as no bound profile corresponds to a
  positive energy. In that case :galacticus-class:`darkMatterProfileScaleRadiusJohnson2021` falls back to using the scale
  radius of the primary progenitor unchanged. The rate at which this happens depends on both the mass resolution and the
  power spectrum, and so is a bias in precisely the regime in which the model is calibrated. It is therefore counted here
  and reported at the end of the run.

  These counters live in this module, rather than in the
  :galacticus-class:`darkMatterProfileScaleRadiusJohnson2021` class itself, for two reasons. First, they are properties of
  the run as a whole, not of any one object - the class is deep-copied per thread, and the counts must accumulate across all
  copies. Second, and decisively, the report must be made via the *static* ``outputFileClose`` hook, which requires a plain
  module: functionClass implementations are compiled as submodules, whose procedures are not accessible to the generated
  static caller. The dynamic ``outputFileCloseEventGlobal`` hook is not usable for this purpose either, as every scale radius
  object is destroyed (and so detached) before the output file is closed, leaving no hooks attached at the point the event
  fires.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  implicit none
  private
  public :: johnson2021EnergyModelApplied, johnson2021EnergyModelFailed

  ! Counts of the nodes to which the energy model was applied, and of those for which it failed. Deliberately not
  ! threadprivate: they accumulate across all threads, and are updated atomically by the caller.
  integer(c_size_t) :: johnson2021EnergyModelApplied=0_c_size_t, johnson2021EnergyModelFailed=0_c_size_t

contains

  !![
  <outputFileClose function="Johnson2021_Report_Statistics"/>
  !!]
  subroutine Johnson2021_Report_Statistics()
    !!{RST
    Report the rate at which the :cite:t:`johnson_random_2021` scale radius energy model failed with a positive total energy.
    !!}
    use :: Display           , only : displayMessage, displayMagenta, displayReset, verbosityLevelSilent, &
         &                            verbosityLevelStandard
    use :: ISO_Varying_String, only : varying_string, assignment(=), operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type(varying_string) :: message

    ! The energy model was never applied to any node. Report this rather than passing over it in silence: "applied but never
    ! failed" must not look identical to "never applied at all". Note that no rate can be formed here - computing one would
    ! divide by zero, which this build traps.
    if (johnson2021EnergyModelApplied <= 0_c_size_t) then
       message='The Johnson et al. (2021) scale radius energy model was constructed but never applied to any node, so every scale radius came from the fall-back relation.'
       call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message,verbosityLevelSilent)
       return
    end if
    message='The Johnson et al. (2021) scale radius energy model failed (the total energy was positive, so the scale radius of the primary progenitor was used unchanged) for '
    message=message//johnson2021EnergyModelFailed//' of the '//johnson2021EnergyModelApplied//' nodes to which it was applied ('
    message=message//int(1.0d4*dble(johnson2021EnergyModelFailed)/dble(johnson2021EnergyModelApplied))//'e-4). This rate depends on both mass resolution and power spectrum, so treat it as a systematic when calibrating this model.'
    if (johnson2021EnergyModelFailed > 0_c_size_t) then
       ! Failures occurred - report at `verbosityLevelSilent` so that the message is shown at any verbosity. This is a bias in
       ! precisely the regime in which this model is calibrated, so it must not be possible to miss it.
       call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message,verbosityLevelSilent)
    else
       ! No failures - report a single line at the default verbosity, so that a clean run is positively confirmed rather than
       ! merely silent.
       call displayMessage(message,verbosityLevelStandard)
    end if
    return
  end subroutine Johnson2021_Report_Statistics

end module Johnson2021_Statistics
