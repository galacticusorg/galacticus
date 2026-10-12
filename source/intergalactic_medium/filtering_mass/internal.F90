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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a filtering mass class which returns the filtering mass solved for along with an internal model of the :term:`IGM`
  state.
  !!}

  use :: Intergalactic_Medium_State, only : intergalacticMediumState, intergalacticMediumStateClass, intergalacticMediumStateInternal

  !![
  <intergalacticMediumFilteringMass name="intergalacticMediumFilteringMassInternal" docformat="rst">
   <description>
   A filtering mass class which returns the filtering mass solved for self-consistently with the state of the :term:`IGM` by
   the :galacticus-class:`universeOperatorIntergalacticMediumStateEvolve` universe operator. The :term:`IGM` state must be
   :galacticus-class:`intergalacticMediumStateInternal`. The filtering mass is interpolated linearly in time between the
   times at which it has been solved for. Since the state is solved for one step at a time as the universe evolves, the
   filtering mass may be requested up to one step beyond the latest time solved for, over which it is extrapolated linearly. The
   fraction of baryons accreted is given by the :cite:t:`gnedin_effect_2000` relation.

   Since the state is not known until it has been solved for, this class can not be used by classes which require the filtering
   mass at all times before the universe is evolved (such as
   :galacticus-class:`criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy`).
   </description>
  </intergalacticMediumFilteringMass>
  !!]
  type, extends(intergalacticMediumFilteringMassClass) :: intergalacticMediumFilteringMassInternal
     !!{RST
     A filtering mass class which returns the filtering mass solved for along with an internal model of the :term:`IGM` state.
     !!}
     private
     class(intergalacticMediumStateClass), pointer :: intergalacticMediumState_ => null()
   contains
     final     ::                              internalDestructor
     procedure :: massFiltering             => internalMassFiltering
     procedure :: massFilteringRateOfChange => internalMassFilteringRateOfChange
  end type intergalacticMediumFilteringMassInternal

  interface intergalacticMediumFilteringMassInternal
     !!{RST
     Constructors for the :galacticus-class:`intergalacticMediumFilteringMassInternal` filtering mass class.
     !!}
     module procedure internalConstructorParameters
     module procedure internalConstructorInternal
  end interface intergalacticMediumFilteringMassInternal

contains

  function internalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`intergalacticMediumFilteringMassInternal` filtering mass class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (intergalacticMediumFilteringMassInternal)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(intergalacticMediumStateClass           ), pointer       :: intergalacticMediumState_

    !![
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    !!]
    self=intergalacticMediumFilteringMassInternal(intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="intergalacticMediumState_"/>
    !!]
    return
  end function internalConstructorParameters

  function internalConstructorInternal(intergalacticMediumState_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`intergalacticMediumFilteringMassInternal` filtering mass class.
    !!}
    use :: Display           , only : displayGreen, displayReset
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    type (intergalacticMediumFilteringMassInternal)                        :: self
    class(intergalacticMediumStateClass           ), intent(in   ), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="*intergalacticMediumState_"/>
    !!]

    select type (intergalacticMediumState_)
    type is (intergalacticMediumStateInternal)
       ! This is the required class.
    class default
       call Error_Report(                                                                                                          &
            &            '[intergalacticMediumFilteringMass]=internal requires [intergalacticMediumState]=internal, but found ['// &
            &            char(self%intergalacticMediumState_%objectType())//']'//char(10)                                       // &
            &            displayGreen()//'HELP:'//displayReset()                                                                // &
            &            ' set [intergalacticMediumState]=internal and [universeOperator]=intergalacticMediumStateEvolve to'    // &
            &            ' solve for the filtering mass along with the IGM state, or choose another'                            // &
            &            ' [intergalacticMediumFilteringMass], such as `gnedin2000`'                                            // &
            &            {introspection:location}                                                                                  &
            &           )
    end select
    return
  end function internalConstructorInternal

  subroutine internalDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`intergalacticMediumFilteringMassInternal` filtering mass class.
    !!}
    implicit none
    type(intergalacticMediumFilteringMassInternal), intent(inout) :: self

    !![
    <objectDestructor name="self%intergalacticMediumState_"/>
    !!]
    return
  end subroutine internalDestructor

  double precision function internalMassFiltering(self,time)
    !!{RST
    Return the filtering mass at the given ``time``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (intergalacticMediumFilteringMassInternal), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    select type (intergalacticMediumState_ => self%intergalacticMediumState_)
    type is (intergalacticMediumStateInternal)
       internalMassFiltering=intergalacticMediumState_%massFiltering(time)
    class default
       internalMassFiltering=0.0d0
       call Error_Report('[intergalacticMediumState] must be `internal`'//{introspection:location})
    end select
    return
  end function internalMassFiltering

  double precision function internalMassFilteringRateOfChange(self,time)
    !!{RST
    Return the rate of change of the filtering mass at the given ``time``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (intergalacticMediumFilteringMassInternal), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    select type (intergalacticMediumState_ => self%intergalacticMediumState_)
    type is (intergalacticMediumStateInternal)
       internalMassFilteringRateOfChange=intergalacticMediumState_%massFilteringRateOfChange(time)
    class default
       internalMassFilteringRateOfChange=0.0d0
       call Error_Report('[intergalacticMediumState] must be `internal`'//{introspection:location})
    end select
    return
  end function internalMassFilteringRateOfChange
