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
  Implements a node operator class that accumulates the ram pressure history of orbiting satellites.
  !!}

  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForceClass

  !![
  <nodeOperator name="nodeOperatorRamPressureHistory" docformat="rst">
   <description>
   A node operator class that accumulates the history of the ram pressure experienced by a satellite along its orbit, as needed
   by models of ram pressure stripping in which the outcome depends on the duration of the encounter as well as on its peak
   strength :cite:`koppen_ram_2018`. Two quantities are tracked in the orbiting satellite component: the time integral of the ram
   pressure, :math:`I`, and its running maximum, :math:`P`. These are evolved as

   .. math::

      \dot{I} = p - \epsilon_\mathrm{decay} I / T_\mathrm{orb},

   and

   .. math::

      \dot{P} = - \epsilon_\mathrm{decay} P / T_\mathrm{orb},

   with :math:`P` reset to :math:`p` after any ODE step at the end of which :math:`p > P`. Here :math:`p` is the instantaneous ram
   pressure, given by the :galacticus-class:`hotHaloRamPressureForceClass`, :math:`T_\mathrm{orb}` is the satellite's orbital
   period, and :math:`\epsilon_\mathrm{decay}=`\ ``[efficiencyDecay]``. The decay terms make both quantities refer, effectively, to
   the most recent pericentric passage rather than to the entire orbital history, so that a disk which regrows between passages
   can be stripped again. Setting ``[efficiencyDecay]`` to zero instead accumulates both quantities over the whole time since
   infall. Both are reset to zero when the satellite's virial orbit is set at infall.

   The instantaneous ram pressure must be evaluated at the satellite's current position, so a
   :galacticus-class:`hotHaloRamPressureForceClass` such as :galacticus-class:`hotHaloRamPressureForceOrbitalPosition` should be
   used here. One which returns a value fixed along the orbit, such as :galacticus-class:`hotHaloRamPressureForceFont2008` (which
   evaluates the ram pressure at pericenter), would give meaningless results. This operator requires the ``orbiting`` satellite
   component.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorRamPressureHistory
     !!{RST
     A node operator class that accumulates the ram pressure history of orbiting satellites.
     !!}
     private
     class           (hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_ => null()
     class           (darkMatterHaloScaleClass    ), pointer :: darkMatterHaloScale_     => null()
     double precision                                        :: efficiencyDecay
   contains
     final     ::                                  ramPressureHistoryDestructor
     procedure :: differentialEvolution         => ramPressureHistoryDifferentialEvolution
     procedure :: differentialEvolutionPostStep => ramPressureHistoryDifferentialEvolutionPostStep
  end type nodeOperatorRamPressureHistory

  interface nodeOperatorRamPressureHistory
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorRamPressureHistory` node operator class.
     !!}
     module procedure ramPressureHistoryConstructorParameters
     module procedure ramPressureHistoryConstructorInternal
  end interface nodeOperatorRamPressureHistory

contains

  function ramPressureHistoryConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorRamPressureHistory` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorRamPressureHistory)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (hotHaloRamPressureForceClass  ), pointer       :: hotHaloRamPressureForce_
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    double precision                                                :: efficiencyDecay

    !![
    <inputParameter docformat="rst">
      <name>efficiencyDecay</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      Efficiency of the decay of the time integral and running maximum of the ram pressure, in units of the inverse orbital period.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"     name="darkMatterHaloScale_"     source="parameters"/>
    !!]
    self=nodeOperatorRamPressureHistory(efficiencyDecay,hotHaloRamPressureForce_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloRamPressureForce_"/>
    <objectDestructor name="darkMatterHaloScale_"    />
    !!]
    return
  end function ramPressureHistoryConstructorParameters

  function ramPressureHistoryConstructorInternal(efficiencyDecay,hotHaloRamPressureForce_,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorRamPressureHistory` node operator class.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none
    type            (nodeOperatorRamPressureHistory)                        :: self
    class           (hotHaloRamPressureForceClass  ), intent(in   ), target :: hotHaloRamPressureForce_
    class           (darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                , intent(in   )         :: efficiencyDecay
    !![
    <constructorAssign variables="efficiencyDecay, *hotHaloRamPressureForce_, *darkMatterHaloScale_"/>
    !!]

    if (efficiencyDecay < 0.0d0) call Error_Report('[efficiencyDecay] must be non-negative'//{introspection:location})
    if     (                                                                                                                         &
         &  .not.(                                                                                                                   &
         &        defaultSatelliteComponent%ramPressureTimeIntegratedIsGettable()                                                    &
         &   .and.                                                                                                                   &
         &        defaultSatelliteComponent%ramPressureMaximumIsGettable       ()                                                    &
         &  )                                                                                                                        &
         & ) call Error_Report(                                                                                                      &
         &                     'this node operator requires a satellite component that provides the gettable properties '//         &
         &                     '"ramPressureTimeIntegrated" and "ramPressureMaximum"'                                    //char(10)// &
         &                     '   (the "orbiting" satellite component provides these)'                                  //           &
         &                     {introspection:location}                                                                              &
         &                    )
    return
  end function ramPressureHistoryConstructorInternal

  subroutine ramPressureHistoryDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorRamPressureHistory` node operator class.
    !!}
    implicit none
    type(nodeOperatorRamPressureHistory), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloRamPressureForce_"/>
    <objectDestructor name="self%darkMatterHaloScale_"    />
    !!]
    return
  end subroutine ramPressureHistoryDestructor

  subroutine ramPressureHistoryDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Accumulate the time integral of the ram pressure, and decay both it and the running maximum of the ram pressure.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Satellite_Orbits, only : Satellite_Orbit_Period
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorRamPressureHistory), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentSatellite        )               , pointer :: satellite
    double precision                                                         :: periodOrbital
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    satellite => node%satellite()
    ! Do not compute rates at zero radius, where the orbital period is undefined.
    if (Vector_Magnitude(satellite%position()) <= 0.0d0) return
    call satellite%ramPressureTimeIntegratedRate(+self%hotHaloRamPressureForce_%force(node))
    if (self%efficiencyDecay <= 0.0d0) return
    periodOrbital=Satellite_Orbit_Period(node,self%darkMatterHaloScale_)
    call satellite%ramPressureTimeIntegratedRate(                                       &
         &                                       -self     %efficiencyDecay             &
         &                                       *satellite%ramPressureTimeIntegrated() &
         &                                       /          periodOrbital               &
         &                                      )
    call satellite%ramPressureMaximumRate       (                                       &
         &                                       -self     %efficiencyDecay             &
         &                                       *satellite%ramPressureMaximum       () &
         &                                       /          periodOrbital               &
         &                                      )
    return
  end subroutine ramPressureHistoryDifferentialEvolution

  subroutine ramPressureHistoryDifferentialEvolutionPostStep(self,node,status)
    !!{RST
    Update the running maximum of the ram pressure after each accepted ODE step.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Interface_GSL   , only : GSL_Success           , GSL_Continue
    implicit none
    class           (nodeOperatorRamPressureHistory), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    integer                                         , intent(inout) :: status
    class           (nodeComponentSatellite        ), pointer       :: satellite
    double precision                                                :: ramPressure

    if (.not.node%isSatellite()) return
    satellite   => node                         %satellite(    )
    ramPressure =  self%hotHaloRamPressureForce_%force    (node)
    if (ramPressure > satellite%ramPressureMaximum()) then
       call satellite%ramPressureMaximumSet(ramPressure)
       ! Signal that the node state was changed, without requiring the ODE evolution to be reset.
       if (status == GSL_Success) status=GSL_Continue
    end if
    return
  end subroutine ramPressureHistoryDifferentialEvolutionPostStep
