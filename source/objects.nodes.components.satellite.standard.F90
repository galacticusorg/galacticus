!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module of satellite orbit tree node methods.
module Node_Component_Satellite_Standard
  !% Implements the standard satellite component.
  use Galacticus_Nodes
  use Kepler_Orbits
  implicit none
  private
  public :: Node_Component_Satellite_Standard_Scale_Set          , Node_Component_Satellite_Standard_Create         , &
       &    Node_Component_Satellite_Standard_Rate_Compute       , Node_Component_Satellite_Standard_Initialize     , &
       &    Node_Component_Satellite_Standard_Halo_Formation_Task, Node_Component_Satellite_Standard_Tree_Initialize

  !# <component>
  !#  <class>satellite</class>
  !#  <name>standard</name>
  !#  <isDefault>yes</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mergeTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Standard_Merge_Time</getFunction>
  !#     <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
  !#   </property>
  !#   <property>
  !#     <name>timeOfMerging</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Standard_Time_Of_Merging</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>boundMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>selfBasicComponent%mass()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Bound mass of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>virialOrbit</name>
  !#     <type>keplerOrbit</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="set:get" />
  !#     <output condition="[[satelliteOutputVirialOrbit]]" comment="Virial orbital parameters of the satellite node."/>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.satellite.standard.bound_functions.inc</functions>
  !# </component>

  ! Record of whether the module has been initialized.
  logical            :: moduleInitialized                   =.false.

  ! Option indicating whether or not satellite virial orbital parameters will be stored.
  logical            :: satelliteOrbitStoreOrbitalParameters

  ! Option indicating whether or not to reset satellite orbits on halo formation events.
  logical            :: satelliteOrbitResetOnHaloFormation

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical, parameter :: acceptUnboundOrbits                 =.false.

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Satellite_Standard_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
   subroutine Node_Component_Satellite_Standard_Initialize()
     !% Initializes the standard satellite orbit component module.
     use Input_Parameters
     implicit none
     type(nodeComponentSatelliteStandard) :: satelliteComponent

     ! Test whether module is already initialize.
     !$omp critical (Node_Component_Satellite_Standard_Initialize)
     if (satelliteComponent%standardIsActive().and..not.moduleInitialized) then
        ! Determine if satellite orbits are to be stored.
        !@ <inputParameter>
        !@   <name>satelliteOrbitStoreOrbitalParameters</name>
        !@   <defaultValue>true</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     Specifies whether satellite virial orbital parameters should be stored (otherwise they are computed
        !@     again---possibly at random---each time they are requested).
        !@   </description>
        !@   <type>boolean</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('satelliteOrbitStoreOrbitalParameters',satelliteOrbitStoreOrbitalParameters,defaultValue=.true.)
        ! Determine if satellite orbits are to be reset on halo formation events.
        !@ <inputParameter>
        !@   <name>satelliteOrbitResetOnHaloFormation</name>
        !@   <defaultValue>false</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     Specifies whether satellite virial orbital parameters should be reset on halo formation events.
        !@   </description>
        !@   <type>boolean</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('satelliteOrbitResetOnHaloFormation',satelliteOrbitResetOnHaloFormation,defaultValue=.false.)
        ! Specify the function to use for setting virial orbits.
        call satelliteComponent%virialOrbitSetFunction(Node_Component_Satellite_Standard_Virial_Orbit_Set)
        call satelliteComponent%virialOrbitFunction   (Node_Component_Satellite_Standard_Virial_Orbit    )
        ! Record that the module is now initialized.
        moduleInitialized=.true.
     end if
     !$omp end critical (Node_Component_Satellite_Standard_Initialize)

     return
   end subroutine Node_Component_Satellite_Standard_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Satellite_Standard_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Satellite_Standard_Tree_Initialize(thisNode)
    !% Initialize the standard satellite component.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (thisNode%isSatellite()) call Node_Component_Satellite_Standard_Create(thisNode)
    return
  end subroutine Node_Component_Satellite_Standard_Tree_Initialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Satellite_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Satellite_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the time until satellite merging rate of change.
    use Dark_Matter_Halos_Mass_Loss_Rates
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    logical                                 , intent(inout)          :: interrupt
    procedure       (                      ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentSatellite)               , pointer :: satelliteComponent
    double precision                                                 :: massLossRate

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Ensure that it is of the standard class.
    select type (satelliteComponent)
    class is (nodeComponentSatelliteStandard)
       if (thisNode%isSatellite()) then
          massLossRate=Dark_Matter_Halos_Mass_Loss_Rate(thisNode)
          call satelliteComponent%mergeTimeRate(-1.0d0      )
          call satelliteComponent%boundMassRate(massLossRate)
       end if
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Rate_Compute

  function Node_Component_Satellite_Standard_Virial_Orbit(self)
    !% Return the orbit of the satellite at the virial radius.
    use Virial_Orbits
    implicit none
    type (keplerOrbit                   )                :: Node_Component_Satellite_Standard_Virial_Orbit
    class(nodeComponentSatelliteStandard), intent(inout) :: self
    type (treeNode                      ), pointer       :: hostNode                                      , selfNode

    selfNode => self%host()
    if (selfNode%isSatellite().or..not.selfNode%isPrimaryProgenitor().and.associated(selfNode%parent)) then
       if (satelliteOrbitStoreOrbitalParameters) then
          Node_Component_Satellite_Standard_Virial_Orbit=self%virialOrbitValue()
       else
          hostNode => selfNode%parent
          Node_Component_Satellite_Standard_Virial_Orbit=Virial_Orbital_Parameters(selfNode,hostNode,acceptUnboundOrbits)
       end if
    else
       call Node_Component_Satellite_Standard_Virial_Orbit%reset()
    end if
    return
  end function Node_Component_Satellite_Standard_Virial_Orbit

  subroutine Node_Component_Satellite_Standard_Virial_Orbit_Set(self,thisOrbit)
    !% Set the orbit of the satellite at the virial radius.
    use Satellite_Merging_Timescales
    implicit none
    class           (nodeComponentSatellite         ), intent(inout) :: self
    type            (keplerOrbit                    ), intent(in   ) :: thisOrbit
    type            (treeNode                       ), pointer       :: selfNode
    class           (satelliteMergingTimescalesClass), pointer       :: satelliteMergingTimescalesDefault
    double precision                                                 :: mergeTime
    type            (keplerOrbit                    )                :: virialOrbit

    select type (self)
    class is (nodeComponentSatelliteStandard)
       ! Ensure the orbit is defined.
       call thisOrbit%assertIsDefined()
       ! Get the node.
       selfNode => self%host()
       ! Update the stored time until merging to reflect the new orbit.
       virialOrbit=thisOrbit
       satelliteMergingTimescalesDefault => satelliteMergingTimescales()
       mergeTime=satelliteMergingTimescalesDefault%timeUntilMerging(selfNode,virialOrbit)
       if (mergeTime >= 0.0d0) call self%mergeTimeSet(mergeTime)
       ! Store the orbit.
       call self%virialOrbitSetValue(thisOrbit)
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Virial_Orbit_Set

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Satellite_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Satellite_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    class           (nodeComponentSatellite)               , pointer :: satelliteComponent
    class           (nodeComponentBasic    )               , pointer :: thisBasicComponent
    double precision                        , parameter              :: massScaleFractional=1.0d-6, timeScale=1.0d-3

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Ensure that it is of the standard class.
    select type (satelliteComponent)
    class is (nodeComponentSatelliteStandard)
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Set scale for time.
       call satelliteComponent%mergeTimeScale(timeScale                                    )
       ! Set scale for bound mass.
       call satelliteComponent%boundMassScale(massScaleFractional*thisBasicComponent%mass())
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Scale_Set

  !# <haloFormationTask>
  !#  <unitName>Node_Component_Satellite_Standard_Halo_Formation_Task</unitName>
  !# </haloFormationTask>
  subroutine Node_Component_Satellite_Standard_Halo_Formation_Task(thisNode)
    !% Reset the orbits of satellite galaxies on halo formation events.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(treeNode)               , pointer :: satelliteNode

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%standardIsActive()) return

    ! Return immediately if orbits are not to be reset.
    if (.not.satelliteOrbitResetOnHaloFormation) return

    ! Loop over all satellites.
    satelliteNode => thisNode%firstSatellite
    do while (associated(satelliteNode))
       ! Create a new orbit for this satellite.
       call Node_Component_Satellite_Standard_Create(satelliteNode)
       satelliteNode => satelliteNode%sibling
    end do

    return
  end subroutine Node_Component_Satellite_Standard_Halo_Formation_Task

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Satellite_Standard_Create</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Satellite_Standard_Create</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Satellite_Standard_Create(thisNode)
    !% Create a satellite orbit component and assign a time until merging and a bound mass equal initially to the total halo mass.
    use Virial_Orbits
    use Satellite_Merging_Timescales
    implicit none
    type            (treeNode                       ), intent(inout), pointer :: thisNode
    type            (treeNode                       )               , pointer :: hostNode
    class           (nodeComponentSatellite         )               , pointer :: satelliteComponent
    class           (nodeComponentBasic             )               , pointer :: basicComponent
    class           (satelliteMergingTimescalesClass)               , pointer :: satelliteMergingTimescalesDefault
    logical                                                                   :: isNewSatellite
    double precision                                                          :: mergeTime
    type            (keplerOrbit                    )                         :: thisOrbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%standardIsActive()) return

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Determine if the satellite component exists already.
    isNewSatellite=.false.
    select type (satelliteComponent)
    type is (nodeComponentSatellite)
       isNewSatellite=.true.
    end select

    ! If this is a new satellite, create the component and set the bound mass.
    if (isNewSatellite) then
       satelliteComponent => thisNode%satellite(autoCreate=.true.)
       select type (satelliteComponent)
       class is (nodeComponentSatelliteStandard)
          ! Set the bound mass of the satellite.
          basicComponent => thisNode%basic()
          call satelliteComponent%boundMassSet(basicComponent%mass())
       end select
    end if

    select type (satelliteComponent)
    class is (nodeComponentSatelliteStandard)
       ! Ensure the module has been initialized.
       call Node_Component_Satellite_Standard_Initialize()
       ! Get an orbit for this satellite.
       hostNode => thisNode%parent
       thisOrbit=Virial_Orbital_Parameters(thisNode,hostNode,acceptUnboundOrbits)
       ! Store the orbit if necessary.
       if (satelliteOrbitStoreOrbitalParameters) call satelliteComponent%virialOrbitSet(thisOrbit)
       ! Compute and store a time until merging.
       satelliteMergingTimescalesDefault => satelliteMergingTimescales()
       mergeTime=satelliteMergingTimescalesDefault%timeUntilMerging(thisNode,thisOrbit)
       if (mergeTime >= 0.0d0) call satelliteComponent%mergeTimeSet(mergeTime)
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Create

end module Node_Component_Satellite_Standard
