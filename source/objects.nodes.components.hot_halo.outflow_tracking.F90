!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements an extension of the standard hot halo node component which tracks the metals arriving from
!% outflows.

module Node_Component_Hot_Halo_Outflow_Tracking
  !% Implements an extension of the standard hot halo node component which tracks the metals arriving from               
  !% outflows.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute, Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>outflowTracking</name>
  !#  <extends>
  !#   <class>hotHalo</class>
  !#   <name>standard</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>trackedOutflowMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass in the hot phase of the hot halo arrived via direct outflow."/>
  !#   </property>
  !#   <property>
  !#     <name>trackedOutflowAbundances</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the hot phase of the hot halo arrived via direct outflow."/>
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="massRemovalRate" function="Node_Component_Hot_Halo_Outflow_Tracking_Mass_Removal_Rate" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.hot_halo.outflow_tracking.bound_functions.inc</functions>
  !# </component>

contains

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                    )           , intent(inout), pointer :: node
    logical                                                  , intent(in   )          :: odeConverged
    logical                                                  , intent(inout)          :: interrupt
    procedure       (interruptTask               )           , intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo        )                          , pointer :: hotHalo
    class           (darkMatterHaloScaleClass    )                          , pointer :: darkMatterHaloScale_
    double precision                                                                  :: massReturnRate
    type            (abundances                  )                                    :: abundancesReturnRate
    !GCC$ attributes unused :: interrupt, interruptProcedure, odeConverged
    
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Act only if this hot halo is of our class.
    select type (hotHalo)
    class is (nodeComponentHotHaloOutflowTracking)
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       ! Add the rate of abundances return from the outflowed component.
       massReturnRate      =hotHaloOutflowReturnRate*hotHalo%outflowedMass      ()/darkMatterHaloScale_%dynamicalTimescale(node)
       abundancesReturnRate=hotHaloOutflowReturnRate*hotHalo%outflowedAbundances()/darkMatterHaloScale_%dynamicalTimescale(node)
       call hotHalo%trackedOutflowMassRate      (      massReturnRate)
       call hotHalo%trackedOutflowAbundancesRate(abundancesReturnRate)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set(node)
    !% Set scales for properties of {\tt node}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                      , parameter              :: scaleMassRelative=1.0d-3
    double precision                                               :: massVirial

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the standard class.
    select type (hotHalo)
    class is (nodeComponentHotHaloOutflowTracking)
       ! Get the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial=basic%mass()
       ! Set a scale for the tracked abundances.
       call hotHalo%trackedOutflowMassScale      (               massVirial*scaleMassRelative)
       call hotHalo%trackedOutflowAbundancesScale(unitAbundances*massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set

end module Node_Component_Hot_Halo_Outflow_Tracking
