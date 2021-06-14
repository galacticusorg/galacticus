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

!% Contains a module which implements a preset position component.

module Node_Component_Position_Trace_Dark_Matter
  !% Implements a preset position component.
  use :: Dark_Matter_Halo_Scales       , only : darkMatterHaloScale                       , darkMatterHaloScaleClass
  use :: Satellite_Oprhan_Distributions, only : satelliteOrphanDistributionTraceDarkMatter
  implicit none
  private
  public :: Node_Component_Position_Trace_Dark_Matter_Initialize         , Node_Component_Position_Trace_Dark_Matter_Thread_Initialize, &
       &    Node_Component_Position_Trace_Dark_Matter_Thread_Uninitialize, Node_Component_Position_Trace_Dark_Matter_Update

  !# <component>
  !#  <class>position</class>
  !#  <name>traceDarkMatter</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>position</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Objects used by this component.
  class(darkMatterHaloScaleClass                  ), pointer :: darkMatterHaloScale_
  type (satelliteOrphanDistributionTraceDarkMatter), pointer :: satelliteOrphanDistribution_
  !$omp threadprivate(darkMatterHaloScale_,satelliteOrphanDistribution_)

contains

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Position_Trace_Dark_Matter_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Position_Trace_Dark_Matter_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultPositionComponent%traceDarkMatterIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       allocate(satelliteOrphanDistribution_)
       !# <referenceConstruct object="satelliteOrphanDistribution_" constructor="satelliteOrphanDistributionTraceDarkMatter(darkMatterHaloScale_)"/>
    end if
    return
  end subroutine Node_Component_Position_Trace_Dark_Matter_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Position_Trace_Dark_Matter_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Position_Trace_Dark_Matter_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none

    if (defaultPositionComponent%traceDarkMatterIsActive()) then
       !# <objectDestructor name="darkMatterHaloScale_"        />
       !# <objectDestructor name="satelliteOrphanDistribution_"/>
    end if
    return
  end subroutine Node_Component_Position_Trace_Dark_Matter_Thread_Uninitialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Position_Trace_Dark_Matter_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Position_Trace_Dark_Matter_Initialize(node)
    !% Initialize the position of {\normalfont \ttfamily node}.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentPosition, nodeComponentPositionTraceDarkMatter, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentPosition)               , pointer :: position

    if (.not.defaultPositionComponent%traceDarkMatterIsActive()) return
    position => node%position(autoCreate=.true.)
    select type (position)
    class is (nodeComponentPositionTraceDarkMatter)
       call Node_Component_Position_Trace_Dark_Matter_Assign(position,checkSatelliteStatus=.true.)
    class default
       call Galacticus_Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine Node_Component_Position_Trace_Dark_Matter_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Position_Trace_Dark_Matter_Update</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Position_Trace_Dark_Matter_Update</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Position_Trace_Dark_Matter_Update(node)
    !% Initialize the position of {\normalfont \ttfamily node}.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentPosition, nodeComponentPositionTraceDarkMatter, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentPosition), pointer       :: position

    if (.not.defaultPositionComponent%traceDarkMatterIsActive()) return
    position => node%position()
    select type (position)
    class is (nodeComponentPositionTraceDarkMatter)
       call Node_Component_Position_Trace_Dark_Matter_Assign(position,checkSatelliteStatus=.false.)
    class default
       call Galacticus_Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine Node_Component_Position_Trace_Dark_Matter_Update

  subroutine Node_Component_Position_Trace_Dark_Matter_Assign(self,checkSatelliteStatus)
    !% Assign a position to a satellite.
    use :: Galacticus_Nodes, only : nodeComponentPositionTraceDarkMatter
    implicit none
    class  (nodeComponentPositionTraceDarkMatter), intent(inout) :: self
    logical                                      , intent(in   ) :: checkSatelliteStatus

    if (.not.checkSatelliteStatus .or. self%hostNode%isSatellite()) then
       call self%positionSet(satelliteOrphanDistribution_%position(self%hostNode))
    else
       call self%positionSet([0.0d0,0.0d0,0.0d0])
    end if
    return
  end subroutine Node_Component_Position_Trace_Dark_Matter_Assign

end module Node_Component_Position_Trace_Dark_Matter
