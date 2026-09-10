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
Contains a module which implements a position component in Cartesian coordinates.
!!}

module Node_Component_Position_Cartesian
  !!{RST
  Implements a position component in Cartesian coordinates.
  !!}
  implicit none
  private
  public :: Node_Component_Position_Cartesian_Thread_Initialize, Node_Component_Position_Cartesian_Thread_Uninitialize
  ! The `position` and `velocity` properties below are rank-1, 3-element arrays. As with all bare 3-element
  ! position/velocity arrays in Galacticus, their components are *Cartesian*, in the order (x,y,z) --- the
  ! coordinate system is implicit in the representation and is not carried by the type. Code which needs a
  ! system-aware position should assign these into a `coordinateCartesian` object, which then converts to
  ! other systems automatically on assignment. See the convention documented in the `Coordinates` module.
  !![
  <component docformat="rst">
   <class>position</class>
   <name>cartesian</name>
   <description>
   The three-dimensional position and velocity of a node in Cartesian coordinates, optionally accompanied by a history of its position in six-dimensional phase space---most often used for satellite nodes. Positions and velocities do not evolve for a given node: they are assigned during merger tree construction or by node operators (see :galacticus-class:`nodeOperatorPositionDiscrete` and :galacticus-class:`nodeOperatorPositionTraceDarkMatter`). When output, if a phase space history is available the entry closest to the output time is used, since no simple interpolation can correctly capture satellite orbital dynamics. See :ref:`manual-sec-GalacticusVelocityDefinitions` for important notes on velocity definitions.
   </description>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>position</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" unitsDescription="Mpc" unitsQuantity="Mpc" comment="Position of the node, in Cartesian components (x,y,z) (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="kilo" unitsDescription="km/s" unitsQuantity="km/s" comment="Velocity of the node, in Cartesian components (vx,vy,vz) (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>positionHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentThreadInitializationTask function="Node_Component_Position_Cartesian_Thread_Initialize"/>
  !!]
  subroutine Node_Component_Position_Cartesian_Thread_Initialize(parameters_)
    !!{RST
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent      , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultPositionComponent%cartesianIsActive()) &
         call nodePromotionEvent%attach(thread,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentPositionCartesian')
    return
  end subroutine Node_Component_Position_Cartesian_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask function="Node_Component_Position_Cartesian_Thread_Uninitialize"/>
  !!]
  subroutine Node_Component_Position_Cartesian_Thread_Uninitialize()
    !!{RST
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none

    if (defaultPositionComponent%cartesianIsActive() .and. nodePromotionEvent%isAttached(thread,nodePromotion)) &
         & call nodePromotionEvent%detach(thread,nodePromotion)
    return
  end subroutine Node_Component_Position_Cartesian_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !!{RST
    Ensure that ``node`` is ready for promotion to its parent. In this case, update the position of ``node`` to that of the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentPositionCartesian, treeNode
    implicit none
    class(*                    ), intent(inout)          :: self
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentPosition)               , pointer :: positionParent, position
    !$GLC attributes unused :: self
    
    position       => node       %position()
    positionParent => node%parent%position()
    select type (positionParent)
    class is (nodeComponentPositionCartesian)
       call position%       positionSet(positionParent%position       ())
       call position%       velocitySet(positionParent%velocity       ())
       call position%positionHistorySet(positionParent%positionHistory())
    end select
    return
  end subroutine nodePromotion

end module Node_Component_Position_Cartesian
