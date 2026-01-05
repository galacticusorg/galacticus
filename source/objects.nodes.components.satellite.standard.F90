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

!!{
Contains a module of satellite orbit tree node methods.
!!}
module Node_Component_Satellite_Standard
  !!{
  Implements the standard satellite component.
  !!}
  implicit none
  private
  public :: Node_Component_Satellite_Standard_Initialize         , Node_Component_Satellite_Standard_Inactive, Node_Component_Satellite_Standard_Scale_Set, Node_Component_Satellite_Standard_Thread_Initialize, &
       &    Node_Component_Satellite_Standard_Thread_Uninitialize

  !![
  <component>
   <class>satellite</class>
   <name>standard</name>
   <extends>
     <class>satellite</class>
     <name>mergeTime</name>
   </extends>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>boundMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>selfBasic%mass()</classDefault>
      <output unitsInSI="massSolar" comment="Bound mass of the node."/>
    </property>
    <property>
      <name>virialOrbit</name>
      <type>keplerOrbit</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

  ! Record of whether satellite bound mass is an inactive variable.
  logical :: inactiveBoundMass

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
    <unitName>Node_Component_Satellite_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Initialize(parameters)
    !!{
    Initializes the standard satellite orbit component module.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatelliteStandard
    use :: Input_Parameters, only : inputParameter                , inputParameters
    implicit none
    type(inputParameters               ), intent(inout) :: parameters
    type(nodeComponentSatelliteStandard)                :: satellite
    type(inputParameters               )                :: subParameters

    if (satellite%standardIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentSatellite')
       !![
       <inputParameter>
         <name>inactiveBoundMass</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not the bound mass variable of the standard satellite component is inactive (i.e. does not appear in any ODE being solved).</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Satellite_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Satellite_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the tree node standard satellite module.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Events_Hooks    , only : nodePromotionEvent       , openMPThreadBindingAtLevel, subhaloPromotionEvent, dependencyExact, &
         &                          dependencyDirectionBefore
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyExact), dimension(1)  :: dependenciesSubhaloPromotion
    !$GLC attributes unused :: parameters

    if (defaultSatelliteComponent%standardIsActive()) then
       dependenciesSubhaloPromotion(1)=dependencyExact(dependencyDirectionBefore,'mergerTreeNodeEvolver')
       call    nodePromotionEvent%attach(thread,   nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSatelliteStandard'                                          )
       call subhaloPromotionEvent%attach(thread,subhaloPromotion,openMPThreadBindingAtLevel,label='nodeComponentSatelliteStandard',dependencies=dependenciesSubhaloPromotion)
    end if
    return
  end subroutine Node_Component_Satellite_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Satellite_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the tree node standard satellite module.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Events_Hooks    , only : nodePromotionEvent       , subhaloPromotionEvent
    implicit none
    
    if (defaultSatelliteComponent%standardIsActive()) then
       if (   nodePromotionEvent%isAttached(thread,   nodePromotion)) call    nodePromotionEvent%detach(thread,   nodePromotion)
       if (subhaloPromotionEvent%isAttached(thread,subhaloPromotion)) call subhaloPromotionEvent%detach(thread,subhaloPromotion)
    end if
    return
  end subroutine Node_Component_Satellite_Standard_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply copy any preexisting satellite orbit
    from the parent.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentSatellite, nodeComponentSatelliteStandard
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite, satelliteParent
    !$GLC attributes unused :: self

    satelliteParent => node%parent%satellite()
    select type (satelliteParent)
    type is (nodeComponentSatelliteStandard)
       satellite => node%satellite()
       select type (satellite)
       type is (nodeComponentSatellite)
          ! This is as expected - nothing to do.
       class default
          call Error_Report('multiple satellite components defined on branch'//{introspection:location})
       end select
       call node%parent%satelliteMove(node,overwrite=.true.)
    end select
    return
  end subroutine nodePromotion

  subroutine subhaloPromotion(self,node,nodePromotion)
    !!{
    Remove the satellite component from the subhalo about to be promoted to an isolated halo (which should have no satellite component).    
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node, nodePromotion
     !$GLC attributes unused :: self, nodePromotion
    
    call node%satelliteRemove(1)
    return
  end subroutine subhaloPromotion
  
  !![
  <inactiveSetTask>
    <unitName>Node_Component_Satellite_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentSatelliteStandard, treeNode
    implicit none
    type (treeNode              ), intent(inout), pointer :: node
    class(nodeComponentSatellite)               , pointer :: satellite
    
    ! Get the satellite component.
    satellite => node%satellite()
    ! Check if an standard satellite component exists.
    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       if (inactiveBoundMass) call satellite%boundMassInactive()
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Inactive
  
  !![
  <scaleSetTask>
    <unitName>Node_Component_Satellite_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, nodeComponentSatelliteStandard, treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentSatellite)               , pointer :: satellite
    class           (nodeComponentBasic    )               , pointer :: basic
    double precision                        , parameter              :: massScaleFractional=1.0d-6

    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the standard class.
    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       ! Get the basic component.
       basic => node%basic()
       ! Set scale for bound mass.
       call satellite%boundMassScale(massScaleFractional*basic%mass())
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Scale_Set
 
end module Node_Component_Satellite_Standard
