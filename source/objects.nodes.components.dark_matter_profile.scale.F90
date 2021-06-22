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

!!{
Contains a module which implements a dark matter profile method that provides a scale radius.
!!}

module Node_Component_Dark_Matter_Profile_Scale
  !!{
  Implements a dark matter profile method that provides a scale radius.
  !!}
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Scale_Set        , Node_Component_Dark_Matter_Profile_Scale_Plausibility       , &
       &    Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize, Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize

  !![
  <component>
   <class>darkMatterProfile</class>
   <name>scale</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>scale</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
      <classDefault>-1.0d0</classDefault>
    </property>
    <property>
      <name>scaleGrowthRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultDarkMatterProfileComponent%scaleIsActive()) &
         & call nodePromotionEvent%attach(defaultDarkMatterProfileComponent,nodePromotion,openMPThreadBindingAtLevel,label="nodeComponentDarkMatterProfileScale")
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize()
    !!{
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scaleIsActive()) &
         & call nodePromotionEvent%detach(defaultDarkMatterProfileComponent,nodePromotion)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Plausibility</unitName>
   <after>Node_Component_Basic_Standard_Plausibility</after>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility(node)
    !!{
    Determines whether the dark matter profile is physically plausible for radius solving tasks.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type   (treeNode                      ), intent(inout) :: node
    class  (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    ! Return immediately if already non-plausible.
    if (.not.node%isPhysicallyPlausible.and..not.node%isSolvable) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScale)
       if (darkMatterProfile%scale() <= 0.0d0) then
          node%isPhysicallyPlausible=.false.
          node%isSolvable           =.false.
       end if
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the growth rate of {\normalfont \ttfamily node}
    to be that of its parent.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    class(*                             ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfileParent, darkMatterProfile
    class(nodeComponentBasic            ), pointer       :: basicParent            , basic
    !$GLC attributes unused :: self    

    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basic                   => node       %basic            ()
    basicParent             => node%parent%basic            ()
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the scale radius and its growth rate to that of the parent node.
    call darkMatterProfile%scaleSet          (darkMatterProfileParent%scale          ())
    call darkMatterProfile%scaleGrowthRateSet(darkMatterProfileParent%scaleGrowthRate())
    return
  end subroutine nodePromotion
  
  !![
  <scaleSetTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScale)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(darkMatterProfile%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale
