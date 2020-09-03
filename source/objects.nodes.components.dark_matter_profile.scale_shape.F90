!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements a dark matter profile method that provides a scale radius and a shape parameter.

module Node_Component_Dark_Matter_Profile_Scale_Shape
  !% Implements a dark matter profile method that provides a scale radius and a shape parameter.
  use :: Dark_Matter_Profiles_Shape, only : darkMatterProfileShapeClass
  use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfileScaleShape, treeNode
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute, Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set   , Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize     , &
       &    Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Init , Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Uninit

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scaleShape</name>
  !#  <isDefault>false</isDefault>
  !#  <extends>
  !#   <class>darkMatterProfile</class>
  !#   <name>scale</name>
  !#  </extends>
  !#  <properties>
  !#   <property>
  !#     <name>shape</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="Shape parameter of the dark matter profile."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>shapeGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

  ! Objects used by this component.
  class(darkMatterProfileShapeClass), pointer :: darkMatterProfileShape_
  !$omp threadprivate(darkMatterProfileShape_)

  ! Flag indicating whether scale radius and shape data should be output when full merger trees are output.
  logical :: mergerTreeStructureOutputDarkMatterProfileShape

  ! Queriable dark matter profile object.
  type(nodeComponentDarkMatterProfileScaleShape) :: darkMatterProfile

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize(parameters_)
    !% Initializes the ``scale'' implementation of the dark matter halo profile component.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>mergerTreeStructureOutputDarkMatterProfileShape</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Determines whether or not dark matter halo shape parameter is included in outputs of merger trees.</description>
    !#   <source>parameters_</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    ! Bind the shape get function.
    call darkMatterProfile%shapeFunction(Node_Component_Dark_Matter_Profile_Scale_Shape_Shape)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Init</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Init(parameters_)
    !% Initializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultDarkMatterProfileComponent%scaleShapeIsActive()) then
       !# <objectBuilder class="darkMatterProfileShape" name="darkMatterProfileShape_" source="parameters_"/>
       call nodePromotionEvent%attach(defaultDarkMatterProfileComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentDarkMatterProfileScaleShape')
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Init

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Uninit</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Uninit()
    !% Uninitializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scaleShapeIsActive()) then
       !# <objectDestructor name="darkMatterProfileShape_"/>
       call nodePromotionEvent%detach(defaultDarkMatterProfileComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Thread_Uninit

  double precision function Node_Component_Dark_Matter_Profile_Scale_Shape_Shape(self)
    !% Return the shape parameter in the dark matter halo profile.
    implicit none
    class(nodeComponentDarkMatterProfileScaleShape), intent(inout) :: self
    type (treeNode                                ), pointer       :: selfNode

    ! Return the shape parameter, setting it if it has not yet been set.
    if (self%shapeValue() < 0.0d0) then
       ! Get the host halo.
       selfNode                => self                  %host()
       ! Set the shape parameter of the halo.
       call self%shapeSet(darkMatterProfileShape_%shape(selfNode))
    end if
    Node_Component_Dark_Matter_Profile_Scale_Shape_Shape=self%shapeValue()
    return
  end function Node_Component_Dark_Matter_Profile_Scale_Shape_Shape

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the rate of change of the scale radius.
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScaleShape, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type     (treeNode                      ), intent(inout), pointer :: node
    logical                                  , intent(inout)          :: interrupt
    procedure(                              ), intent(inout), pointer :: interruptProcedure
    integer                                  , intent(in   )          :: propertyType
    class    (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultDarkMatterProfileComponent%scaleShapeIsActive()) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale+shape class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScaleShape)
       call darkMatterProfile%shapeRate(darkMatterProfile%shapeGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize(node)
    !% Initialize the scale radius of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfileParent, darkMatterProfile
    class           (nodeComponentBasic            )               , pointer :: basicParent            , basic
    double precision                                                         :: deltaTime              , shape

    if (defaultDarkMatterProfileComponent%scaleShapeIsActive()) then
       ! Get the dark matter profile component, creating it if necessary.
       darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
       ! Get the shape parameter - this will initialize the shape if necessary.
       shape                          = darkMatterProfile%shape()
       ! Check if this node is the primary progenitor.
       if (node%isPrimaryProgenitor()) then
          ! It is, so compute the shape parameter growth rate.
          ! Now compute the growth rate.
          basic       => node       %basic()
          basicParent => node%parent%basic()
          deltaTime=basicParent%time()-basic%time()
          if (deltaTime > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile(autoCreate=.true.)
             call darkMatterProfile%shapeGrowthRateSet(                                               &
                  &                                                 (                                 &
                  &                                                   darkMatterProfileParent%shape() &
                  &                                                  -darkMatterProfile      %shape() &
                  &                                                 )                                 &
                  &                                                 /deltaTime                        &
                  &                                                )
          else
             call darkMatterProfile%shapeGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set shape parameter growth rates to zero.
          call    darkMatterProfile%shapeGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the growth rate of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScaleShape, treeNode
    implicit none
    class(*                             ), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfileParent, darkMatterProfile
    class(nodeComponentBasic            )               , pointer :: basicParent            , basic
    !$GLC attributes unused :: self
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basic                   => node       %basic            ()
    basicParent             => node%parent%basic            ()
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the shape parameter and its growth rate to that of the parent node.
    call darkMatterProfile%shapeSet          (darkMatterProfileParent%shape          ())
    call darkMatterProfile%shapeGrowthRateSet(darkMatterProfileParent%shapeGrowthRate())
    return
  end subroutine nodePromotion

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale+shape class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScale)
          ! Set scale for the scale radius.
       call darkMatterProfile%shapeScale(darkMatterProfile%shape())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale_Shape
