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

!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Node_Component_Dark_Matter_Profile_Scale_Preset
  !% Implements a dark matter profile method that provides a scale radius.
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute , Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set    , Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Init    , &
       &    Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Uninit

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scalePreset</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>scale</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
  !#   </property>
  !#   <property>
  !#     <name>scaleGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Init</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Init(parameters_)
    !% Initializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultDarkMatterProfileComponent%scalePresetIsActive()) &
         & call nodePromotionEvent%attach(defaultDarkMatterProfileComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentDarkMatterProfileScalePreset')
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Init

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Uninit</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Uninit()
    !% Uninitializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scalePresetIsActive()) &
         & call nodePromotionEvent%detach(defaultDarkMatterProfileComponent,nodePromotion)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Thread_Uninit

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the rate of change of the scale radius.
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScalePreset, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type            (treeNode                      ), intent(inout)          :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (                              ), intent(inout), pointer :: interruptProcedure
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultDarkMatterProfileComponent%scalePresetIsActive()) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScalePreset)
       call darkMatterProfile%scaleRate(darkMatterProfile%scaleGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize(node)
    !% Initialize the scale radius of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfileParent, darkMatterProfile
    class           (nodeComponentBasic            )               , pointer :: basicParent            , basic
    double precision                                                         :: deltaTime

    if (defaultDarkMatterProfileComponent%scalePresetIsActive()) then
        ! Get the dark matter profile component.
       darkMatterProfile => node%darkMatterProfile()
       ! Check if this node is the primary progenitor.
       if (node%isPrimaryProgenitor()) then
          ! It is, so compute the scale radius growth rate.
          basic   => node       %basic()
          basicParent => node%parent%basic()
          deltaTime     =basicParent%time()-basic%time()
          if (deltaTime > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile()
             call darkMatterProfile%scaleGrowthRateSet(                                   &
                  &                                     (                                 &
                  &                                       darkMatterProfileParent%scale() &
                  &                                      -darkMatterProfile      %scale() &
                  &                                     )                                 &
                  &                                    /deltaTime                         &
                  &                                   )
          else
             call darkMatterProfile%scaleGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set scale radius growth rate to zero.
          call    darkMatterProfile%scaleGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the growth rate of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentDarkMatterProfile, treeNode
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
    ! Adjust the scale radius to that of the parent node.
    call darkMatterProfile%scaleSet          (darkMatterProfileParent%scale          ())
    ! Adjust the growth rate to that of the parent node.
    call darkMatterProfile%scaleGrowthRateSet(darkMatterProfileParent%scaleGrowthRate())
    return
  end subroutine nodePromotion

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScalePreset, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScalePreset)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(darkMatterProfile%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale_Preset
