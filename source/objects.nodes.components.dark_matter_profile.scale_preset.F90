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

!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Node_Component_Dark_Matter_Profile_Scale_Preset
  !% Implements a dark matter profile method that provides a scale radius.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute, Node_Component_Dark_Matter_Profile_Scale_Preset_Promote        , &
       &    Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set   , Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scalePreset</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>scale</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
  !#   </property>
  !#   <property>
  !#     <name>scaleGrowthRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the scale radius.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    logical                                         , intent(inout)          :: interrupt
    procedure       (                              ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => thisNode%darkMatterProfile()
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
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    class           (nodeComponentDarkMatterProfile)               , pointer :: parentDarkMatterProfile, thisDarkMatterProfile
    class           (nodeComponentBasic            )               , pointer :: parentBasic            , thisBasic
    double precision                                                         :: deltaTime

    if (defaultDarkMatterProfileComponent%scalePresetIsActive()) then
        ! Get the dark matter profile component.
       thisDarkMatterProfile => thisNode%darkMatterProfile()
       ! Check if this node is the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          ! It is, so compute the scale radius growth rate.
          thisBasic   => thisNode       %basic()
          parentBasic => thisNode%parent%basic()
          deltaTime     =parentBasic%time()-thisBasic%time()
          if (deltaTime > 0.0d0) then
             parentDarkMatterProfile => thisNode%parent%darkMatterProfile()
             call thisDarkMatterProfile%scaleGrowthRateSet(                                   &
                  &                                         (                                 &
                  &                                           parentDarkMatterProfile%scale() &
                  &                                          -thisDarkMatterProfile  %scale() &
                  &                                         )                                 &
                  &                                        /deltaTime                         &
                  &                                       )
          else
             call thisDarkMatterProfile%scaleGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set scale radius growth rate to zero.
          call    thisDarkMatterProfile%scaleGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Tree_Initialize

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the growth rate of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile)               , pointer :: parentDarkMatterProfile, thisDarkMatterProfile
    class(nodeComponentBasic            )               , pointer :: parentBasic            , thisBasic

    ! Get the dark matter profile component.
    thisDarkMatterProfile => thisNode%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (thisDarkMatterProfile)
       class is (nodeComponentDarkMatterProfileScalePreset)
       parentDarkMatterProfile => thisNode%parent%darkMatterProfile()
       thisBasic               => thisNode       %basic            ()
       parentBasic             => thisNode%parent%basic            ()
       if (thisBasic%time() /= parentBasic%time()) call Galacticus_Error_Report('Node_Component_Dark_Matter_Profile_Scale_Preset_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the scale radius to that of the parent node.
       call thisDarkMatterProfile%scaleSet          (parentDarkMatterProfile%scale          ())
       ! Adjust the growth rate to that of the parent node.
       call thisDarkMatterProfile%scaleGrowthRateSet(parentDarkMatterProfile%scaleGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Promote

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => thisNode%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScalePreset)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(darkMatterProfile%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Preset_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale_Preset
