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

!% Contains a module which extends the standard implementation of basic component to track the
!% Bertschinger mass.

module Node_Component_Basic_Standard_Extended
  !% Extends the standard implementation of basic component to track the Bertschinger mass.
  use Virial_Density_Contrast
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Basic_Standard_Extended_Initialize, Node_Component_Basic_Standard_Extended_Node_Merger , &
       &    Node_Component_Basic_Standard_Extended_Promote   , Node_Component_Basic_Standard_Extended_Rate_Compute, &
       &    Node_Component_Basic_Standard_Extended_Scale_Set

  !# <component>
  !#  <class>basic</class>
  !#  <name>standardExtended</name>
  !#  <extends>
  !#   <class>basic</class>
  !#   <name>standard</name>
  !#  </extends>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massBertschinger</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" /> 
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>accretionRateBertschinger</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>radiusTurnaround</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>radiusTurnaroundGrowthRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Options controlling spherical collapse model to use.
  integer                                        :: nodeComponentBasicExtendedSphericalCollapseType
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeLambda         =0
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeDE             =1
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998=2

  ! Initialization state.
  logical                                          :: moduleInitialized                                  =.false.

  ! Virial density contrast object.
  logical                                          :: virialDensityContrastInitialized                   =.false.
  class  (virialDensityContrastClass), allocatable :: virialDensityContrast_
  !$omp threadprivate(virialDensityContrast_,virialDensityContrastInitialized)


contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Standard_Extended_Initialize(node)
    !% Set the mass accretion rate for {\tt node}.
    use Dark_Matter_Profile_Mass_Definitions
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type            (treeNode                  ), intent(inout), pointer :: node
    type            (treeNode                  )               , pointer :: thisNode              , parentNode    , &
         &                                                                  childNode
    class           (nodeComponentBasic        )               , pointer :: basic                 , childBasic    , &
         &                                                                  thisBasic             , parentBasic
    double precision                                                     :: massUnresolved        , deltaTime     , &
         &                                                                  progenitorMassTotal   , radiusVirial
    type            (varying_string            )                         :: nodeComponentBasicExtendedSphericalCollapseTypeText

    ! Determine spherical collapse model to use.
    if (.not.moduleInitialized) then
       !$omp critical(basicExtendedInitialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>nodeComponentBasicExtendedSphericalCollapseType</name>
          !@   <defaultValue>matterLambda</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The type of spherical collapse model to assume in the extended basic node component class.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('nodeComponentBasicExtendedSphericalCollapseType',nodeComponentBasicExtendedSphericalCollapseTypeText,defaultValue='matterLambda')
          select case (char(nodeComponentBasicExtendedSphericalCollapseTypeText))
          case ('matterLambda'    )
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeLambda
          case ('matterDarkEnergy')
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeDE
          case ('bryanNorman')
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998
          end select
          moduleInitialized=.true.
       end if
       !$omp end critical(basicExtendedInitialize)
    end if
    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard extended class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Get required objects.
       if (.not.virialDensityContrastInitialized) then
          select case (nodeComponentBasicExtendedSphericalCollapseType)
          case (nodeComponentBasicExtendedSphericalCollapseTypeLambda         )
             allocate(virialDensityContrastSphericalCollapseMatterLambda :: virialDensityContrast_)
             select type (virialDensityContrast_)
             type is (virialDensityContrastSphericalCollapseMatterLambda)
                virialDensityContrast_=virialDensityContrastSphericalCollapseMatterLambda()
             end select
          case (nodeComponentBasicExtendedSphericalCollapseTypeDE             )
             allocate(virialDensityContrastSphericalCollapseMatterDE     :: virialDensityContrast_)
             select type (virialDensityContrast_)
             type is (virialDensityContrastSphericalCollapseMatterDE    )
                virialDensityContrast_=virialDensityContrastSphericalCollapseMatterDE    ()
             end select
          case (nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998)
             allocate(virialDensityContrastBryanNorman1998                :: virialDensityContrast_)
             select type (virialDensityContrast_)
             type is (virialDensityContrastBryanNorman1998              )
                virialDensityContrast_=virialDensityContrastBryanNorman1998              ()
             end select
          end select
          virialDensityContrastInitialized=.true.
       end if
       ! Ensure Bertschinger mass is set for all nodes in this tree.
       if (basic%massBertschinger() <= 0.0d0) then
          thisNode => node
          do while (associated(thisNode%parent))
             thisNode => thisNode%parent
          end do
          do while (associated(thisNode))
             thisBasic => thisNode%basic()
             ! Set the Bertschinger mass.
             if (thisBasic%massBertschinger() <= 0.0d0) then
                call thisBasic%massBertschingerSet(               &
                     &  Dark_Matter_Profile_Mass_Definition(      &
                     &   thisNode                               , &
                     &   virialDensityContrast_%densityContrast(  &
                     &    thisBasic%mass()                      , &
                     &    thisBasic%time()                        &
                     &                                         ), &
                     &   radius=radiusVirial                      &
                     &                                     )      &
                     &                                   )
                call thisBasic%radiusTurnaroundSet(virialDensityContrast_%turnAroundOverVirialRadii(time=thisBasic%time())*radiusVirial)
             end if
             ! Move to next node.
             call thisNode%walkTree(thisNode)
          end do
       end if
        ! Determine if this node has a descendent.
       if (.not.associated(node%parent)) then
          ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
          ! progenitor, if it has one.
          childNode => node%firstChild
          if (associated(childNode)) then
             ! Get the basic component of the child node.
             childBasic => childNode%basic()
             ! Ensure the child has a mass growth rate computed.
             call Node_Component_Basic_Standard_Extended_Initialize(childNode)
             ! Get the growth rate of the child.
             call basic% accretionRateBertschingerSet(childBasic%             accretionRate())
             call basic%radiusTurnaroundGrowthRateSet(childBasic%radiusTurnaroundGrowthRate())
          else
             ! Parentless node has no child - set a zero growth rate.
             call basic%accretionRateBertschingerSet (0.0d0                                  )
             call basic%radiusTurnaroundGrowthRateSet(0.0d0                                  )
          end if
       else
          ! Get the parent node.
          parentNode => node%parent
          ! Get the basic component of the parent node.
          parentBasic => parentNode%basic()
          ! Compute the unresolved mass.
          massUnresolved=Node_Component_Basic_Standard_Extended_Unresolved_Mass(parentNode)
          if (massUnresolved > 0.0d0) then
             ! Positive mass growth - assume this occurs entirely in the main progenitor.
             if (node%isPrimaryProgenitor()) then
                ! Main progenitor - compute required growth rate.
                deltaTime=parentBasic%time()-basic%time()
                if (deltaTime > 0.0d0) then
                   call basic% accretionRateBertschingerSet(massUnresolved/deltaTime)
                   call basic%radiusTurnaroundGrowthRateSet((parentBasic%radiusTurnaround()-basic%radiusTurnaround())/deltaTime)
                end if
             else
                ! Non-main progenitor - assume zero growth rate.
                call basic%             accretionRateSet(0.0d0)
                call basic%radiusTurnaroundGrowthRateSet(0.0d0)
             end if
          else
             ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
             ! Compute the total mass in progenitors.
             progenitorMassTotal=parentBasic%massBertschinger()-massUnresolved
             ! Compute the time available for accretion.
             deltaTime=parentBasic%time()-basic%time()
             ! Compute mass growth rate.
             if (deltaTime > 0.0d0) then
                call basic% accretionRateBertschingerSet((massUnresolved/deltaTime)*(basic%massBertschinger()/progenitorMassTotal))
                call basic%radiusTurnaroundGrowthRateSet((parentBasic%radiusTurnaround()-basic%radiusTurnaround())/deltaTime)
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Initialize

  double precision function Node_Component_Basic_Standard_Extended_Unresolved_Mass(node)
    !% Return the unresolved mass for {\tt node}.
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: child
    class(nodeComponentBasic)               , pointer :: childBasic, basic

    ! Get the basic component.
    basic => node%basic()
    ! Initialize the unresolved mass to the mass of the current node's basic component.
    Node_Component_Basic_Standard_Extended_Unresolved_Mass=basic%massBertschinger()
    ! Remove the mass of all child nodes.
    child => node%firstChild
    do while (associated(child))
       childBasic                                             => child%basic()
       Node_Component_Basic_Standard_Extended_Unresolved_Mass =  Node_Component_Basic_Standard_Extended_Unresolved_Mass-childBasic%massBertschinger()
       child                                                  => child%sibling
    end do
    return
  end function Node_Component_Basic_Standard_Extended_Unresolved_Mass


  !# <rateComputeTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Basic_Standard_Extended_Rate_Compute(node,interrupt,interruptProcedure)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    implicit none
    type     (treeNode          ), intent(inout), pointer :: node
    logical                      , intent(inout)          :: interrupt
    procedure(                  ), intent(inout), pointer :: interruptProcedure
    class    (nodeComponentBasic)               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Mass rate of change is set to the accretion rate.
       call basic%massBertschingerRate(basic%accretionRateBertschinger ())
       ! Set the radius growth rate.
       call basic%radiusTurnaroundRate(basic%radiusTurnaroundGrowthRate())
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Basic_Standard_Extended_Scale_Set(node)
    !% Set scales for properties in the standard implementation of the basic component.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    double precision                    , parameter              :: scaleLengthRelative=1.0d-6
    double precision                    , parameter              :: scaleMassRelative  =1.0d-6
    class           (nodeComponentBasic)               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Set scale for mass.
       call basic%massBertschingerScale(basic%massBertschinger()*scaleMassRelative  )
       ! Set scale for raadius.
       call basic%radiusTurnaroundScale(basic%radiusTurnaround()*scaleLengthRelative)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Scale_Set

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Basic_Standard_Extended_Node_Merger(node)
    !% Switch off accretion of new mass onto this node once it becomes a satellite.
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    class(nodeComponentBasic)               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Shut down mass accretion onto the halo and growth now that it is a satellite.
       call basic%accretionRateBertschingerSet (0.0d0)
       call basic%radiusTurnAroundGrowthRateSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Node_Merger

   !# <nodePromotionTask>
   !#  <unitName>Node_Component_Basic_Standard_Extended_Promote</unitName>
   !# </nodePromotionTask>
   subroutine Node_Component_Basic_Standard_Extended_Promote(node)
     !% Ensure that {\tt node} is ready for promotion to its parent. In this case, we simply
     !% update the mass of {\tt node} to be that of its parent.
     use Galacticus_Error
     implicit none
     type (treeNode          ), intent(inout), pointer :: node
     type (treeNode          )               , pointer :: parentNode
     class(nodeComponentBasic)               , pointer :: parentBasic, basic

     ! Get the basic component.
     basic => node%basic()
     ! Ensure that it is of the standard class.
     select type (basic)
     class is (nodeComponentBasicStandardExtended)
        ! Get the parent node and its basic component.
        parentNode  => node      %parent
        parentBasic => parentNode%basic ()
        ! Ensure the two halos exist at the same time.
        if (basic%time() /= parentBasic%time())                                               &
             & call Galacticus_Error_Report(                                                  &
             &                              'Node_Component_Basic_Standard_Extended_Promote', &
             &                              'node has not been evolved to its parent'         &
             &                             )
        ! Adjust the mass to that of the parent node.
        call basic%massBertschingerSet          (parentBasic%massBertschinger          ())
        ! Adjust the accretion rate to that of the parent node.
        call basic%accretionRateBertschingerSet (parentBasic%accretionRateBertschinger ())
        ! Adjust the turnaround radius to that of the parent node.
        call basic%radiusTurnaroundSet          (parentBasic%radiusTurnaround          ())
        ! Adjust the accretion rate to that of the parent node.
        call basic%radiusTurnaroundGrowthRateSet(parentBasic%radiusTurnaroundGrowthRate())
     end select
     return
   end subroutine Node_Component_Basic_Standard_Extended_Promote

end module Node_Component_Basic_Standard_Extended
