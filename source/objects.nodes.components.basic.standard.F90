!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module with the standard implementation of basic tree node methods.

module Node_Component_Basic_Standard
  !% The standard implementation of basic tree node methods.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Basic_Standard_Rate_Compute   , Node_Component_Basic_Standard_Scale_Set     , &
       &    Node_Component_Basic_Standard_Tree_Initialize, Node_Component_Basic_Standard_Stop_Accretion, &
       &    Node_Component_Basic_Standard_Promote

  !# <component>
  !#  <class>basic</class>
  !#  <name>standard</name>
  !#  <isDefault>yes</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Total mass of the node, assuming univeral baryon fraction."/>
  !#   </property>
  !#   <property>
  !#     <name>time</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>timeLastIsolated</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">BasicStandardTimeLastIsolated</getFunction>
  !#     <output unitsInSI="gigaYear" comment="Time at which node was last an isolated halo."/>
  !#   </property>
  !#   <property>
  !#     <name>accretionRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.basic.standard.bound_functions.inc</functions>
  !# </component>

contains

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Basic_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Basic_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    implicit none
    type     (treeNode          ), intent(inout), pointer :: thisNode            
    logical                      , intent(inout)          :: interrupt           
    procedure(                  ), intent(inout), pointer :: interruptProcedure  
    class    (nodeComponentBasic)               , pointer :: basicComponent      
    
    ! Get the basic component.                                                                          
    basicComponent => thisNode%basic()
    ! Ensure that it is of the standard class.
    select type (basicComponent)
    class is (nodeComponentBasicStandard)
       ! Mass rate of change is set to the accretion rate.
       call basicComponent%massRate(basicComponent%accretionRate())
       ! Time rate of change is unity, by definition.
       call basicComponent%timeRate(1.0d0                         )
    end select
    return
  end subroutine Node_Component_Basic_Standard_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Basic_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Basic_Standard_Scale_Set(thisNode)
    !% Set scales for properties in the standard implementation of the basic component.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode                  
    double precision                    , parameter              :: timeScale        =1.0d-3  
    double precision                    , parameter              :: scaleMassRelative=1.0d-6  
    class           (nodeComponentBasic)               , pointer :: basicComponent            
    
    ! Get the basic component.                                                                                       
    basicComponent => thisNode%basic()
    ! Ensure that it is of the standard class.
    select type (basicComponent)
    class is (nodeComponentBasicStandard)
       ! Set scale for time.
       call basicComponent%timeScale(timeScale                              )
       ! Set scale for mass.
       call basicComponent%massScale(basicComponent%mass()*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Basic_Standard_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Standard_Tree_Initialize(thisNode)
    !% Set the mass accretion rate for {\tt thisNode}.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode                                                        
    type            (treeNode          )               , pointer :: childNode          , parentNode                                 
    class           (nodeComponentBasic)               , pointer :: childBasicComponent, parentBasicComponent, thisBasicComponent   
    double precision                                             :: deltaTime          , massUnresolved      , progenitorMassTotal  
    
    ! Get the basic component.                                                                                                                             
    thisBasicComponent => thisNode%basic()
    ! Ensure that it is of the standard class.
    select type (thisBasicComponent)
    class is (nodeComponentBasicStandard)
       ! Set the last isolated time to the current time at the farthest point along the future of this branch.
       parentNode => thisNode
       do while (associated(parentNode%parent).and.parentNode%isPrimaryProgenitor())
          parentNode => parentNode%parent
       end do
       parentBasicComponent => parentNode%basic()
       call thisBasicComponent%timeLastIsolatedSet(parentBasicComponent%time())
       ! Determine if this node has a descendent.
       if (.not.associated(thisNode%parent)) then
          ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
          ! progenitor, if it has one.
          childNode => thisNode%firstChild
          if (associated(childNode)) then
             ! Get the basic component of the child node.
             childBasicComponent => childNode%basic()
             ! Ensure the child has a mass growth rate computed.
             call Node_Component_Basic_Standard_Tree_Initialize(childNode)
             ! Get the growth rate of the child.
             call thisBasicComponent%accretionRateSet(childBasicComponent%accretionRate())
          else
             ! Parentless node has no child - set a zero growth rate.
             call thisBasicComponent%accretionRateSet(0.0d0                              )
          end if
       else
          ! Get the parent node.
          parentNode => thisNode%parent
          ! Get the basic component of the parent node.
          parentBasicComponent => parentNode%basic()
          ! Compute the unresolved mass.
          massUnresolved=Node_Component_Basic_Standard_Unresolved_Mass(parentNode)
          if (massUnresolved > 0.0d0) then
             ! Positive mass growth - assume this occurs entirely in the main progenitor.
             if (thisNode%isPrimaryProgenitor()) then
                ! Main progenitor - compute required growth rate.
                deltaTime=parentBasicComponent%time()-thisBasicComponent%time()
                if (deltaTime > 0.0d0) call thisBasicComponent%accretionRateSet(massUnresolved/deltaTime)
             else
                ! Non-main progenitor - assume zero growth rate.
                call thisBasicComponent%accretionRateSet(0.0d0)
             end if
          else
             ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
             ! Compute the total mass in progenitors.
             progenitorMassTotal=parentBasicComponent%mass()-massUnresolved
             ! Compute the time available for accretion.
             deltaTime=parentBasicComponent%time()-thisBasicComponent%time()
             ! Compute mass growth rate.
             if (deltaTime > 0.0d0) call thisBasicComponent%accretionRateSet((massUnresolved/deltaTime)*(thisBasicComponent%mass()/progenitorMassTotal))
          end if
       end if

    end select
    return
  end subroutine Node_Component_Basic_Standard_Tree_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Basic_Standard_Stop_Accretion</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Basic_Standard_Stop_Accretion(thisNode)
    !% Switch off accretion of new mass onto this node once it becomes a satellite.
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode            
    class(nodeComponentBasic)               , pointer :: thisBasicComponent  
    
    ! Get the basic component.                                                                      
    thisBasicComponent => thisNode%basic()
    ! Ensure that it is of the standard class.
    select type (thisBasicComponent)
    class is (nodeComponentBasicStandard)
       ! Shut down mass accretion onto the halo now that it is a satellite.
       call thisBasicComponent%accretionRateSet   (0.0d0                    )
       ! Record the time at which the node became a satellite - used for computing halo scales etc.
       call thisBasicComponent%timeLastIsolatedSet(thisBasicComponent%time())
    end select
    return
  end subroutine Node_Component_Basic_Standard_Stop_Accretion

   !# <nodePromotionTask>
   !#  <unitName>Node_Component_Basic_Standard_Promote</unitName>
   !# </nodePromotionTask>
   subroutine Node_Component_Basic_Standard_Promote(thisNode)
     !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the mass of {\tt thisNode}
     !% to be that of its parent.
     use Galacticus_Error
     implicit none
     type (treeNode          ), intent(inout), pointer :: thisNode                                  
     type (treeNode          )               , pointer :: parentNode                                
     class(nodeComponentBasic)               , pointer :: parentBasicComponent, thisBasicComponent  
     
     ! Get the basic component.                                                                                            
     thisBasicComponent => thisNode%basic()
     ! Ensure that it is of the standard class.
     select type (thisBasicComponent)
     class is (nodeComponentBasicStandard)
        ! Get the parent node and its basic component.
        parentNode           => thisNode  %parent
        parentBasicComponent => parentNode%basic()
        ! Ensure the two halos exist at the same time.
        if (thisBasicComponent%time() /= parentBasicComponent%time()) call Galacticus_Error_Report('Node_Component_Basic_Standard_Promote','thisNode&
             & has not been evolved to its parent')
        ! Adjust the mass to that of the parent node.
        call thisBasicComponent%massSet         (parentBasicComponent%mass         ())
        ! Adjust the accretion rate to that of the parent node.
        call thisBasicComponent%accretionRateSet(parentBasicComponent%accretionRate())
     end select
     return
   end subroutine Node_Component_Basic_Standard_Promote

  double precision function Node_Component_Basic_Standard_Unresolved_Mass(thisNode)
    !% Return the unresolved mass for {\tt thisNode}.
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode                                 
    type (treeNode          )               , pointer :: childNode                                
    class(nodeComponentBasic)               , pointer :: childBasicComponent, thisBasicComponent  
    
    ! Get the basic component.                                                                                           
    thisBasicComponent => thisNode%basic()
    ! Initialize the unresolved mass to the mass of the current node's basic component.
    Node_Component_Basic_Standard_Unresolved_Mass=thisBasicComponent%mass()
    ! Remove the mass of all child nodes.
    childNode => thisNode%firstChild
    do while (associated(childNode))
       childBasicComponent                           => childNode%basic()
       Node_Component_Basic_Standard_Unresolved_Mass =  Node_Component_Basic_Standard_Unresolved_Mass-childBasicComponent%mass()
       childNode                                     => childNode%sibling
    end do
    return
  end function Node_Component_Basic_Standard_Unresolved_Mass

end module Node_Component_Basic_Standard
