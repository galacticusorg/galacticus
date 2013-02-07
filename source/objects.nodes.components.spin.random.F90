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

!% Contains a module of spin tree node methods.

module Node_Component_Spin_Random
  !% Implement random spin tree node method.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Spin_Random_Initialize  , Node_Component_Spin_Random_Initialize_Spins, &
       &    Node_Component_Spin_Random_Promote
  
  !# <component>
  !#  <class>spin</class>
  !#  <name>random</name>
  !#  <isDefault>yes</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>spin</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="0.0d0" comment="Spin parameter of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>spinGrowthRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>yes</isVirtual>
  !#     <getFunction bindsTo="component">SpinRandomSpinGrowthRate</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.spin.random.bound_functions.inc</functions>
  !# </component>

  ! Record of whether the module has been initialized.
  logical          :: moduleInitialized=.false.

  ! The factor by which the mass of a node must increase before its spin parameter is re-chosen.
  double precision :: randomSpinResetMassFactor

contains

  subroutine Node_Component_Spin_Random_Initialize()
    !% Initializes the random spin component module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    
    ! Test whether module is already initialize.
    !$omp critical (Node_Component_Spin_Random_Initialize)
    if (.not.moduleInitialized) then
       
       !@ <inputParameter>
       !@   <name>randomSpinResetMassFactor</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which a node must increase in mass before its spin parameter is reset.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('randomSpinResetMassFactor',randomSpinResetMassFactor,defaultValue=2.0d0)
       ! Record that the module is now initialized.
       moduleInitialized=.true.     
    end if
    !$omp end critical (Node_Component_Spin_Random_Initialize)
    return
  end subroutine Node_Component_Spin_Random_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Spin_Random_Initialize_Spins</unitName>
  !#  <sortName>spin</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Spin_Random_Initialize_Spins(thisNode)
    !% Initialize the spin of {\tt thisNode}.
    use Cosmological_Parameters
    use Halo_Spin_Distributions
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    type (treeNode          ), pointer                :: relatedNode
    class(nodeComponentSpin ), pointer                :: thisSpinComponent,relatedSpinComponent
    class(nodeComponentBasic), pointer                :: relatedBasicComponent
    double precision                                  :: previousSetMass,previousSetSpin

    ! Check if we are the default method.
    if (defaultSpinComponent%randomIsActive()) then
       ! Ensure the module is initialized.
       call Node_Component_Spin_Random_Initialize()
       ! Get the basic component.
       thisSpinComponent => thisNode%spin()
       ! Ensure that the spin has not yet been assigned for this node.
       select type (thisSpinComponent)
       type is (nodeComponentSpin)
          ! Walk the tree back along primary children to the earliest such progenitor.
          relatedNode => thisNode
          do while (associated(relatedNode%firstChild))
             relatedNode => relatedNode%firstChild
          end do
          ! Walk forward through the branch, assigning spins. If the mass of the halo exceeds that of the halo for which we last
          ! selected a spin by a given factor, then select a new spin from the distribution. Otherwise, use the previously
          ! assigned spin.
          relatedSpinComponent  => relatedNode%spin (autoCreate=.true.)
          relatedBasicComponent => relatedNode%basic(                 )
          previousSetSpin=Halo_Spin_Distribution_Sample(relatedNode)
          previousSetMass=relatedBasicComponent%mass()
          call relatedSpinComponent%spinSet(previousSetSpin)
          do while (relatedNode%isPrimaryProgenitor())
             relatedNode           => relatedNode%parent
             relatedBasicComponent => relatedNode%basic()
             if (relatedBasicComponent%mass() > randomSpinResetMassFactor*previousSetMass) then
                previousSetSpin=Halo_Spin_Distribution_Sample(relatedNode)
                previousSetMass=relatedBasicComponent%mass()
             end if
             relatedSpinComponent => relatedNode%spin(autoCreate=.true.)
             call relatedSpinComponent%spinSet(previousSetSpin)
          end do
       end select
    end if
    return
  end subroutine Node_Component_Spin_Random_Initialize_Spins

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Spin_Random_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Spin_Random_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the spin of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    type (treeNode          ), pointer                :: parentNode
    class(nodeComponentSpin ), pointer                :: thisSpinComponent ,parentSpinComponent
    class(nodeComponentBasic), pointer                :: thisBasicComponent,parentBasicComponent

    ! Ensure that the spin component is of the random class.
    thisSpinComponent => thisNode%spin()
    select type (thisSpinComponent)
    class is (nodeComponentSpinRandom)
       parentNode           => thisNode  %parent
       thisBasicComponent   => thisNode  %basic()
       parentBasicComponent => parentNode%basic()
       if (thisBasicComponent%time() /= parentBasicComponent%time()) call Galacticus_Error_Report('Node_Component_Spin_Random_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the mass to that of the parent node.
       parentSpinComponent => parentNode%spin()
       call thisSpinComponent%spinSet(parentSpinComponent%spin())
    end select
    return
  end subroutine Node_Component_Spin_Random_Promote

end module Node_Component_Spin_Random
