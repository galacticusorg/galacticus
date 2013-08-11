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

!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Node_Component_Dark_Matter_Profile_Scale
  !% Implements a dark matter profile method that provides a scale radius.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Rate_Compute, Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_Promote     , Node_Component_Dark_Matter_Profile_Scale_Scale_Set      , &
       &    Node_Component_Dark_Matter_Profile_Scale_Tree_Output , Node_Component_Dark_Matter_Profile_Scale_Plausibility   , &
       &    Node_Component_Dark_Matter_Profile_Scale_Initialize

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scale</name>
  !#  <isDefault>yes</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>scale</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
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

  ! Parameters of the method.
  double precision                                      :: darkMatterProfileMaximumConcentration                  , darkMatterProfileMinimumConcentration

  ! Flag indicating whether scale radius data should be output when full merger trees are output.
  logical                                               :: mergerTreeStructureOutputDarkMatterProfileScale

  ! Record of whether the module has been initialized.
  logical                                               :: moduleInitialized                              =.false.

  ! Queriable dark matter profile object.
  type            (nodeComponentDarkMatterProfileScale) :: darkMatterProfile

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize()
    !% Initializes the ``scale'' implementation of the dark matter halo profile component.
    use Input_Parameters
    implicit none

    ! Check if this implementation is selected.
    !$omp critical (Node_Component_Dark_Matter_Profile_Scale_Initialize)
    if (.not.moduleInitialized) then
       ! Get parameters.
       !@ <inputParameter>
       !@   <name>darkMatterProfileMinimumConcentration</name>
       !@   <defaultValue>4</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum concentration allowed for dark matter profiles.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMinimumConcentration',darkMatterProfileMinimumConcentration,defaultValue=4.0d0)
        !@ <inputParameter>
       !@   <name>darkMatterProfileMaximumConcentration</name>
       !@   <defaultValue>100</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum concentration allowed for dark matter profiles.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMaximumConcentration',darkMatterProfileMaximumConcentration,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputDarkMatterProfileScale</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not dark matter halo scale radius is included in outputs of merger trees.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputDarkMatterProfileScale',mergerTreeStructureOutputDarkMatterProfileScale,defaultValue=.false.)
       ! Bind the scale get function.
       call darkMatterProfile%scaleFunction(Node_Component_Dark_Matter_Profile_Scale_Scale)
       ! Record that the module is now initialize.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Dark_Matter_Profile_Scale_Initialize)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize

  double precision function Node_Component_Dark_Matter_Profile_Scale_Scale(self)
    !% Return the scale radius in the dark matter halo profile.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (nodeComponentDarkMatterProfileScale), intent(inout) :: self
    type            (treeNode                           ), pointer       :: selfNode
    double precision                                                     :: scaleLengthMaximum, scaleLengthMinimum

    selfNode => self%host()
    scaleLengthMaximum=Dark_Matter_Halo_Virial_Radius(selfNode)/darkMatterProfileMinimumConcentration
    scaleLengthMinimum=Dark_Matter_Halo_Virial_Radius(selfNode)/darkMatterProfileMaximumConcentration
    Node_Component_Dark_Matter_Profile_Scale_Scale= &
         & min(                                     &
         &     scaleLengthMaximum    ,              &
         &     max(                                 &
         &         scaleLengthMinimum,              &
         &         self%scaleValue()                &
         &        )                                 &
         &    )
    return
  end function Node_Component_Dark_Matter_Profile_Scale_Scale

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the scale radius.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    logical                                         , intent(inout)          :: interrupt
    procedure       (                              ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , growthRate

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)
       ! Find the concentration of this halo.
       concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/thisDarkMatterProfileComponent%scale()
       ! Find the growth rate and limit to ensure minimum and maximum concentrations are not exceeded.
       growthRate=thisDarkMatterProfileComponent%scaleGrowthRate()
       if (concentration <= darkMatterProfileMinimumConcentration) growthRate=min(growthRate,Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/darkMatterProfileMinimumConcentration)
       if (concentration >= darkMatterProfileMaximumConcentration) growthRate=max(growthRate,Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/darkMatterProfileMaximumConcentration)
       call thisDarkMatterProfileComponent%scaleRate(growthRate)
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Rate_Compute

  subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize_Scale(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    use Dark_Matter_Profiles_Concentrations
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    class           (nodeComponentDarkMatterProfile)               , pointer :: newDarkMatterProfileComponent, thisDarkMatterProfileComponent
    double precision                                                         :: concentration

    ! Ensure that the module is initialized.
    call Node_Component_Dark_Matter_Profile_Scale_Initialize()
    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    select type (thisDarkMatterProfileComponent)
    type is (nodeComponentDarkMatterProfile)
       newDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)
       ! Set the scale radius of the halo.
       concentration=max(Dark_Matter_Profile_Concentration(thisNode),darkMatterProfileMinimumConcentration)
       call newDarkMatterProfileComponent%scaleSet(Dark_Matter_Halo_Virial_Radius(thisNode)/concentration)
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize_Scale

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the dark matter profile is physically plausible for radius solving tasks.
    implicit none
    type   (treeNode                      ), intent(inout), pointer :: thisNode
    logical                                , intent(inout)          :: galaxyIsPhysicallyPlausible
    class  (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)
       if (thisDarkMatterProfileComponent%scale() <= 0.0d0) galaxyIsPhysicallyPlausible=.false.
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    class           (nodeComponentDarkMatterProfile)               , pointer :: parentDarkMatterProfileComponent, thisDarkMatterProfileComponent
    class           (nodeComponentBasic            )               , pointer :: parentBasicComponent            , thisBasicComponent
    double precision                                                         :: deltaTime

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       ! Ensure that current node has its scale set.
       call Node_Component_Dark_Matter_Profile_Scale_Initialize_Scale(thisNode)
       ! Get the dark matter profile component.
       thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
       ! Check if this node is the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          ! It is, so compute the scale radius growth rate.
          ! First ensure that parent node has scale radius set.
          call Node_Component_Dark_Matter_Profile_Scale_Initialize_Scale(thisNode%parent)
          ! Now compute the growth rate.
          thisBasicComponent   => thisNode       %basic()
          parentBasicComponent => thisNode%parent%basic()
          deltaTime=parentBasicComponent%time()-thisBasicComponent%time()
          if (deltaTime > 0.0d0) then
             parentDarkMatterProfileComponent => thisNode%parent%darkMatterProfile()
             call thisDarkMatterProfileComponent%scaleGrowthRateSet(                                           &
                  &                                                 (                                          &
                  &                                                   parentDarkMatterProfileComponent%scale() &
                  &                                                  -thisDarkMatterProfileComponent  %scale() &
                  &                                                 )                                          &
                  &                                                 /deltaTime                                 &
                  &                                                )
          else
             call thisDarkMatterProfileComponent%scaleGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set scale radius growth rate to zero.
          call    thisDarkMatterProfileComponent%scaleGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the growth rate of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile)               , pointer :: parentDarkMatterProfileComponent, thisDarkMatterProfileComponent
    class(nodeComponentBasic            )               , pointer :: parentBasicComponent            , thisBasicComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)
       parentDarkMatterProfileComponent => thisNode%parent%darkMatterProfile()
       thisBasicComponent               => thisNode       %basic            ()
       parentBasicComponent             => thisNode%parent%basic            ()
       if (thisBasicComponent%time() /= parentBasicComponent%time()) call Galacticus_Error_Report('Node_Component_Dark_Matter_Profile_Scale_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the scale radius to that of the parent node.
       call thisDarkMatterProfileComponent%scaleSet          (parentDarkMatterProfileComponent%scale          ())
       ! Adjust the growth rate to that of the parent node.
       call thisDarkMatterProfileComponent%scaleGrowthRateSet(parentDarkMatterProfileComponent%scaleGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Promote

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)
       ! Set scale for the scale radius.
       call thisDarkMatterProfileComponent%scaleScale(thisDarkMatterProfileComponent%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set

  !# <mergerTreeStructureOutputTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Tree_Output</unitName>
  !# </mergerTreeStructureOutputTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Output(baseNode,nodeProperty,treeGroup)
    !% Write the scale radius property to a full merger tree output.
    use IO_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode                      )              , intent(in   ), pointer :: baseNode
    double precision                                , dimension(:), intent(inout)          :: nodeProperty
    type            (hdf5Object                    )              , intent(inout)          :: treeGroup
    type            (treeNode                      )                             , pointer :: thisNode
    integer                                                                                :: nodeCount
    class           (nodeComponentDarkMatterProfile)                             , pointer :: baseDarkMatterProfileComponent, thisDarkMatterProfileComponent
    type            (hdf5Object                    )                                       :: nodeDataset

    ! Check if scale radius is to be included in merger tree outputs.
    if (mergerTreeStructureOutputDarkMatterProfileScale) then
       ! Get the dark matter profile component.
       baseDarkMatterProfileComponent => baseNode%darkMatterProfile()
       ! Ensure it is of the scale class.
       select type (baseDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)
          ! Extract node scale radius and output to file.
          nodeCount=0
          thisNode => baseNode
          do while (associated(thisNode))
             thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=thisDarkMatterProfileComponent%scale()
             call thisNode%walkTree(thisNode)
          end do
          call treeGroup  %writeDataset  (nodeProperty,'darkMatterScaleRadius','Scale radius of the dark matter profile [Mpc].',datasetReturned=nodeDataset)
          call nodeDataset%writeAttribute(megaParsec,"unitsInSI")
          call nodeDataset%close         ()
       end select
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Output

end module Node_Component_Dark_Matter_Profile_Scale
