!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module that implements a very simple disk component.

module Node_Component_Disk_Very_Simple_Size
  !% Implements a very simple disk component.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility, Node_Component_Disk_Very_Simple_Size_Radius_Solver, &
       &    Node_Component_Disk_Very_Simple_Size_Initialize

  !# <component>
  !#  <class>disk</class>
  !#  <name>verySimpleSize</name>
  !#  <extends>
  !#   <class>disk</class>
  !#   <name>verySimple</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaparsec" comment="Radial scale length in the disk."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction>Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity of the disk."/>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.disk.very_simple.size.bound_functions.inc</functions>
  !# </component>

  ! Parameters controlling the physical implementation.
  double precision :: diskMassToleranceAbsolute

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized        =.false.

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Size_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Disk_Very_Simple_Size_Initialize()
    !% Initializes the tree node exponential disk methods module.
    use Input_Parameters
    implicit none

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Disk_Very_Simple_Size_Initialize)
    if (defaultDiskComponent%verySimpleSizeIsActive().and..not.moduleInitialized) then
       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>diskMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the disk is physically plausible.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskMassToleranceAbsolute',diskMassToleranceAbsolute,defaultValue=1.0d-6)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Disk_Very_Simple_Size_Initialize)
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Initialize
  
  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode         ), intent(inout) :: thisNode
    logical                            , intent(inout) :: galaxyIsPhysicallyPlausible
    class           (nodeComponentDisk), pointer       :: thisDisk

    ! Return immediately if our method is not selected.
    if (.not.defaultDiskComponent%verySimpleSizeIsActive()) return

     ! Determine the plausibility of the current disk.
     thisDisk => thisNode%disk()
     select type (thisDisk)
     class is (nodeComponentDiskVerySimpleSize)
        galaxyIsPhysicallyPlausible=(thisDisk%massStellar()+thisDisk%massGas() >= -diskMassToleranceAbsolute)
     end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Size_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver(thisNode,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get&
       &,Radius_Set,Velocity_Get,Velocity_Set)
    !% Interface for the size solver algorithm.
    use Dark_Matter_Halo_Spins
    implicit none
    type            (treeNode                                       ), intent(inout)          :: thisNode
    logical                                                          , intent(  out)          :: componentActive
    logical                                                          , intent(in   )          :: specificAngularMomentumRequired
    double precision                                                 , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Disk_Very_Simple_Size_Radius    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Disk_Very_Simple_Size_Radius_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentDisk                              )               , pointer :: thisDisk
    class           (nodeComponentBasic                             )               , pointer :: thisBasic

    ! Determine if thisNode has an active disk component supported by this module.
    componentActive =  .false.
    thisDisk        => thisNode%disk()
    select type (thisDisk)
    class is (nodeComponentDiskVerySimpleSize)
       componentActive        =  .true.
       if (specificAngularMomentumRequired) then
          thisBasic              => thisNode%basic()
          specificAngularMomentum=  Dark_Matter_Halo_Angular_Momentum(thisNode)/thisBasic%mass()
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Disk_Very_Simple_Size_Radius
       Radius_Set   => Node_Component_Disk_Very_Simple_Size_Radius_Set
       Velocity_Get => Node_Component_Disk_Very_Simple_Size_Velocity
       Velocity_Set => Node_Component_Disk_Very_Simple_Size_Velocity_Set
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver

  double precision function Node_Component_Disk_Very_Simple_Size_Radius(thisNode)
    !% Return the radius of the disk used in structure solvers.
    implicit none
    type (treeNode         ), intent(inout) :: thisNode
    class(nodeComponentDisk), pointer       :: thisDisk

    thisDisk => thisNode%disk()
    Node_Component_Disk_Very_Simple_Size_Radius=thisDisk%radius()
    return
  end function Node_Component_Disk_Very_Simple_Size_Radius

  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Set(thisNode,radius)
    !% Set the radius of the disk used in structure solvers.
    implicit none
    type            (treeNode         ), intent(inout) :: thisNode
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: thisDisk

    thisDisk => thisNode%disk()
    call thisDisk%radiusSet(radius)
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Set

  double precision function Node_Component_Disk_Very_Simple_Size_Velocity(thisNode)
    !% Return the circular velocity of the disk.
    implicit none
    type (treeNode         ), intent(inout) :: thisNode
    class(nodeComponentDisk), pointer       :: thisDisk

    thisDisk => thisNode%disk()
    Node_Component_Disk_Very_Simple_Size_Velocity=thisDisk%velocity()
    return
  end function Node_Component_Disk_Very_Simple_Size_Velocity

  subroutine Node_Component_Disk_Very_Simple_Size_Velocity_Set(thisNode,velocity)
    !% Set the circular velocity of the disk.
    implicit none
    type            (treeNode         ), intent(inout) :: thisNode
    double precision                   , intent(in   ) :: velocity
    class           (nodeComponentDisk), pointer       :: thisDisk

    thisDisk => thisNode%disk()
    call thisDisk%velocitySet(velocity)
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Velocity_Set

end module Node_Component_Disk_Very_Simple_Size
