!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module that implements a very simple disk component.
!!}

module Node_Component_Disk_Very_Simple_Size
  !!{
  Implements a very simple disk component.
  !!}
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility, Node_Component_Disk_Very_Simple_Size_Radius_Solver    , &
       &    Node_Component_Disk_Very_Simple_Size_Initialize                , Node_Component_Disk_Very_Simple_Size_Thread_Initialize, &
       &    Node_Component_Disk_Very_Simple_Size_State_Store               , Node_Component_Disk_Very_Simple_Size_State_Retrieve   , &
       &    Node_Component_Disk_Very_Simple_Size_Thread_Uninitialize

  !![
  <component>
   <class>disk</class>
   <name>verySimpleSize</name>
   <extends>
    <class>disk</class>
    <name>verySimple</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>radius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="megaparsec" comment="Radial scale length in the disk."/>
    </property>
    <property>
      <name>halfMassRadius</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius</getFunction>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="kilo" comment="Circular velocity of the disk."/>
    </property>
   </properties>
   <bindings>
    <binding method="massDistribution" function="Node_Component_Disk_Very_Simple_Size_Mass_Distribution"/>
   </bindings>
   <functions>objects.nodes.components.disk.very_simple.size.bound_functions.inc</functions>
  </component>
  !!]

  ! Parameters controlling the physical implementation.
  double precision :: toleranceAbsoluteMass

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_Initialize(parameters)
    !!{
    Initializes the tree node exponential disk methods module.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    if (defaultDiskComponent%verySimpleSizeIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>toleranceAbsoluteMass</name>
         <defaultValue>1.0d-6</defaultValue>
         <description>The mass tolerance used to judge whether the disk is physically plausible.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Initialize
  
  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_Thread_Initialize(parameters)
    !!{
    Initializes the tree node very simple size disk module.
    !!}
    use :: Error                                    , only : Error_Report
    use :: Galacticus_Nodes                         , only : defaultDiskComponent
    use :: Input_Parameters                         , only : inputParameter             , inputParameters
    use :: Mass_Distributions                       , only : massDistributionCylindrical
    use :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_   , massDistributionGas_
    use :: Galactic_Structure_Options               , only : componentTypeDisk          , massTypeStellar     , massTypeGaseous
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    if (defaultDiskComponent%verySimpleSizeIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       !![
       <objectBuilder class="massDistribution" parameterName="massDistributionDisk" name="massDistributionStellar_" source="subParameters" threadPrivate="yes">
        <default>
         <massDistributionDisk value="exponentialDisk">
          <dimensionless value="true"/>
         </massDistributionDisk>
        </default>
       </objectBuilder>
       !!]
      ! Validate the disk mass distribution.
       select type (massDistributionStellar_)
       class is (massDistributionCylindrical)
          ! The disk mass distribution must have cylindrical symmetry. So, this is acceptable.        
       class default
          call Error_Report('only cylindrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       if (.not.massDistributionStellar_%isDimensionless()) call Error_Report('disk mass distribution must be dimensionless'//{introspection:location})
       ! Duplicate the dimensionless mass distribution to use for the gas component, and set component and mass types in both.
       !$omp critical(diskVerySimpleSizeDeepCopy)
       allocate(massDistributionGas_,mold=massDistributionStellar_)
       !![
       <deepCopyReset variables="massDistributionStellar_"/>
       <deepCopy source="massDistributionStellar_" destination="massDistributionGas_"/>
       <deepCopyFinalize variables="massDistributionGas_"/>  
       !!]
       !$omp end critical(diskVerySimpleSizeDeepCopy)
       call massDistributionStellar_%setTypes(componentTypeDisk,massTypeStellar)
       call massDistributionGas_    %setTypes(componentTypeDisk,massTypeGaseous)
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_Thread_Uninitialize()
    !!{
    Uninitializes the tree node standard merging statistics module.
    !!}
    use :: Galacticus_Nodes                         , only : defaultDiskComponent
    use :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_, massDistributionGas_
    implicit none

    if (defaultDiskComponent%verySimpleSizeIsActive()) then
       !![
       <objectDestructor name="massDistributionStellar_"/>
       <objectDestructor name="massDistributionGas_"    />
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Thread_Uninitialize

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility</unitName>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility(node)
    !!{
    Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskVerySimpleSize, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    ! Return immediately if our method is not selected.
    if (.not.defaultDiskComponent%verySimpleSizeIsActive()) return

     ! Determine the plausibility of the current disk.
     disk => node%disk()
     select type (disk)
     class is (nodeComponentDiskVerySimpleSize)
        if (disk%massStellar()+disk%massGas() < -toleranceAbsoluteMass) node%isPhysicallyPlausible=.false.
     end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver_Plausibility

  !![
  <radiusSolverTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_Radius_Solver</unitName>
  </radiusSolverTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver(node,componentActive,component,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get&
       &,Radius_Set,Velocity_Get,Velocity_Set)
    !!{
    Interface for the size solver algorithm.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic          , nodeComponentDisk, nodeComponentDiskVerySimpleSize, treeNode, &
         &                                    nodeComponentSpin
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, componentTypeDisk
    implicit none
    type            (treeNode                                       ), intent(inout)          :: node
    logical                                                          , intent(  out)          :: componentActive
    type            (enumerationComponentTypeType                   ), intent(  out)          :: component
    logical                                                          , intent(in   )          :: specificAngularMomentumRequired
    double precision                                                 , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Disk_Very_Simple_Size_Radius    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Disk_Very_Simple_Size_Radius_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentDisk                              )               , pointer :: disk
    class           (nodeComponentBasic                             )               , pointer :: basic
    class           (nodeComponentSpin                              )               , pointer :: spin

    ! Determine if node has an active disk component supported by this module.
    componentActive         =  .false.
    component               =  componentTypeDisk
    specificAngularMomentum =  0.0d0
    disk                    => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimpleSize)
       componentActive=  .true.
       if (specificAngularMomentumRequired) then
          basic                   =>  node               %basic()
          spin                    =>  node               %spin ()
          specificAngularMomentum =  +spin%angularMomentum     () &
               &                     /basic              %mass ()
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Disk_Very_Simple_Size_Radius
       Radius_Set   => Node_Component_Disk_Very_Simple_Size_Radius_Set
       Velocity_Get => Node_Component_Disk_Very_Simple_Size_Velocity
       Velocity_Set => Node_Component_Disk_Very_Simple_Size_Velocity_Set
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Solver

  double precision function Node_Component_Disk_Very_Simple_Size_Radius(node)
    !!{
    Return the radius of the disk used in structure solvers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Very_Simple_Size_Radius=disk%radius()
    return
  end function Node_Component_Disk_Very_Simple_Size_Radius

  subroutine Node_Component_Disk_Very_Simple_Size_Radius_Set(node,radius)
    !!{
    Set the radius of the disk used in structure solvers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%radiusSet(radius)
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Radius_Set

  double precision function Node_Component_Disk_Very_Simple_Size_Velocity(node)
    !!{
    Return the circular velocity of the disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Very_Simple_Size_Velocity=disk%velocity()
    return
  end function Node_Component_Disk_Very_Simple_Size_Velocity

  subroutine Node_Component_Disk_Very_Simple_Size_Velocity_Set(node,velocity)
    !!{
    Set the circular velocity of the disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: velocity
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%velocitySet(velocity)
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_Velocity_Set

  !![
  <stateStoreTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the tabulation state to file.
    !!}
    use            :: Display                                  , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                            , only : c_ptr                   , c_size_t
    use            :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_, massDistributionGas_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionStellar_ massDistributionGas_"/>
    !!]
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Disk_Very_Simple_Size_State_Retrieve</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Size_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display                                  , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                            , only : c_ptr                   , c_size_t
    use            :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_, massDistributionGas_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionStellar_ massDistributionGas_"/>
    !!]
    return
  end subroutine Node_Component_Disk_Very_Simple_Size_State_Retrieve

end module Node_Component_Disk_Very_Simple_Size
