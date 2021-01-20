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

!% Contains a module with the standard implementation of N-body component method.

module Node_Component_NBody_Generic
  !% An implementation of the N-body component which supports generic properties.
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: Node_Component_NBody_Generic_Thread_Initialize, Node_Component_NBody_Generic_Initialize         , &
       &    Node_Component_NBody_Generic_Output_Names     , Node_Component_NBody_Generic_Output_Count       , &
       &    Node_Component_NBody_Generic_Output           , Node_Component_NBody_Generic_Thread_Uninitialize

  !# <component>
  !#  <class>nBody</class>
  !#  <name>generic</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>reals</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>integers</name>
  !#     <type>longInteger</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#    <binding method="addRealProperty" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>integer</type>
  !#      <rank>0</rank>
  !#      <self pass="true" intent="inout" />
  !#      <argument>character(len=*), intent(in   ) :: propertyName</argument>
  !#     </interface>
  !#    </binding>
  !#    <binding method="addIntegerProperty" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>integer</type>
  !#      <rank>0</rank>
  !#      <self pass="true" intent="inout" />
  !#      <argument>character(len=*), intent(in   ) :: propertyName</argument>
  !#     </interface>
  !#    </binding>
  !#    <binding method="setRealProperty" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>void</type>
  !#      <self pass="true" intent="inout" />
  !#      <argument>integer         , intent(in   ) :: propertyIndex</argument>
  !#      <argument>double precision, intent(in   ) :: propertyValue</argument>
  !#     </interface>
  !#    </binding>
  !#    <binding method="setIntegerProperty" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>void</type>
  !#      <self pass="true" intent="inout" />
  !#      <argument>integer           , intent(in   ) :: propertyIndex</argument>
  !#      <argument>integer(kind_int8), intent(in   ) :: propertyValue</argument>
  !#     </interface>
  !#    </binding>
  !#  </bindings>
  !# </component>
  
  ! Property names.
  type   (varying_string), allocatable, dimension(:) :: propertyNamesReal        , propertyNamesInteger

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_NBody_Generic_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_NBody_Generic_Initialize(parameters_)
    !% Initializes the generic N-body component module.
    use :: Galacticus_Nodes, only : defaultNBodyComponent, nodeComponentNBodyGeneric
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters             ), intent(inout) :: parameters_
    type(nodeComponentNBodyGeneric   )                :: nbodyComponent
    !$GLC attributes unused :: parameters_

    ! Initialize the module if necessary.
    if (defaultNBodyComponent%genericIsActive()) then
       call nbodyComponent%addRealPropertyFunction   (Node_Component_NBody_Generic_Add_Real_Property   )
       call nbodyComponent%addIntegerPropertyFunction(Node_Component_NBody_Generic_Add_Integer_Property)
       call nbodyComponent%setRealPropertyFunction   (Node_Component_NBody_Generic_Set_Real_Property   )
       call nbodyComponent%setIntegerPropertyFunction(Node_Component_NBody_Generic_Set_Integer_Property)
    end if
    return
  end subroutine Node_Component_NBody_Generic_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_NBody_Generic_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_NBody_Generic_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent   , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultNBodyComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultNBodyComponent%genericIsActive()) &
         call nodePromotionEvent%attach(defaultNBodyComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentNBodyGeneric')
    return
  end subroutine Node_Component_NBody_Generic_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_NBody_Generic_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_NBody_Generic_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultNBodyComponent
    implicit none

    if (defaultNBodyComponent%genericIsActive()) &
         & call nodePromotionEvent%detach(defaultNBodyComponent,nodePromotion)
    return
  end subroutine Node_Component_NBody_Generic_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the
    !% properties of {\normalfont \ttfamily node} to be those of its parent.
    use :: Galacticus_Nodes, only : nodeComponentNBody, nodeComponentNBodyGeneric, treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), target  :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentNBody)               , pointer :: nBodyParent, nBody
    !$GLC attributes unused :: self
    
    nBody       => node      %nBody ()
    nodeParent  => node      %parent
    nBodyParent => nodeParent%nBody ()
    ! Copy properties from the parent.
    call nBody%   realsSet(nBodyParent%reals   ())
    call nBody%integersSet(nBodyParent%integers())
    return
  end subroutine nodePromotion

  function Node_Component_NBody_Generic_Add_Real_Property(self,propertyName) result(propertyIndex)
    !% Add a new real-valued property to this component, and return the index of the property.
    use :: Galacticus_Nodes  , only : nodeComponentNBodyGeneric
    use :: ISO_Varying_String, only : operator(==)             , assignment(=)
    integer                                                         :: propertyIndex
    class    (nodeComponentNBodyGeneric), intent(inout)             :: self
    character(len=*                    ), intent(in   )             :: propertyName
    type     (varying_string           ), allocatable, dimension(:) :: propertyNamesTmp
    logical                                                         :: nameExists
    integer                                                         :: i
    !$GLC attributes unused :: self

    ! Check for prior existance of the property.
    !$omp critical(nbodyGenericAccess)
    nameExists=.false.
    if (allocated(propertyNamesReal)) then
       do i=1,size(propertyNamesReal)
          if (propertyNamesReal(i) == propertyName) then
             nameExists=.true.
             propertyIndex=i
             exit
          end if
       end do
    end if
    if (.not.nameExists) then
       ! Property must be added.
       if (allocated(propertyNamesReal)) then
          call move_alloc(propertyNamesReal,propertyNamesTmp)
          allocate(propertyNamesReal(size(propertyNamesTmp)+1))
          propertyNamesReal(1:size(propertyNamesTmp))=propertyNamesTmp
          deallocate(propertyNamesTmp)
       else
          allocate(propertyNamesReal(1))
       end if
       propertyIndex=size(propertyNamesReal)
       propertyNamesReal(propertyIndex)=propertyName
    end if
    !$omp end critical(nbodyGenericAccess)
    return
  end function Node_Component_NBody_Generic_Add_Real_Property

  function Node_Component_NBody_Generic_Add_Integer_Property(self,propertyName) result(propertyIndex)
    !% Add a new real-valued property to this component, and return the index of the property.
    use :: Galacticus_Nodes  , only : nodeComponentNBodyGeneric
    use :: ISO_Varying_String, only : operator(==)             , assignment(=)
    implicit none
    integer                                                         :: propertyIndex
    class    (nodeComponentNBodyGeneric), intent(inout)             :: self
    character(len=*                    ), intent(in   )             :: propertyName
    type     (varying_string           ), allocatable, dimension(:) :: propertyNamesTmp
    logical                                                         :: nameExists
    integer                                                         :: i
    !$GLC attributes unused :: self

    ! Check for prior existance of the property.
    !$omp critical(nbodyGenericAccess)
    nameExists=.false.
    if (allocated(propertyNamesInteger)) then
       do i=1,size(propertyNamesInteger)
          if (propertyNamesInteger(i) == propertyName) then
             nameExists=.true.
             propertyIndex=i
             exit
          end if
       end do
    end if
    if (.not.nameExists) then
       ! Property must be added.
       if (allocated(propertyNamesInteger)) then
          call move_alloc(propertyNamesInteger,propertyNamesTmp)
          allocate(propertyNamesInteger(size(propertyNamesTmp)+1))
          propertyNamesInteger(1:size(propertyNamesTmp))=propertyNamesTmp
          deallocate(propertyNamesTmp)
       else
          allocate(propertyNamesInteger(1))
       end if
       propertyIndex=size(propertyNamesInteger)
       propertyNamesInteger(propertyIndex)=propertyName
    end if
    !$omp end critical(nbodyGenericAccess)
    return
  end function Node_Component_NBody_Generic_Add_Integer_Property

  subroutine Node_Component_NBody_Generic_Set_Real_Property(self,propertyIndex,propertyValue)
    !% Set the value of the indexed real property.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentNBodyGeneric
    implicit none
    class           (nodeComponentNBodyGeneric), intent(inout)               :: self
    integer                                    , intent(in   )               :: propertyIndex
    double precision                           , intent(in   )               :: propertyValue
    double precision                           , allocatable  , dimension(:) :: propertyValuesCurrent, propertyValuesNew

    if (.not.allocated(propertyNamesReal).or.propertyIndex > size(propertyNamesReal)) call Galacticus_Error_Report('property index is out of range'//{introspection:location})
    allocate(propertyValuesNew(size(propertyNamesReal)))
    propertyValuesCurrent           =self         %reals()
    if (size(propertyValuesCurrent) > 0) propertyValuesNew(1:size(propertyValuesCurrent))=propertyValuesCurrent
    propertyValuesNew(propertyIndex)=propertyValue
    call self%realsSet(propertyValuesNew)
    return
  end subroutine Node_Component_NBody_Generic_Set_Real_Property

  subroutine Node_Component_NBody_Generic_Set_Integer_Property(self,propertyIndex,propertyValue)
    !% Set the value of the indexed real property.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentNBodyGeneric
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class           (nodeComponentNBodyGeneric), intent(inout)               :: self
    integer                                    , intent(in   )               :: propertyIndex
    integer         (kind_int8                ), intent(in   )               :: propertyValue
    integer         (kind_int8                ), allocatable  , dimension(:) :: propertyValuesCurrent, propertyValuesNew

    if (.not.allocated(propertyNamesInteger).or.propertyIndex > size(propertyNamesInteger)) call Galacticus_Error_Report('property index is out of range'//{introspection:location})
    allocate(propertyValuesNew(size(propertyNamesReal)))
    propertyValuesCurrent           =self         %integers()
    if (size(propertyValuesCurrent) > 0) propertyValuesNew(1:size(propertyValuesCurrent))=propertyValuesCurrent
    propertyValuesNew(propertyIndex)=propertyValue
    call self%integersSet(propertyValuesNew)
    return
  end subroutine Node_Component_NBody_Generic_Set_Integer_Property

  !# <mergerTreeOutputNames>
  !#  <unitName>Node_Component_NBody_Generic_Output_Names</unitName>
  !#  <sortName>Node_Component_NBody_Generic_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Node_Component_NBody_Generic_Output_Names(node,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI ,doubleProperty,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use :: Galacticus_Nodes  , only : treeNode
    use :: ISO_Varying_String, only : char
    use :: String_Handling   , only : String_Upper_Case_First, char
    implicit none
    type            (treeNode)              , intent(inout) :: node
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                 :: i
    !$GLC attributes unused :: node, time


    !$omp critical(nbodyGenericAccess)
    if (allocated(propertyNamesInteger)) then
       do i=1,size(propertyNamesInteger)
          integerProperty                         =integerProperty+1
          integerPropertyNames   (integerProperty)='nBody'//String_Upper_Case_First(char(propertyNamesInteger(i)))
          integerPropertyComments(integerProperty)=''
          integerPropertyUnitsSI (integerProperty)=0.0d0
       end do
    end if
    if (allocated(propertyNamesReal)) then
       do i=1,size(propertyNamesReal)
          doubleProperty                        =doubleProperty+1
          doublePropertyNames   (doubleProperty)='nBody'//String_Upper_Case_First(char(propertyNamesReal(i)))
          doublePropertyComments(doubleProperty)=''
          doublePropertyUnitsSI (doubleProperty)=0.0d0
       end do
    end if
    !$omp end critical(nbodyGenericAccess)
    return
  end subroutine Node_Component_NBody_Generic_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Node_Component_NBody_Generic_Output_Count</unitName>
  !#  <sortName>Node_Component_NBody_Generic_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Node_Component_NBody_Generic_Output_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !$GLC attributes unused :: node, time

    !$omp critical(nbodyGenericAccess)
    if (allocated(propertyNamesInteger)) integerPropertyCount=integerPropertyCount+size(propertyNamesInteger)
    if (allocated(propertyNamesReal   ))  doublePropertyCount= doublePropertyCount+size(propertyNamesReal   )
    !$omp end critical(nbodyGenericAccess)
    return
  end subroutine Node_Component_NBody_Generic_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Node_Component_NBody_Generic_Output</unitName>
  !#  <sortName>Node_Component_NBody_Generic_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Node_Component_NBody_Generic_Output(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance)
    !% Store black hole properties in the \glc\ output file buffers.
    use :: Galacticus_Nodes, only : nodeComponentNBody, treeNode
    use :: Kind_Numbers    , only : kind_int8
    use :: Multi_Counters  , only : multiCounter
    implicit none
    double precision                    , intent(in   )               :: time
    type            (treeNode          ), intent(inout)               :: node
    integer                             , intent(inout)               :: doubleBufferCount          , doubleProperty, integerBufferCount, &
         &                                                               integerProperty
    integer         (kind=kind_int8    ), intent(inout)               :: integerBuffer         (:,:)
    double precision                    , intent(inout)               :: doubleBuffer          (:,:)
    type            (multiCounter      ), intent(inout )              :: instance
    class           (nodeComponentNBody)               , pointer      :: nBody
    integer         (kind=kind_int8    ), allocatable  , dimension(:) :: propertyValuesInteger
    double precision                    , allocatable  , dimension(:) :: propertyValuesReal
    integer                                                           :: i
    !$GLC attributes unused :: time, instance

    nBody => node%nBody()
    !$omp critical(nbodyGenericAccess)
    if (allocated(propertyNamesInteger)) then
       propertyValuesInteger=nBody%integers()
       do i=1,size(propertyNamesInteger)
          integerProperty=integerProperty+1
          if (i > size(propertyValuesInteger)) then
             integerBuffer(integerBufferCount,integerProperty)=0_kind_int8
          else
             integerBuffer(integerBufferCount,integerProperty)=propertyValuesInteger(i)
          end if
       end do
    end if
    if (allocated(propertyNamesReal)) then
       propertyValuesReal=nBody%reals()
       do i=1,size(propertyNamesReal)
          doubleProperty=doubleProperty+1
          if (i > size(propertyValuesReal)) then
             doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
          else
             doubleBuffer(doubleBufferCount,doubleProperty)=propertyValuesReal(i)
          end if
       end do
    end if
    !$omp end critical(nbodyGenericAccess)
    return
  end subroutine Node_Component_NBody_Generic_Output

end module Node_Component_NBody_Generic
