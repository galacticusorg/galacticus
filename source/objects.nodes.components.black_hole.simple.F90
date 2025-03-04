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
Contains a module which implements the simple black hole node component.
!!}

module Node_Component_Black_Hole_Simple
  !!{
  Implements the simple black hole node component.
  !!}
  use :: Black_Hole_Binary_Mergers , only : blackHoleBinaryMergerClass
  implicit none
  private
  public :: Node_Component_Black_Hole_Simple_Thread_Uninitialize, Node_Component_Black_Hole_Simple_Thread_Initialize, &
       &    Node_Component_Black_Hole_Simple_State_Store        , Node_Component_Black_Hole_Simple_State_Restore    , &
       &    Node_Component_Black_Hole_Simple_Scale_Set 

  !![
  <component>
   <class>blackHole</class>
   <name>simple</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of the black hole."/>
    </property>
   </properties>
   <bindings>
     <binding method="massDistribution" function="Node_Component_Black_Hole_Simple_Mass_Distribution" bindsTo="component"/>
     <binding method="massBaryonic"     function="Node_Component_Black_Hole_Simple_Mass_Baryonic"     bindsTo="component"/>
   </bindings>
   <functions>objects.nodes.components.black_hole.simple.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(blackHoleBinaryMergerClass), pointer :: blackHoleBinaryMerger_
  !$omp threadprivate(blackHoleBinaryMerger_)

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Black_Hole_Simple_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Thread_Initialize(parameters)
    !!{
    Initializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent     , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters

    if (defaultBlackHoleComponent%simpleIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentBlackHole')
       !![
       <objectBuilder class="blackHoleBinaryMerger" name="blackHoleBinaryMerger_" source="subParameters"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(thread,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentBlackHoleSimple',dependencies=dependencies)
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Black_Hole_Simple_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Thread_Uninitialize()
    !!{
    Uninitializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none

    if (defaultBlackHoleComponent%simpleIsActive()) then
       !![
       <objectDestructor name="blackHoleBinaryMerger_"/>
       !!]
       if (satelliteMergerEvent%isAttached(thread,satelliteMerger)) call satelliteMergerEvent%detach(thread,satelliteMerger)
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Thread_Uninitialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Black_Hole_Simple_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole   , nodeComponentBlackHoleSimple, nodeComponentSpheroid, treeNode, &
         &                          defaultBlackHoleComponent
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    double precision                        , parameter              :: massScaleAbsolute=1.0d+0, massScaleRelative=1.0d-3

    ! Check if we are the default method.
    if (.not.defaultBlackHoleComponent%simpleIsActive()) return
    ! Get the black hole component.
    blackHole => node%blackHole()
    ! Ensure that it is of the simple class.
    select type (blackHole)
    class is (nodeComponentBlackHoleSimple)
       ! Get the spheroid component.
       spheroid => node%spheroid()
       ! Set scale for mass.
       call blackHole%massScale(                                                         &
             &                  max(                                                     &
            &                           massScaleRelative*spheroid %massStellar      (), &
            &                       max(                                                 &
            &                                                       massScaleAbsolute  , &
            &                                             blackHole%mass             ()  &
            &                          )                                                 &
            &                      )                                                     &
            &                  )
    end select
    return
  end subroutine Node_Component_Black_Hole_Simple_Scale_Set

  subroutine satelliteMerger(self,node)
    !!{
    Merge (instantaneously) any simple black hole associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    class           (*                     ), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: nodeHost
    class           (nodeComponentBlackHole), pointer       :: blackHoleHost   , blackHole
    double precision                                        :: massBlackHoleNew, spinBlackHoleNew
    !$GLC attributes unused :: self
    
    ! Find the node to merge with.
    nodeHost      => node    %mergesWith(                 )
    ! Get the black holes.
    blackHole     => node    %blackHole (autoCreate=.true.)
    blackHoleHost => nodeHost%blackHole (autoCreate=.true.)
    ! Compute the effects of the merger.
    call blackHoleBinaryMerger_%merge(blackHole    %mass(), &
         &                            blackHoleHost%mass(), &
         &                            0.0d0               , &
         &                            0.0d0               , &
         &                            massBlackHoleNew    , &
         &                            spinBlackHoleNew      &
         &                           )
    ! Move the black hole to the host.
    call blackHoleHost%massSet(massBlackHoleNew)
    call blackHole    %massSet(           0.0d0)
    return
  end subroutine satelliteMerger

  !![
  <stateStoreTask>
   <unitName>Node_Component_Black_Hole_Simple_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentBlackHole -> simple',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="blackHoleBinaryMerger_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_Simple_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Black_Hole_Simple_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentBlackHole -> simple',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="blackHoleBinaryMerger_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_Simple_State_Restore

end module Node_Component_Black_Hole_Simple
