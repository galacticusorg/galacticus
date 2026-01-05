!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module implementing a scalar angular momentum component for dark matter halos.
!!}

module Node_Component_Halo_Angular_Momentum_Scalar
  !!{
  Implement a scalar spin component for tree nodes.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Halo_Angular_Momentum_Scalar_Thread_Initialize, Node_Component_Halo_Angular_Momentum_Scalar_Thread_Uninitialize, &
       &    Node_Component_Halo_Angular_Momentum_Scalar_Scale_Set        , Node_Component_Halo_Angular_Momentum_Scalar_State_Store        , &
       &    Node_Component_Halo_Angular_Momentum_Scalar_State_Restore

  !![
  <component>
   <class>spin</class>
   <name>scalar</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum magnitude of the DMO halo."/>
    </property>
    <property>
      <name>angularMomentumGrowthRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  !$omp threadprivate(darkMatterHaloScale_)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Scalar_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Scalar_Thread_Initialize(parameters)
    !!{
    Initializes the halo scalar angular momentum module.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (defaultSpinComponent%scalarIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
       !!]
     end if
     return
  end subroutine Node_Component_Halo_Angular_Momentum_Scalar_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Scalar_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Scalar_Thread_Uninitialize()
    !!{
    Uninitializes the halo scalar angular momentum module.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%scalarIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Scalar_Thread_Uninitialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Scalar_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Scalar_Scale_Set(node)
    !!{
    Set scales for properties in the scalar implementation of the spin component.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , nodeComponentSpinScalar, treeNode, defaultSpinComponent
    implicit none
    type            (treeNode         ), intent(inout), pointer   :: node
    class           (nodeComponentSpin)               , pointer   :: spin
    double precision                                  , parameter :: spinMinimum=1.0d-6

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%scalarIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the scalar class.
    select type (spin)
    class is (nodeComponentSpinScalar)
       ! Set scale for spin.
       call spin%angularMomentumScale(max(spin%angularMomentum(),spinMinimum*Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterHaloScale_,useBullockDefinition=.true.)))
    end select
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Scalar_Scale_Set

  !![
  <stateStoreTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Scalar_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Scalar_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSpin -> scalar',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Scalar_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Scalar_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Scalar_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSpin -> scalar',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Scalar_State_Restore

end module Node_Component_Halo_Angular_Momentum_Scalar
