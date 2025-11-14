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
Contains a module which implements the a vector angular momentum component for halos.
!!}

module Node_Component_Halo_Angular_Momentum_Vector
  !!{
  Implements the vector spin component.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize, Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize, &
       &    Node_Component_Halo_Angular_Momentum_Vector_Scale_Set        , Node_Component_Halo_Angular_Momentum_Vector_Initialize         , &
       &    Node_Component_Halo_Angular_Momentum_Vector_State_Store      , Node_Component_Halo_Angular_Momentum_Vector_State_Restore
  
  !![
  <component>
   <class>spin</class>
   <name>vector</name>
   <extends>
    <class>spin</class>
    <name>scalar</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>angularMomentumVector</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum vector of the DMO halo."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>angularMomentumVectorGrowthRate</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>angularMomentumGrowthRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" />
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  !$omp threadprivate(darkMatterHaloScale_)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Initialize(parameters)
    !!{
    Initializes the tree node vector halo angular momentum methods module.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent, nodeComponentSpinVector
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters        ), intent(inout) :: parameters
    type(nodeComponentSpinVector)                :: spinVectorComponent
    !$GLC attributes unused :: parameters

    if (defaultSpinComponent%vectorIsActive())                                                      &
         & call spinVectorComponent%angularMomentumGrowthRateFunction(angularMomentumGrowthRateGet)
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize(parameters)
    !!{
    Initializes the halo vector angular momentum module.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    
    if (defaultSpinComponent%vectorIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
       !!]
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize()
    !!{
    Uninitializes the tree node preset spin module.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%vectorIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize
 
  !![
  <scaleSetTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Scale_Set(node)
    !!{
    Set scales for properties in the preset implementation of the spin component.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , nodeComponentSpinVector, treeNode, defaultSpinComponent
    implicit none
    type            (treeNode         ), intent(inout), pointer   :: node
    class           (nodeComponentSpin)               , pointer   :: spin
    double precision                                  , parameter :: spinMinimum=1.0d-6

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vectorIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the preset class.
    select type (spin)
    class is (nodeComponentSpinVector)
       ! Set scale for spin.
       call spin%angularMomentumVectorScale([1.0d0,1.0d0,1.0d0]*max(spin%angularMomentum(),spinMinimum*Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterHaloScale_,useBullockDefinition=.true.)))
    end select
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Scale_Set

  double precision function angularMomentumGrowthRateGet(self)
    !!{
    Compute the rate of growth of the magnitude of halo angular momentum.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpinVector
    class(nodeComponentSpinVector), intent(inout) :: self

    if     (                                                      &
         &   all(self%angularMomentumVectorGrowthRate() == 0.0d0) &
         &  .or.                                                  &
         &   all(self%angularMomentumVector          () == 0.0d0) &
         & ) then
       angularMomentumGrowthRateGet=+0.0d0
    else
       angularMomentumGrowthRateGet=+sum(                                        &
            &                            +self%angularMomentumVector          () &
            &                            *self%angularMomentumVectorGrowthRate() &
            &                           )                                        &
            &                       /     self%angularMomentum                ()
    end if
    return
  end function angularMomentumGrowthRateGet

  !![
  <stateStoreTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSpin -> vector',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSpin -> vector',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_State_Restore

end module Node_Component_Halo_Angular_Momentum_Vector
