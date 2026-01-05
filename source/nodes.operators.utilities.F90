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
Contains a module of globally-accessible functions supporting the \refClass{nodeOperatorClass} class.
!!}

module Nodes_Operators_Utilities
  !!{
  Provides globally-accessible functions supporting the \refClass{nodeOperatorClass} class.
  !!}
  private
  public :: nodeOperatorConstruct       , nodeOperatorDestruct    , nodeOperatorDeepCopy  , nodeOperatorDeepCopyReset              , &
       &    nodeOperatorDeepCopyFinalize, nodeOperatorStateRestore, nodeOperatorStateStore, nodeOperatorPredeterminedSolveAnalytics

  ! Module-scope pointer to our task object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: nodeOperator__
  !$omp threadprivate(nodeOperator__)

contains

  !![
  <functionGlobal>
   <unitName>nodeOperatorConstruct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout), target  :: parameters   </arguments>
   <arguments>class(*              ), intent(  out), pointer :: nodeOperator_</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorConstruct(parameters,nodeOperator_)
    !!{
    Build a {\normalfont \ttfamily nodeOperator} object from a given parameter set. This is a globally-callable function
    to allow us to subvert the class/module hierarchy.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter   , inputParameters
    use :: Nodes_Operators , only : nodeOperatorClass, nodeOperator
    implicit none
    type (inputParameters), intent(inout), target  :: parameters
    class(*              ), intent(  out), pointer :: nodeOperator_
    type (inputParameters)               , pointer :: parametersCurrent

    parametersCurrent => parameters
    do while (.not.parametersCurrent%isPresent('nodeOperator').and.associated(parametersCurrent%parent))
       parametersCurrent => parametersCurrent%parent
    end do
    if (.not.parametersCurrent%isPresent('nodeOperator')) parametersCurrent => parameters
    nodeOperator__ => nodeOperator(parametersCurrent)
    select type (nodeOperator__)
    class is (nodeOperatorClass) 
       !![
       <referenceCountIncrement object="nodeOperator__"/>
       !!]
       call nodeOperator__%autoHook()
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    nodeOperator_ => nodeOperator__
    return
  end subroutine nodeOperatorConstruct

  !![
  <functionGlobal>
   <unitName>nodeOperatorPredeterminedSolveAnalytics</unitName>
   <type>void</type>
   <module>Galacticus_Nodes, only : treeNode</module>
   <arguments>class           (*       ), intent(inout) :: nodeOperator_</arguments>
   <arguments>type            (treeNode), intent(inout) :: node         </arguments>
   <arguments>double precision          , intent(in   ) :: time         </arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorPredeterminedSolveAnalytics(nodeOperator_,node,time)
    !!{
    Evaluate analytic ppre-determined roperties using a {\normalfont \ttfamily nodeOperator} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error             , only : Error_Report
    use :: Nodes_Operators   , only : nodeOperatorClass
    use :: Galacticus_Nodes  , only : treeNode
    use :: ISO_Varying_String, only : char
    implicit none
    class           (*       ), intent(inout) :: nodeOperator_
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    
    select type (nodeOperator_)
    class is (nodeOperatorClass)
       call nodeOperator_%predeterminedSolveAnalytics(node,time)
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine nodeOperatorPredeterminedSolveAnalytics
  
  !![
  <functionGlobal>
   <unitName>nodeOperatorDestruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: nodeOperator_</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorDestruct(nodeOperator_)
    !!{
    Destruct a {\normalfont \ttfamily taskEvolveForests} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error          , only : Error_Report
    use :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class(*), intent(inout), pointer :: nodeOperator_

    nodeOperator__ => nodeOperator_
    select type (nodeOperator__)
    class is (nodeOperatorClass)
       !![
       <objectDestructor name="nodeOperator__"/>
       !!]
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine nodeOperatorDestruct

  !![
  <functionGlobal>
   <unitName>nodeOperatorStateRestore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self</arguments>
   <arguments>integer          , intent(in   ) :: stateFile</arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile</arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
   </functionGlobal>
  !!]
  subroutine nodeOperatorStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding  , only : c_ptr            , c_size_t
    use            :: Error          , only : Error_Report
    use            :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (nodeOperatorClass)
       call self%stateRestore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine nodeOperatorStateRestore

  !![
  <functionGlobal>
   <unitName>nodeOperatorStateStore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self</arguments>
   <arguments>integer          , intent(in   ) :: stateFile</arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile</arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorStateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding  , only : c_ptr            , c_size_t
    use            :: Error          , only : Error_Report
    use            :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (nodeOperatorClass)
       call self%stateStore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine nodeOperatorStateStore

  !![
  <functionGlobal>
   <unitName>nodeOperatorDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorDeepCopyReset(self)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error          , only : Error_Report
    use :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (nodeOperatorClass)
       call self%deepCopyReset()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine nodeOperatorDeepCopyReset
  
  !![
  <functionGlobal>
   <unitName>nodeOperatorDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorDeepCopyFinalize(self)
    !!{
    Finalize a deep copy of galactic structure objects.
    !!}
    use :: Error          , only : Error_Report
    use :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (nodeOperatorClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine nodeOperatorDeepCopyFinalize
  
  !![
  <functionGlobal>
   <unitName>nodeOperatorDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine nodeOperatorDeepCopy(self,destination)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error          , only : Error_Report
    use :: Nodes_Operators, only : nodeOperatorClass
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (nodeOperatorClass)
       select type (destination)
       class is (nodeOperatorClass)
          call self%deepCopy(destination)
       class default
          call Error_Report("unexpected class"//{introspection:location})
       end select
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine nodeOperatorDeepCopy

end module Nodes_Operators_Utilities
