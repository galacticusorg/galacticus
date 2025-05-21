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
  Implements a non-evolving class for evolving merger trees.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Merger_Tree_Initialization, only : mergerTreeInitializorClass
  
  !![
  <mergerTreeEvolver name="mergerTreeEvolverNonEvolving">
   <description>A non-evolving merger tree evolver.</description>
  </mergerTreeEvolver>
  !!]
  type, extends(mergerTreeEvolverClass) :: mergerTreeEvolverNonEvolving
     !!{
     Implementation of a non-evolving  merger tree evolver.
     !!}
     private
     class  (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class  (mergerTreeInitializorClass), pointer :: mergerTreeInitializor_ => null()
     logical                                      :: pruneTree
   contains
     final     ::           nonEvolvingDestructor
     procedure :: evolve => nonEvolvingEvolve
  end type mergerTreeEvolverNonEvolving

  interface mergerTreeEvolverNonEvolving
     !!{
     Constructors for the \refClass{mergerTreeEvolverNonEvolving} merger tree evolver.
     !!}
     module procedure nonEvolvingConstructorParameters
     module procedure nonEvolvingConstructorInternal
  end interface mergerTreeEvolverNonEvolving

contains

  function nonEvolvingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolverNonEvolving} merger tree evolver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeEvolverNonEvolving)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    logical                                              :: pruneTree
    class  (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class  (mergerTreeInitializorClass  ), pointer       :: mergerTreeInitializor_

    !![
    <inputParameter>
      <name>pruneTree</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, prune the tree to the evolve-to-time after each evolution.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"     name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="mergerTreeInitializor"  name="mergerTreeInitializor_" source="parameters"/>
    !!]
    self=mergerTreeEvolverNonEvolving(pruneTree,cosmologyFunctions_,mergerTreeInitializor_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="mergerTreeInitializor_"/>
    !!]
    return
  end function nonEvolvingConstructorParameters

  function nonEvolvingConstructorInternal(pruneTree,cosmologyFunctions_,mergerTreeInitializor_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolverNonEvolving} merger tree evolver class.
    !!}
    implicit none
    type   (mergerTreeEvolverNonEvolving)                        :: self
    class  (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    class  (mergerTreeInitializorClass  ), intent(in   ), target :: mergerTreeInitializor_
    logical                              , intent(in   )         :: pruneTree
    !![
    <constructorAssign variables="pruneTree, *cosmologyFunctions_, *mergerTreeInitializor_"/>
    !!]

    return
  end function nonEvolvingConstructorInternal

  subroutine nonEvolvingDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolverNonEvolving} merger tree evolver class.
    !!}
    implicit none
    type(mergerTreeEvolverNonEvolving), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeInitializor_"/>
    <objectDestructor name="self%cosmologyFunctions_"   />
    !!]
    return
  end subroutine nonEvolvingDestructor

  subroutine nonEvolvingEvolve(self,tree,timeEnd,treeDidEvolve,suspendTree,deadlockReporting,systemClockMaximum,initializationLock,status)
    !!{
    Evolves all properties of a merger tree to the specified time.
    !!}
    use    :: Error                , only : errorStatusSuccess
    use    :: Merger_Tree_Operators, only : mergerTreeOperatorPruneByTime
    !$ use :: OMP_Lib              , only : OMP_Set_Lock                 , OMP_Unset_Lock, omp_lock_kind
    implicit none
    class           (mergerTreeEvolverNonEvolving )                   , intent(inout) :: self
    integer                                        , optional         , intent(  out) :: status
    type            (mergerTree                   )          , target , intent(inout) :: tree
    double precision                                                  , intent(in   ) :: timeEnd
    logical                                                           , intent(  out) :: treeDidEvolve     , suspendTree
    logical                                                           , intent(in   ) :: deadlockReporting
    integer         (kind_int8                    ), optional         , intent(in   ) :: systemClockMaximum
    integer         (omp_lock_kind                ), optional         , intent(inout) :: initializationLock
    type            (mergerTree                   )          , pointer                :: currentTree
    type            (mergerTreeOperatorPruneByTime)                                   :: pruner
    
    !$GLC attributes unused :: self, deadlockReporting, systemClockMaximum

    if (present(status)) status=errorStatusSuccess
    if (self%pruneTree) pruner=mergerTreeOperatorPruneByTime(timeEnd,0.0d0,huge(0.0d0),self%cosmologyFunctions_)
    suspendTree   =  .false.
    treeDidEvolve =  .true.
    currentTree   => tree
    do while (associated(currentTree))
       if (associated(currentTree%nodeBase)) then
          !$ if (present(initializationLock)) call OMP_Set_Lock  (initializationLock)
          call self%mergerTreeInitializor_%initialize(currentTree,timeEnd)
          !$ if (present(initializationLock)) call OMP_Unset_Lock(initializationLock)
          if (self%pruneTree) call pruner%operatePreEvolution(currentTree)
       end if
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine nonEvolvingEvolve
