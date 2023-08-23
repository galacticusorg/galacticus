!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a merger tree build controller class which limits tree building to a single step.
!!}
  
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass 

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerSingleStep">
   <description>A merger tree build controller class which limits tree building to a single step.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerSingleStep
     !!{
     A merger tree build controller class which limits tree building to a single step.
     !!}
     private
     class           (cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_        => null()
     class           (mergerTreeBuildControllerClass), pointer :: mergerTreeBuildController_ => null()
     class           (linearGrowthClass             ), pointer :: linearGrowth_              => null()
     class           (criticalOverdensityClass      ), pointer :: criticalOverdensity_       => null()
     double precision                                          :: redshiftStep                        , criticalOverdensityStep
   contains
     final     ::                               singleStepDestructor
     procedure :: control                    => singleStepControl
     procedure :: timeMinimum                => singleStepTimeMinimum
     procedure :: timeMaximum                => singleStepTimeMaximum
     procedure :: branchingProbabilityObject => singleStepBranchingProbabilityObject
     procedure :: nodesInserted              => singleStepNodesInserted
  end type mergerTreeBuildControllerSingleStep

  interface mergerTreeBuildControllerSingleStep
     !!{
     Constructors for the ``singleStep'' merger tree build controller class.
     !!}
     module procedure singleStepConstructorParameters
     module procedure singleStepConstructorInternal
  end interface mergerTreeBuildControllerSingleStep

contains

  function singleStepConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``singleStep'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildControllerSingleStep)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                  ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass           ), pointer       :: criticalOverdensity_ 
    class           (mergerTreeBuildControllerClass     ), pointer       :: mergerTreeBuildController_
    double precision                                                     :: redshiftStep              , criticalOverdensityStep, &
         &                                                                  timeStep
    
    !![
    <inputParameter>
      <name>redshiftStep</name>
      <source>parameters</source>
      <description>The redshift to which to take a single step.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="linearGrowth"              name="linearGrowth_"              source="parameters"/>
    <objectBuilder class="criticalOverdensity"       name="criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="mergerTreeBuildController" name="mergerTreeBuildController_" source="parameters"/>
    !!]
    timeStep               = cosmologyFunctions_ %cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftStep))
    criticalOverdensityStep=+criticalOverdensity_%value     (                                                    timeStep ) &
            &               /linearGrowth_       %value     (                                                    timeStep )
    self                   = mergerTreeBuildControllerSingleStep(criticalOverdensityStep,cosmologyFunctions_,criticalOverdensity_,linearGrowth_,mergerTreeBuildController_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="linearGrowth_"             />
    <objectDestructor name="criticalOverdensity_"      />
    <objectDestructor name="mergerTreeBuildController_"/>
    !!]
    return
  end function singleStepConstructorParameters

  function singleStepConstructorInternal(criticalOverdensityStep,cosmologyFunctions_,criticalOverdensity_,linearGrowth_,mergerTreeBuildController_) result(self)
    !!{
    Internal constructor for the ``singleStep'' merger tree build controller class .
    !!}
    implicit none
    type            (mergerTreeBuildControllerSingleStep)                        :: self
    double precision                                     , intent(in   )         :: criticalOverdensityStep
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                  ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (mergerTreeBuildControllerClass     ), intent(in   ), target :: mergerTreeBuildController_
    !![
    <constructorAssign variables="criticalOverdensityStep, *cosmologyFunctions_, *criticalOverdensity_, *linearGrowth_, *mergerTreeBuildController_"/>
    !!]

    self%redshiftStep=self%cosmologyFunctions_ %redshiftFromExpansionFactor(                          &
         &            self%cosmologyFunctions_ %expansionFactor             (                         &
         &            self%criticalOverdensity_%timeOfCollapse               (                        &
         &                                                                    criticalOverdensityStep &
         &                                                                   )                        &
         &                                                                  )                         &
         &                                                                 )
    return
  end function singleStepConstructorInternal

  subroutine singleStepDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily singleStep} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerSingleStep), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"       />
    <objectDestructor name="self%linearGrowth_"             />
    <objectDestructor name="self%criticalOverdensity_"      />
    <objectDestructor name="self%mergerTreeBuildController_"/>
    !!]
    return
  end subroutine singleStepDestructor

  logical function singleStepControl(self,node,treeWalker_)
    !!{
    Apply control to merger tree building.
    !!}
    implicit none
    class(mergerTreeBuildControllerSingleStep), intent(inout)           :: self
    type (treeNode                           ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass              ), intent(inout), optional :: treeWalker_
    !$GLC attributes unused :: self, treeWalker_

    singleStepControl=.not.associated(node%parent)
    return
  end function singleStepControl

  double precision function singleStepTimeMinimum(self,node,massBranch,criticalOverdensityBranch)
    !!{
    Return the maximum allowed time for this node.
    !!}
    implicit none
    class           (mergerTreeBuildControllerSingleStep), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: massBranch, criticalOverdensityBranch

    singleStepTimeMinimum=self%criticalOverdensityStep
    return
  end function singleStepTimeMinimum
  
  double precision function singleStepTimeMaximum(self,node,massBranch,criticalOverdensityBranch)
    !!{
    Return the maximum allowed time for this node.
    !!}
    implicit none
    class           (mergerTreeBuildControllerSingleStep), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: massBranch, criticalOverdensityBranch

    singleStepTimeMaximum=self%criticalOverdensityStep
    return
  end function singleStepTimeMaximum
  
  function singleStepBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerSingleStep), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBuildController_%branchingProbabilityObject(node)
    return
  end function singleStepBranchingProbabilityObject

  subroutine singleStepNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class(mergerTreeBuildControllerSingleStep), intent(inout)           :: self
    type (treeNode                           ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type (treeNode                           ), intent(inout), optional :: nodeProgenitor2

    call self%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeProgenitor1,nodeProgenitor2)
    return
  end subroutine singleStepNodesInserted
