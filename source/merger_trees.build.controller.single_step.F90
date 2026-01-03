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
Implements a merger tree build controller class which limits tree building to a single step.
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
     logical                                                   :: haltAfterStep
     integer                                                   :: primaryLabelID                      , secondaryLabelID       , &
          &                                                       smoothLabelID                       , singleLabelID
   contains
     final     ::                               singleStepDestructor
     procedure :: control                    => singleStepControl
     procedure :: timeMinimum                => singleStepTimeMinimum
     procedure :: timeMaximum                => singleStepTimeMaximum
     procedure :: controlTimeMaximum         => singleStepControlTimeMaximum
     procedure :: branchingProbabilityObject => singleStepBranchingProbabilityObject
     procedure :: nodesInserted              => singleStepNodesInserted
  end type mergerTreeBuildControllerSingleStep

  interface mergerTreeBuildControllerSingleStep
     !!{
     Constructors for the \refClass{mergerTreeBuildControllerSingleStep} merger tree build controller class.
     !!}
     module procedure singleStepConstructorParameters
     module procedure singleStepConstructorInternal
  end interface mergerTreeBuildControllerSingleStep

contains

  function singleStepConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildControllerSingleStep} merger tree build controller class which takes a parameter set as input.
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
    logical                                                              :: haltAfterStep
 
    !![
    <inputParameter>
      <name>redshiftStep</name>
      <source>parameters</source>
      <description>The redshift to which to take a single step.</description>
    </inputParameter>
    <inputParameter>
      <name>haltAfterStep</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, cease building the tree after the first step. Otherwise, continue to build with no further step limitations.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="linearGrowth"              name="linearGrowth_"              source="parameters"/>
    <objectBuilder class="criticalOverdensity"       name="criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="mergerTreeBuildController" name="mergerTreeBuildController_" source="parameters"/>
    !!]
    timeStep               = cosmologyFunctions_ %cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftStep))
    criticalOverdensityStep=+criticalOverdensity_%value     (                                                    timeStep ) &
            &               /linearGrowth_       %value     (                                                    timeStep )
    self                   = mergerTreeBuildControllerSingleStep(criticalOverdensityStep,haltAfterStep,cosmologyFunctions_,criticalOverdensity_,linearGrowth_,mergerTreeBuildController_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="linearGrowth_"             />
    <objectDestructor name="criticalOverdensity_"      />
    <objectDestructor name="mergerTreeBuildController_"/>
    !!]
    return
  end function singleStepConstructorParameters

  function singleStepConstructorInternal(criticalOverdensityStep,haltAfterStep,cosmologyFunctions_,criticalOverdensity_,linearGrowth_,mergerTreeBuildController_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildControllerSingleStep} merger tree build controller class .
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type            (mergerTreeBuildControllerSingleStep)                        :: self
    double precision                                     , intent(in   )         :: criticalOverdensityStep
    logical                                              , intent(in   )         :: haltAfterStep
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                  ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (mergerTreeBuildControllerClass     ), intent(in   ), target :: mergerTreeBuildController_
    !![
    <constructorAssign variables="criticalOverdensityStep, haltAfterStep, *cosmologyFunctions_, *criticalOverdensity_, *linearGrowth_, *mergerTreeBuildController_"/>
    !!]

    self%redshiftStep    =self%cosmologyFunctions_ %redshiftFromExpansionFactor(                          &
         &                self%cosmologyFunctions_ %expansionFactor             (                         &
         &                self%criticalOverdensity_%timeOfCollapse               (                        &
         &                                                                        criticalOverdensityStep &
         &                                                                       )                        &
         &                                                                      )                         &
         &                                                                     )
    self%  primaryLabelID=nodeLabelRegister('progenitorFirst' ,'Identifies progenitors that were sampled first from the progenitor mass distribution.'                           )
    self%secondaryLabelID=nodeLabelRegister('progenitorSecond','Identifies progenitors that were sampled second from the progenitor mass distribution.'                          )
    self%   smoothLabelID=nodeLabelRegister('progenitorSmooth','Identifies progenitors resulting from smooth/sub-resolution accretion.'                                          )
    self%   singleLabelID=nodeLabelRegister('progenitorSingle','Identifies progenitors on a single branch all of which were sampled first from the progenitor mass distribution.')
    return
  end function singleStepConstructorInternal

  subroutine singleStepDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildControllerSingleStep} merger tree build controller class.
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
    use :: Nodes_Labels, only : nodeLabelSet
    implicit none
    class(mergerTreeBuildControllerSingleStep), intent(inout)           :: self
    type (treeNode                           ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass              ), intent(inout), optional :: treeWalker_
    
    ! Mark the root halo as being on the single branch.
    if (associated(node) .and. .not.associated(node%parent)) &
         & call nodeLabelSet(self%singleLabelID,node)
    ! First call the decorated controller.
    singleStepControl=self%mergerTreeBuildController_%control(node,treeWalker_)
    ! If we are to halt after the first step, then override the decorated controller as necessary.
    if (self%haltAfterStep .and. associated(node) .and. associated(node%parent)) singleStepControl=.false.
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

    if (associated(node%parent)) then
       singleStepTimeMinimum=self%mergerTreeBuildController_%timeMinimum(node,massBranch,criticalOverdensityBranch)
    else
       singleStepTimeMinimum=self%criticalOverdensityStep
    end if
    return
  end function singleStepTimeMinimum
  
  double precision function singleStepTimeMaximum(self,node,massBranch,criticalOverdensityBranch,timeReference,insertNode)
    !!{
    Return the maximum allowed time for this node.
    !!}
    implicit none
    class           (mergerTreeBuildControllerSingleStep), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: massBranch   , criticalOverdensityBranch, &
         &                                                                  timeReference
    logical                                              , intent(  out) :: insertNode

    if (associated(node%parent)) then
       singleStepTimeMaximum=self%mergerTreeBuildController_%timeMaximum(node,massBranch,criticalOverdensityBranch,timeReference,insertNode)
    else
       singleStepTimeMaximum=self%criticalOverdensityStep
    end if
    insertNode=.false.
    return
  end function singleStepTimeMaximum
  
  logical function singleStepControlTimeMaximum(self,node,massBranch,criticalOverdensityBranch,nodeIndex)
    !!{
    Control when the maximum time is reached.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: Kind_Numbers    , only : kind_int8
    use :: Nodes_Labels    , only : nodeLabelSet      , nodeLabelIsPresent
    implicit none
    class           (mergerTreeBuildControllerSingleStep), intent(inout)         :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: massBranch, criticalOverdensityBranch
    integer         (kind=kind_int8                     ), intent(inout)         :: nodeIndex
    type            (treeNode                           ), pointer               :: nodeNew
    class           (nodeComponentBasic                 ), pointer               :: basic

    ! After the first step, allow our decorated controller to act on this event.
    if (associated(node%parent)) then
       singleStepControlTimeMaximum=self%mergerTreeBuildController_%controlTimeMaximum(node,massBranch,criticalOverdensityBranch,nodeIndex)
       return
    end if
    ! The single step time has been reached. If a progenitor was inserted, we have nothing more to do here.
    if (associated(node%firstChild)) then
       singleStepControlTimeMaximum=.false.
       return
    end if
    ! No progenitors were inserted after the first step. We must insert one here to avoid getting stuck in an infinite loop.
    nodeIndex            =  nodeIndex+1
    nodeNew              => treeNode(nodeIndex,node%hostTree)
    basic                => nodeNew%basic(autoCreate=.true.)
    nodeNew  %parent     => node
    node     %firstChild => nodeNew
    nodeNew  %sibling    => null()
    call basic%massSet(massBranch               )
    call basic%timeSet(criticalOverdensityBranch)
    call nodeLabelSet(self%smoothLabelID,nodeNew)
    if (nodeLabelIsPresent(self%singleLabelID,node)) &
         & call nodeLabelSet(self%singleLabelID,nodeNew)
    ! Return false indicating that the current node is finished, so building should continue from its progenitor nodes.
    singleStepControlTimeMaximum=.false.
    return
  end function singleStepControlTimeMaximum
  
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

  subroutine singleStepNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    use :: Nodes_Labels, only : nodeLabelSet, nodeLabelIsPresent
    implicit none
    class  (mergerTreeBuildControllerSingleStep), intent(inout)           :: self
    type   (treeNode                           ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                           ), intent(inout), optional :: nodeProgenitor2
    logical                                     , intent(in   ), optional :: didBranch
    !![
    <optionalArgument name="didBranch" defaultsTo=".false."/>
    !!]
    
    call self%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeProgenitor1,nodeProgenitor2)
    if (didBranch_) then
       call        nodeLabelSet(self%  primaryLabelID,nodeProgenitor1)
       if (present(nodeProgenitor2))                                   &
            & call nodeLabelSet(self%secondaryLabelID,nodeProgenitor2)
    else
       call        nodeLabelSet(self%   smoothLabelID,nodeProgenitor1)  
    end if
    if (nodeLabelIsPresent(self%singleLabelID,nodeCurrent)) &
         & call nodeLabelSet(self%singleLabelID,nodeProgenitor1)    
    return
  end subroutine singleStepNodesInserted
