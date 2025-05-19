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
Implements a merger tree build controller class that forces tree time steps to exactly coincidence with output times.
!!}

  use :: Cosmological_Density_Field, only : criticalOverdensityClass, cosmologicalMassVarianceClass
  use :: Output_Times              , only : outputTimesClass

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerOutputTimeSnap">
   <description>A merger tree build controller class that forces tree time steps to exactly coincidence with output times.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerOutputTimeSnap
     !!{
     A merger tree build controller class that forces tree time steps to exactly coincidence with output times.
     !!}
     private
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     class(criticalOverdensityClass           ), pointer :: criticalOverdensity_            => null()
     class(cosmologicalMassVarianceClass      ), pointer :: cosmologicalMassVariance_       => null()
     class(outputTimesClass                   ), pointer :: outputTimes_                    => null()
   contains
     final     ::                               outputTimeSnapDestructor
     procedure :: control                    => outputTimeSnapControl
     procedure :: branchingProbabilityObject => outputTimeSnapBranchingProbabilityObject
     procedure :: nodesInserted              => outputTimeSnapNodesInserted
     procedure :: timeMaximum                => outputTimeSnapTimeMaximum
  end type mergerTreeBuildControllerOutputTimeSnap

  interface mergerTreeBuildControllerOutputTimeSnap
     !!{
     Constructors for the ``outputTimeSnap'' merger tree build controller class.
     !!}
     module procedure outputTimeSnapConstructorParameters
     module procedure outputTimeSnapConstructorInternal
  end interface mergerTreeBuildControllerOutputTimeSnap

contains

  function outputTimeSnapConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``outputTimeSnap'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeBuildControllerOutputTimeSnap)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(mergerTreeBranchingProbabilityClass    ), pointer       :: mergerTreeBranchingProbability_
    class(cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass               ), pointer       :: criticalOverdensity_ 
    class(outputTimesClass                       ), pointer       :: outputTimes_

    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="outputTimes"                    name="outputTimes_"                    source="parameters"/>
    !!]
    self=mergerTreeBuildControllerOutputTimeSnap(mergerTreeBranchingProbability_,outputTimes_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="outputTimes_"                   />
    !!]
    return
  end function outputTimeSnapConstructorParameters

  function outputTimeSnapConstructorInternal(mergerTreeBranchingProbability_,outputTimes_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the ``outputTimeSnap'' merger tree build controller class .
    !!}
    implicit none
    type (mergerTreeBuildControllerOutputTimeSnap)                        :: self
    class(mergerTreeBranchingProbabilityClass    ), intent(in   ), target :: mergerTreeBranchingProbability_
    class(cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_
    class(criticalOverdensityClass               ), intent(in   ), target :: criticalOverdensity_
    class(outputTimesClass                       ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="*mergerTreeBranchingProbability_, *criticalOverdensity_, *cosmologicalMassVariance_, *outputTimes_"/>
    !!]
    
    return
  end function outputTimeSnapConstructorInternal

  subroutine outputTimeSnapDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily outputTimeSnap} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerOutputTimeSnap), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%outputTimes_"                   />
    !!]
    return
  end subroutine outputTimeSnapDestructor

  logical function outputTimeSnapControl(self,node,treeWalker_) result(control)
    !!{
    Apply control to merger tree building.
    !!}
    implicit none
    class(mergerTreeBuildControllerOutputTimeSnap), intent(inout)           :: self
    type (treeNode                               ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass                  ), intent(inout), optional :: treeWalker_
    !$GLC attributes unused :: self, node, treeWalker_

    ! Always return true as we never want to halt tree building.
    control=.true.
    return
  end function outputTimeSnapControl

  function outputTimeSnapBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass    ), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerOutputTimeSnap), intent(inout) :: self
    type (treeNode                               ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function outputTimeSnapBranchingProbabilityObject

  subroutine outputTimeSnapNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerOutputTimeSnap), intent(inout)           :: self
    type   (treeNode                               ), intent(inout)           :: nodeCurrent     , nodeProgenitor1
    type   (treeNode                               ), intent(inout), optional :: nodeProgenitor2
    logical                                         , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine outputTimeSnapNodesInserted

  double precision function outputTimeSnapTimeMaximum(self,node,massBranch,criticalOverdensityBranch,timeReference,insertNode) result(timeMaximum)
    !!{
    Return the maximum allowed time for this node.
    !!}
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (mergerTreeBuildControllerOutputTimeSnap), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    double precision                                         , intent(in   ) :: massBranch   , criticalOverdensityBranch, &
         &                                                                      timeReference
    logical                                                  , intent(  out) :: insertNode
    double precision                                                         :: time         , timeNext
    integer         (c_size_t                               )                :: indexOutput

    insertNode =.true.
    time       =self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=criticalOverdensityBranch,mass       =massBranch ,node=node)
    timeNext   =self%outputTimes_        %timePrevious  (timeCurrent        =time                     ,indexOutput=indexOutput          )
    timeMaximum=huge(0.0d0)
    ! Return if there is no earlier time.
    if (timeNext <= 0.0d0) return
    ! Check for times very close to an output time - if we are at such a time most likely we should be precisely at this time, and
    ! this is just a rounding error. In such cases move to the next earlier time (if such exists).
    if (Values_Agree(time,timeNext,relTol=1.0d-8)) then
       ! If there is no earlier time, return (with huge time maximum).
       if (indexOutput == 1_c_size_t) return
       ! Set the time to the next earliest time.
       timeNext=self%outputTimes_%time(indexOutput-1_c_size_t)
    end if
    ! Compute the w parameter corresponding to this time.
    timeMaximum=+self%criticalOverdensity_     %value       (time=timeNext     ,mass=massBranch,node=node) &
         &      *self%cosmologicalMassVariance_%rootVariance(time=timeReference,mass=massBranch          ) &
         &      /self%cosmologicalMassVariance_%rootVariance(time=timeNext     ,mass=massBranch          )
    ! If the maximum time has already been reached, return an infinite maximum time instead.
    if (timeMaximum <= criticalOverdensityBranch) timeMaximum=huge(0.0d0)
    return
  end function outputTimeSnapTimeMaximum
