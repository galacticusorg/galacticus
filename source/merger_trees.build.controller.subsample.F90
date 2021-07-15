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

!!{
Contains a module which implements a merger tree build controller class which performs subsampling of branches.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerSubsample">
   <description>A merger tree build controller class which performs subsampling of branches.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerSubsample
     !!{     
     A merger tree build controller class which performs subsampling of branches. A branch of mass $M$ will be retained with
     probability
     \begin{equation}
       P(M) = \left\{ \begin{array}{ll} 1 & \hbox{if } M \ge M_0 \\ P_0 (M/M_0)^\alpha & \hbox{if } M &lt; M_0, \end{array} \right.
     \end{equation}
     where $M_0=${\normalfont \ttfamily [massThreshold]}, $P_0=${\normalfont \ttfamily [subsamplingRateAtThreshold]} and
     $\alpha=${\normalfont \ttfamily [exponent]}, otherwise being pruned. Node weights are adjusted to account for this pruning.     
     !!}
     private
     double precision :: massThreshold, subsamplingRateAtThreshold, &
          &              exponent
  contains
     procedure :: control => controlSubsample
  end type mergerTreeBuildControllerSubsample

  interface mergerTreeBuildControllerSubsample
     !!{
     Constructors for the ``subsample'' merger tree build controller class.
     !!}
     module procedure subsampleConstructorParameters
     module procedure subsampleConstructorInternal
  end interface mergerTreeBuildControllerSubsample

contains

  function subsampleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``subsample'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildControllerSubsample)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: massThreshold, subsamplingRateAtThreshold, &
         &                                                                 exponent
    
    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The mass threshold, $M_0$, below which subsampling is applied.</description>
    </inputParameter>
    <inputParameter>
      <name>subsamplingRateAtThreshold</name>
      <source>parameters</source>
      <description>The subsampling rate at the mass treshold, $P_0$.</description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <description>The exponent, $\alpha$, of mass in the subsampling probability, i.e. $P(M) = (M/M_0)^\alpha$ for $M &lt; M_0$.</description>
    </inputParameter>
    !!]
    self=mergerTreeBuildControllerSubsample(massThreshold,subsamplingRateAtThreshold,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function subsampleConstructorParameters

  function subsampleConstructorInternal(massThreshold,subsamplingRateAtThreshold,exponent) result(self)
    !!{
    Internal constructor for the ``subsample'' merger tree build controller class.
    !!}
    implicit none
    type            (mergerTreeBuildControllerSubsample)                :: self
    double precision                                    , intent(in   ) :: massThreshold, subsamplingRateAtThreshold, &
         &                                                                 exponent
    !![
    <constructorAssign variables="massThreshold, subsamplingRateAtThreshold, exponent"/>
    !!]
    
    return
  end function subsampleConstructorInternal

  logical function controlSubsample(self,node,treeWalker_)
    !!{
    Subsample branches of a tree under construction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (mergerTreeBuildControllerSubsample), intent(inout)          :: self    
    type            (treeNode                          ), intent(inout), pointer :: node
    class           (mergerTreeWalkerClass             ), intent(inout)          :: treeWalker_
    type            (treeNode                          )               , pointer :: nodeNext       , nodeChild
    class           (nodeComponentBasic                )               , pointer :: basic
    double precision                                                             :: rateSubsampling

    controlSubsample=.true.
    ! Root node is not eligible for pruning.
    if (.not.associated(node%parent)) return
    ! Set the subsampling weight for this node to equal that of its parent.
    call node%subsamplingWeightSet(node%parent%subsamplingWeight())
    ! Primary progenitors are not eligible for pruning.
    if (node%isPrimaryProgenitor()) return
    ! Nodes above the mass threshold are not eligible for pruning.
    basic => node%basic()
    if (basic%mass() >= self%massThreshold) return
    ! Compute subsampling rate, perform sampling.
    rateSubsampling=+self%subsamplingRateAtThreshold &
         &          *(                               &
         &            +basic%mass         ()         &
         &            /self %massThreshold           &
         &           )**self%exponent
    if (node%hostTree%randomNumberGenerator_%uniformSample() < rateSubsampling) then
       ! Node is to be kept - increase its weight to account for corresponding branches which will have been lost.
       call node%subsamplingWeightSet(node%subsamplingWeight()/rateSubsampling)
    else
       ! Prune the node.
       !! Get the next node to walk to in the tree.
       controlSubsample=treeWalker_%next(nodeNext)
       !! Decouple the node from the tree.
       nodeChild => node%parent%firstChild
       do while (.not.associated(nodeChild%sibling,node))
          nodeChild => nodeChild%sibling
       end do
       nodeChild%sibling => node%sibling
       ! Destroy and deallocate the node.
       call node%destroy()
       deallocate(node)
       ! Set the current node to the next node in the tree walk.
       node => nodeNext
    end if
    return
  end function controlSubsample
