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
Implements an ISM mass output analysis property extractor class.
!!}

  use :: Output_Times, only : outputTimes, outputTimesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDescendants">
   <description>
    A node property extractor which extracts the index of the node containing the galaxy to which each current galaxy will
    belong at the next output time (i.e. the \gls{forwardDescendant}). To clarify, this will be the index of the node into
    which the galaxy descends, or the index of a node with which it merges prior to the next output time (and if that node
    merges with another, the index will be of that node and so on).
  
    Note that, to operate correctly, information about which node a given node may merge with (and when this merger will
    happen) must be available. This is typically available in merger trees read from file (i.e. using the ``{\normalfont
    \ttfamily read}'' \refClass{mergerTreeConstructorClass}) providing {\normalfont \ttfamily [presetMergerNodes]} and {\normalfont
    \ttfamily [presetMergerTimes]} are both set to {\normalfont \ttfamily true}. When using randomly assigned satellite orbits
    and merger times, information on when merging occurs does not exist until a node becomes a satellite. Thus, if the node
    becomes a satellite after the current output, but before the next output, there is no way to know which node it will belong
    to at the next output (in such cases, the fallback assumption is no merging).
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorDescendants
     !!{
     A node property extractor descendant indices.
     !!}
     private
     class(outputTimesClass), pointer :: outputTimes_ => null()
   contains
     final     ::                descendantsDestructor
     procedure :: extract     => descendantsExtract
     procedure :: name        => descendantsName
     procedure :: description => descendantsDescription
  end type nodePropertyExtractorDescendants

  interface nodePropertyExtractorDescendants
     !!{
     Constructors for the {\normalfont \ttfamily descendants} output analysis class.
     !!}
     module procedure descendantsConstructorParameters
     module procedure descendantsConstructorInternal
  end interface nodePropertyExtractorDescendants

contains

  function descendantsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily descendants} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorDescendants)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(outputTimesClass                ), pointer       :: outputTimes_

    !![
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDescendants(outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function descendantsConstructorParameters

  function descendantsConstructorInternal(outputTimes_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily descendants} node property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorDescendants)                        :: self
    class(outputTimesClass                ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="*outputTimes_"/>
    !!]

    return
  end function descendantsConstructorInternal

  subroutine descendantsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily descendants} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDescendants), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine descendantsDestructor

  function descendantsExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily descendants} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    implicit none
    integer         (kind_int8                       )                          :: descendantsExtract
    class           (nodePropertyExtractorDescendants), intent(inout)           :: self
    type            (treeNode                        ), intent(inout), target   :: node
    double precision                                  , intent(in   )           :: time
    type            (multiCounter                    ), intent(inout), optional :: instance
    type            (treeNode                        ), pointer                 :: nodeDescendant
    class           (nodeComponentBasic              ), pointer                 :: basic
    class           (nodeComponentSatellite          ), pointer                 :: satellite
    double precision                                                            :: outputTimeNext
    logical                                                                     :: foundDescendant
    !$GLC attributes unused :: self, instance

    satellite       => node%satellite            (    )
    outputTimeNext  =  self%outputTimes_%timeNext(time)
    foundDescendant =  .false.
    if (outputTimeNext < 0.0d0) then
       ! There is no next output time.
       descendantsExtract=node%index()
       foundDescendant=.true.
    else if (node%isSatellite()) then
       ! Node is a satellite, so its node index will remain unchanged.
       if (satellite%timeOfMerging() > outputTimeNext) then
          ! Satellite will not have merged prior to the next output time, so retains its own index.
          descendantsExtract=node%index()
          foundDescendant=.true.
       else
          ! Satellite will merge prior to the next output time - find the node it merges with.
          nodeDescendant => node%mergesWith()
       end if
    else
       ! Node is not a satellite, so set the initial descendant to itself.
       nodeDescendant => node
    end if
    ! Check if we still need to find the descendant.
    if (.not.foundDescendant) then
       ! No descendant has yet been found, so trace forward in time until we find one.
       ! Continue until the tree base is reached, or the next output time is reached.
       do while (.not.foundDescendant)
          ! Get the satellite component.
          satellite => nodeDescendant%satellite()
          basic     => nodeDescendant%basic    ()
          ! If the next output time has been surpassed, then we are finished.
          if (basic%time() >= outputTimeNext) then
             foundDescendant=.true.
          else
             ! Test whether this node is the primary progenitor.
             if (nodeDescendant%isPrimaryProgenitor()) then
                ! It is, so simply move to the parent node.
                nodeDescendant => nodeDescendant%parent
             else
                ! It is not, so it becomes a satellite. Test whether it has a merge target associated with it.
                if (associated(nodeDescendant%mergeTarget)) then
                   ! It does. If merging occurs before the next output time, jump to that node. Otherwise we are finished.
                   if (satellite%timeOfMerging() <= outputTimeNext) then
                      nodeDescendant => nodeDescendant%mergeTarget
                   else
                      foundDescendant=.true.
                   end if
                else
                   ! We no longer can tell if this node will exist as a separate entity at the next output time. Assume that it
                   ! will, and therefore we are finished.
                   foundDescendant=.true.
                end if
             end if
          end if
       end do
       ! If the descendant exists after the next output time, then we have gone one step too far - back up a step.
       if (basic%time() > outputTimeNext .and. associated(nodeDescendant%firstChild)) nodeDescendant => nodeDescendant%firstChild
       ! Store the descendant index.
       descendantsExtract=nodeDescendant%index()
    end if
    return
  end function descendantsExtract


  function descendantsName(self)
    !!{
    Return the name of the descendants property.
    !!}
    implicit none
    type (varying_string                  )                :: descendantsName
    class(nodePropertyExtractorDescendants), intent(inout) :: self
    !$GLC attributes unused :: self

    descendantsName=var_str('descendantIndex')
    return
  end function descendantsName

  function descendantsDescription(self)
    !!{
    Return a description of the descendants property.
    !!}
    implicit none
    type (varying_string                  )                :: descendantsDescription
    class(nodePropertyExtractorDescendants), intent(inout) :: self
    !$GLC attributes unused :: self

    descendantsDescription=var_str('ID of the node which this node will have descended into by the next timestep.')
    return
  end function descendantsDescription
