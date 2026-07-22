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

!!{RST
Implements a node evolution trajectory event property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTrajectoryEvent" docformat="rst">
   <description>
   Extracts the identifier of the event which caused the current record of a node's evolutionary trajectory to be made by the
   :galacticus-class:`nodeOperatorOutputTrajectory` node operator. This is meaningful only in trajectory output---in a regular
   output it simply reports whichever event caused the most recent trajectory record to be made on the current thread (or zero if
   no trajectory record has been made at all). The mapping from identifier to event name is given in the description of the
   ``trajectoryEvent`` dataset.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorTrajectoryEvent
     !!{RST
     A node evolution trajectory event property extractor.
     !!}
     private
   contains
     procedure :: extract     => trajectoryEventExtract
     procedure :: name        => trajectoryEventName
     procedure :: description => trajectoryEventDescription
  end type nodePropertyExtractorTrajectoryEvent

  interface nodePropertyExtractorTrajectoryEvent
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorTrajectoryEvent` property extractor class.
     !!}
     module procedure trajectoryEventConstructorParameters
  end interface nodePropertyExtractorTrajectoryEvent

contains

  function trajectoryEventConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorTrajectoryEvent` property extractor class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorTrajectoryEvent)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorTrajectoryEvent()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function trajectoryEventConstructorParameters

  function trajectoryEventExtract(self,node,time,instance)
    !!{RST
    Implement a ``trajectoryEvent`` node property extractor.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventCurrent
    implicit none
    integer         (kind_int8                           )                          :: trajectoryEventExtract
    class           (nodePropertyExtractorTrajectoryEvent), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: time
    type            (multiCounter                        ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, node, instance, time

    trajectoryEventExtract=int(nodeTrajectoryEventCurrent,kind=kind_int8)
    return
  end function trajectoryEventExtract

  function trajectoryEventName(self)
    !!{RST
    Return the name of the trajectory event property.
    !!}
    implicit none
    type (varying_string                      )                :: trajectoryEventName
    class(nodePropertyExtractorTrajectoryEvent), intent(inout) :: self
    !$GLC attributes unused :: self

    trajectoryEventName=var_str('trajectoryEvent')
    return
  end function trajectoryEventName

  function trajectoryEventDescription(self)
    !!{RST
    Return a description of the trajectory event property.
    !!}
    use :: Node_Trajectory_Events, only : enumerationNodeTrajectoryEventDecode, nodeTrajectoryEventMin, nodeTrajectoryEventMax
    use :: ISO_Varying_String    , only : varying_string                      , var_str               , assignment(=)         , &
         &                                operator(//)
    use :: String_Handling       , only : operator(//)
    implicit none
    type   (varying_string                      )                :: trajectoryEventDescription
    class  (nodePropertyExtractorTrajectoryEvent), intent(inout) :: self
    integer                                                      :: event
    !$GLC attributes unused :: self

    trajectoryEventDescription=var_str('Identifier of the event at which this trajectory record was made:')
    do event=nodeTrajectoryEventMin,nodeTrajectoryEventMax
       trajectoryEventDescription=trajectoryEventDescription//' '//event//'='//enumerationNodeTrajectoryEventDecode(event)//';'
    end do
    return
  end function trajectoryEventDescription
