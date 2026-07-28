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
Contains a module which enumerates the events at which node evolution trajectories may be recorded.
!!}

module Node_Trajectory_Events
  !!{RST
  Enumerates the events at which node evolution trajectories may be recorded by the
  :galacticus-class:`nodeOperatorOutputTrajectory` node operator. The identifier of the event which caused each record to be made
  is available in the output via the :galacticus-class:`nodePropertyExtractorTrajectoryEvent` property extractor.
  !!}
  private

  !![
  <enumeration docformat="rst">
   <name>nodeTrajectoryEvent</name>
   <description>
   Enumeration of the events at which the evolutionary trajectory of a node may be recorded.
   </description>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="nodeInitialize"               description="Recorded when the node is initialized prior to evolution."       />
   <entry label="differentialEvolutionPostStep" description="Recorded after each accepted step of the ODE solver."           />
   <entry label="differentialEvolutionPost"    description="Recorded after the node has been evolved to the end of its step."/>
   <entry label="nodePromote"                  description="Recorded immediately before the node is promoted to its parent." />
   <entry label="nodesMerge"                   description="Recorded when the node merges with another node."                />
   <entry label="galaxiesMerge"                description="Recorded when the galaxy in the node merges with another galaxy."/>
  </enumeration>
  !!]

  ! The identifier of the event responsible for the trajectory record currently being made. This is set by the
  ! :galacticus-class:`nodeOperatorOutputTrajectory` node operator immediately before it triggers output of a node, and is read by
  ! the :galacticus-class:`nodePropertyExtractorTrajectoryEvent` property extractor during that output. Since the extractor is
  ! always called synchronously, on the same thread, from within the output triggered by the operator, a thread-private module
  ! variable suffices here - and avoids having to store a meta-property on every node.
  integer, public :: nodeTrajectoryEventCurrent=0
  !$omp threadprivate(nodeTrajectoryEventCurrent)

end module Node_Trajectory_Events
