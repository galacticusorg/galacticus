!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a simple time-stepping criterion for merger tree evolution.

module Merger_Tree_Timesteps_Simple
  !% Implements a simple time-stepping criterion for merger tree evolution.
  implicit none
  private
  public :: Merger_Tree_Timestep_Simple

  ! Variable inidicating if module is initialized.
  logical          :: timestepSimpleInitialized=.false.

  ! Parameters controlling the size of timesteps.
  double precision :: timestepSimpleAbsolute           , timestepSimpleRelative

contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_Simple</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_Simple(node,timeStep,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Determine a suitable timestep for {\normalfont \ttfamily node} using the simple method. This simply selects the smaller of {\normalfont \ttfamily
    !% timestepSimpleAbsolute} and {\normalfont \ttfamily timestepSimpleRelative}$H^{-1}(t)$.
    use Galacticus_Nodes
    use Input_Parameters
    use Cosmology_Functions
    use Evolve_To_Time_Reports
    use ISO_Varying_String
    implicit none
    type            (treeNode               ), intent(inout)          , pointer :: node
    procedure       (                       ), intent(inout)          , pointer :: End_Of_Timestep_Task
    double precision                         , intent(inout)                    :: timeStep
    logical                                  , intent(in   )                    :: report
    type            (treeNode               ), intent(inout), optional, pointer :: lockNode
    type            (varying_string         ), intent(inout), optional          :: lockType
    class           (nodeComponentBasic     )                         , pointer :: basic
    class           (cosmologyFunctionsClass)                         , pointer :: cosmologyFunctions_
    double precision                                                            :: expansionFactor      , expansionTimescale, &
         &                                                                         ourTimeStep          , time
    !GCC$ attributes unused :: End_Of_Timestep_Task
    
    if (.not.timestepSimpleInitialized) then
       !$omp critical (timestepSimpleInitialize)
       if (.not.timestepSimpleInitialized) then
          !# <inputParameter>
          !#   <name>timestepSimpleRelative</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.1d0</defaultValue>
          !#   <description>The maximum allowed relative change in time for a single step in the evolution of a node.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>timestepSimpleAbsolute</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>1.0d0</defaultValue>
          !#   <description>The maximum allowed absolute change in time (in Gyr) for a single step in the evolution of a node.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          timestepSimpleInitialized=.true.
       end if
       !$omp end critical (timestepSimpleInitialize)
    end if

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()

    ! Get current cosmic time.
    basic => node%basic()
    time=basic%time()
    
    ! Find current expansion timescale.
    if (timestepSimpleRelative > 0.0d0) then
       expansionFactor   =      cosmologyFunctions_%expansionFactor(           time)
       expansionTimescale=1.0d0/cosmologyFunctions_%expansionRate  (expansionFactor)
       ourTimeStep=min(timestepSimpleRelative*expansionTimescale,timestepSimpleAbsolute)
    else
       ourTimeStep=                                              timestepSimpleAbsolute
    end if
       
    ! Set return value if our timestep is smaller than current one.
    if (ourTimeStep < timeStep) then
       if (present(lockNode)) lockNode => node
       if (present(lockType)) lockType =  "simple"
       timeStep=ourTimeStep
    end if

    if (report) call Evolve_To_Time_Report("simple: ",timeStep)
    return
  end subroutine Merger_Tree_Timestep_Simple

end module Merger_Tree_Timesteps_Simple
