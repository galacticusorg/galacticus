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

!% Contains a module which implements a filter which passes only halos that have always been isolated.

  !# <galacticFilter name="galacticFilterHaloAlwaysIsolated">
  !#  <description>A filter which passes only halos that have always been isolated.</description>
  !# </galacticFilter>
  type, extends(galacticFilterClass) :: galacticFilterHaloAlwaysIsolated
     !% A galactic filter class which passes only halos that have always been isolated.
     private
   contains
     procedure :: passes => haloAlwaysIsolatedPasses
  end type galacticFilterHaloAlwaysIsolated

  interface galacticFilterHaloAlwaysIsolated
     !% Constructors for the ``haloAlwaysIsolated'' galactic filter class.
     module procedure haloAlwaysIsolatedConstructorParameters
  end interface galacticFilterHaloAlwaysIsolated

contains

  function haloAlwaysIsolatedConstructorParameters(parameters) result(self)
    !% Constructor for the ``haloAlwaysIsolated'' galactic filter class which takes a parameter set as input.
    use :: Galacticus_Error, only : Galacticus_Component_List        , Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterHaloAlwaysIsolated)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    if (.not.defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumIsGettable())                                                           &
         & call Galacticus_Error_Report                                                                                                         &
         &      (                                                                                                                               &
         &       'requires a merging statistics component that provides a gettable "nodeHierarchyLevelMaximum" property.'                    // &
         &       Galacticus_Component_List(                                                                                                     &
         &                                 'mergingStatistics'                                                                                , &
         &                                  defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumAttributeMatch(requireGettable=.true.)   &
         &                                )                                                                                                  // &
         &       {introspection:location}                                                                                                       &
         &      )
    self=galacticFilterHaloAlwaysIsolated()
    return
  end function haloAlwaysIsolatedConstructorParameters

  logical function haloAlwaysIsolatedPasses(self,node)
    !% Implement a galactic filter which passes only isolated halos.
    use :: Galacticus_Nodes, only : nodeComponentMergingStatistics
    implicit none
    class(galacticFilterHaloAlwaysIsolated), intent(inout)         :: self
    type (treeNode                        ), intent(inout), target :: node
    class(nodeComponentMergingStatistics  ), pointer               :: mergingStatistics
    !$GLC attributes unused :: self

    mergingStatistics        => node             %mergingStatistics        (autoCreate=.true.)
    haloAlwaysIsolatedPasses =  mergingStatistics%nodeHierarchyLevelMaximum(                 ) == 0
    return
  end function haloAlwaysIsolatedPasses
