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

  !+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a hot halo ram pressure history class which returns the history accumulated by the
  :galacticus-class:`nodeOperatorRamPressureHistory` node operator.
  !!}

  !![
  <hotHaloRamPressureHistory name="hotHaloRamPressureHistoryAccumulated" docformat="rst">
   <description>
   A hot halo ram pressure history class which returns the running maximum, :math:`P`, and time integral, :math:`I`, of the ram
   pressure accumulated along the satellite's orbit by the :galacticus-class:`nodeOperatorRamPressureHistory` node operator,
   which must be included in the model (and which itself requires a satellite component, such as the ``orbiting`` component,
   which tracks the satellite's position and velocity). Any decay of :math:`P` and :math:`I` is controlled by that operator.
   </description>
  </hotHaloRamPressureHistory>
  !!]
  type, extends(hotHaloRamPressureHistoryClass) :: hotHaloRamPressureHistoryAccumulated
     !!{RST
     Implementation of a hot halo ram pressure history class which returns the history accumulated by the
     :galacticus-class:`nodeOperatorRamPressureHistory` node operator.
     !!}
     private
     integer :: ramPressureTimeIntegratedID, ramPressureMaximumID
   contains
     procedure :: forceMaximum        => accumulatedForceMaximum
     procedure :: forceTimeIntegrated => accumulatedForceTimeIntegrated
  end type hotHaloRamPressureHistoryAccumulated

  interface hotHaloRamPressureHistoryAccumulated
     !!{RST
     Constructors for the :galacticus-class:`hotHaloRamPressureHistoryAccumulated` hot halo ram pressure history class.
     !!}
     module procedure accumulatedConstructorParameters
     module procedure accumulatedConstructorInternal
  end interface hotHaloRamPressureHistoryAccumulated

contains

  function accumulatedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`hotHaloRamPressureHistoryAccumulated` hot halo ram pressure history class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(hotHaloRamPressureHistoryAccumulated)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=hotHaloRamPressureHistoryAccumulated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function accumulatedConstructorParameters

  function accumulatedConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`hotHaloRamPressureHistoryAccumulated` hot halo ram pressure history class.
    !!}
    implicit none
    type(hotHaloRamPressureHistoryAccumulated) :: self

    !![
    <addMetaProperty component="satellite" name="ramPressureTimeIntegrated" id="self%ramPressureTimeIntegratedID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="satellite" name="ramPressureMaximum"        id="self%ramPressureMaximumID"        isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function accumulatedConstructorInternal

  double precision function accumulatedForceMaximum(self,node)
    !!{RST
    Return the running maximum of the ram pressure force accumulated by the :galacticus-class:`nodeOperatorRamPressureHistory`
    node operator.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(hotHaloRamPressureHistoryAccumulated), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentSatellite              ), pointer       :: satellite

    if (node%isSatellite()) then
       satellite               => node     %satellite                (                        )
       accumulatedForceMaximum =  satellite%floatRank0MetaPropertyGet(self%ramPressureMaximumID)
    else
       accumulatedForceMaximum =  0.0d0
    end if
    return
  end function accumulatedForceMaximum

  double precision function accumulatedForceTimeIntegrated(self,node)
    !!{RST
    Return the time integral of the ram pressure force accumulated by the :galacticus-class:`nodeOperatorRamPressureHistory` node
    operator.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(hotHaloRamPressureHistoryAccumulated), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentSatellite              ), pointer       :: satellite

    if (node%isSatellite()) then
       satellite                      => node     %satellite                (                               )
       accumulatedForceTimeIntegrated =  satellite%floatRank0MetaPropertyGet(self%ramPressureTimeIntegratedID)
    else
       accumulatedForceTimeIntegrated =  0.0d0
    end if
    return
  end function accumulatedForceTimeIntegrated
