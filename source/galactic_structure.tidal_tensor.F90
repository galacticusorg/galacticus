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

!+    Contributions to this file made by: Andrew Benson.

!% Contains a module which implements calculations of gravitational tidal tensors.

module Galactic_Structure_Tidal_Tensors
  !% Implements calculations of the gravitational tidal tensors.
  use :: Kind_Numbers, only : kind_int8
  implicit none
  private
  public :: Galactic_Structure_Tidal_Tensor

  ! Module scope variables used in mapping over components.
  integer                                        :: componentTypeShared    , massTypeShared
  double precision                , dimension(3) :: positionCartesianShared
  !$omp threadprivate(massTypeShared,componentTypeShared,positionCartesianShared)

contains

  function Galactic_Structure_Tidal_Tensor(node,positionCartesian,componentType,massType)
    !% Compute the gravitational tidal tensor at a given position.
    use :: Galactic_Structure_Options, only : componentTypeAll               , massTypeAll
    use :: Galacticus_Nodes          , only : optimizeForTidalTensorSummation, reductionSummation, treeNode
    use :: Tensors                   , only : tensorRank2Dimension3Symmetric , assignment(=)
    !# <include directive="tidalTensorTask" type="moduleUse">
    include 'galactic_structure.tidal_tensor.tasks.modules.inc'
    !# </include>
    implicit none
    type            (tensorRank2Dimension3Symmetric)                              :: Galactic_Structure_Tidal_Tensor
    type            (treeNode                      ), intent(inout)               :: node
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    integer                                         , intent(in   ), optional     :: componentType                    , massType
    procedure       (Component_Tidal_Tensor        ), pointer                     :: componentTidalTensorFunction
    type            (tensorRank2Dimension3Symmetric)                              :: componentTidalTensor

    ! Copy the position to module-scope.
    positionCartesianShared=positionCartesian
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Initialize pointer to function that supplies the tidal tensor for all components.
    componentTidalTensorFunction => Component_Tidal_Tensor
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared     =massType
    else
       massTypeShared     =massTypeAll
    end if
    ! Compute the tidal tensor.
    Galactic_Structure_Tidal_Tensor=node%mapTensorR2D3(componentTidalTensorFunction,reductionSummation,optimizeFor=optimizeForTidalTensorSummation)
    !# <include directive="tidalTensorTask" type="functionCall" functionType="function" returnParameter="componentTidalTensor">
    !#  <functionArgs>node,positionCartesianShared,componentTypeShared,massTypeShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Tidal_Tensor=Galactic_Structure_Tidal_Tensor+componentTidalTensor</onReturn>
    include 'galactic_structure.tidal_tensor.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Tidal_Tensor

  function Component_Tidal_Tensor(component)
    !% Function returning the tidal tensor in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    use :: Tensors         , only : tensorRank2Dimension3Symmetric
    implicit none
    class(nodeComponent                 ), intent(inout) :: component
    type (tensorRank2Dimension3Symmetric)                :: Component_Tidal_Tensor

    Component_Tidal_Tensor=component%tidalTensor(positionCartesianShared,componentTypeShared,massTypeShared)
    return
  end function Component_Tidal_Tensor

end module Galactic_Structure_Tidal_Tensors












