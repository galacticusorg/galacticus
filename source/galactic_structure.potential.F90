!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of gravitationl potential.

module Galactic_Structure_Potentials
  !% Implements calculations of the gravitational potential.
  use :: Kind_Numbers, only : kind_int8
  implicit none
  private
  public :: Galactic_Structure_Potential, Galactic_Structure_Potential_Standard_Reset

  ! Module scope variables used in mapping over components.
  integer                          :: componentTypeShared, massTypeShared, statusShared
  double precision                 :: radiusShared
  !$omp threadprivate(massTypeShared,componentTypeShared,radiusShared,statusShared)

  ! Precomputed values.
  integer         (kind=kind_int8) :: lastUniqueID           =-1_kind_int8
  logical                          :: potentialOffsetComputed=.false.
  double precision                 :: potentialOffset
    !$omp threadprivate(lastUniqueID,potentialOffsetComputed,potentialOffset)

contains

  double precision function Galactic_Structure_Potential(node,radius,componentType,massType,status)
    !% Solve for the gravitational potential at a given radius. Assumes the galactic structure has already been computed.
    use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale          , darkMatterHaloScaleClass
    use :: Galactic_Structure_Options, only : componentTypeAll             , massTypeAll             , structureErrorCodeSuccess
    use :: Galacticus_Nodes          , only : optimizeForPotentialSummation, reductionSummation      , treeNode
    !# <include directive="potentialTask" type="moduleUse">
    include 'galactic_structure.potential.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode                ), intent(inout)                    :: node
    integer                                   , intent(in   ), optional          :: componentType             , massType
    double precision                          , intent(in   )                    :: radius
    integer                                   , intent(  out), optional          :: status
    procedure       (Component_Potential     )                         , pointer :: componentPotentialFunction
    class           (darkMatterHaloScaleClass)                         , pointer :: darkMatterHaloScale_
    double precision                                                             :: componentPotential

    ! Initialize status.
    if (present(status)) status=structureErrorCodeSuccess
    ! Initialize pointer to function that supplies the potential for all components.
    componentPotentialFunction => Component_Potential
    ! Reset calculations if this is a new node.
    if (node%uniqueID() /= lastUniqueID) call Galactic_Structure_Potential_Standard_Reset(node)
    ! Evaluate the potential at the halo virial radius.
    if (.not.potentialOffsetComputed) then
       componentTypeShared  =  componentTypeAll
       massTypeShared       =  massTypeAll
       darkMatterHaloScale_ => darkMatterHaloScale()
       radiusShared         =  darkMatterHaloScale_%virialRadius(node)
       statusShared         =  structureErrorCodeSuccess
       Galactic_Structure_Potential=node%mapDouble0(componentPotentialFunction,reductionSummation,optimizeFor=optimizeForPotentialSummation)
       if (statusShared /= structureErrorCodeSuccess) status=statusShared
       !# <include directive="potentialTask" type="functionCall" functionType="function" returnParameter="componentPotential">
       !#  <functionArgs>node,radiusShared,componentTypeShared,massTypeShared,status</functionArgs>
       !#  <onReturn>Galactic_Structure_Potential=Galactic_Structure_Potential+componentPotential</onReturn>
       include 'galactic_structure.potential.tasks.inc'
       !# </include>
       potentialOffset        =-Galactic_Structure_Potential-darkMatterHaloScale_%virialVelocity(node)**2
       potentialOffsetComputed=.true.
    end if
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
    ! Store the radius.
    radiusShared=radius
    ! Compute the potential offset such that the total gravitational potential at the virial radius is -V^2 where V is the virial
    ! velocity.
    statusShared=structureErrorCodeSuccess
    Galactic_Structure_Potential=+node%mapDouble0(componentPotentialFunction,reductionSummation,optimizeFor=optimizeForPotentialSummation) &
         &                       +potentialOffset
    if (statusShared /= structureErrorCodeSuccess) status=statusShared
    include 'galactic_structure.potential.tasks.inc'
    return
  end function Galactic_Structure_Potential

  double precision function Component_Potential(component)
    !% Unary function returning the potential in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Potential=component%potential(radiusShared,componentTypeShared,massTypeShared,statusShared)
    return
  end function Component_Potential

  !# <calculationResetTask>
  !# <unitName>Galactic_Structure_Potential_Standard_Reset</unitName>
  !# </calculationResetTask>
  subroutine Galactic_Structure_Potential_Standard_Reset(node)
    !% Reset calculations for galactic structure potentials.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(in   ) :: node

    potentialOffsetComputed=.false.
    lastUniqueID           =node%uniqueID()
    return
  end subroutine Galactic_Structure_Potential_Standard_Reset

end module Galactic_Structure_Potentials












