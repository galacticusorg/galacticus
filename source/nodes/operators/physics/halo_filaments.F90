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
  Implements a node operator which assigns each halo the cosmic web filaments which feed it.
  !!}

  use :: Dark_Matter_Halo_Filaments, only : darkMatterHaloFilamentsClass

  !![
  <nodeOperator name="nodeOperatorHaloFilaments" docformat="rst">
   <description>
   A node operator which assigns each halo the cosmic web filaments along which it accretes, as given by the
   :galacticus-class:`darkMatterHaloFilamentsClass` object. Filaments are assigned to all nodes when the tree is initialized.

   This operator is optional: filaments are also assigned on demand, by :galacticus-class:`virialOrbitFilamentary`, to any
   branch which does not yet have them when an orbit is required---for example, for orbits assigned while a tree is being
   constructed. The filaments of a branch do not depend on which of these occurs first. Including this operator ensures that
   every halo is assigned filaments, whether or not it acquires any satellites.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloFilaments
     !!{RST
     A node operator which assigns each halo the cosmic web filaments which feed it.
     !!}
     private
     class(darkMatterHaloFilamentsClass), pointer :: darkMatterHaloFilaments_ => null()
   contains
     final     ::                       haloFilamentsDestructor
     procedure :: nodeTreeInitialize => haloFilamentsNodeTreeInitialize
  end type nodeOperatorHaloFilaments

  interface nodeOperatorHaloFilaments
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class.
     !!}
     module procedure haloFilamentsConstructorParameters
     module procedure haloFilamentsConstructorInternal
  end interface nodeOperatorHaloFilaments

contains

  function haloFilamentsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorHaloFilaments   )                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterHaloFilamentsClass), pointer       :: darkMatterHaloFilaments_

    !![
    <objectBuilder class="darkMatterHaloFilaments" name="darkMatterHaloFilaments_" source="parameters"/>
    !!]
    self=nodeOperatorHaloFilaments(darkMatterHaloFilaments_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloFilaments_"/>
    !!]
    return
  end function haloFilamentsConstructorParameters

  function haloFilamentsConstructorInternal(darkMatterHaloFilaments_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class.
    !!}
    implicit none
    type (nodeOperatorHaloFilaments   )                        :: self
    class(darkMatterHaloFilamentsClass), intent(in   ), target :: darkMatterHaloFilaments_
    !![
    <constructorAssign variables="*darkMatterHaloFilaments_"/>
    !!]

    return
  end function haloFilamentsConstructorInternal

  subroutine haloFilamentsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloFilaments), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloFilaments_"/>
    !!]
    return
  end subroutine haloFilamentsDestructor

  subroutine haloFilamentsNodeTreeInitialize(self,node)
    !!{RST
    Assign filaments to the branch of a node, if it does not already have them.
    !!}
    implicit none
    class           (nodeOperatorHaloFilaments), intent(inout), target         :: self
    type            (treeNode                 ), intent(inout), target         :: node
    double precision                           , allocatable  , dimension(:,:) :: axes

    axes=self%darkMatterHaloFilaments_%axes(node)
    return
  end subroutine haloFilamentsNodeTreeInitialize
