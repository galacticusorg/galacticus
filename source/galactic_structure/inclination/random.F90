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
  Implements a galactic inclination class in which each galaxy is assigned a random orientation.
  !!}

  !![
  <galacticInclination name="galacticInclinationRandom" docformat="rst">
   <description>
   A galactic inclination class in which each galaxy is assigned an isotropically distributed orientation---uniform
   in :math:`\cos i`---which it then keeps.

   The angle is *not* drawn here. It is drawn once per galaxy by
   :galacticus-class:`nodeOperatorInclinationRandom`, which stores it in a meta-property of the ``basic`` component;
   this class only reads it back. Requiring that operator is what makes the orientation stable: every consumer sees
   the same angle for a given galaxy, the angle is the same at every output, and it survives the promotion of a node
   up its branch, since meta-properties move with the component.

   Drawing the angle at the point of use instead would fail all three: separate consumers would disagree, a galaxy
   would be differently oriented at each output, and consuming randoms during output would perturb every subsequent
   draw in the model.
   </description>
  </galacticInclination>
  !!]
  type, extends(galacticInclinationClass) :: galacticInclinationRandom
     !!{RST
     A galactic inclination class in which each galaxy is assigned a random orientation.
     !!}
     private
     integer :: inclinationID
   contains
     procedure :: inclination => randomInclination
  end type galacticInclinationRandom

  interface galacticInclinationRandom
     !!{RST
     Constructors for the :galacticus-class:`galacticInclinationRandom` galactic inclination class.
     !!}
     module procedure randomConstructorParameters
     module procedure randomConstructorInternal
  end interface galacticInclinationRandom

contains

  function randomConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticInclinationRandom` galactic inclination class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticInclinationRandom)                :: self
    type(inputParameters          ), intent(inout) :: parameters

    self=galacticInclinationRandom()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function randomConstructorParameters

  function randomConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticInclinationRandom` galactic inclination class.
    !!}
    implicit none
    type(galacticInclinationRandom) :: self

    ! This class only reads the meta-property; `nodeOperatorInclinationRandom` creates and sets it.
    !![
    <addMetaProperty component="basic" name="inclination" id="self%inclinationID" isEvolvable="no"/>
    !!]
    return
  end function randomConstructorInternal

  double precision function randomInclination(self,node) result(inclination)
    !!{RST
    Return the stored random inclination of the galaxy.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticInclinationRandom), intent(inout)         :: self
    type (treeNode                 ), intent(inout), target :: node
    class(nodeComponentBasic       ), pointer               :: basic

    basic       => node %basic                    (                  )
    inclination =  basic%floatRank0MetaPropertyGet(self%inclinationID)
    return
  end function randomInclination
