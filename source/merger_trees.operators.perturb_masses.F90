!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% Contains a module which implements a merger tree operator which perturbs halo masses by some error model.

  use Statistics_NBody_Halo_Mass_Errors
  use Statistics_Distributions
  
  !# <mergerTreeOperator name="mergerTreeOperatorPerturbMasses" defaultThreadPrivate="yes">
  !#  <description>
  !#   A merger tree operator which perturbs halo masses by some error model.
  !# </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPerturbMasses
     !% A merger tree operator class perturbs halo masses by some error model.
     private
     class(nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_
     type (distributionNormal     )          :: standardNormal
   contains
     final     ::             perturbMassesDestructor
     procedure :: operate  => perturbMassesOperate
  end type mergerTreeOperatorPerturbMasses
  
  interface mergerTreeOperatorPerturbMasses
     !% Constructors for the mass perturbing merger tree operator class.
     module procedure perturbMassesConstructorParameters
     module procedure perturbMassesConstructorInternal
  end interface mergerTreeOperatorPerturbMasses

contains

  function perturbMassesConstructorParameters(parameters) result(self)
    !% Constructor for the mass perturbing merger tree operator class which takes a
    !% parameter set as input.
    use Input_Parameters2
    implicit none
    type (mergerTreeOperatorPerturbMasses)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(nbodyHaloMassErrorClass        ), pointer       :: nbodyHaloMassError_
    
    !# <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    self=mergerTreeOperatorPerturbMasses(nbodyHaloMassError_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function perturbMassesConstructorParameters

  function perturbMassesConstructorInternal(nbodyHaloMassError_) result(self)
    !% Constructor for the mass perturbing merger tree operator class which takes a
    !% parameter set as input.
    use Input_Parameters2
    implicit none
    type (mergerTreeOperatorPerturbMasses)                        :: self
    class(nbodyHaloMassErrorClass        ), intent(in   ), target :: nbodyHaloMassError_
    !# <constructorAssign variables="*nbodyHaloMassError_"/>

    ! Build a standard normal distribution.
    self%standardNormal=distributionNormal(0.0d0,1.0d0)
    return
  end function perturbMassesConstructorInternal

  subroutine perturbMassesDestructor(self)
    !% Destructor for the mass perturbing merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPerturbMasses), intent(inout) :: self
    
    !# <objectDestructor name="self%nbodyHaloMassError_"/>
    return
  end subroutine perturbMassesDestructor
  
  subroutine perturbMassesOperate(self,tree)
    !% Perform a mass perturbing operation on a merger tree.
    implicit none
    class           (mergerTreeOperatorPerturbMasses), intent(inout)         :: self
    type            (mergerTree                     ), intent(inout), target :: tree
    type            (treeNode                       ), pointer               :: node
    class           (nodeComponentBasic             ), pointer               :: basic
    type            (mergerTree                     ), pointer               :: treeCurrent
    double precision                                                         :: rootVarianceMassFractional, perturbationMassFractional
    
    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Get root node of the tree.
       node => treeCurrent%baseNode
       ! Walk the tree.
       do while (associated(node))
          ! Determine the fractional error in the mass of this node.
          rootVarianceMassFractional=+self%nbodyHaloMassError_%errorFractional(node                                                   )
          ! Draw a fractional error at random from a normal distribution with this root-variance.
          perturbationMassFractional=+rootVarianceMassFractional                                                                        &
               &                     *self%standardNormal     %sample         (randomNumberGenerator=treeCurrent%randomNumberGenerator)
          ! Perturb the mass of the halo.
          basic => node%basic()
          call basic%massSet(                                        &
               &             +basic%mass(                          ) &
               &             *exp       (perturbationMassFractional) &
               &            )
          ! Walk to the next node in the tree.
          node => node%walkTree()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine perturbMassesOperate

