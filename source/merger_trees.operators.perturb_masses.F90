!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implements a merger tree operator which perturbs halo masses by some error model.
  !!}

  use :: Statistics_Distributions         , only : distributionFunction1DNormal
  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorPerturbMasses">
   <description>
    A merger tree operator which perturbs halo masses by some error model.
  </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPerturbMasses
     !!{
     A merger tree operator class perturbs halo masses by some error model.
     !!}
     private
     class(nbodyHaloMassErrorClass     ), pointer :: nbodyHaloMassError_ => null()
     type (distributionFunction1DNormal)          :: standardNormal
   contains
     final     ::                         perturbMassesDestructor
     procedure :: operatePreEvolution  => perturbMassesOperatePreEvolution
  end type mergerTreeOperatorPerturbMasses

  interface mergerTreeOperatorPerturbMasses
     !!{
     Constructors for the mass perturbing merger tree operator class.
     !!}
     module procedure perturbMassesConstructorParameters
     module procedure perturbMassesConstructorInternal
  end interface mergerTreeOperatorPerturbMasses

contains

  function perturbMassesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the mass perturbing merger tree operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeOperatorPerturbMasses)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(nbodyHaloMassErrorClass        ), pointer       :: nbodyHaloMassError_

    !![
    <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    !!]
    self=mergerTreeOperatorPerturbMasses(nbodyHaloMassError_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyHaloMassError_"/>
    !!]
    return
  end function perturbMassesConstructorParameters

  function perturbMassesConstructorInternal(nbodyHaloMassError_) result(self)
    !!{
    Constructor for the mass perturbing merger tree operator class which takes a
    parameter set as input.
    !!}
    implicit none
    type (mergerTreeOperatorPerturbMasses)                         :: self
    class(nbodyHaloMassErrorClass        ), intent(in   ), target  :: nbodyHaloMassError_
    !![
    <constructorAssign variables="*nbodyHaloMassError_"/>
    !!]
    
    self%standardNormal=distributionFunction1DNormal(0.0d0,1.0d0)
    return
  end function perturbMassesConstructorInternal

  subroutine perturbMassesDestructor(self)
    !!{
    Destructor for the mass perturbing merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorPerturbMasses), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_"/>
    !!]
    return
  end subroutine perturbMassesDestructor

  subroutine perturbMassesOperatePreEvolution(self,tree)
    !!{
    Perform a mass perturbing operation on a merger tree. Perturbations are applied to each branch of the tree, and are
    independent of perturbations in all other branches. Within each branch, the perturbation to each node mass is drawn from a
    log-normal distribution with variance and correlation specified by the selected N-body statistics class.
    !!}
    use            :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: Linear_Algebra     , only : assignment(=)                , matrix            , operator(*), vector
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class           (mergerTreeOperatorPerturbMasses), intent(inout), target         :: self
    type            (mergerTree                     ), intent(inout), target         :: tree
    type            (treeNode                       ), pointer                       :: node                      , nodeChild1                  , &
         &                                                                              nodeChild2
    class           (nodeComponentBasic             ), pointer                       :: basic
    double precision                                 , allocatable  , dimension(:  ) :: deviates                  , perturbations
    double precision                                 , allocatable  , dimension(:,:) :: covariance
    type            (mergerTreeWalkerIsolatedNodes  )                                :: treeWalker
    type            (matrix                         )                                :: cholesky
    type            (vector                         )                                :: deviateVector             , perturbationVector
    integer         (c_size_t                       )                                :: nodeCount                 , iNode1                      , &
         &                                                                              iNode2
    double precision                                                                 :: rootVarianceMassFractional1, rootVarianceMassFractional2, &
         &                                                                              correlation

    ! Iterate over nodes.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! Identify nodes which are at the end (i.e. latest time) of their branch.
       if (.not.associated(node%parent).or..not.node%isPrimaryProgenitor()) then
          ! Count nodes in this branch.
          nodeCount  =  0_c_size_t
          nodeChild1 => node
          do while (associated(nodeChild1))
             nodeCount  =  nodeCount           +1_c_size_t
             nodeChild1 => nodeChild1%firstChild
          end do
          ! Allocate storage for the correlation matrix.
          allocate(covariance   (nodeCount,nodeCount))
          allocate(deviates     (nodeCount          ))
          allocate(perturbations(nodeCount          ))
          ! Walk the branch, populating the covariance matrix.
          iNode1     =  0_c_size_t
          nodeChild1 => node
          do while (associated(nodeChild1))
             iNode1=iNode1+1_c_size_t
             ! Determine the fractional error in the mass of this node.
             rootVarianceMassFractional1=+self%nbodyHaloMassError_%errorFractional(nodeChild1)
             ! Generate a random deviate.
             deviates(iNode1)=self%standardNormal%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
             ! Walk the remainder of the branch to build the covariance matrix.
             iNode2     =  iNode1-1_c_size_t
             nodeChild2 => nodeChild1
             do while (associated(nodeChild2))
                iNode2=iNode2+1_c_size_t
                ! Determine the fractional error in the mass of this node.
                rootVarianceMassFractional2=+self%nbodyHaloMassError_%errorFractional(           nodeChild2)
                ! Determine the correlation
                correlation                =+self%nbodyHaloMassError_%correlation    (nodeChild1,nodeChild2)
                ! Populate the covariance matrix.
                covariance       (iNode1,iNode2)=+rootVarianceMassFractional1 &
                     &                           *rootVarianceMassFractional2 &
                     &                           *correlation
                covariance       (iNode2,iNode1)=                             &
                     & covariance(iNode1,iNode2)
                nodeChild2 => nodeChild2%firstChild
             end do
             nodeChild1 => nodeChild1%firstChild
          end do
          ! Get Cholesky decomposition of the covariance matrix.
          cholesky=covariance
          call cholesky%choleskyDecomposition()
          ! Construct realization of fractional mass perturbations.
          deviateVector     =deviates
          perturbationVector=cholesky          *deviateVector
          perturbations     =perturbationVector
          ! Perturb the masses of the halos.
          iNode1     =  0_c_size_t
          nodeChild1 => node
          do while (associated(nodeChild1))
             basic => nodeChild1%basic()
             call basic%massSet(                                   &
                  &             +basic%mass(                     ) &
                  &             *exp       (perturbations(iNode1)) &
                  &            )
             nodeChild1 => nodeChild1%firstChild
          end do
          deallocate(covariance   )
          deallocate(deviates     )
          deallocate(perturbations)
       end if
    end do
    return
  end subroutine perturbMassesOperatePreEvolution

