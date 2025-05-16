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
  Implements a node operator class that initializes halo angular momenta using spins drawn at random from a distribution.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Halo_Spin_Distributions, only : haloSpinDistributionClass

  !![
  <nodeOperator name="nodeOperatorHaloAngularMomentumRandom">
   <description>A node operator class that initializes halo spins to random values drawn from a distribution.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumRandom
     !!{
     A node operator class that initializes halo spins to random values drawn from a distribution.
     !!}
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (haloSpinDistributionClass), pointer :: haloSpinDistribution_ => null()
     double precision                                     :: factorReset
   contains
     final     ::                   haloAngularMomentumRandomDestructor
     procedure :: nodeInitialize => haloAngularMomentumRandomNodeInitialize
  end type nodeOperatorHaloAngularMomentumRandom
  
  interface nodeOperatorHaloAngularMomentumRandom
     !!{
     Constructors for the \refClass{nodeOperatorHaloAngularMomentumRandom} node operator class.
     !!}
     module procedure haloAngularMomentumRandomConstructorParameters
     module procedure haloAngularMomentumRandomConstructorInternal
  end interface nodeOperatorHaloAngularMomentumRandom
  
contains
  
  function haloAngularMomentumRandomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAngularMomentumRandom} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloAngularMomentumRandom)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_
    class           (haloSpinDistributionClass            ), pointer       :: haloSpinDistribution_
    double precision                                                       :: factorReset
     
    !![
    <inputParameter>
      <name>factorReset</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The factor by which a node must increase in mass before its spin parameter is reset.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=nodeOperatorHaloAngularMomentumRandom(factorReset,haloSpinDistribution_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="haloSpinDistribution_"/>
    !!]
    return
  end function haloAngularMomentumRandomConstructorParameters

  function haloAngularMomentumRandomConstructorInternal(factorReset,haloSpinDistribution_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHaloAngularMomentumRandom} node operator class.
    !!}
    implicit none
    type            (nodeOperatorHaloAngularMomentumRandom)                        :: self
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    class           (haloSpinDistributionClass            ), intent(in   ), target :: haloSpinDistribution_
    double precision                                       , intent(in   )         :: factorReset
    !![
    <constructorAssign variables="factorReset, *haloSpinDistribution_, *darkMatterHaloScale_"/>
    !!]

    return
  end function haloAngularMomentumRandomConstructorInternal

  subroutine haloAngularMomentumRandomDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorHaloAngularMomentumRandom} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAngularMomentumRandom), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%haloSpinDistribution_"/>
    !!]
    return
  end subroutine haloAngularMomentumRandomDestructor

  subroutine haloAngularMomentumRandomNodeInitialize(self,node)
    !!{
    Assign a randomly-drawn spin to a node.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentBasic                     , nodeComponentSpin
    implicit none
    class           (nodeOperatorHaloAngularMomentumRandom), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    type            (treeNode                             )               , pointer :: nodeProgenitor
    class           (nodeComponentSpin                    )               , pointer :: spinProgenitor , spin
    class           (nodeComponentBasic                   )               , pointer :: basicProgenitor
    double precision                                                                :: massPrevious   , spinPrevious, &
         &                                                                             angularMomentum

    ! Ensure that the spin has not yet been assigned for this node.
    spin => node%spin()
    select type (spin)
    type is (nodeComponentSpin)
       ! Walk the tree back along primary children to the earliest such progenitor.
       nodeProgenitor => node
       do while (associated(nodeProgenitor%firstChild))
          nodeProgenitor => nodeProgenitor%firstChild
       end do
       ! Walk forward through the branch, assigning spins/angular momenta. If the mass of the halo exceeds that of the halo for
       ! which we last selected a spin by a given factor, then select a new spin from the distribution. Otherwise, use the
       ! previously assigned spin. If available, the vector angular momentum is also set. As this model makes no prediction for
       ! the direction of the angular momentum vector, it is assumed to be always aligned with the first axis.
       spinProgenitor  => nodeProgenitor                       %spin  (autoCreate=.true.        )
       basicProgenitor => nodeProgenitor                       %basic (                         )
       spinPrevious    =  self           %haloSpinDistribution_%sample(           nodeProgenitor)
       massPrevious    =  basicProgenitor                      %mass  (                         )
       angularMomentum =  spinPrevious*Dark_Matter_Halo_Angular_Momentum_Scale(nodeProgenitor,self%darkMatterHaloScale_)
       call spinProgenitor%angularMomentumSet(spinPrevious*Dark_Matter_Halo_Angular_Momentum_Scale(nodeProgenitor,self%darkMatterHaloScale_))
       call spinProgenitor%angularMomentumSet(angularMomentum)
       if (spinProgenitor%angularMomentumVectorIsSettable()) &
            & call spinProgenitor%angularMomentumVectorSet([angularMomentum,0.0d0,0.0d0])
       do while (nodeProgenitor%isPrimaryProgenitor())
          nodeProgenitor  => nodeProgenitor%parent
          basicProgenitor => nodeProgenitor%basic ()
          if (basicProgenitor%mass() > self%factorReset*massPrevious) then
             spinPrevious=self           %haloSpinDistribution_%sample(nodeProgenitor)
             massPrevious=basicProgenitor                      %mass  (              )
          end if
          spinProgenitor => nodeProgenitor%spin(autoCreate=.true.)
          angularMomentum =  spinPrevious*Dark_Matter_Halo_Angular_Momentum_Scale(nodeProgenitor,self%darkMatterHaloScale_)
          call       spinProgenitor%angularMomentumSet       ( angularMomentum             )
          if (spinProgenitor%angularMomentumVectorIsSettable()) &
               & call spinProgenitor%angularMomentumVectorSet([angularMomentum,0.0d0,0.0d0])
       end do
    end select
    return
  end subroutine haloAngularMomentumRandomNodeInitialize
  
