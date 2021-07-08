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

  !!{
  Implements a node operator class that initializes halo spins to random values drawn from a distribution.
  !!}

  use :: Halo_Spin_Distributions, only : haloSpinDistributionClass

  !![
  <nodeOperator name="nodeOperatorHaloSpinRandom">
   <description>A node operator class that initializes halo spins to random values drawn from a distribution.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloSpinRandom
     !!{
     A node operator class that initializes halo spins to random values drawn from a distribution.
     !!}
     private
     class           (haloSpinDistributionClass), pointer :: haloSpinDistribution_ => null()
     double precision                                     :: factorReset
   contains
     final     ::                   haloSpinRandomDestructor
     procedure :: nodeInitialize => haloSpinRandomNodeInitialize
  end type nodeOperatorHaloSpinRandom
  
  interface nodeOperatorHaloSpinRandom
     !!{
     Constructors for the {\normalfont \ttfamily haloSpinRandom} node operator class.
     !!}
     module procedure haloSpinRandomConstructorParameters
     module procedure haloSpinRandomConstructorInternal
  end interface nodeOperatorHaloSpinRandom
  
contains
  
  function haloSpinRandomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloSpinRandom} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloSpinRandom)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (haloSpinDistributionClass ), pointer       :: haloSpinDistribution_
    double precision                                            :: factorReset
     
    !![
    <inputParameter>
      <name>factorReset</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The factor by which a node must increase in mass before its spin parameter is reset.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    !!]
    self=nodeOperatorHaloSpinRandom(factorReset,haloSpinDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"/>
    !!]
    return
  end function haloSpinRandomConstructorParameters

  function haloSpinRandomConstructorInternal(factorReset,haloSpinDistribution_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily haloSpinRandom} node operator class.
    !!}
    implicit none
    type            (nodeOperatorHaloSpinRandom)                        :: self
    class           (haloSpinDistributionClass ), intent(in   ), target :: haloSpinDistribution_
    double precision                            , intent(in   )         :: factorReset
    !![
    <constructorAssign variables="factorReset, *haloSpinDistribution_"/>
    !!]

    return
  end function haloSpinRandomConstructorInternal

  subroutine haloSpinRandomDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily haloSpinRandom} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloSpinRandom), intent(inout) :: self

    !![
    <objectDestructor name="self%haloSpinDistribution_"/>
    !!]
    return
  end subroutine haloSpinRandomDestructor

  subroutine haloSpinRandomNodeInitialize(self,node)
    !!{
    Assign a randomly-drawn spin to a node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    implicit none
    class           (nodeOperatorHaloSpinRandom), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), target  :: node
    type            (treeNode                  )               , pointer :: nodeProgenitor
    class           (nodeComponentSpin         )               , pointer :: spinProgenitor , spin
    class           (nodeComponentBasic        )               , pointer :: basicProgenitor
    double precision                                                     :: massPrevious   , spinPrevious

    ! Ensure that the spin has not yet been assigned for this node.
    spin => node%spin()
    select type (spin)
    type is (nodeComponentSpin)
       ! Walk the tree back along primary children to the earliest such progenitor.
       nodeProgenitor => node
       do while (associated(nodeProgenitor%firstChild))
          nodeProgenitor => nodeProgenitor%firstChild
       end do
       ! Walk forward through the branch, assigning spins. If the mass of the halo exceeds that of the halo for which we last
       ! selected a spin by a given factor, then select a new spin from the distribution. Otherwise, use the previously
       ! assigned spin.
       spinProgenitor  => nodeProgenitor                       %spin  (autoCreate=.true.        )
       basicProgenitor => nodeProgenitor                       %basic (                         )
       spinPrevious    =  self           %haloSpinDistribution_%sample(           nodeProgenitor)
       massPrevious    =  basicProgenitor                      %mass  (                         )
       call spinProgenitor%spinSet(spinPrevious)
       do while (nodeProgenitor%isPrimaryProgenitor())
          nodeProgenitor  => nodeProgenitor%parent
          basicProgenitor => nodeProgenitor%basic ()
          if (basicProgenitor%mass() > self%factorReset*massPrevious) then
             spinPrevious=self           %haloSpinDistribution_%sample(nodeProgenitor)
             massPrevious=basicProgenitor                      %mass  (              )
          end if
          spinProgenitor => nodeProgenitor%spin(autoCreate=.true.)
          call spinProgenitor%spinSet(spinPrevious)
       end do
    end select
    return
  end subroutine haloSpinRandomNodeInitialize
  
