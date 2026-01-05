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

  !!{
  Implements a node operator class that initializes halo angular momenta using a random walk in angular momentum.
  !!}

  use :: Halo_Spin_Distributions, only : haloSpinDistributionClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorHaloAngularMomentumRandomWalk">
   <description>
    A node operator class that initializes halo spins using a random walk in angular momentum. The three components of the
    angular momentum vector are treated as independent Wiener processes with time-dependent variance. Specifically, each
    component of the angular momentum vector obeys:
    \begin{equation}
     J_\mathrm{i}(t_2) = J_\mathrm{i}(t_1) + \left[ \sigma^2 \left( J_\mathrm{v}^2(t_2) - J_\mathrm{v}^2(t_1) \right) \right]^{1/2} N(0,1)
    \end{equation}
    where $J_\mathrm{v}(t) = M_\mathrm{v}(t) V_\mathrm{v}(t) R_\mathrm{v}(t)$ is the characteristic virial angular momentum,
    $M_\mathrm{v}(t)$, $V_\mathrm{v}(t)$, and $R_\mathrm{v}(t)$ are the virial mass, velocity, and radius respectively,
    $\sigma^2$ represents the variance in angular momentum per unit increase in $J_\mathrm{v}^2$, and $N(0,1)$ is a random
    variable distributed as a standard normal.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumRandomWalk
     !!{
     A node operator class that initializes halo spins using a random walk in angular momentum.
     !!}
     private
     class           (haloSpinDistributionClass), pointer :: haloSpinDistribution_           => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_            => null()
     double precision                                     :: angularMomentumVarianceSpecific
   contains
     final     ::                   haloAngularMomentumRandomWalkDestructor
     procedure :: nodeInitialize => haloAngularMomentumRandomWalkNodeInitialize
  end type nodeOperatorHaloAngularMomentumRandomWalk
  
  interface nodeOperatorHaloAngularMomentumRandomWalk
     !!{
     Constructors for the \refClass{nodeOperatorHaloAngularMomentumRandomWalk} node operator class.
     !!}
     module procedure haloAngularMomentumRandomWalkConstructorParameters
     module procedure haloAngularMomentumRandomWalkConstructorInternal
  end interface nodeOperatorHaloAngularMomentumRandomWalk
  
contains
  
  function haloAngularMomentumRandomWalkConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAngularMomentumRandomWalk} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloAngularMomentumRandomWalk)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (haloSpinDistributionClass                ), pointer       :: haloSpinDistribution_
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    double precision                                                           :: angularMomentumVarianceSpecific
     
    !![
    <inputParameter>
      <name>angularMomentumVarianceSpecific</name>
      <description>The variance in the difference in the angular momentum of a halo per unit mass growth.</description>
      <source>parameters</source>
      <defaultValue>0.0029d0</defaultValue>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=nodeOperatorHaloAngularMomentumRandomWalk(angularMomentumVarianceSpecific,haloSpinDistribution_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function haloAngularMomentumRandomWalkConstructorParameters

  function haloAngularMomentumRandomWalkConstructorInternal(angularMomentumVarianceSpecific,haloSpinDistribution_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHaloAngularMomentumRandomWalk} node operator class.
    !!}
    implicit none
    type            (nodeOperatorHaloAngularMomentumRandomWalk)                        :: self
    class           (haloSpinDistributionClass                ), intent(in   ), target :: haloSpinDistribution_
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                           , intent(in   )         :: angularMomentumVarianceSpecific
    !![
    <constructorAssign variables="angularMomentumVarianceSpecific, *haloSpinDistribution_, *darkMatterHaloScale_"/>
    !!]

    return
  end function haloAngularMomentumRandomWalkConstructorInternal

  subroutine haloAngularMomentumRandomWalkDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorHaloAngularMomentumRandomWalk} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAngularMomentumRandomWalk), intent(inout) :: self

    !![
    <objectDestructor name="self%haloSpinDistribution_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine haloAngularMomentumRandomWalkDestructor

  subroutine haloAngularMomentumRandomWalkNodeInitialize(self,node)
    !!{
    Assign a randomly-drawn spin to a node.
    !!}
    use :: Galacticus_Nodes      , only : nodeComponentBasic                     , nodeComponentSpin
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Vectors               , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorHaloAngularMomentumRandomWalk), intent(inout), target       :: self
    type            (treeNode                                 ), intent(inout), target       :: node
    type            (treeNode                                 )               , pointer      :: nodeProgenitor
    class           (nodeComponentSpin                        )               , pointer      :: spinProgenitor         , spin
    class           (nodeComponentBasic                       )               , pointer      :: basicProgenitor
    double precision                                                          , dimension(3) :: angularMomentumVector
    double precision                                                                         :: angularMomentumPrevious, angularMomentumCurrent, &
         &                                                                                      angularMomentumScalar
    integer                                                                                  :: i
    
    ! Ensure that the spin has not yet been assigned for this node.
    spin => node%spin()
    select type (spin)
    type is (nodeComponentSpin)
       ! Walk the tree back along primary children to the earliest such progenitor.
       nodeProgenitor => node
       do while (associated(nodeProgenitor%firstChild))
          nodeProgenitor => nodeProgenitor%firstChild
       end do
       ! Select a angular momentum for the initial halo using the spin distribution function.
       basicProgenitor       =>  nodeProgenitor                      %basic (                         )
       spinProgenitor        =>  nodeProgenitor                      %spin  (autoCreate=.true.        )
       angularMomentumScalar =  +self          %haloSpinDistribution_%sample(           nodeProgenitor)            &
            &                   *Dark_Matter_Halo_Angular_Momentum_Scale(nodeProgenitor,self%darkMatterHaloScale_)
       call spinProgenitor%angularMomentumSet(angularMomentumScalar)
       ! Compute the initial angular momentum vector. We choose this to be aligned along the x-axis. As we only care about the
       ! magnitude of the angular momentum any choice of initial vector direction is equivalent.
       angularMomentumVector  =[                       &
            &                   angularMomentumScalar, &
            &                   0.0d0                , &
            &                   0.0d0                  &
            &                  ]
       ! Compute the characteristic angular momentum of this halo.
       angularMomentumPrevious=+self           %darkMatterHaloScale_%radiusVirial  (nodeProgenitor) &
            &                  *self           %darkMatterHaloScale_%velocityVirial(nodeProgenitor) &
            &                  *basicProgenitor                     %mass          (              )
       ! Walk up through descendants.
       do while (associated(nodeProgenitor))
          ! Set the spin of the current halo from the current angular momentum vector.
          basicProgenitor => nodeProgenitor%basic(                 )
          spinProgenitor  => nodeProgenitor%spin (autoCreate=.true.)
          ! Compute the characteristic angular momentum of this halo.
          angularMomentumCurrent=+self           %darkMatterHaloScale_%radiusVirial  (nodeProgenitor) &
               &                 *self           %darkMatterHaloScale_%velocityVirial(nodeProgenitor) &
               &                 *basicProgenitor                     %mass          (              )
          ! Perform the random walk in each dimension for the angular momentum vector.
          do i=1,3
             angularMomentumVector(i)=+angularMomentumVector(i)                                              &
                  &                   +sqrt(                                                                 &
                  &                         +self%angularMomentumVarianceSpecific                            &
                  &                         *max(                                                            &
                  &                              +angularMomentumCurrent **2                                 &
                  &                              -angularMomentumPrevious**2,                                &
                  &                              +0.0d0                                                      &
                  &                             )                                                            &
                  &                        )                                                                 &
                  &                   *nodeProgenitor%hostTree%randomNumberGenerator_%standardNormalSample()
          end do
          ! Store the scalar angular momentum.
          call        spinProgenitor%angularMomentumSet      (Vector_Magnitude(angularMomentumVector))
          if (spinProgenitor%angularMomentumVectorIsSettable())                                        &
               & call spinProgenitor%angularMomentumVectorSet(                 angularMomentumVector )
          ! Store the current characteristic angular momentum.
          angularMomentumPrevious=max(angularMomentumPrevious,angularMomentumCurrent)
          ! Move to the next descendant halo.
          if (nodeProgenitor%isPrimaryProgenitor()) then
             nodeProgenitor  => nodeProgenitor%parent
          else
             nodeProgenitor => null()
          end if
       end do
    end select
    return
  end subroutine haloAngularMomentumRandomWalkNodeInitialize

