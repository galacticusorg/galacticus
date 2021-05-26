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

  !% Implements a node operator class that initializes halo spins using a random walk in angular momentum.

  use :: Halo_Spin_Distributions , only : haloSpinDistributionClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass

  !# <nodeOperator name="nodeOperatorHaloSpinRandomWalk">
  !#  <description>
  !#   A node operator class that initializes halo spins using a random walk in angular momentum. The three components of the
  !#   angular momentum vector are treated as independent Wiener processes with time-dependent variance. Specifically, each
  !#   component of the angular momentum vector obeys:
  !#   \begin{equation}
  !#    J_\mathrm{i}(t_2) = J_\mathrm{i}(t_1) + \left[ \sigma^2 \left( J_\mathrm{v}^2(t_2) - J_\mathrm{v}^2(t_1) \right) \right]^{1/2} N(0,1)
  !#   \end{equation}
  !#   where $J_\mathrm{v}(t) = M_\mathrm{v}(t) V_\mathrm{v}(t) R_\mathrm{v}(t)$ is the characteristic virial angular momentum,
  !#   $M_\mathrm{v}(t)$, $V_\mathrm{v}(t)$, and $R_\mathrm{v}(t)$ are the virial mass, velocity, and radius respectively,
  !#   $\sigma^2$ represents the variance in angular momentum per unit increase in $J_\mathrm{v}^2$, and $N(0,1)$ is a random
  !#   variable distributed as a standard normal.
  !#
  !#   The spin is then found from the magnitude of the total angular momentum, $|J| = \sum_{i=1}^3 J_\mathrm{i}^2$.
  !#  </description>
  !# </nodeOperator>
  type, extends(nodeOperatorClass) :: nodeOperatorHaloSpinRandomWalk
     !% A node operator class that initializes halo spins using a random walk in angular momentum.
     private
     class           (haloSpinDistributionClass), pointer :: haloSpinDistribution_ => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     double precision                                     :: spinVarianceSpecific
   contains
     final     ::                   haloSpinRandomWalkDestructor
     procedure :: nodeInitialize => haloSpinRandomWalkNodeInitialize
  end type nodeOperatorHaloSpinRandomWalk
  
  interface nodeOperatorHaloSpinRandomWalk
     !% Constructors for the {\normalfont \ttfamily haloSpinRandomWalk} node operator class.
     module procedure haloSpinRandomWalkConstructorParameters
     module procedure haloSpinRandomWalkConstructorInternal
  end interface nodeOperatorHaloSpinRandomWalk
  
contains
  
  function haloSpinRandomWalkConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily haloSpinRandomWalk} node operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloSpinRandomWalk)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (haloSpinDistributionClass     ), pointer       :: haloSpinDistribution_
    class           (darkMatterProfileDMOClass     ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    double precision                                                :: spinVarianceSpecific
     
    !# <inputParameter>
    !#   <name>spinVarianceSpecific</name>
    !#   <description>The variance in the difference in the angular momentum of a halo per unit mass growth.</description>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0029d0</defaultValue>
    !# </inputParameter>
    !# <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    self=nodeOperatorHaloSpinRandomWalk(spinVarianceSpecific,haloSpinDistribution_,darkMatterProfileDMO_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="haloSpinDistribution_"/>
    !# <objectDestructor name="darkMatterProfileDMO_"/>
    !# <objectDestructor name="darkMatterHaloScale_" />
    return
  end function haloSpinRandomWalkConstructorParameters

  function haloSpinRandomWalkConstructorInternal(spinVarianceSpecific,haloSpinDistribution_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily haloSpinRandomWalk} node operator class.
    implicit none
    type            (nodeOperatorHaloSpinRandomWalk)                        :: self
    class           (haloSpinDistributionClass     ), intent(in   ), target :: haloSpinDistribution_
    class           (darkMatterProfileDMOClass     ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                , intent(in   )         :: spinVarianceSpecific
    !# <constructorAssign variables="spinVarianceSpecific, *haloSpinDistribution_, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>

    return
  end function haloSpinRandomWalkConstructorInternal

  subroutine haloSpinRandomWalkDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloSpinRandomWalk} node operator class.
    implicit none
    type(nodeOperatorHaloSpinRandomWalk), intent(inout) :: self

    !# <objectDestructor name="self%haloSpinDistribution_"/>
    !# <objectDestructor name="self%darkMatterProfileDMO_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine haloSpinRandomWalkDestructor

  subroutine haloSpinRandomWalkNodeInitialize(self,node)
    !% Assign a randomly-drawn spin to a node.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum, Dark_Matter_Halo_Spin
    implicit none
    class           (nodeOperatorHaloSpinRandomWalk), intent(inout)               :: self
    type            (treeNode                      ), intent(inout), target       :: node
    type            (treeNode                      )               , pointer      :: nodeProgenitor
    class           (nodeComponentSpin             )               , pointer      :: spinProgenitor         , spin
    class           (nodeComponentBasic            )               , pointer      :: basicProgenitor
    double precision                                               , dimension(3) :: angularMomentumVector
    double precision                                                              :: angularMomentumPrevious, angularMomentumCurrent, &
         &                                                                           spinScalar
    integer                                                                       :: i
    
    ! Ensure that the spin has not yet been assigned for this node.
    spin => node%spin()
    select type (spin)
    type is (nodeComponentSpin)
       ! Walk the tree back along primary children to the earliest such progenitor.
       nodeProgenitor => node
       do while (associated(nodeProgenitor%firstChild))
          nodeProgenitor => nodeProgenitor%firstChild
       end do
       ! Select a spin for the initial halo from a distribution function.
       basicProgenitor => nodeProgenitor                      %basic (                         )
       spinProgenitor  => nodeProgenitor                      %spin  (autoCreate=.true.        )
       spinScalar      =  self          %haloSpinDistribution_%sample(           nodeProgenitor)
       call spinProgenitor%spinSet(spinScalar)
       ! Compute the initial angular momentum vector. We choose this to be aligned along the x-axis. As we only care about the
       ! magnitude of the spin any choice of initial vector direction is equivalent.
       angularMomentumVector  =[                                                                              &
            &                   Dark_Matter_Halo_Angular_Momentum(nodeProgenitor,self%darkMatterProfileDMO_), &
            &                   0.0d0                                                                       , &
            &                   0.0d0                                                                         &
            &                  ]
       ! Compute the characteristic angular momentum of this halo.
       angularMomentumPrevious=+self           %darkMatterHaloScale_%virialRadius  (nodeProgenitor) &
            &                  *self           %darkMatterHaloScale_%virialVelocity(nodeProgenitor) &
            &                  *basicProgenitor                     %mass          (              )
       ! Walk up through descendents.
       do while (associated(nodeProgenitor))
          ! Set the spin of the current halo from the current angular momentum vector.
          basicProgenitor => nodeProgenitor%basic(                 )
          spinProgenitor  => nodeProgenitor%spin (autoCreate=.true.)
          ! Compute the characteristic angular momentum of this halo.
          angularMomentumCurrent=+self           %darkMatterHaloScale_%virialRadius  (nodeProgenitor) &
               &                 *self           %darkMatterHaloScale_%virialVelocity(nodeProgenitor) &
               &                 *basicProgenitor                     %mass          (              )
          ! Perform the random walk in each dimension for the angular momentum vector.
          do i=1,3
             angularMomentumVector(i)=+angularMomentumVector(i)                                              &
                  &                   +sqrt(                                                                 &
                  &                         +self%spinVarianceSpecific                                       &
                  &                         *(                                                               &
                  &                           +angularMomentumCurrent**2                                     &
                  &                           -angularMomentumPrevious**2                                    &
                  &                          )                                                               &
                  &                        )                                                                 &
                  &                   *nodeProgenitor%hostTree%randomNumberGenerator_%standardNormalSample()
          end do
          ! Compute and store the scalar spin.
          call spinProgenitor%spinSet(Dark_Matter_Halo_Spin(nodeProgenitor,sqrt(sum(angularMomentumVector**2)),self%darkMatterProfileDMO_))
          ! Store the current characteristic angular momentum.
          angularMomentumPrevious=angularMomentumCurrent
          ! Move to the next descendent halo.
          if (nodeProgenitor%isPrimaryProgenitor()) then
             nodeProgenitor  => nodeProgenitor%parent
          else
             nodeProgenitor => null()
          end if
       end do
    end select
    return
  end subroutine haloSpinRandomWalkNodeInitialize

