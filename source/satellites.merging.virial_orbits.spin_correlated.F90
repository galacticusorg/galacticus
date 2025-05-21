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
  An implementation of virial orbits which assumes that the orbital angular momentum is correlated with the spin vector of the
  host halo.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <virialOrbit name="virialOrbitSpinCorrelated">
   <description>A virial orbit class which assigns infall directions and tangential velocities directions to an orbit return by another virial orbit class such that the orbital angular momentum of the satellite is correlated with the spin vector of the host. Specifically, the angle, $\theta$, between the orbital angular momentum of the satellite and the spin of the host is assumed to be distributed such that $P(\cos \theta) = (1 + \alpha | \lambda | \cos \theta)/2$, where $|\lambda|$ is the magnitude of the host halo spin, and $\alpha$ is a parameter.</description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitSpinCorrelated
     !!{
     A virial orbit class which assumes that the orbital angular momentum is correlated with the spin vector of the host halo.
     !!}
     private
     double precision                                     :: alpha
     class           (virialOrbitClass         ), pointer :: virialOrbit_          => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
   contains
     final     ::                                    spinCorrelatedDestructor
     procedure :: orbit                           => spinCorrelatedOrbit
     procedure :: densityContrastDefinition       => spinCorrelatedDensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => spinCorrelatedVelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => spinCorrelatedVelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => spinCorrelatedAngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => spinCorrelatedAngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => spinCorrelatedVelocityTotalRootMeanSquared
     procedure :: energyMean                      => spinCorrelatedEnergyMean
     procedure :: isAngularlyResolved             => spinCorrelatedIsAngularlyResolved
  end type virialOrbitSpinCorrelated

  interface virialOrbitSpinCorrelated
     !!{
     Constructors for the \refClass{virialOrbitSpinCorrelated} virial orbit class.
     !!}
     module procedure spinCorrelatedConstructorParameters
     module procedure spinCorrelatedConstructorInternal
  end interface virialOrbitSpinCorrelated

contains

  function spinCorrelatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitSpinCorrelated} satellite virial orbit class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitSpinCorrelated)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (virialOrbitClass         ), pointer       :: virialOrbit_
    class           (darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass), pointer       :: darkMatterProfileDMO_
    double precision                                           :: alpha

    !![
    <inputParameter>
      <name>alpha</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\alpha$ which expresses the strength of the correlation between satellite orbital angular momentum and the spin of the host halo.</description>
    </inputParameter>
    <objectBuilder class="virialOrbit"          name="virialOrbit_"          source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=virialOrbitSpinCorrelated(alpha,virialOrbit_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"        />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function spinCorrelatedConstructorParameters

  function spinCorrelatedConstructorInternal(alpha,virialOrbit_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitSpinCorrelated} virial orbits class.
    !!}
    use            :: Error               , only : Component_List      , Error_Report
    use            :: Galacticus_Nodes    , only : defaultSpinComponent
    implicit none
    type            (virialOrbitSpinCorrelated)                        :: self
    double precision                           , intent(in   )         :: alpha
    class           (virialOrbitClass         ), intent(in   ), target :: virialOrbit_
    class           (darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="alpha, *virialOrbit_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    if (.not.defaultSpinComponent%angularMomentumVectorIsGettable())                                                             &
            & call Error_Report                                                                                                  &
            &      (                                                                                                             &
            &       'spin-correlated orbits require that the angularMomentumVector property of the spin component be gettable'// &
            &       Component_List(                                                                                              &
            &                      'spin'                                                                                     ,  &
            &                       defaultSpinComponent%angularMomentumVectorAttributeMatch(requireGettable=.true.)             &
            &                     )                                                                                           // &
            &       {introspection:location}                                                                                     &
            &      )
    return
  end function spinCorrelatedConstructorInternal

  subroutine spinCorrelatedDestructor(self)
    !!{
    Destructor for the \refClass{virialOrbitSpinCorrelated} virial orbits class.
    !!}
    implicit none
    type(virialOrbitSpinCorrelated), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"         />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine spinCorrelatedDestructor

  function spinCorrelatedOrbit(self,node,host,acceptUnboundOrbits)
    !!{
    Return spinCorrelated orbital parameters for a satellite.
    !!}
    use :: Coordinates             , only : assignment(=)                          , coordinateCartesian
    use :: Dark_Matter_Halo_Spins  , only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentSpin                      , treeNode
    use :: Numerical_Constants_Math, only : Pi
    use :: Vectors                 , only : Vector_Magnitude                       , Vector_Product
    implicit none
    type            (keplerOrbit              )                         :: spinCorrelatedOrbit
    class           (virialOrbitSpinCorrelated), intent(inout), target  :: self
    type            (treeNode                 ), intent(inout)          :: host                     , node
    logical                                    , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentSpin        )               , pointer :: spinHost
    double precision                           , dimension(3)           :: spinVector               , position                , &
         &                                                                 velocity                 , angularMomentum
    integer                                    , parameter              :: trialMaximum       =10000
    double precision                                                    :: probabilityRetain        , spinMagnitude           , &
         &                                                                 cosineTheta              , angularMomentumMagnitude
    logical                                                             :: retain
    integer                                                             :: trial
    type            (coordinateCartesian      )                         :: coordinates

    ! Get the underlying orbit.
    spinCorrelatedOrbit=self%virialOrbit_%orbit(node,host,acceptUnboundOrbits)
    ! Generate orbital positions until we find one that is acceptable.
    retain=.false.
    trial =0
    do while (.not.retain)
       ! Abort if too many trials have occurred.
       trial=trial+1
       if (trial > trialMaximum) call Error_Report('unable to find acceptable orbit'//{introspection:location})
       ! Sample from an isotropic distribution.
       call spinCorrelatedOrbit%phiSet    (     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()       )
       call spinCorrelatedOrbit%thetaSet  (acos(2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0))
       call spinCorrelatedOrbit%epsilonSet(     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()       )
       ! Compute the cosine of the angle between the angular momentum of this orbit and the spin of the host halo.
       spinHost  =>  host                   %spin  ()
       spinVector               =  +spinHost%angularMomentumVector() &
            &                      /Dark_Matter_Halo_Angular_Momentum_Scale(host,self%darkMatterHaloScale_,self%darkMatterProfileDMO_)
       coordinates              =   spinCorrelatedOrbit%position  ()
       position                 =   coordinates
       coordinates              =   spinCorrelatedOrbit%velocity  ()
       velocity                 =   coordinates
       spinMagnitude            =   Vector_Magnitude(spinVector                   )
       angularMomentum          =   Vector_Product  (position     ,velocity       )
       angularMomentumMagnitude =   Vector_Magnitude(              angularMomentum)
       if (spinMagnitude <= 0.0d0.or.angularMomentumMagnitude <= 0.0d0) then
          cosineTheta=0.0d0
       else
          cosineTheta  =+Dot_Product(spinVector   ,angularMomentum         ) &
               &        /            spinMagnitude                           &
               &        /                          angularMomentumMagnitude
       end if
       ! Decide whether to keep this orbit using rejection sampling.
       probabilityRetain=+(1.0d0+self%alpha*spinMagnitude*cosineTheta) &
            &            /(1.0d0+self%alpha*spinMagnitude            )
       retain           =node%hostTree%randomNumberGenerator_%uniformSample() <= probabilityRetain
    end do
    return
  end function spinCorrelatedOrbit

  function spinCorrelatedDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: spinCorrelatedDensityContrastDefinition
    class(virialOrbitSpinCorrelated ), intent(inout) :: self

    spinCorrelatedDensityContrastDefinition => self%virialOrbit_%densityContrastDefinition()
    return
  end function spinCorrelatedDensityContrastDefinition

  double precision function spinCorrelatedVelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    implicit none
    class(virialOrbitSpinCorrelated), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node, host

    spinCorrelatedVelocityTangentialMagnitudeMean=self%virialOrbit_%velocityTangentialMagnitudeMean(node,host)
    return
  end function spinCorrelatedVelocityTangentialMagnitudeMean

  function spinCorrelatedVelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , treeNode
    implicit none
    double precision                           , dimension(3)  :: spinCorrelatedVelocityTangentialVectorMean
    class           (virialOrbitSpinCorrelated), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node                                       , host
    class           (nodeComponentSpin        ), pointer       :: spinHost

    spinHost                                   =>  host                   %spin            (         )
    spinCorrelatedVelocityTangentialVectorMean =  +self    %alpha                                                                                     &
         &                                        *self    %velocityTangentialMagnitudeMean(node,host)                                                &
         &                                        *spinHost%angularMomentumVector          (         )                                                &
         &                                        /Dark_Matter_Halo_Angular_Momentum_Scale(host,self%darkMatterHaloScale_,self%darkMatterProfileDMO_) &
         &                                        /3.0d0
    return
  end function spinCorrelatedVelocityTangentialVectorMean

  double precision function spinCorrelatedAngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    implicit none
    class(virialOrbitSpinCorrelated), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node, host

    spinCorrelatedAngularMomentumMagnitudeMean=self%virialOrbit_%angularMomentumMagnitudeMean(node,host)
    return
  end function spinCorrelatedAngularMomentumMagnitudeMean

  function spinCorrelatedAngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector angular momentum.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , treeNode
    implicit none
    double precision                           , dimension(3)  :: spinCorrelatedAngularMomentumVectorMean
    class           (virialOrbitSpinCorrelated), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node                                   , host
    class           (nodeComponentSpin        ), pointer       :: spinHost

    spinHost                                =>  host    %spin                        (         )
    spinCorrelatedAngularMomentumVectorMean =  +self    %alpha                                                                                     &
         &                                     *self    %angularMomentumMagnitudeMean(node,host)                                                   &
         &                                     *spinHost%angularMomentumVector       (         )                                                   &
         &                                     /Dark_Matter_Halo_Angular_Momentum_Scale(host,self%darkMatterHaloScale_,self%darkMatterProfileDMO_) &
         &                                     /3.0d0
    return
  end function spinCorrelatedAngularMomentumVectorMean

  double precision function spinCorrelatedVelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the root mean squared of the total velocity.
    !!}
    implicit none
    class(virialOrbitSpinCorrelated), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node, host

    spinCorrelatedVelocityTotalRootMeanSquared=self%virialOrbit_%velocityTotalRootMeanSquared(node,host)
    return
  end function spinCorrelatedVelocityTotalRootMeanSquared

  double precision function spinCorrelatedEnergyMean(self,node,host)
    !!{
    Return the mean of the total energy.
    !!}
    implicit none
    class(virialOrbitSpinCorrelated), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node, host

    spinCorrelatedEnergyMean=self%virialOrbit_%energyMean(node,host)
    return
  end function spinCorrelatedEnergyMean

  logical function spinCorrelatedIsAngularlyResolved(self)
    !!{
    Return true indicating that orbits are angularly-resolved.
    !!}
    implicit none
    class(virialOrbitSpinCorrelated), intent(inout) :: self
    !$GLC attributes unused :: self

    spinCorrelatedIsAngularlyResolved=.true.
    return
  end function spinCorrelatedIsAngularlyResolved
  
