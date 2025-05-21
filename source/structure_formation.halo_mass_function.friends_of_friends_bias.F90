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
Contains a  module which  implements a  dark matter  halo mass function  class which  modifies another  mass function  to mimic
systematic errors arising in the friends-of-friends halo finding algorithm.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <haloMassFunction name="haloMassFunctionFofBias">
   <description>
    The halo mass function is computed by modifying another halo mass function to mimic systematic
    errors arising in the friends-of-friends halo finding algorithm. Specifically, a
    systematic shift in mass motivated by the results of percolation theory \cite[their eqn. B11]{more_overdensity_2011}
    is applied. In particular, $M_\mathrm{particle}=${\normalfont \ttfamily [massParticle]} is the mass of the particle in the simulation
    to which the friends-of-friends algorithm was applied.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionFofBias
     !!{
     A halo mass function class which modifies another mass function to mimic systematic errors arising in the
     friends-of-friends halo finding algorithm.
     !!}
     private
     double precision                                     :: massParticle                     , massInfiniteToMassSharpEdge, &
          &                                                  linkingLength
     logical                                              :: linkingLengthIsComoving
     class           (haloMassFunctionClass    ), pointer :: massFunctionIntrinsic   => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_    => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_   => null()
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_     => null()
    contains
     final     ::                 fofBiasDestructor
     procedure :: differential => fofBiasDifferential
  end type haloMassFunctionFofBias

  interface haloMassFunctionFofBias
     !!{
     Constructors for the \refClass{haloMassFunctionFofBias} halo mass function class.
     !!}
     module procedure fofBiasConstructorParameters
     module procedure fofBiasConstructorInternal
  end interface haloMassFunctionFofBias

contains

  function fofBiasConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionFofBias} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionFofBias  )                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (haloMassFunctionClass    ), pointer       :: massFunctionIntrinsic
    class           (darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass), pointer       :: darkMatterProfileDMO_
    class           (cosmologyParametersClass ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    double precision                                           :: massParticle           , massInfiniteToMassSharpEdge, &
         &                                                        linkingLength
    logical                                                    :: linkingLengthIsComoving

    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>Parameter $M_\mathrm{particle}$ appearing in model for friends-of-friends errors in the halo mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>massInfiniteToMassSharpEdge</name>
      <source>parameters</source>
      <defaultValue>0.98d0</defaultValue>
      <defaultSource>\cite[estimate based on comments in text]{more_overdensity_2011}</defaultSource>
      <description>The ratio of the friends-of-friends mass in the limit of infinite number of particles to the mass of the halo enclosed within a sharp-edged sphere bounding an isodensity surface equal to the critical density for percolation.</description>
    </inputParameter>
    <inputParameter>
      <name>linkingLength</name>
      <source>parameters</source>
      <description>The linking length (in physical Mpc) used in the friends-of-friends algorithm.</description>
    </inputParameter>
    <inputParameter>
      <name>linkingLengthIsComoving</name>
      <source>parameters</source>
      <description>Specifies whether or not the given linking length is in comoving units.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="haloMassFunction"     name="massFunctionIntrinsic" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=haloMassFunctionFofBias(massFunctionIntrinsic,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,massParticle,linkingLength,linkingLengthIsComoving,massInfiniteToMassSharpEdge)
    !![
    <inputParametersValidate source="parameters"/>
   <objectDestructor name="cosmologyParameters_" />
   <objectDestructor name="cosmologyFunctions_"  />
   <objectDestructor name="massFunctionIntrinsic"/>
   <objectDestructor name="darkMatterHaloScale_" />
   <objectDestructor name="darkMatterProfileDMO_"/>
   !!]
   return
  end function fofBiasConstructorParameters

  function fofBiasConstructorInternal(massFunctionIntrinsic,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,massParticle,linkingLength,linkingLengthIsComoving,massInfiniteToMassSharpEdge) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionFofBias} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionFofBias  )                        :: self
    class           (haloMassFunctionClass    ), target, intent(in   ) :: massFunctionIntrinsic
    class           (cosmologyParametersClass ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass  ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass ), target, intent(in   ) :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass), target, intent(in   ) :: darkMatterProfileDMO_
    double precision                                   , intent(in   ) :: massParticle           , massInfiniteToMassSharpEdge, &
         &                                                                linkingLength
    logical                                            , intent(in   ) :: linkingLengthIsComoving
    !![
    <constructorAssign variables="*massFunctionIntrinsic, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfileDMO_, massParticle, linkingLength, linkingLengthIsComoving, massInfiniteToMassSharpEdge"/>
    !!]

    return
  end function fofBiasConstructorInternal

  subroutine fofBiasDestructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionFofBias} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionFofBias), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunctionIntrinsic"/>
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine fofBiasDestructor

  double precision function fofBiasDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Coordinates             , only : coordinateSpherical  , assignment(=)
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentBasic   , treeNode
    use :: Mass_Distributions      , only : massDistributionClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (haloMassFunctionFofBias), intent(inout), target   :: self
    double precision                         , intent(in   )           :: time                                                                , &
         &                                                                mass
    type            (treeNode               ), intent(inout), optional :: node
    double precision                         , parameter               :: correctionCoefficient                                 =  0.220000d+0 ! More et al. (2011; eqn. B11).
    double precision                         , parameter               :: fofFractionalAccuracyExponent                         =  1.330000d+0 ! More et al. (2011; nu).
    double precision                         , parameter               :: numberDensityPercolation                              =  0.652960d+0 ! More et al. (2011; n_c).
    double precision                         , parameter               :: solutionAccuracy                                      =  1.000000d-6
    double precision                         , parameter               :: correctionFactorMaximum                               =  2.000000d+0
    integer                                  , parameter               :: iterationCountMaximum                                 =1000
    type            (treeNode               ), pointer                 :: nodeWork
    class           (nodeComponentBasic     ), pointer                 :: basic
    class           (massDistributionClass  ), pointer                 :: massDistribution_
    type            (coordinateSpherical    )                          :: coordinates
    integer                                                            :: iterationCount
    double precision                                                   :: fractionalAccuracyParameter                                         , &
         &                                                                densityProfileLogarithmicSlope                                      , &
         &                                                                gradientFractionalAccuracyParameterMass                             , &
         &                                                                gradientMassLogarithmicDensity                                      , &
         &                                                                gradientMassLogarithmicProbabilityCritical                          , &
         &                                                                gradientRadiusHaloMass                                              , &
         &                                                                massHaloInfinite                                                    , &
         &                                                                correctionFactor                                                    , &
         &                                                                massJacobian                                                        , &
         &                                                                gradientGradientMassLogarithmicProbabilityCriticalMass              , &
         &                                                                gradientProbabilityCriticalDensity                                  , &
         &                                                                radiusHalo                                                          , &
         &                                                                numberDensityCritical                                               , &
         &                                                                linkingLength                                                       , &
         &                                                                accuracyIteration                                                   , &
         &                                                                massHaloInfinitePrevious

    ! Compute linking length in physical coordinates.
    if (self%linkingLengthIsComoving) then
       ! Linking length is in comoving units, convert it.
       linkingLength=self%linkingLength*self%cosmologyFunctions_%expansionFactor(time)
    else
       ! Linking length is in physical units, just use it.
       linkingLength=self%linkingLength
    end if
    ! Create a work node.
    nodeWork => treeNode      (                 )
    basic    => nodeWork%basic(autoCreate=.true.)
    ! Initialize iterations.
    iterationCount          =0
    massHaloInfinitePrevious=mass
    accuracyIteration       =huge(0.0d0)
    do while (                                           &
         &     accuracyIteration > solutionAccuracy      &
         &    .and.                                      &
         &     iterationCount    < iterationCountMaximum &
         &   )
       ! Set the properties of the work node.
       call basic%massSet(massHaloInfinitePrevious)
       call basic%timeSet(time                    )
       ! Get the radius of the halo.
       radiusHalo                             =+self%darkMatterHaloScale_%radiusVirial                       (nodeWork)
       gradientRadiusHaloMass                 =+self%darkMatterHaloScale_%radiusVirialGradientLogarithmicMass(nodeWork) &
            &                                  *radiusHalo                                                              &
            &                                  /massHaloInfinitePrevious
       ! Compute friends-of-friends fractional accuracy parameter (Lsize from More et al. 2011).
       fractionalAccuracyParameter            =+2.0d0         &
            &                                  *radiusHalo    &
            &                                  /linkingLength
       gradientFractionalAccuracyParameterMass=+2.0d0                  &
            &                                  *gradientRadiusHaloMass &
            &                                  /linkingLength
       ! Compute negative a logarithmic slope of the density profile at the outer edge of the halo.
       coordinates                    =  [radiusHalo,0.0d0,0.0d0]
       massDistribution_              =>  self             %darkMatterProfileDMO_%get                  (nodeWork                      )
       densityProfileLogarithmicSlope =  +massDistribution_                      %densityGradientRadial(coordinates,logarithmic=.true.)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       ! Compute absolute value of the rate of change of logarithmic mass at halo outer edge as the
       ! percolation critical probability changes.
       numberDensityCritical                                 =+numberDensityPercolation          &
            &                                                 /linkingLength                 **3
       gradientProbabilityCriticalDensity                    =+Pi                                &
            &                                                 /6.0d0                             &
            &                                                 *linkingLength                 **3 &
            &                                                 /self%massParticle                 &
            &                                                 *exp(                              &
            &                                                      -Pi                           &
            &                                                      /6.0d0                        &
            &                                                      *numberDensityPercolation     &
            &                                                     )
       gradientMassLogarithmicDensity                        =+4.0d0                             &
            &                                                 *Pi                                &
            &                                                 *radiusHalo                    **3 &
            &                                                 /massHaloInfinitePrevious          &
            &                                                 /densityProfileLogarithmicSlope
       gradientMassLogarithmicProbabilityCritical            =+abs(                                    &
            &                                                      +gradientMassLogarithmicDensity     &
            &                                                      /gradientProbabilityCriticalDensity &
            &                                                     )
       gradientGradientMassLogarithmicProbabilityCriticalMass=+4.0d0                                 &
            &                                                 *Pi                                    &
            &                                                 *(                                     &
            &                                                   +3.0d0                               &
            &                                                   *gradientRadiusHaloMass              &
            &                                                   *radiusHalo                      **2 &
            &                                                   /massHaloInfinitePrevious            &
            &                                                   -radiusHalo                      **3 &
            &                                                   /massHaloInfinitePrevious        **2 &
            &                                                  )                                     &
            &                                                 /densityProfileLogarithmicSlope        &
            &                                                 /gradientProbabilityCriticalDensity
       ! Handle the abs() function in the definition of gradientMassLogarithmicProbabilityCritical - requires us to flip the sign of
       ! the derivative if gradientMassLogarithmicDensity (and, therefore, gradientMassLogarithmicProbabilityCritical) is negative.
       if (gradientMassLogarithmicDensity < 0.0d0)                    &
            & gradientGradientMassLogarithmicProbabilityCriticalMass= &
            & -gradientGradientMassLogarithmicProbabilityCriticalMass
       ! Compute the mass the halo would have under the friends-of-friends algorithm, given infinite resolution. This correction
       ! utilizes the percolation theory-motivated results of More et al. (2011; their equation B11). The minus sign below arises
       ! because More et al. define their α as the negative logarithmic slope of the halo density. We do not allow the
       ! correction factor to exceed some maximum value (set to unity above) since this would correspond to halos far below the
       ! resolution limit of the simulation and so is not realistic.
       correctionFactor=min(                                                                                    &
            &               +1.00                                                                               &
            &               -correctionCoefficient                                                              &
            &               *densityProfileLogarithmicSlope                                                     &
            &               /fractionalAccuracyParameter               **(1.0d0/fofFractionalAccuracyExponent)  &
            &               *gradientMassLogarithmicProbabilityCritical                                       , &
            &               +correctionFactorMaximum                                                            &
            &             )
       massHaloInfinite=+mass                             &
            &           /self%massInfiniteToMassSharpEdge &
            &           /correctionFactor
       ! Update current iteration accuracy.
       iterationCount   =+iterationCount                                 &
            &            +1
       accuracyIteration=+abs(massHaloInfinite-massHaloInfinitePrevious) &
            &            /   (massHaloInfinite+massHaloInfinitePrevious) &
            &            /0.5d0
       massHaloInfinitePrevious=massHaloInfinite
    end do
    if (iterationCount >= iterationCountMaximum) call Error_Report('failed to converge after maximum iterations'//{introspection:location})
    ! Compute the Jacobian of the transform. We currently ignore
    ! d/dm(densityProfileLogarithmicSlope) as we have no straightforward way to evaluate
    ! this. It should be a relatively minor correction as this term is likely small for
    ! realistic halos. If necessary it could be computed numerically though.
    massJacobian=+massHaloInfinite                                                                                      &
         &       /mass                                                                                                  &
         &       -mass                                                                                                  &
         &       /self%massInfiniteToMassSharpEdge                                                                      &
         &       /correctionFactor                                        ** 2                                          &
         &       *correctionCoefficient                                                                                 &
         &       *densityProfileLogarithmicSlope                                                                        &
         &       *(                                                                                                     &
         &         -gradientGradientMassLogarithmicProbabilityCriticalMass                                              &
         &         /fractionalAccuracyParameter                           **(     +1.0d0/fofFractionalAccuracyExponent) &
         &         +gradientMassLogarithmicProbabilityCritical                                                          &
         &         /fofFractionalAccuracyExponent                                                                       &
         &         /fractionalAccuracyParameter                           **(1.0d0+1.0d0/fofFractionalAccuracyExponent) &
         &         *gradientFractionalAccuracyParameterMass                                                             &
         &        )
    ! Clean up our work node.
    call nodeWork%destroy()
    deallocate(nodeWork)
    ! Compute the mass function from the intrinsic mass function evaluated at the halo mass
    ! corresponding to infinite resolution, and modified by the Jacobian of the mass
    ! transformation.
    fofBiasDifferential=+self%massFunctionIntrinsic%differential(time,massHaloInfinite,node=node) &
         &              /massJacobian
    return
  end function fofBiasDifferential
