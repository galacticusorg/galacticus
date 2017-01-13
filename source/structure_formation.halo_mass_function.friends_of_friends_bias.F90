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

!% Contains a  module which  implements a  dark matter  halo mass function  class which  modifies another  mass function  to mimic
!% systematic errors arising in the friends-of-friends halo finding algorithm.

  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  
  !# <haloMassFunction name="haloMassFunctionFofBias">
  !#  <description>
  !#   The halo mass function is computed by modifying another halo mass function to mimic systematic
  !#   errors arising in the friends-of-friends halo finding algorithm. Specifically, a
  !#   systematic shift in mass motivated by the results of percolation theory \cite[their eqn. B11]{more_overdensity_2011}
  !#   is applied. In particular, $M_{\rm particle}=${\normalfont \ttfamily [massParticle]} is the mass of the particle in the simulation
  !#   to which the friends-of-friends algorithm was applied.
  !#  </description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionFofBias
     !% A halo mass function class which modifies another mass function to mimic systematic errors arising in the
     !% friends-of-friends halo finding algorithm.
     private
     double precision                                    :: massParticle           , massInfiniteToMassSharpEdge, &
          &                                                 linkingLength
     logical                                             :: linkingLengthIsComoving
     class           (haloMassFunctionClass   ), pointer :: massFunctionIntrinsic
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
     class           (darkMatterProfileClass  ), pointer :: darkMatterProfile_
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
    contains
     final     ::                 fofBiasDestructor
     procedure :: differential => fofBiasDifferential
  end type haloMassFunctionFofBias

  interface haloMassFunctionFofBias
     !% Constructors for the {\normalfont \ttfamily fofBias} halo mass function class.
     module procedure fofBiasConstructorParameters
     module procedure fofBiasConstructorInternal
  end interface haloMassFunctionFofBias

contains

  function fofBiasConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily fofBias} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionFofBias)                :: fofBiasConstructorParameters
    type(inputParameters        ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>fofBiasConstructorParameters%massParticle</variable>
    !#   <description>Parameter $M_{\rm particle}$ appearing in model for friends-of-friends errors in the halo mass function.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massInfiniteToMassSharpEdge</name>
    !#   <source>parameters</source>
    !#   <variable>fofBiasConstructorParameters%massInfiniteToMassSharpEdge</variable>
    !#   <defaultValue>0.98d0</defaultValue>
    !#   <defaultSource>\cite[estimate based on comments in text]{more_overdensity_2011}</defaultSource>
    !#   <description>The ratio of the friends-of-friends mass in the limit of infinite number of particles to the mass of the halo enclosed within a sharp-edged sphere bounding an isodensity surface equal to the critical density for percolation.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>linkingLength</name>
    !#   <source>parameters</source>
    !#   <variable>fofBiasConstructorParameters%linkingLength</variable>
    !#   <description>The linking length (in physical Mpc) used in the friends-of-friends algorithm.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>linkingLengthIsComoving</name>
    !#   <source>parameters</source>
    !#   <variable>fofBiasConstructorParameters%linkingLengthIsComoving</variable>
    !#   <description>Specifies whether or not the given linking length is in comoving units.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="fofBiasConstructorParameters%cosmologyParameters_"  source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="fofBiasConstructorParameters%cosmologyFunctions_"   source="parameters"/>
    !# <objectBuilder class="haloMassFunction"    name="fofBiasConstructorParameters%massFunctionIntrinsic" source="parameters"/>
    fofBiasConstructorParameters%darkMatterHaloScale_ => darkMatterHaloScale()
    fofBiasConstructorParameters%darkMatterProfile_   => darkMatterProfile  ()
   return
  end function fofBiasConstructorParameters

  function fofBiasConstructorInternal(massFunctionIntrinsic,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfile_,massParticle,linkingLength,massInfiniteToMassSharpEdge)
    !% Internal constructor for the {\normalfont \ttfamily fofBias} halo mass function class.
    implicit none
    type            (haloMassFunctionFofBias )                        :: fofBiasConstructorInternal
    class           (haloMassFunctionClass   ), target, intent(in   ) :: massFunctionIntrinsic
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass), target, intent(in   ) :: darkMatterHaloScale_
    class           (darkMatterProfileClass  ), target, intent(in   ) :: darkMatterProfile_
    double precision                                  , intent(in   ) :: massParticle                , massInfiniteToMassSharpEdge   , &
         &                                                               linkingLength
    
    fofBiasConstructorInternal%massFunctionIntrinsic       => massFunctionIntrinsic
    fofBiasConstructorInternal%cosmologyParameters_        => cosmologyParameters_
    fofBiasConstructorInternal%cosmologyFunctions_         => cosmologyFunctions_
    fofBiasConstructorInternal%darkMatterHaloScale_        => darkMatterHaloScale_
    fofBiasConstructorInternal%darkMatterProfile_          => darkMatterProfile_
    fofBiasConstructorInternal%massParticle                =  massParticle
    fofBiasConstructorInternal%linkingLength               =  linkingLength
    fofBiasConstructorInternal%massInfiniteToMassSharpEdge =  massInfiniteToMassSharpEdge
    return
  end function fofBiasConstructorInternal
  
  subroutine fofBiasDestructor(self)
    !% Destructor for the {\normalfont \ttfamily fofBias} halo mass function class.
    implicit none
    type(haloMassFunctionFofBias), intent(inout) :: self

    !# <objectDestructor name="self%massFunctionIntrinsic" />
    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    return
  end subroutine fofBiasDestructor

  double precision function fofBiasDifferential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    use Galacticus_Error
    use Galacticus_Nodes
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionFofBias), intent(inout) :: self
    double precision                         , intent(in   ) :: time                                                                , &
         &                                                      mass
    double precision                         , parameter     :: correctionCoefficient                                 =  0.220000d+0 ! More et al. (2011; eqn. B11).
    double precision                         , parameter     :: fofFractionalAccuracyExponent                         =  1.330000d+0 ! More et al. (2011; nu).
    double precision                         , parameter     :: numberDensityPercolation                              =  0.652960d+0 ! More et al. (2011; n_c).
    double precision                         , parameter     :: solutionAccuracy                                      =  1.000000d-6
    double precision                         , parameter     :: correctionFactorMaximum                               =  2.000000d+0
    integer                                  , parameter     :: iterationCountMaximum                                 =1000
    type            (treeNode               ), pointer       :: node
    class           (nodeComponentBasic     ), pointer       :: basic
    integer                                                  :: iterationCount
    double precision                                         :: fractionalAccuracyParameter                                         , &
         &                                                      densityProfileLogarithmicSlope                                      , &
         &                                                      gradientFractionalAccuracyParameterMass                             , &
         &                                                      gradientMassLogarithmicDensity                                      , &
         &                                                      gradientMassLogarithmicProbabilityCritical                          , &
         &                                                      gradientRadiusHaloMass                                              , &
         &                                                      massHaloInfinite                                                    , &
         &                                                      correctionFactor                                                    , &
         &                                                      massJacobian                                                        , &
         &                                                      gradientGradientMassLogarithmicProbabilityCriticalMass              , &
         &                                                      gradientProbabilityCriticalDensity                                  , &
         &                                                      radiusHalo                                                          , &
         &                                                      numberDensityCritical                                               , &
         &                                                      linkingLength                                                       , &
         &                                                      accuracyIteration                                                   , &
         &                                                      massHaloInfinitePrevious
    
    ! Compute linking length in physical coordinates.
    if (self%linkingLengthIsComoving) then
       ! Linking length is in comoving units, convert it.
       linkingLength=self%linkingLength*self%cosmologyFunctions_%expansionFactor(time)
    else
       ! Linking length is in physical units, just use it.
       linkingLength=self%linkingLength
    end if
    ! Create a work node.
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
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
       radiusHalo                             =+self%darkMatterHaloScale_%virialRadius                       (node)
       gradientRadiusHaloMass                 =+self%darkMatterHaloScale_%virialRadiusGradientLogarithmicMass(node) &
            &                                  *radiusHalo                                                          &
            &                                  /massHaloInfinitePrevious
       ! Compute friends-of-friends fractional accuracy parameter (Lsize from More et al. 2011).
       fractionalAccuracyParameter            =+2.0d0         &
            &                                  *radiusHalo    &
            &                                  /linkingLength
       gradientFractionalAccuracyParameterMass=+2.0d0                  &
            &                                  *gradientRadiusHaloMass &
            &                                  /linkingLength    
       ! Compute negative a logarithmic slope of the density profile at the outer edge of the halo.
       densityProfileLogarithmicSlope=+self%darkMatterProfile_ %densityLogSlope(node,radiusHalo)
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
       ! because More et al. define their alpha as the negative logarithmic slope of the halo density. We do not allow the
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
    if (iterationCount >= iterationCountMaximum) call Galacticus_Error_Report('fofBiasDifferential','failed to converge after maximum iterations')    
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
    call node%destroy()
    deallocate(node)
    ! Compute the mass function from the intrinsic mass function evaluated at the halo mass
    ! corresponding to infinte resolution, and modified by the Jacobian of the mass
    ! transformation.
    fofBiasDifferential=+self%massFunctionIntrinsic%differential(time,massHaloInfinite)&
         &                /massJacobian
    return
  end function fofBiasDifferential
