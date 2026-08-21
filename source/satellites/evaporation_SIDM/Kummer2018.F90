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

  !+    Contributions to this file made by: Niusha Ahvazi, Andrew Benson, Claude.

  !!{RST
  Implementation of a satellite evaporation due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Numerical_Interpolation , only : interpolator2D
  use :: Numerical_Ranges        , only : rangeLattice

  !![
  <satelliteEvaporationSIDM name="satelliteEvaporationSIDMKummer2018" docformat="rst">
   <description>
   A satellite evaporation due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
   </description>
  </satelliteEvaporationSIDM>
  !!]
  type, extends(satelliteEvaporationSIDMClass) :: satelliteEvaporationSIDMKummer2018
     !!{RST
     Implementation of a satellite evaporation due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
     !!}
     private
     class           (cosmologyParametersClass ), pointer     :: cosmologyParameters_  => null()
     class           (darkMatterParticleClass  ), pointer     :: darkMatterParticle_   => null()
     class           (darkMatterHaloScaleClass ), pointer     :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer     :: darkMatterProfileDMO_ => null()
     type            (interpolator2D           ), allocatable :: evaporationFactor_
     double precision                                         :: xMaximum                       , velocityMinimum   , &
          &                                                      velocityMaximum                , fractionDarkMatter
     type            (rangeLattice             )              :: latticeX                       , latticeVelocity
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the evaporation factor." method="tabulate"          />
       <method description="Return the evaporation factor."   method="evaporationFactor" />
     </methods>
     !!]
     final     ::                      kummer2018Destructor
     procedure :: massLossRate      => kummer2018MassLossRate
     procedure :: tabulate          => kummer2018Tabulate
     procedure :: evaporationFactor => kummer2018EvaporationFactor
  end type satelliteEvaporationSIDMKummer2018

  interface satelliteEvaporationSIDMKummer2018
     !!{RST
     Constructors for the :galacticus-class:`satelliteEvaporationSIDMKummer2018` satellite evaporation due to dark matter self-interactions class.
     !!}
     module procedure kummer2018ConstructorParameters
     module procedure kummer2018ConstructorInternal
  end interface satelliteEvaporationSIDMKummer2018

contains

  function kummer2018ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteEvaporationSIDMKummer2018` satellite evaporation due to dark matter self-interactions class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteEvaporationSIDMKummer2018)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class(darkMatterParticleClass           ), pointer       :: darkMatterParticle_
    class(darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
  
    !![
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteEvaporationSIDMKummer2018(cosmologyParameters_,darkMatterParticle_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterParticle_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function kummer2018ConstructorParameters

  function kummer2018ConstructorInternal(cosmologyParameters_,darkMatterParticle_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteEvaporationSIDMKummer2018` satellite evaporation due to dark matter self-interactions class.
    !!}
    implicit none
    type (satelliteEvaporationSIDMKummer2018)                        :: self
    class(cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    class(darkMatterParticleClass           ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*cosmologyParameters_, *darkMatterParticle_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    ! Initialize the maximum tabulated x to an unphysical value. This will force tabulation on the first attempt to evaluate the
    ! evaporation factor.
    self%xMaximum       =-     1.0d0
    self%velocityMinimum=+huge(0.0d0)
    self%velocityMaximum=-huge(0.0d0)
    ! Evaluate the universal dark matter fraction.
    self%fractionDarkMatter=+(                                         & 
         &                    +self%cosmologyParameters_%OmegaMatter() &
         &                    -self%cosmologyParameters_%OmegaBaryon() &
         &                   )                                         &
         &                  /  self%cosmologyParameters_%OmegaMatter()
    return
  end function kummer2018ConstructorInternal

  subroutine kummer2018Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteEvaporationSIDMKummer2018` satellite evaporation due to dark matter self-interactions class.
    !!}
    implicit none
    type(satelliteEvaporationSIDMKummer2018), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterParticle_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine kummer2018Destructor

  double precision function kummer2018MassLossRate(self,node)
    !!{RST
    Return a evaporation for satellites due to dark matter self-interactions using the formulation of :cite:t:`kummer_effective_2018`.
    !!}
    use :: Coordinates                     , only : coordinateSpherical                        , coordinateCartesian        , assignment(=)
    use :: Galactic_Structure_Options      , only : coordinateSystemCartesian                  , radiusLarge                , massTypeDark
    use :: Galacticus_Nodes                , only : nodeComponentSatellite                     , nodeComponentBasic
    use :: Mass_Distributions              , only : massDistributionClass                      , kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal             , megaParsec                 , gigaYear     , massSolar
    use :: Vectors                         , only : Vector_Magnitude
    use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Numerical_Constants_Prefixes    , only : centi                                      , milli                      , kilo
    implicit none
    class           (satelliteEvaporationSIDMKummer2018), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    class           (nodeComponentSatellite            ), pointer       :: satellite
    class           (nodeComponentBasic                ), pointer       :: basic
    type            (treeNode                          ), pointer       :: nodeHost
    class           (massDistributionClass             ), pointer       :: massDistribution_     , massDistributionHost_      , &
         &                                                                 massDistributionDark
    class           (kinematicsDistributionClass       ), pointer       :: kinematics_           , kinematicsHost_
    double precision                                    , dimension(3)  :: position              , velocity
    double precision                                                    :: radiusOrbital         , speedOrbital               , &
         &                                                                 densityHost           , rateScattering             , &
         &                                                                 massBoundary          , radiusBoundary             , &
         &                                                                 velocityEscape        , speedRelative              , &
         &                                                                 velocityDispersionHost, velocityDispersionSatellite, &
         &                                                                 x                     , radiusHalfMass             , &
         &                                                                 velocityDispersion    , potentialEscape            , &
         &                                                                 speedHalfMass         , rateScatteringNormalization
    type            (coordinateSpherical                )                :: coordinates          , coordinatesHost            , &
         &                                                                  coordinatesBoundary  , coordinatesHalfMass
    type            (coordinateCartesian                )                :: coordinatesCartesian
    
    ! Set zero mass loss rate by default.
    kummer2018MassLossRate=0.0d0
    ! Evaluate satellite and host properties.
    nodeHost               =>  node                 %mergesWith      (                    )
    satellite              =>  node                 %satellite       (                    )
    massDistributionHost_  =>  nodeHost             %massDistribution(                    )
    position               =   satellite            %position        (                    )
    velocity               =   satellite            %velocity        (                    )
    radiusOrbital          =   Vector_Magnitude                      (position            )
    speedOrbital           =   Vector_Magnitude                      (velocity            )
    coordinatesCartesian   =                                          position
    densityHost            =   massDistributionHost_%density         (coordinatesCartesian)
    !![
    <objectDestructor name="massDistributionHost_"/>
    !!]
    ! If the dark matter particle is not self-interacting there is no scattering, so we can return immediately. The
    ! normalization of the scattering rate itself can not be evaluated here, as the cross section must be evaluated at the speed
    ! with which host particles arrive at the half-mass radius of the subhalo, which is not yet known.
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Self-interacting, so scattering can occur.
    class default
       return
    end select
    ! Find the escape velocity from the half-mass radius of the subhalo. This is equal to the potential difference between the
    ! half-mass radius and outer boundary of the subhalo, plus the potential difference from the outer boundary to infinity (for
    ! which we can treat the subhalo as a point mass).
    !
    ! Note that Kummer et al. (2018) suggest computing an effective escape velocity which is the escape velocity averaged over the
    ! mass profile of the subhalo, M⁻¹∫ρ(r) v(r) d³r. Evaluating that integral would be computationally expensive. Instead, we use
    ! the escape velocity at the half-mass radius of the subhalo. For a scale-free Hernquist profile, M⁻¹∫ρ(r) v(r) d³r =
    ! 8√2/15=0.754, while v(r½=1+√2)=√(2/(2+√2))=0.765 - so the two choices result in very similar effective escape velocities.
    basic        =>    node      %basic             ()
    massBoundary =  +min(                               &
         &              satellite%boundMass         (), &
         &              basic    %     mass         ()  &
         &             )                                &
         &          *   self     %fractionDarkMatter
    if (massBoundary > 0.0d0) then
       massDistribution_    => node%massDistribution(                           )
       massDistributionDark => node%massDistribution(massType=      massTypeDark)
       if (      massBoundary > massDistributionDark%massEnclosedBySphere(self%darkMatterHaloScale_%radiusVirial(node))) then
          radiusBoundary=self                %darkMatterHaloScale_%radiusVirial       (     node        )
       else
          radiusBoundary=massDistributionDark                     %radiusEnclosingMass(mass=massBoundary)
       end if
       if (0.5d0*massBoundary > massDistributionDark%massEnclosedBySphere(self%darkMatterHaloScale_%radiusVirial(node))) then
          radiusHalfMass=self                %darkMatterHaloScale_%radiusVirial       (     node        )
       else
          radiusHalfMass=massDistributionDark                     %radiusEnclosingMass(mass=0.5d0*massBoundary)
       end if
       if (radiusBoundary < 0.5d0*radiusLarge) then
          coordinatesBoundary=[radiusBoundary,0.0d0,0.0d0]
          coordinatesHalfMass=[radiusHalfMass,0.0d0,0.0d0]
          potentialEscape    =+massDistribution_%potentialDifference(coordinatesBoundary,coordinatesHalfMass) &
               &              +gravitationalConstant_internal                                                 &
               &              *massDistribution_%massEnclosedBySphere(radiusBoundary)                         &
               &              /radiusBoundary
          if (potentialEscape > 0.0d0) then
             velocityEscape=sqrt(2.0d0*potentialEscape)
          else
             velocityEscape=0.0d0
          end if
       else
          velocityEscape=0.0d0
       end if
       !![
       <objectDestructor name="massDistribution_"   />
       <objectDestructor name="massDistributionDark"/>
       !!]
       ! Find the combined velocity dispersion of satellite and host.
       massDistribution_           =>  self                 %darkMatterProfileDMO_%get                   (node                                                       )
       massDistributionHost_       =>  self                 %darkMatterProfileDMO_%get                   (nodeHost                                                   ) 
       kinematics_                 =>  massDistribution_                          %kinematicsDistribution(                                                           )
       kinematicsHost_             =>  massDistributionHost_                      %kinematicsDistribution(                                                           )
       coordinates                 =  [radiusHalfMass,0.0d0,0.0d0]
       coordinatesHost             =  [radiusOrbital ,0.0d0,0.0d0]
       velocityDispersionHost      =  +kinematicsHost_                            %velocityDispersion1D  (coordinatesHost,massDistributionHost_,massDistributionHost_)
       velocityDispersionSatellite =  +kinematics_                                %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
       velocityDispersion          =  +sqrt(                                &
            &                               +velocityDispersionHost     **2 &
            &                               +velocityDispersionSatellite**2 &
            &                              )
       !![
       <objectDestructor name="massDistribution_"    />
       <objectDestructor name="massDistributionHost_"/>
       <objectDestructor name="kinematics_"          />
       <objectDestructor name="kinematicsHost_"      />
       !!]
       ! Find the effective speed of the subhalo relative to the host - this is the speed of host particles in the rest frame of
       ! the subhalo far from the subhalo, broadened by the velocity dispersion as suggested in equation (A4) of Kummer et
       ! al. (2018). Host particles reaching the half-mass radius of the subhalo are further boosted by the energy they gain in
       ! falling in, arriving with speed √(speedRelative²+velocityEscape²).
       speedRelative=+sqrt(                       &
            &              +speedOrbital      **2 &
            &              +velocityDispersion**2 &
            &             )
       ! Evaluate the ratio, x, of the escape velocity from the subhalo to that relative speed. Note that x is defined in terms
       ! of the relative speed far from the subhalo, and *not* the speed at the point of scattering: the infall boost is already
       ! accounted for within the definitions of the critical scattering angle and of the evaporation factor, both of which are
       ! functions of this x. To see this, note that a subhalo particle recoiling with speed √(speedRelative²+velocityEscape²)
       ! sin(θ/2) escapes the subhalo only for θ > π-θ_c, while the scattered host particle escapes only for θ < θ_c, where
       ! θ_c = cos⁻¹([x²-1]/[x²+1]) is the critical angle tabulated in `kummer2018Tabulate`. Net evaporation therefore requires
       ! π-θ_c < θ < θ_c - the range integrated over in that tabulation - which is non-empty only for x < 1. For a constant
       ! cross section the tabulated integral is χ_e = (1-x²)/(1+x²), the closed form given by Kummer et al. (2018), which fixes
       ! the sense of x: small x is a weakly bound subhalo moving quickly through its host, from which every scattering
       ! evaporates a particle.
       x            =+velocityEscape &
            &        /speedRelative
       ! Evaporation occurs for x<1.
       if (x < 1.0d0) then
          ! Find the speed with which host particles arrive at the half-mass radius of the subhalo: the effective relative speed,
          ! boosted by the energy they gain in falling in. This is the typical speed of the particles which actually scatter, and
          ! so is the speed at which the cross section - which may be velocity dependent - is evaluated, both here and in the
          ! tabulation of the evaporation factor.
          speedHalfMass=+sqrt(                   &
               &              +speedRelative **2 &
               &              +velocityEscape**2 &
               &             )
          ! Compute the normalization of the scattering rate in units such that when multiplied by a velocity in km s⁻¹, and a
          ! density in units of M⊙ Mpc⁻³, we get a rate in units of Gyr⁻¹.
          rateScatteringNormalization=0.0d0
          select type (darkMatterParticle_ => self%darkMatterParticle_)
          class is (darkMatterParticleSelfInteractingDarkMatter)
             rateScatteringNormalization=+darkMatterParticle_%crossSectionSelfInteraction(speedHalfMass) &
                  &                      *centi    **2/milli                                             & ! Convert cross-section from cm² g⁻¹ to m² kg⁻¹.
                  &                      *kilo                                                           & ! Convert velocity from km s⁻¹ to m s⁻¹.
                  &                      *massSolar   /megaParsec**3                                     & ! Convert density from M⊙ Mpc⁻³ to kg m⁻³.
                  &                      *gigaYear                                                         ! Convert rate from s⁻¹ to Gyr⁻¹.
          end select
          ! Evaluate the scattering rate and mass loss rate. Note that the flux of host particles through the subhalo is set by
          ! the speed with which the subhalo moves through its host, not by the speed at the point of scattering.
          rateScattering        =+speedOrbital                             &
               &                 *densityHost                              &
               &                 *rateScatteringNormalization
          kummer2018MassLossRate=-massBoundary                             &
               &                 *rateScattering                           &
               &                 *self%evaporationFactor(x,speedHalfMass)
       end if
    end if
    return
  end function kummer2018MassLossRate

  subroutine kummer2018Tabulate(self,xMaximum,velocityMinimum,velocityMaximum)
    !!{RST
    Tabulate the evaporation factor, :math:`\chi_\mathrm{d}`.
    !!}
    use :: Dark_Matter_Particles   , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Range_Pinned                               , gridSchemePerUnit, gridSchemePerDecade
    use :: Numerical_Constants_Math, only : Pi
    use :: Error                   , only : Error_Report
    implicit none
    class           (satelliteEvaporationSIDMKummer2018         ), intent(inout)                 :: self
    double precision                                             , intent(in   )                 :: xMaximum                  , velocityMinimum   , &
         &                                                                                          velocityMaximum
    class           (darkMatterParticleSelfInteractingDarkMatter), pointer                       :: darkMatterParticleSIDM_
    double precision                                             , allocatable  , dimension(:,:) :: evaporationFactor_
    double precision                                             , allocatable  , dimension(:  ) :: x                         , velocity
    integer                                                      , parameter                     :: countPerUnit           =10, countPerDex    =10
    type            (integrator                                 )                                :: integrator_
    double precision                                                                             :: thetaCritical
    integer                                                                                      :: i                         , j                 , &
         &                                                                                          countX                    , countVelocity

    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Tabulate the χ_d factor.
       ! Pin both axes to absolute lattices, so that the tabulation - and every value interpolated from it - depends only on
       ! which anchored intervals the request fell in, rather than on the requested x and speed themselves. The x axis is
       ! linear and runs from zero, so the *whole* range [0,xMaximum] is the target - `limitMinimum` only clamps a range which
       ! would otherwise extend below zero, it does not make one reach down to it - and its lower edge is clamped there; the
       ! speed axis
       ! is logarithmic. Taking the union with the lattices already in use guarantees the tabulated range never shrinks.
       !
       ! The table is rebuilt rather than extended. Each entry depends only on its own two abscissae, so entries could in
       ! principle be carried over, but the computed values are held only inside the `interpolator2D` built from them and are
       ! not recoverable from it; doing so would mean storing the raw array alongside. Left as a possible follow-up.
       self%latticeX       =Range_Pinned(                                                  &
            &                                           [0.0d0,xMaximum]                 , &
            &                                           countPerUnit                     , &
            &                                           gridSchemePerUnit                , &
            &                            marginOffset  =0.0d0                            , &
            &                            limitMinimum  =0.0d0                            , &
            &                            latticeCurrent=self%latticeX                      &
            &                           )
       self%latticeVelocity=Range_Pinned(                                                  &
            &                                           [velocityMinimum,velocityMaximum], &
            &                                           countPerDex                      , &
            &                                           gridSchemePerDecade              , &
            &                            latticeCurrent=self%latticeVelocity               &
            &                           )
       self%xMaximum       =self%latticeX       %maximum()
       self%velocityMaximum=self%latticeVelocity%maximum()
       self%velocityMinimum=self%latticeVelocity%minimum()
       countX       =self%latticeX       %count
       countVelocity=self%latticeVelocity%count
       allocate(x                 (countX              ))
       allocate(velocity          (       countVelocity))
       allocate(evaporationFactor_(countX,countVelocity))
       x           =self%latticeX       %values()
       velocity    =self%latticeVelocity%values()
       darkMatterParticleSIDM_ => darkMatterParticle_
       integrator_             =  integrator(integrandEvaporationFactor,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
       do i=1,countX
          thetaCritical=acos(                  & ! Equation (2) from Kummer et al. (2018).
               &             +(+x(i)**2-1.0d0) &
               &             /(+x(i)**2+1.0d0) &
               &            )
          do j=1,countVelocity
             evaporationFactor_(i,j)=+integrator_        %integrate                  (Pi-thetaCritical,thetaCritical) &
                  &                  /darkMatterParticle_%crossSectionSelfInteraction(velocity(j)                   )
          end do
       end do
       if (allocated(self%evaporationFactor_)) deallocate(self%evaporationFactor_)
       allocate(self%evaporationFactor_)
       self%evaporationFactor_=interpolator2D(x,velocity,evaporationFactor_)
    class default
       call Error_Report('this class expects an SIDM particle'//{introspection:location})
    end select
    return

  contains
    
    double precision function integrandEvaporationFactor(theta)
      !!{RST
      The integrand used to evaluate :math:`\chi_\mathrm{e}` from :cite:t:`kummer_effective_2018`.
      !!}
      implicit none
      double precision, intent(in   ) :: theta

      integrandEvaporationFactor=+darkMatterParticleSIDM_%crossSectionSelfInteractionDifferential(theta,velocity(j))
      return
    end function integrandEvaporationFactor

  end subroutine kummer2018Tabulate
  
  double precision function kummer2018EvaporationFactor(self,x,velocityRelative)
    !!{RST
    Evaluate the evaporation factor from :cite:t:`kummer_effective_2018`.
    !!}
    implicit none
    class           (satelliteEvaporationSIDMKummer2018), intent(inout) :: self
    double precision                                    , intent(in   ) :: x
    double precision                                    , intent(in   ) :: velocityRelative
    double precision                                    , parameter     :: xCritical       =1.0d0

    if     (                                         &
         &   x                > self%xMaximum        &
         &  .or.                                     &
         &   velocityRelative > self%velocityMaximum &
         &  .or.                                     &
         &   velocityRelative < self%velocityMinimum &
         & ) then
       if (x < xCritical)                                                          &
            & call self%tabulate(                                                  &
            &                    x+1.0d0                                         , &
            &                    min(0.5d0*velocityRelative,self%velocityMinimum), &
            &                    max(2.0d0*velocityRelative,self%velocityMaximum)  &
            &                   )
    end if
    if (x < xCritical) then
       kummer2018EvaporationFactor=self%evaporationFactor_%interpolate(x,velocityRelative)
    else
       kummer2018EvaporationFactor=0.0d0
    end if
    return
  end function kummer2018EvaporationFactor
