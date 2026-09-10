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
  Implementation of a satellite deceleration due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Numerical_Interpolation , only : interpolator2D
  use :: Numerical_Ranges        , only : rangeLattice

  !![
  <satelliteDecelerationSIDM name="satelliteDecelerationSIDMKummer2018" docformat="rst">
   <description>
   A satellite deceleration due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
   </description>
  </satelliteDecelerationSIDM>
  !!]
  type, extends(satelliteDecelerationSIDMClass) :: satelliteDecelerationSIDMKummer2018
     !!{RST
     Implementation of a satellite deceleration due to dark matter self-interactions using the model of :cite:t:`kummer_effective_2018`.
     !!}
     private
     class           (cosmologyParametersClass ), pointer     :: cosmologyParameters_  => null()
     class           (darkMatterParticleClass  ), pointer     :: darkMatterParticle_   => null()
     class           (darkMatterHaloScaleClass ), pointer     :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer     :: darkMatterProfileDMO_ => null()
     type            (interpolator2D           ), allocatable :: decelerationFactor_
     double precision                                         :: xMaximum                       , velocityMinimum   , &
          &                                                      velocityMaximum                , fractionDarkMatter
     type            (rangeLattice             )              :: latticeX                       , latticeVelocity
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the deceleration factor."                  method="tabulate"           />
       <method description="Return the dynamical-friction deceleration factor." method="decelerationFactor" />
     </methods>
     !!]
     final     ::                       kummer2018Destructor
     procedure :: acceleration       => kummer2018Acceleration
     procedure :: tabulate           => kummer2018Tabulate
     procedure :: decelerationFactor => kummer2018DecelerationFactor
  end type satelliteDecelerationSIDMKummer2018

  interface satelliteDecelerationSIDMKummer2018
     !!{RST
     Constructors for the :galacticus-class:`satelliteDecelerationSIDMKummer2018` satellite deceleration due to dark matter self-interactions class.
     !!}
     module procedure kummer2018ConstructorParameters
     module procedure kummer2018ConstructorInternal
  end interface satelliteDecelerationSIDMKummer2018

contains

  function kummer2018ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteDecelerationSIDMKummer2018` satellite deceleration due to dark matter self-interactions class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteDecelerationSIDMKummer2018)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    class(darkMatterParticleClass            ), pointer       :: darkMatterParticle_
    class(darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
  
    !![
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteDecelerationSIDMKummer2018(cosmologyParameters_,darkMatterParticle_,darkMatterHaloScale_,darkMatterProfileDMO_)
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
    Internal constructor for the :galacticus-class:`satelliteDecelerationSIDMKummer2018` satellite deceleration due to dark matter self-interactions class.
    !!}
    implicit none
    type (satelliteDecelerationSIDMKummer2018)                        :: self
    class(cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class(darkMatterParticleClass            ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*cosmologyParameters_, *darkMatterParticle_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    ! Initialize the maximum tabulated x to an unphysical value. This will force tabulation on the first attempt to evaluate the
    ! deceleration factor.
    self%xMaximum=-     1.0d0
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
    Destructor for the :galacticus-class:`satelliteDecelerationSIDMKummer2018` satellite deceleration due to dark matter self-interactions class.
    !!}
    implicit none
    type(satelliteDecelerationSIDMKummer2018), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterParticle_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine kummer2018Destructor

  function kummer2018Acceleration(self,node)
    !!{RST
    Return a deceleration for satellites due to dark matter self-interactions using the formulation of :cite:t:`kummer_effective_2018`.
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
    double precision                                     , dimension(3)  :: kummer2018Acceleration
    class           (satelliteDecelerationSIDMKummer2018), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentSatellite             ), pointer       :: satellite
    class           (nodeComponentBasic                 ), pointer       :: basic
    type            (treeNode                           ), pointer       :: nodeHost
    class           (massDistributionClass              ), pointer       :: massDistribution_          , massDistributionHost_ , &
         &                                                                  massDistributionDark
    class           (kinematicsDistributionClass        ), pointer       :: kinematics_                , kinematicsHost_
    double precision                                     , dimension(3)  :: position                   , velocity
    double precision                                                     :: radiusOrbital              , speedOrbital          , &
         &                                                                  densityHost                , rateScattering        , &
         &                                                                  massBoundary               , radiusBoundary        , &
         &                                                                  velocityEscape             , velocityDispersionHost, &
         &                                                                  velocityDispersionSatellite, radiusHalfMass        , &
         &                                                                  velocityDispersion         , dispersionFactor      , &
         &                                                                  potentialEscape            , speedHalfMass         , &
         &                                                                  rateScatteringNormalization, x
    type            (coordinateSpherical                )                :: coordinates                , coordinatesHost       , &
         &                                                                  coordinatesBoundary        , coordinatesHalfMass
    type            (coordinateCartesian                )                :: coordinatesCartesian
    
    ! Set zero acceleration by default.
    kummer2018Acceleration=0.0d0
    ! Evaluate satellite and host properties.
    nodeHost                     =>  node                               %mergesWith(                                           )
    satellite                    =>  node                               %satellite (                                           )
    massDistributionHost_        =>  nodeHost                           %massDistribution(                                     )
    position                     =   satellite                          %position  (                                           )
    velocity                     =   satellite                          %velocity  (                                           )
    radiusOrbital                =   Vector_Magnitude                              (         position                          )
    speedOrbital                 =   Vector_Magnitude                              (         velocity                          )
    coordinatesCartesian         =                                                           position
    densityHost                  =   massDistributionHost_%density                 (coordinatesCartesian)
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
       massDistribution_    => node%massDistribution(                     )
       massDistributionDark => node%massDistribution(massType=massTypeDark)
       if (      massBoundary > massDistributionDark%massEnclosedBySphere(self%darkMatterHaloScale_%radiusVirial(node))) then
          radiusBoundary=self                %darkMatterHaloScale_%radiusVirial       (           node        )
       else
          radiusBoundary=massDistributionDark                     %radiusEnclosingMass(mass=      massBoundary)
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
       ! Evaluate the ratio, x, of the escape velocity from the subhalo to the speed of the subhalo relative to the host. Note
       ! that the relative speed here is the speed far from the subhalo - the boost which host particles acquire as they fall
       ! in to the half-mass radius, arriving with speed √(speedOrbital²+velocityEscape²), is already accounted for within the
       ! definitions of the critical scattering angle and of the deceleration factor, both of which are functions of this x. To
       ! see this, note that a particle leaving the scattering with speed √(speedOrbital²+velocityEscape²) cos(θ/2) escapes the
       ! subhalo only for θ < θ_c = cos⁻¹([x²-1]/[x²+1]), which is precisely the θ_c tabulated in `kummer2018Tabulate`.
       ! Expanding that tabulated integral at large x for a constant cross section reproduces the asymptotic form
       ! 1-8/3x²+16/5x⁴ used in `kummer2018DecelerationFactor`, which fixes the sense of x: large x is a deeply bound subhalo
       ! from which nothing escapes, so that momentum transfer is complete and χ_d → 1.
       x            =+velocityEscape &
            &        /speedOrbital
       ! Find the combined velocity dispersion of satellite and host, and evaluate the correction factor given in Appendix A of
       ! Kummer et al. (2018).
       massDistribution_           =>  self                 %darkMatterProfileDMO_%get                   (node           )
       massDistributionHost_       =>  self                 %darkMatterProfileDMO_%get                   (nodeHost       )       
       kinematics_                 =>  massDistribution_                          %kinematicsDistribution(               )
       kinematicsHost_             =>  massDistributionHost_                      %kinematicsDistribution(               )
       coordinates                 =  [radiusHalfMass,0.0d0,0.0d0]
       coordinatesHost             =  [radiusOrbital ,0.0d0,0.0d0]
       velocityDispersionHost      =  +kinematicsHost_                            %velocityDispersion1D  (coordinatesHost,massDistributionHost_,massDistributionHost_)
       velocityDispersionSatellite =  +kinematics_                                %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
       velocityDispersion          =  +sqrt(                                &
            &                               +velocityDispersionHost     **2 &
            &                               +velocityDispersionSatellite**2 &
            &                              )
       dispersionFactor            =  +1.0d0                  &
            &                         /(                      &
            &                           +1.0d0                &
            &                           +(                    &
            &                             +velocityDispersion &
            &                             /speedOrbital       &
            &                            )**3                 &
            &                          )
       !![
       <objectDestructor name="massDistribution_"    />
       <objectDestructor name="massDistributionHost_"/>
       <objectDestructor name="kinematics_"          />
       <objectDestructor name="kinematicsHost_"      />
       !!]
       ! Find the speed with which host particles arrive at the half-mass radius of the subhalo: the speed of the subhalo
       ! relative to the host, broadened by the velocity dispersion and boosted by the energy the particles gain in falling in.
       ! This is the typical speed of the particles which actually scatter, and so is the speed at which the cross section -
       ! which may be velocity dependent - is evaluated, both here and in the tabulation of the deceleration factor. Note that
       ! this is distinct from the speed which enters x and the dispersion correction factor above: those are functions of the
       ! orbital speed alone, with the infall boost and the dispersion accounted for within their own definitions.
       speedHalfMass                =  +sqrt(                       &
            &                                +speedOrbital      **2 &
            &                                +velocityEscape    **2 &
            &                                +velocityDispersion**2 &
            &                               )
       ! Compute the normalization of the scattering rate in units such that when multiplied by a velocity in km s⁻¹, and a
       ! density in units of M⊙ Mpc⁻³, we get a rate in units of Gyr⁻¹.
       rateScatteringNormalization  =  +0.0d0
       select type (darkMatterParticle_ => self%darkMatterParticle_)
       class is (darkMatterParticleSelfInteractingDarkMatter)
          rateScatteringNormalization=+darkMatterParticle_%crossSectionSelfInteraction(speedHalfMass) &
               &                      *centi    **2/milli                                            & ! Convert cross-section from cm² g⁻¹ to m² kg⁻¹.
               &                      *kilo                                                          & ! Convert velocity from km s⁻¹ to m s⁻¹.
               &                      *massSolar   /megaParsec**3                                    & ! Convert density from M⊙ Mpc⁻³ to kg m⁻³.
               &                      *gigaYear                                                        ! Convert rate from s⁻¹ to Gyr⁻¹.
       end select
       ! Evaluate the scattering rate and acceleration. Note that the flux of host particles through the subhalo is set by the
       ! speed with which the subhalo moves through its host, not by the speed at the point of scattering.
       rateScattering               =  +speedOrbital                             &
            &                          *densityHost                              &
            &                          *rateScatteringNormalization
       kummer2018Acceleration       =  -velocity                                 &
            &                          *rateScattering                           &
            &                          *self%decelerationFactor(x,speedHalfMass) &
            &                          *dispersionFactor
    end if
    return
  end function kummer2018Acceleration

  subroutine kummer2018Tabulate(self,xMaximum,velocityMinimum,velocityMaximum)
    !!{RST
    Tabulate the deceleration factor, :math:`\chi_\mathrm{d}`.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Range_Pinned                               , gridSchemePerUnit, gridSchemePerDecade
    implicit none
    class           (satelliteDecelerationSIDMKummer2018        ), intent(inout)                 :: self
    double precision                                             , intent(in   )                 :: xMaximum                  , velocityMinimum   , &
         &                                                                                          velocityMaximum
    class           (darkMatterParticleSelfInteractingDarkMatter), pointer                       :: darkMatterParticleSIDM_
    double precision                                             , allocatable  , dimension(:,:) :: decelerationFactor_
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
       allocate(x                  (countX       ))
       allocate(velocity                  (       countVelocity))
       allocate(decelerationFactor_(countX,countVelocity))
       x           =self%latticeX       %values()
       velocity    =self%latticeVelocity%values()
       darkMatterParticleSIDM_ => darkMatterParticle_
       integrator_             =  integrator(integrandDecelerationFactor,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
       do i=1,countX
          thetaCritical=+acos(                  & ! Equation (2) from Kummer et al. (2018).
               &              +(+x(i)**2-1.0d0) &
               &              /(+x(i)**2+1.0d0) &
               &             )
          do j=1,countVelocity             
             decelerationFactor_(i,j)=+1.0d0                                                                &
               &                      -integrator_        %integrate                  (0.0d0,thetaCritical) &
               &                      /darkMatterParticle_%crossSectionSelfInteraction(velocity(j)               ) &
               &                      *sqrt(2.0d0)
          end do
       end do
       if (allocated(self%decelerationFactor_)) deallocate(self%decelerationFactor_)
       allocate(self%decelerationFactor_)
       self%decelerationFactor_=interpolator2D(x,velocity,decelerationFactor_)
    end select
    return

  contains
    
    double precision function integrandDecelerationFactor(theta)
      !!{RST
      The integrand used to evaluate :math:`\chi_\mathrm{d}` from :cite:t:`kummer_effective_2018`.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: theta

      integrandDecelerationFactor=+                                   cos(0.5d0*theta)                              &
           &                      *sqrt(1.0d0-x(i)**2+(1.0d0+x(i)**2)*cos(      theta))                             &
           &                      *(                                                                                &
           &                        +darkMatterParticleSIDM_%crossSectionSelfInteractionDifferential(  +theta,velocity(j)) &
           &                        +darkMatterParticleSIDM_%crossSectionSelfInteractionDifferential(Pi-theta,velocity(j)) &
           &                       )
      return
    end function integrandDecelerationFactor

  end subroutine kummer2018Tabulate
  
  double precision function kummer2018DecelerationFactor(self,x,velocityRelative)
    !!{RST
    Compute the deceleration factor, :math:`\chi_\mathrm{d}`, as defined by :cite:t:`kummer_effective_2018`.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSIDMVelocityDependent, darkMatterParticleSelfInteractingDarkMatterConstant
    use :: Error                , only : Error_Report
    implicit none
    class           (satelliteDecelerationSIDMKummer2018), intent(inout) :: self
    double precision                                     , intent(in   ) :: x
    double precision                                     , intent(in   ) :: velocityRelative
    double precision                                     , parameter     :: xCritical       =100.0d0
    double precision                                                     :: q                       , velocityCharacteristic

    if (x < xCritical) then
       if     (                                                                    &
            &   x                > self%xMaximum                                   &
            &  .or.                                                                &
            &   velocityRelative > self%velocityMaximum                            &
            &  .or.                                                                &
            &   velocityRelative < self%velocityMinimum                            &
            & )                                                                    &
            & call self%tabulate(                                                  &
            &                    x+1.0d0                                         , &
            &                    min(0.5d0*velocityRelative,self%velocityMinimum), &
            &                    max(2.0d0*velocityRelative,self%velocityMaximum)  &
            &                   )
       kummer2018DecelerationFactor=self%decelerationFactor_%interpolate(x,velocityRelative)
    else
       ! Large-x (x > xCritical) regime: use the analytic asymptotic form of the deceleration factor. The form is
       ! cross-section-model specific, so dispatch on the dark matter particle class.
       select type (darkMatterParticle_ => self%darkMatterParticle_)
       type is (darkMatterParticleSIDMVelocityDependent)
          ! Velocity-dependent cross section: read the characteristic velocity scale, v_c, from the particle.
          velocityCharacteristic      = darkMatterParticle_%velocityCharacteristic()
          q                           =+velocityRelative       &
               &                       /velocityCharacteristic
          kummer2018DecelerationFactor=+1.0d0                                &
               &                       -sqrt(2.0d0     )                     &
               &                       *    (1.0d0+q**2)                     &
               &                       *(                                    &
               &                         +2.0d0*sqrt(2.0d0)     / 3.0d0/x**2 &
               &                         -4.0d0*sqrt(2.0d0)     / 5.0d0/x**4 &
               &                         -8.0d0*sqrt(2.0d0)*q**2/15.0d0/x**4 &
               &                        )
       type is (darkMatterParticleSelfInteractingDarkMatterConstant)
          ! Constant cross section.
          kummer2018DecelerationFactor=+ 1.0d0            &
               &                       - 8.0d0/3.0d0/x**2 &
               &                       +16.0d0/5.0d0/x**4
       class default
          kummer2018DecelerationFactor=0.0d0
          call Error_Report('Kummer (2018) large-x deceleration factor is only implemented for the velocity-dependent and constant self-interacting dark matter particle classes'//{introspection:location})
       end select
    end if
    return
  end function kummer2018DecelerationFactor

