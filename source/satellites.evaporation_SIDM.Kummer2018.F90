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
  Implementation of a satellite evaporation due to dark matter self-interactions using the model of
  \cite{kummer_effective_2018}.
  !!}

  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Numerical_Interpolation , only : interpolator

  !![
  <satelliteEvaporationSIDM name="satelliteEvaporationSIDMKummer2018">
   <description>A satellite evaporation due to dark matter self-interactions using the model of \cite{kummer_effective_2018}.</description>
  </satelliteEvaporationSIDM>
  !!]
  type, extends(satelliteEvaporationSIDMClass) :: satelliteEvaporationSIDMKummer2018
     !!{
     Implementation of a satellite evaporation due to dark matter self-interactions using the model of
     \cite{kummer_effective_2018}.
     !!}
     private
     class           (darkMatterParticleClass  ), pointer     :: darkMatterParticle_         => null()
     class           (darkMatterProfileDMOClass), pointer     :: darkMatterProfileDMO_       => null()
     type            (interpolator             ), allocatable :: evaporationFactor
     double precision                                         :: rateScatteringNormalization          , xMaximum
   contains
     !![
     <methods>
       <method description="Tabulate the evaporation factor." method="tabulate" />
     </methods>
     !!]
     final     ::                 kummer2018Destructor
     procedure :: massLossRate => kummer2018MassLossRate
     procedure :: tabulate     => kummer2018Tabulate
  end type satelliteEvaporationSIDMKummer2018

  interface satelliteEvaporationSIDMKummer2018
     !!{
     Constructors for the \refClass{satelliteEvaporationSIDMKummer2018} satellite evaporation due to dark matter self-interactions class.
     !!}
     module procedure kummer2018ConstructorParameters
     module procedure kummer2018ConstructorInternal
  end interface satelliteEvaporationSIDMKummer2018

contains

  function kummer2018ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteEvaporationSIDMKummer2018} satellite evaporation due to dark matter self-interactions class
    which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteEvaporationSIDMKummer2018)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterParticleClass           ), pointer       :: darkMatterParticle_
    class(darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
  
    !![
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteEvaporationSIDMKummer2018(darkMatterParticle_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function kummer2018ConstructorParameters

  function kummer2018ConstructorInternal(darkMatterParticle_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteEvaporationSIDMKummer2018} satellite evaporation due to dark matter self-interactions
    class.
    !!}
    use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Numerical_Constants_Prefixes    , only : centi                                     , milli   , kilo
    use :: Numerical_Constants_Astronomical, only : megaParsec                                , gigaYear, massSolar
    implicit none
    type (satelliteEvaporationSIDMKummer2018)                        :: self
    class(darkMatterParticleClass           ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterParticle_, *darkMatterProfileDMO_"/>
    !!]

    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Compute the normalization of the scattering rate in units such that when multiplied by a velocity in km s⁻¹, and a
       ! density in units of M☉ Mpc⁻³, we get a rate in units of Gyr⁻¹.
       self%rateScatteringNormalization=+darkMatterParticle_%crossSectionSelfInteraction()*centi    **2/milli         & ! Convert cross-section from cm² g⁻¹ to m² kg⁻¹.
            &                           *                                                  kilo                       & ! Convert velocity from km s⁻¹ to m s⁻¹.
            &                           *                                                  massSolar   /megaParsec**3 & ! Convert density from M☉ Mpc⁻³ to kg m⁻³.
            &                           *                                                  gigaYear                     ! Convert rate from s⁻¹ to Gyr⁻¹.
    class default
       ! No scattering.
       self%rateScatteringNormalization=+0.0d0
    end select
    ! Initialize the maximum tabulated x to an unphysical value. This will force tabulation on the first attempt to evaluate the
    ! evaporation factor.
    self%xMaximum=-1.0d0
    return

  end function kummer2018ConstructorInternal

  subroutine kummer2018Destructor(self)
    !!{
    Destructor for the \refClass{satelliteEvaporationSIDMKummer2018} satellite evaporation due to dark matter self-interactions class.
    !!}
    implicit none
    type(satelliteEvaporationSIDMKummer2018), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine kummer2018Destructor

 double precision function kummer2018MassLossRate(self,node)
    !!{
    Return a evaporation for satellites due to dark matter self-interactions using the formulation of \cite{kummer_effective_2018}.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , coordinateCartesian        , assignment(=)
    use :: Galactic_Structure_Options      , only : coordinateSystemCartesian     , radiusLarge
    use :: Galacticus_Nodes                , only : nodeComponentSatellite        , nodeComponentBasic
    use :: Mass_Distributions              , only : massDistributionClass         , kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (satelliteEvaporationSIDMKummer2018), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    class           (nodeComponentSatellite            ), pointer       :: satellite
    class           (nodeComponentBasic                ), pointer       :: basic
    type            (treeNode                          ), pointer       :: nodeHost
    class           (massDistributionClass             ), pointer       :: massDistribution_     , massDistributionHost_
    class           (kinematicsDistributionClass       ), pointer       :: kinematics_           , kinematicsHost_
    double precision                                    , dimension(3)  :: position              , velocity
    double precision                                                    :: radiusOrbital         , speedOrbital               , &
         &                                                                 densityHost           , rateScattering             , &
         &                                                                 massBoundary          , radiusBoundary             , &
         &                                                                 velocityEscape        , speedHalfMass              , &
         &                                                                 velocityDispersionHost, velocityDispersionSatellite, &
         &                                                                 x                     , radiusHalfMass             , &
         &                                                                 velocityDispersion    , potentialEscape
    type            (coordinateSpherical                )                :: coordinates          , coordinatesHost            , &
         &                                                                  coordinatesBoundary   , coordinatesHalfMass
    type            (coordinateCartesian                )                :: coordinatesCartesian
    
    ! Set zero mass loss rate by default.
    kummer2018MassLossRate=0.0d0
    ! If the scattering cross section is zero, we can return immediately.
    if (self%rateScatteringNormalization == 0.0d0) return
    ! Evaluate satellite and host properties.
    nodeHost              =>  node                 %mergesWith      (                    )
    satellite             =>  node                 %satellite       (                    )
    massDistributionHost_ =>  nodeHost             %massDistribution(                    )
    position              =   satellite            %position        (                    )
    velocity              =   satellite            %velocity        (                    )
    radiusOrbital         =   Vector_Magnitude                      (            position)
    speedOrbital          =   Vector_Magnitude                      (            velocity)
    coordinatesCartesian  =                                                      position
    densityHost           =   massDistributionHost_%density         (coordinatesCartesian)
    !![
    <objectDestructor name="massDistributionHost_"/>
    !!]
    ! Find the escape velocity from the half-mass radius of the subhalo. This is equal to the potential difference between the
    ! half-mass radius and outer boundary of the subhalo, plus the potential difference from the outer boundary to infinity (for
    ! which we can treat the subhalo as a point mass).
    !
    ! Note that Kummer et al. (2018) suggest computing an effective escape velocity which is the escape velocity averaged over the
    ! mass profile of the subhalo, M⁻¹∫ρ(r) v(r) d³r. Evaluating that integral would be computationally expensive. Instead, we use
    ! the escape velocity at the half-mass radius of the subhalo. For a scale-free Hernquist profile, M⁻¹∫ρ(r) v(r) d³r =
    ! 8√2/15=0.754, while v(r½=1+√2)=√(2/(2+√2))=0.765 - so the two choices result in very similar effective escape velocities.
    basic        =>    node     %basic    ()
    massBoundary = min(                       &
         &             satellite%boundMass(), &
         &             basic    %     mass()  &
         &            )
    if (massBoundary > 0.0d0) then
       massDistribution_ => node             %massDistribution   (                       )
       radiusBoundary    =  massDistribution_%radiusEnclosingMass(mass=      massBoundary)
       radiusHalfMass    =  massDistribution_%radiusEnclosingMass(mass=0.5d0*massBoundary)
       if (radiusBoundary < 0.5d0*radiusLarge) then
          coordinatesBoundary=[radiusBoundary,0.0d0,0.0d0]
          coordinatesHalfMass=[radiusHalfMass,0.0d0,0.0d0]
          potentialEscape    =+massDistribution_%potentialDifference(coordinatesBoundary,coordinatesHalfMass) &
               &              +gravitationalConstant_internal                                                 &
               &              *massBoundary                                                                   &
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
       <objectDestructor name="massDistribution_"/>
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
       ! Get the speed of a host particle at the half-mass radius of the subhalo - this is the sum of the kinetic energy or host
       ! particles in the rest-frame of the subhalo, plus the energy they gain by falling in to the half-mass radius of the
       ! subhalo. We include the correction factor of the velocity dispersion as suggested in equation (A4) of Kummer et
       ! al. (2018).
       speedHalfMass=+sqrt(                       &
            &              +speedOrbital      **2 &
            &              +velocityDispersion**2 &
            &              +velocityEscape    **2 &
            &             )
       x            =+      velocityEscape        &
            &        /      speedOrbital
       ! Evaporation occurs for x<1.
       if (x < 1.0d0) then
          call self%tabulate(1.0d0)
          ! Evaluate the scattering rate and mass loss rate.
          rateScattering               =  +     speedOrbital                     &
               &                          *     densityHost                      &
               &                          *self%rateScatteringNormalization
          kummer2018MassLossRate       =  -     massBoundary                     &
               &                          *     rateScattering                   &
               &                          *self%evaporationFactor%interpolate(x)
       end if
    end if
    return
  end function kummer2018MassLossRate

  subroutine kummer2018Tabulate(self,xMaximum)
    !!{
    Tabulate the evaporation factor, $\chi_\mathrm{d}$.
    !!}
    use :: Dark_Matter_Particles   , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range                                 , rangeTypeLinear
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (satelliteEvaporationSIDMKummer2018         ), intent(inout)               :: self
    double precision                                             , intent(in   )               :: xMaximum
    class           (darkMatterParticleSelfInteractingDarkMatter), pointer                     :: darkMatterParticleSIDM_
    double precision                                             , allocatable  , dimension(:) :: x                         , evaporationFactor
    integer                                                      , parameter                   :: countPerUnit           =10
    type            (integrator                                 )                              :: integrator_
    double precision                                                                           :: thetaCritical
    integer                                                                                    :: i                         , countX

    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Tabulate the χ_d factor.
       self%xMaximum=xMaximum
       countX       =int(xMaximum*dble(countPerUnit))+1
       allocate(x                (countX))
       allocate(evaporationFactor(countX))
       x                       =  Make_Range(0.0d0,xMaximum,countX,rangeTypeLinear)
       darkMatterParticleSIDM_ => darkMatterParticle_
       integrator_             =  integrator(integrandEvaporationFactor,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
       do i=1,countX
          thetaCritical        =+acos(                  & ! Equation (2) from Kummer et al. (2018).
               &                      +(+x(i)**2-1.0d0) &
               &                      /(+x(i)**2+1.0d0) &
               &                     )
          evaporationFactor(i)=+integrator_        %integrate                  (Pi-thetaCritical,thetaCritical) &
               &               /darkMatterParticle_%crossSectionSelfInteraction(                              )
       end do
       if (allocated(self%evaporationFactor)) deallocate(self%evaporationFactor)
       allocate(self%evaporationFactor)
       self%evaporationFactor=interpolator(x,evaporationFactor)
    end select
    return

  contains
    
    double precision function integrandEvaporationFactor(theta)
      !!{
      The integrand used to evaluate $\chi_\mathrm{e}$ from \cite[][eqn.~3]{kummer_effective_2018}.
      !!}
      implicit none
      double precision, intent(in   ) :: theta

      integrandEvaporationFactor=+darkMatterParticleSIDM_%crossSectionSelfInteractionDifferential(theta)
      return
    end function integrandEvaporationFactor

  end subroutine kummer2018Tabulate
  
