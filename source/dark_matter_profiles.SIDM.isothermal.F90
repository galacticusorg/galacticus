!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of dark matter halo profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et
  al. (2022), including the effects of a baryonic potential.
  !!}

  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  use            :: Numerical_ODE_Solvers  , only : odeSolver
  
  !![
  <darkMatterProfile name="darkMatterProfileSIDMIsothermal">
    <description>
      Dark matter halo profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et al. (2022). This
      model assumes that the dark matter within the interaction radius, $r_1$, has thermalized and can therefore be described by a
      constant velocity dispersion, $\sigma_0$. Under this assumption the spherical Jeans equation has a solution of the form:
      \begin{equation}
      \rho(r) = \rho_0 \exp\left[-\frac{\phi(r)}{\sigma_0^2}\right],
      \end{equation}
      where $\rho(r)$ is the density $\rho_0$ is the density at $r=0$, and the gravitational potential satisfies (Jiang et al. 2022):
      \begin{equation}
      \nabla^2 \phi(r) = 4 \pi \mathrm{G} \left[ \rho_0 \exp \left( - \frac{\phi(r)}{\sigma_0^2} \right) + \rho_\mathrm{b}(r) \right],
      \end{equation}
      where $\rho_\mathrm{b}(r)$ is the density of the baryonic component. This second-order differential equation is solved using the boundary conditions $\phi(r=0)=0$ and
      $\mathrm{d}\phi/\mathrm{d}r(r=0)=0$. The values of $\rho_0$ and $\sigma_0$ are then found by minimizing a function      
      \begin{equation}
      \delta^2(\rho_0,\sigma_0) = \left[ \frac{\rho(r_1)}{\rho^\prime(r_1)} - 1 \right]^2 + \left[ \frac{M(r_1)}{M^\prime(r_1)} - 1 \right]^2,
      \end{equation}
      where $M(r)$ is the mass contained within radius $r$, and primes indicate the profile prior to SIDM thermalization.
    </description>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileSIDM) :: darkMatterProfileSIDMIsothermal
     !!{
     A dark matter halo profile class implementing profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et al. (2022).
     !!}
     private
     integer         (kind=kind_int8)              :: uniqueIDPrevious
     double precision                              :: velocityDispersionCentral
     logical                                       :: solutionsTabulated
     class           (*             ), pointer     :: galacticStructure_        => null()
     type            (interpolator  ), allocatable :: densityProfile                     , massProfile
   contains
     !![
     <methods>
       <method method="computeSolution"   description="Compute a solution for the isothermal core of a SIDM halo."/>
     </methods>
     !!]
     final     ::                                      sidmIsothermalDestructor
     procedure :: autoHook                          => sidmIsothermalAutoHook
     procedure :: calculationReset                  => sidmIsothermalCalculationReset
     procedure :: density                           => sidmIsothermalDensity
     procedure :: densityLogSlope                   => sidmIsothermalDensityLogSlope
     procedure :: radiusEnclosingDensity            => sidmIsothermalRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => sidmIsothermalRadiusEnclosingMass
     procedure :: radialMoment                      => sidmIsothermalRadialMoment
     procedure :: enclosedMass                      => sidmIsothermalEnclosedMass
     procedure :: potential                         => sidmIsothermalPotential
     procedure :: circularVelocity                  => sidmIsothermalCircularVelocity
     procedure :: circularVelocityMaximum           => sidmIsothermalCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => sidmIsothermalRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => sidmIsothermalRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => sidmIsothermalRotationNormalization
     procedure :: energy                            => sidmIsothermalEnergy
     procedure :: kSpace                            => sidmIsothermalKSpace
     procedure :: freefallRadius                    => sidmIsothermalFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => sidmIsothermalFreefallRadiusIncreaseRate
     procedure :: computeSolution                   => sidmIsothermalComputeSolution
  end type darkMatterProfileSIDMIsothermal

  interface darkMatterProfileSIDMIsothermal
     !!{
     Constructors for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class.
     !!}
     module procedure sidmIsothermalConstructorParameters
     module procedure sidmIsothermalConstructorInternal
  end interface darkMatterProfileSIDMIsothermal

contains

  function sidmIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : galacticStructureConstruct_, galacticStructureDestruct_
    use :: Input_Parameters, only : inputParameter             , inputParameters
    implicit none
    type (darkMatterProfileSIDMIsothermal)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class(darkMatterParticleClass        ), pointer       :: darkMatterParticle_
    class(darkMatterProfileClass         ), pointer       :: darkMatterProfile_
    class(*                              ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfile"    name="darkMatterProfile_"    source="parameters"/>
    !!]
    call galacticStructureConstruct_(parameters,galacticStructure_)
    self=darkMatterProfileSIDMIsothermal(darkMatterProfile_,darkMatterHaloScale_,darkMatterParticle_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="darkMatterProfile_"  />
    !!]
    call galacticStructureDestruct_(self%galacticStructure_)
    return
  end function sidmIsothermalConstructorParameters

  function sidmIsothermalConstructorInternal(darkMatterProfile_,darkMatterHaloScale_,darkMatterParticle_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sidmIsothermal} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type (darkMatterProfileSIDMIsothermal)                        :: self
    class(darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterParticleClass        ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterProfileClass         ), intent(in   ), target :: darkMatterProfile_
    class(*                              ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*darkMatterProfile_, *darkMatterHaloScale_, *darkMatterParticle_, *galacticStructure_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('transfer function expects a self-interacting dark matter particle'//{introspection:location})
    end select
    self%uniqueIDPrevious    =-1_kind_int8
    self%genericLastUniqueID =-1_kind_int8
    self%uniqueIDPreviousSIDM=-1_kind_int8
    return
  end function sidmIsothermalConstructorInternal

  subroutine sidmIsothermalAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileSIDMIsothermal), intent(inout) :: self

    call calculationResetEvent%attach(self,sidmIsothermalCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine sidmIsothermalAutoHook

  subroutine sidmIsothermalDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class.
    !!}
    use :: Functions_Global, only : galacticStructureDestruct_
    use :: Events_Hooks    , only : calculationResetEvent
    implicit none
    type(darkMatterProfileSIDMIsothermal), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfile_"  />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    if (associated(self%galacticStructure_)) call galacticStructureDestruct_(self%galacticStructure_)
    if (calculationResetEvent%isAttached(self,sidmIsothermalCalculationReset)) call calculationResetEvent%detach(self,sidmIsothermalCalculationReset)
    return
  end subroutine sidmIsothermalDestructor

  subroutine sidmIsothermalCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    self%uniqueIDPrevious         =node%uniqueID()
    self%genericLastUniqueID      =node%uniqueID()
    self%uniqueIDPreviousSIDM     =node%uniqueID()
    self%radiusInteractivePrevious=-1.0d0
    self%velocityDispersionCentral=-1.0d0
    if (allocated(self%densityProfile                         )) deallocate(self%densityProfile                         )
    if (allocated(self%massProfile                            )) deallocate(self%massProfile                            )
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine sidmIsothermalCalculationReset

  subroutine sidmIsothermalComputeSolution(self,node)
    !!{
    Compute a solution for the isothermal core of an SIDM halo.
    !!}
    use :: Numerical_ODE_Solvers           , only : odeSolver
    use :: Numerical_Ranges                , only : Make_Range                     , rangeTypeLinear
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Multidimensional_Minimizer      , only : multiDMinimizer
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout)              :: self
    type            (treeNode                       ), intent(inout)              :: node
    integer         (c_size_t                       ), parameter                  :: propertyCount                =2
    integer                                          , parameter                  :: countTable                   =1000
    double precision                                 , parameter                  :: odeToleranceAbsolute         =1.0d-3, odeToleranceRelative     =1.0d-3
    double precision                                 , parameter                  :: fractionRadiusInitial        =1.0d-6
    double precision                                 , dimension(propertyCount+1) :: properties                          , propertyScales
    double precision                                 , dimension(countTable     ) :: radiusTable                         , densityTable                    , &
         &                                                                           massTable
    double precision                                 , dimension(propertyCount  ) :: locationMinimum
    type            (odeSolver                      )                             :: odeSolver_
    type            (multiDMinimizer                )                             :: minimizer_
    integer                                                                       :: i                                   , iteration
    logical                                                                       :: converged
    double precision                                                              :: densityCentral                      , velocityDispersionCentral       , &
         &                                                                           densityInteraction                  , massInteraction                 , &
         &                                                                           radiusInteraction                   , radius                          , &
         &                                                                           velocityDispersionInteraction       , mass                            , &
         &                                                                           density

    ! Find the interaction radius.
    radiusInteraction            =self%radiusInteraction                          (node                  )
    ! Properties of the original density profile at the interaction radius.
    densityInteraction           =self%darkMatterProfile_%density                 (node,radiusInteraction)
    massInteraction              =self%darkMatterProfile_%enclosedMass            (node,radiusInteraction)
    ! Find the velocity dispersion scale.
    velocityDispersionInteraction=sqrt(gravitationalConstantGalacticus*massInteraction/radiusInteraction)
    ! Set ODE solver  scales.
    propertyScales               =[velocityDispersionInteraction**2,velocityDispersionInteraction**2/radiusInteraction,massInteraction]
    ! Construct an ODE solver.
    odeSolver_                   =odeSolver      (propertyCount+1,sidmIsothermalODEs     ,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative,scale=propertyScales)
    ! Construct a minimizer.
    minimizer_                   =multiDMinimizer(propertyCount  ,sidmIsothermalFitMetric                                                                                                   )
    ! Seek the solution.
    call minimizer_%set(x=[0.0d0,1.0d0],stepSize=[1.0d0,1.0d0])
    iteration=0
    converged=.false.
    do while (.not.converged .and. iteration < 100)
       call minimizer_%iterate()
       iteration=iteration+1
       converged=minimizer_%testSize(toleranceAbsolute=1.0d-3)
    end do
    locationMinimum          =minimizer_%x()
    densityCentral           =exp(locationMinimum(1))*densityInteraction
    velocityDispersionCentral=    locationMinimum(2) *velocityDispersionInteraction
    ! Tabulate solutions for density and mass.
    radiusTable=Make_Range(rangeMinimum=0.0d0,rangeMaximum=radiusInteraction,rangeNumber=countTable,rangeType=rangeTypeLinear)
    densityTable(1)=densityCentral
    massTable   (1)=0.0d0
    do i=2,countTable
       radius    =fractionRadiusInitial*radiusInteraction
       properties=0.0d0
       call odeSolver_%solve(radius,radiusTable(i),properties)
       densityTable(i)=+densityCentral                    &
            &          *exp(                              &
            &               -properties(1)                &
            &               /velocityDispersionCentral**2 &
            &              )
       massTable   (i)=+     properties(3)
    end do
    allocate(self%densityProfile)
    allocate(self%   massProfile)
    self%           densityProfile=interpolator(radiusTable,             densityTable)
    self%              massProfile=interpolator(radiusTable,                massTable)
    self%velocityDispersionCentral=                         velocityDispersionCentral
    return
    
  contains
    
    double precision function sidmIsothermalFitMetric(propertiesCentral)
      !!{
      Evaluate the fit metric.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:)               :: propertiesCentral
      double precision               , dimension(propertyCount+1) :: properties
      double precision                                            :: radius
      
      ! Extract current parameters.
      densityCentral           =exp(propertiesCentral(1))*densityInteraction
      velocityDispersionCentral=    propertiesCentral(2) *velocityDispersionInteraction
      ! Solve the ODE to r₁.
      radius    =fractionRadiusInitial*radiusInteraction
      properties=0.0d0
      call odeSolver_%solve(radius,radiusInteraction,properties)
      ! Extract density and mass at r₁.
      density=+densityCentral                    &
           &  *exp(                              &
           &       -properties(1)                &
           &       /velocityDispersionCentral**2 &
           &      )
      mass     =+   properties(3)
      ! Evaluate the fit metric.
      sidmIsothermalFitMetric=+(density/densityInteraction-1.0d0)**2 &
           &                  +(   mass/   massInteraction-1.0d0)**2
      return
    end function sidmIsothermalFitMetric
  
    integer function sidmIsothermalODEs(radius,properties,propertiesRateOfChange)
      !!{
      Define the ODE system to solve for isothermal self-interacting dark matter cores.
      !!}
      use :: Functions_Global                , only : galacticStructureDensity_
      use :: Galactic_Structure_Options      , only : massTypeBaryonic
      use :: Interface_GSL                   , only : GSL_Success
      use :: Numerical_Constants_Math        , only : Pi
      use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
      implicit none
      double precision, intent(in   )               :: radius
      double precision, intent(in   ), dimension(:) :: properties
      double precision, intent(  out), dimension(:) :: propertiesRateOfChange
      double precision                              :: densityDarkMatter     , densityBaryons
      
      densityDarkMatter               =+densityCentral                    &
           &                           *exp(                              &
           &                                -max(properties(1),0.0d0)     &
           &                                /velocityDispersionCentral**2 &
           &                               )
      densityBaryons                  =+galacticStructureDensity_(self%galacticStructure_,node,position=[radius,0.0d0,0.0d0],massType=massTypeBaryonic)
      propertiesRateOfChange       (1)=+properties(2)
      propertiesRateOfChange       (2)=+4.0d0                             &
           &                           *Pi                                &
           &                           *gravitationalConstantGalacticus   &
           &                           *(                                 &
           &                             +densityDarkMatter               &
           &                             +densityBaryons                  &
           &                            )
      if (radius > 0.0d0)                                                 &
           & propertiesRateOfChange(2)=+propertiesRateOfChange(2)         &
           &                           -2.0d0                             &
           &                           *properties            (2)         &
           &                           /radius
      propertiesRateOfChange       (3)=+4.0d0                             &
           &                           *Pi                                &
           &                           *radius**2                         &
           &                           *densityDarkMatter
      sidmIsothermalODEs              = GSL_Success
      return
    end function sidmIsothermalODEs
    
  end subroutine sidmIsothermalComputeSolution
  
  double precision function sidmIsothermalDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    if (radius > self%radiusInteraction(node)) then
       sidmIsothermalDensity=self%darkMatterProfile_%density(node,radius)
    else
       if (node%uniqueID()                /= self%uniqueIDPrevious) call self%calculationReset(node)
       if (self%velocityDispersionCentral <= 0.0d0                ) call self%computeSolution(node)
       sidmIsothermalDensity=self%densityProfile%interpolate(radius)
    end if
    return
  end function sidmIsothermalDensity

  double precision function sidmIsothermalDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    if (radius > self%radiusInteraction(node)) then
       sidmIsothermalDensityLogSlope=self%darkMatterProfile_%densityLogSlope(node,radius)
    else
       if (node%uniqueID()                /= self%uniqueIDPrevious) call self%calculationReset(node)
       if (self%velocityDispersionCentral <= 0.0d0                ) call self%computeSolution(node)
       sidmIsothermalDensityLogSlope=self%densityProfile%derivative(radius)*radius/self%densityProfile%interpolate(radius)
    end if
    return
  end function sidmIsothermalDensityLogSlope

  double precision function sidmIsothermalEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    if (radius > self%radiusInteraction(node)) then
       sidmIsothermalEnclosedMass=self%darkMatterProfile_%enclosedMass(node,radius)
    else
       if (node%uniqueID()                /= self%uniqueIDPrevious) call self%calculationReset(node)
       if (self%velocityDispersionCentral <= 0.0d0                ) call self%computeSolution(node)
       sidmIsothermalEnclosedMass=self%massProfile%interpolate(radius)
    end if
    return
  end function sidmIsothermalEnclosedMass

  double precision function sidmIsothermalRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: density
    
    sidmIsothermalRadiusEnclosingDensity=self%radiusEnclosingDensityNumerical(node,density)
    return
  end function sidmIsothermalRadiusEnclosingDensity

  double precision function sidmIsothermalRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: mass

    sidmIsothermalRadiusEnclosingMass=self%radiusEnclosingMassNumerical(node,mass)
    return
  end function sidmIsothermalRadiusEnclosingMass

  double precision function sidmIsothermalRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum

    sidmIsothermalRadialMoment=self%radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    return
  end function sidmIsothermalRadialMoment

  double precision function sidmIsothermalPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout), target:: self
    type            (treeNode                       ), intent(inout), target   :: node
    double precision                                 , intent(in   )           :: radius
    integer                                          , intent(  out), optional :: status

    if (radius > self%radiusInteraction(node)) then
       sidmIsothermalPotential=self%darkMatterProfile_%potential(node,radius)
    else
       if (node%uniqueID()                /= self%uniqueIDPrevious) call self%calculationReset(node)
       if (self%velocityDispersionCentral <= 0.0d0                ) call self%computeSolution(node)
       sidmIsothermalPotential=self%darkMatterProfile_%potential(node,self%radiusInteraction(node))-self%velocityDispersionCentral**2*log(self%densityProfile%interpolate(radius)/self%densityProfile%interpolate(self%radiusInteraction(node)))
    end if
    return
  end function sidmIsothermalPotential

  double precision function sidmIsothermalCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    sidmIsothermalCircularVelocity=self%circularVelocityNumerical(node,radius)
    return
  end function sidmIsothermalCircularVelocity

  double precision function sidmIsothermalCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmIsothermalCircularVelocityMaximum=self%circularVelocityMaximumNumerical(node)
    return
  end function sidmIsothermalCircularVelocityMaximum

  double precision function sidmIsothermalRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    if (radius > self%radiusInteraction(node)) then
       sidmIsothermalRadialVelocityDispersion=self%darkMatterProfile_%radialVelocityDispersion(node,radius)
    else
       if (node%uniqueID()                /= self%uniqueIDPrevious) call self%calculationReset(node)
       if (self%velocityDispersionCentral <= 0.0d0                ) call self%computeSolution(node)
       sidmIsothermalRadialVelocityDispersion=self%velocityDispersionCentral
    end if
    return
  end function sidmIsothermalRadialVelocityDispersion

  double precision function sidmIsothermalRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout), target :: self
    type            (treeNode                       ), intent(inout)         :: node
    double precision                                 , intent(in   )         :: specificAngularMomentum

    sidmIsothermalRadiusFromSpecificAngularMomentum=self%radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    return
  end function sidmIsothermalRadiusFromSpecificAngularMomentum

  double precision function sidmIsothermalRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmIsothermalRotationNormalization=self%rotationNormalizationNumerical(node)
    return
  end function sidmIsothermalRotationNormalization

  double precision function sidmIsothermalEnergy(self,node)
    !!{
    Return the energy of a sidmIsothermal halo density profile.
    !!}
    implicit none
    class(darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmIsothermalEnergy=self%energyNumerical(node)
    return
  end function sidmIsothermalEnergy

  double precision function sidmIsothermalKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the sidmIsothermal density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout)         :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: waveNumber

    sidmIsothermalKSpace=self%kSpaceNumerical(node,waveNumber)
    return
  end function sidmIsothermalKSpace

  double precision function sidmIsothermalFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the sidmIsothermal density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: time

    sidmIsothermalFreefallRadius=self%freefallRadiusNumerical(node,time)
    return
  end function sidmIsothermalFreefallRadius

  double precision function sidmIsothermalFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the sidmIsothermal density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileSIDMIsothermal), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: time

    sidmIsothermalFreefallRadiusIncreaseRate=self%freefallRadiusIncreaseRateNumerical(node,time)
    return
  end function sidmIsothermalFreefallRadiusIncreaseRate
