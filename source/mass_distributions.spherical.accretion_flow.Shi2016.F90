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
  An accretion flow mass distribution using the framework of \cite{shi_outer_2016}.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Numerical_Interpolation, only : gsl_interp_cspline      , interpolator
  
  ! Note: Throughout this class the following acronyms are used:
  !        * HRTA - "half radius turnaround" - i.e. half of the turnaround radius for a given shell (more precisely, we use the
  !                 ratio of the virial radius to turnaround radius determined by spherical collapse models - this is previously
  !                 1/2 for an Einstein-de Sitter universe, but differs by a small amount for other cosmologies).

  ! Note: Throughout this class different three separate unit systems are used, identified by the variable name suffix:
  !        * Scaled    - these correspond to the scaled, self-similar variables used in Appendix A of Shi (2016) - i.e. column 3 of Table A1, "y" for radius, etc.
  !        * Original  - these correspond to the original variables used in Appendix A of Shi (2016) - i.e. column 2 of Table A1, "R" for radius, etc.
  !        * ScaleHRTA - these correspond to the HRTA unit system - that is, quantities are scaled to Rₕᵣₜₐ(a) and Mₕᵣₜₐ(a).
  !        * Physical  - these correspond to physical units (Mpc, M☉, km/s).
  
  !![
  <massDistribution name="massDistributionShi2016">
    <description>
      A mass distribution for accretion flows using the framework of \cite{shi_outer_2016}.
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionShi2016
     !!{
     A mass distribution for accretion flows using the framework of \cite{shi_outer_2016}.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_                               => null()
     type            (interpolator            ), allocatable               :: interpolatorDensityPhysical                                , interpolatorVelocityPhysical        , &
          &                                                                   interpolatorScaleFactorHalfRadiusTurnaroundScaled          , interpolatorRadiusScaled            , &
          &                                                                   interpolatorRadiusTurnaroundOriginal                       , interpolatorRadiusHRTAOriginal      , &
          &                                                                   interpolatorRadiusComovingInitialOriginal                  , interpolatorMassTurnaroundScaled    , &
          &                                                                   interpolatorMassHRTAOriginal                               , interpolatorMassMultiStreamScaleHRTA
     double precision                          , allocatable, dimension(:) :: radiusScaled                                               , overdensityScaled                   , &
          &                                                                   expansionFactorHRTAScaled                                  , radiusGrowthRateScaled              , &
          &                                                                   radiusComovingInitialOriginal                              , massEnclosedInitialOriginal         , &
          &                                                                   timeTurnaroundScaled                                       , timeHRTAScaled                      , &
          &                                                                   radiusTurnaroundScaled                                     , radiusOrderedOriginal               , &
          &                                                                   massShellOrderedOriginal                                   , densityOrderedOriginal              , &
          &                                                                   massEnclosedOrderedOriginal
     double precision                                                      :: radiusMaximumPhysical                                      , scaleFactorVelocity                 , &
          &                                                                   expansionFactorScaled                                      , radiusMultistreamMinimumScaledHRTA  , &
          &                                                                   radiusMultistreamMaximumScaledHRTA                         , cosmologicalConstantScaled          , &
          &                                                                   radiusTurnaroundNowOriginal                                , radiusHRTANowOriginal               , &
          &                                                                   timeNowScaled                                              , radiusSplashbackTurnaround          , &
          &                                                                   ratioRadiusSplashbackHRTA                                  , radiusSplashbackOriginal            , &
          &                                                                   radiusSplashbackScaled                                     , radiusMinimumPhysical               , &
          &                                                                   radiusVirial                                               , mass                                , &
          &                                                                   massAccretionRate                                          , redshift                            , &
          &                                                                   time                                                       , ratioRadiusTurnaroundVirial
   contains
     !![
     <methods>
       <method description="Solve for the structure of the accretion flow." method="solve" />
     </methods>
     !!]
     procedure :: density               => shi2016Density
     procedure :: densityGradientRadial => shi2016DensityGradientRadial 
     procedure :: solve                 => shi2016Solve
  end type massDistributionShi2016

  interface massDistributionShi2016
     !!{
     Constructors for the \refClass{massDistributionShi2016} mass distribution class.
     !!}
     module procedure massDistributionShi2016ConstructorParameters
     module procedure massDistributionShi2016ConstructorInternal
  end interface massDistributionShi2016

  ! Sub-module scope variables used in root finding.
  double precision :: timeTurnaroundScaled__, radiusTurnaroundScaled__
  !$omp threadprivate(timeTurnaroundScaled__, radiusTurnaroundScaled__)
  
  ! Sub-module scope variables used in ODE solving.
  double precision                                   :: radiusComovingInitialOriginal        , massEnclosedInitialOriginal
  class           (massDistributionShi2016), pointer :: self_
  logical                                            :: noShellCrossing              =.false.
  !$omp threadprivate(self_,massEnclosedInitialOriginal,radiusComovingInitialOriginal,noShellCrossing)

contains

  function massDistributionShi2016ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionShi2016} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionShi2016 )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    type            (varying_string          )                :: componentType
    type            (varying_string          )                :: massType
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                                          :: scaleFactorVelocity , mass                       , &
         &                                                       massAccretionRate   , ratioRadiusTurnaroundVirial, &
         &                                                       radiusVirial        , redshift
    
    !![
    <inputParameter>
      <name>mass</name>
      <description>The mass of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massAccretionRate</name>
      <description>The mass accretion rate of the halo.</description>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
      <name>radiusVirial</name>
      <description>The virial radius of the halo.</description>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
      <name>ratioRadiusTurnaroundVirial</name>
      <description>The ratio of the turnaround to virial radii of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=massDistributionShi2016(mass,massAccretionRate,radiusVirial,ratioRadiusTurnaroundVirial,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),scaleFactorVelocity,cosmologyFunctions_,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function massDistributionShi2016ConstructorParameters

  function massDistributionShi2016ConstructorInternal(mass,massAccretionRate,radiusVirial,ratioRadiusTurnaroundVirial,time,scaleFactorVelocity,cosmologyFunctions_,componentType,massType) result(self)
    !!{
    Internal constructor for ``shi2016'' mass distribution class.
    !!}
    implicit none
    type            (massDistributionShi2016     )                          :: self
    class           (cosmologyFunctionsClass     ), intent(in   ), target   :: cosmologyFunctions_
    double precision                              , intent(in   )           :: mass                , massAccretionRate          , &
         &                                                                     time                , scaleFactorVelocity        , &
         &                                                                     radiusVirial        , ratioRadiusTurnaroundVirial
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="mass, massAccretionRate, radiusVirial, ratioRadiusTurnaroundVirial, time, scaleFactorVelocity, *cosmologyFunctions_,  componentType, massType"/>
    !!]

    self%dimensionless=.false.
    self%redshift     =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    call self%solve()
    return
  end function massDistributionShi2016ConstructorInternal

  subroutine shi2016Destructor(self)
    !!{
    Destructor for the \refClass{massDistributionShi2016} accretion flow mass distribution class.
    !!}
    implicit none
    type(massDistributionShi2016), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine shi2016Destructor

  double precision function shi2016Density(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a accretion flow modeled on the 2-halo correlation function.
    !!}
    implicit none
    class(massDistributionShi2016), intent(inout) :: self
    class(coordinate             ), intent(in   ) :: coordinates
    
    if      (coordinates%rSpherical() > self%radiusMaximumPhysical) then
       ! Beyond the maximum radius for the flow just return the mean matter density.
       density=self%cosmologyFunctions_       %matterDensityEpochal(self       %time        )
    else if (coordinates%rSpherical() < self%radiusMinimumPhysical) then
       density=0.0d0
       call Error_Report('radius is less than minimum tabulated for accretion flow'//{introspection:location})
    else
       density=self%interpolatorDensityPhysical%interpolate        (coordinates%rSpherical())
    end if
    return
  end function shi2016Density

  double precision function shi2016DensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a accretion flow modeled on the 2-halo correlation function.
    !!}
    implicit none
    class  (massDistributionShi2016), intent(inout), target   :: self
    class  (coordinate             ), intent(in   )           :: coordinates
    logical                         , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    if      (coordinates%rSpherical() > self%radiusMaximumPhysical) then
       ! Beyond the maximum radius for the flow just return the mean matter density.
       densityGradientRadial=0.0d0
    else if (coordinates%rSpherical() < self%radiusMinimumPhysical) then
       densityGradientRadial=0.0d0
       call Error_Report('radius is less than minimum tabulated for accretion flow'//{introspection:location})
    else
       densityGradientRadial=self%interpolatorDensityPhysical%derivative(coordinates%rSpherical())
    end if
    return
  end function shi2016DensityGradientRadial
  
  subroutine shi2016Solve(self)
    !!{
    Solve the accretion flow.
    !!}
    use :: Array_Utilities                 , only : Array_Reverse
    use :: Display                         , only : displayCounter           , displayCounterClear          , displayIndent                , displayUnindent, &
          &                                         verbosityLevelWorking
    use :: Elliptic_Integrals              , only : Elliptic_Integral_K      , Elliptic_Integral_Pi
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : var_str
    use :: Numerical_Comparison            , only : Values_Differ
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Ranges                , only : Make_Range               , rangeTypeLogarithmic
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Sorting                         , only : sortIndex
    use :: String_Handling                 , only : operator(//)
    implicit none
    class           (massDistributionShi2016), intent(inout), target       :: self
    double precision                         , allocatable  , dimension(:) :: radiusOriginal                                , densityOriginal                      , &
         &                                                                    velocityOriginal                              , radiusSingleStreamAnalyticPhysical   , &
         &                                                                    velocitySingleStreamAnalyticPhysical          , densitySingleStreamAnalyticPhysical  , &
         &                                                                    expansionFactorHRTAScaled                     , radiusPhysical                       , &
         &                                                                    velocityPhysical                              , densityPhysical
    integer         (c_size_t               ), allocatable  , dimension(:) :: order
    type            (interpolator           ), allocatable                 :: interpolatorMassMultiStreamNewScaleHRTA
    integer                                  , parameter                   :: countRadii                             =1000
    integer                                  , parameter                   :: countCompare                           =  10
    integer                                  , parameter                   :: iterationMaximum                       =  30
    double precision                         , parameter                   :: expansionFactorRelativeInitial         =1.0d-6
    double precision                         , parameter                   :: radiusMultiStreamFractionalSmall       =1.0d-6
    double precision                         , parameter                   :: multistreamToleranceRelative           =5.0d-2
    double precision                                                       :: expansionFactorOriginal                       , timeInitialScaled                    , &
         &                                                                    radiusInitialScaled                           , radiusGrowthRateInitialScaled        , &
         &                                                                    overdensityMinimumScaled                      , overdensityMaximumScaled             , &
         &                                                                    bigA                                          , bigB                                 , &
         &                                                                    bigC                                          , expansionFactorScaledInitial         , &
         &                                                                    massVirialOriginal                            , radiusVirialOriginal                 , &
         &                                                                    growthIndex                                   , radiusScaleHRTA                      , &
         &                                                                    radiusMultistreamMinimumNewScaledHRTA         , radiusMultistreamMaximumNewScaledHRTA, &
         &                                                                    radiusCompareMinimumScaledHRTA                , radiusCompareMaximumScaledHRTA       , &
         &                                                                    massMultiStream                               , massMultiStreamNew                   , &
         &                                                                    changeRelative                                , changeRelativeMaximum                , &
         &                                                                    radiusSplashbackPhysical                      , h
    integer                                                                :: i                                             , j                                    , &
         &                                                                    iTurnaround                                   , iFirstZero                           , &
         &                                                                    iVirial                                       , iteration
    type            (rootFinder             )                              :: finder
    logical                                                                :: firstZeroFound                                , multistreamConverged
    character       (len=12                 )                              :: label
    
    ! The growth index specifies the profile of the initial mass perturbation, δMᵢ/Mᵢ ∝ M^{-1/s}, or, equivalently, the growth rate,
    ! M(t) ∝ aˢ.
    growthIndex=+self                    %massAccretionRate                                                      &
         &      /self%cosmologyFunctions_%expansionRate    (self%cosmologyFunctions_%expansionFactor(self%time)) &
         &      /self                    %mass
    ! Compute the scale-free solution.
    allocate(self%overdensityScaled            (countRadii))
    allocate(self%radiusScaled                 (countRadii))
    allocate(self%radiusGrowthRateScaled       (countRadii))
    allocate(self%expansionFactorHRTAScaled    (countRadii))
    allocate(self%radiusComovingInitialOriginal(countRadii))
    allocate(self%massEnclosedInitialOriginal  (countRadii))
    allocate(self%radiusTurnaroundScaled       (countRadii))
    allocate(self%timeTurnaroundScaled         (countRadii))
    allocate(self%radiusOrderedOriginal        (countRadii))
    allocate(self%massShellOrderedOriginal     (countRadii))
    allocate(self%massEnclosedOrderedOriginal  (countRadii))
    allocate(self%densityOrderedOriginal       (countRadii))
    allocate(self%timeHRTAScaled               (countRadii))
    allocate(     order                        (countRadii))
    allocate(     expansionFactorHRTAScaled    (countRadii))
    allocate(     radiusOriginal               (countRadii))
    ! Create a module-scope pointer to self for use in ODE solver functions.
    self_ => self
    ! Compute the scaled cosmological constant parameter ("w" in the notation of Shi 2016, Table A1).
    self%cosmologicalConstantScaled=+1.0d0/self%cosmologyFunctions_ %OmegaMatterEpochal(self%time)                &
         &                          -1.0d0
    ! Compute expansion factor, and scaled expansion factor ("y" in the notation of Shi 2016, Table A1).
    expansionFactorOriginal        =+self%cosmologyFunctions_       %expansionFactor   (self%time)
    self%expansionFactorScaled     =+self%cosmologicalConstantScaled                              **(1.0d0/3.0d0) &
         &                          *expansionFactorOriginal
    ! Find the scaled time at the current epoch.
    self%timeNowScaled             =+sqrt(1.0d0-self %cosmologyFunctions_%OmegaMatterEpochal(self%time              )) & ! Equation A5 from Shi (2016).
         &                          *           self %cosmologyFunctions_%expansionRate     (expansionFactorOriginal)  &
         &                          *           self%time
    ! Choose an initial epoch and (scaled) radius for the ODE. This is chosen to be an expansion factor much smaller than the
    ! present day such that the perturbations will be small.
    expansionFactorScaledInitial   =+     expansionFactorRelativeInitial &
         &                          *self%expansionFactorScaled
    radiusInitialScaled            =+     expansionFactorScaledInitial
    ! Choose a value of the scaled overdensity, "β" in the notation of Shi (2016), identifying a mass shell. We avoid β=1 because such a shell never collapses.
    overdensityMinimumScaled       =+1.001d0
    ! Choose a maximum overdensity. We use β=10 here as it's more than sufficient to allow most of the accretion stream to be captured,
    overdensityMaximumScaled       =+1.000d1
    ! Build array of overdensities.
    self%overdensityScaled         =Make_Range(overdensityMaximumScaled,overdensityMinimumScaled,countRadii,rangeTypeLogarithmic)
    ! Build a root finder which will be used for finding the time at half the turnaround radius.
    finder=rootFinder(                                            &
         &            rootFunction     =halfRadiusTurnAroundRoot, &
         &            toleranceAbsolute=0.0d+0                  , &
         &            toleranceRelative=1.0d-3                    &
         &           )
    ! Find turnaround radius and mass as a function of time, along with epoch at which half the turnaround radius is
    ! reached. This is all in scale-free units.
    do i=1,countRadii
       ! Find the epoch of turnaround.
       self%radiusTurnaroundScaled(i)=+2.0d0**(2.0d0/3.0d0)                                            & ! Equation A9 from Shi (2016).
            &                         *                  sqrt(      self%overdensityScaled(i)       )  &
            &                         *sin((1.0d0/3.0d0)*asin(1.0d0/self%overdensityScaled(i)**1.5d0))
       bigA                          =+1.0d0                                                           & ! Text after equation of A12 from Shi (2016).
            &                         /self%radiusTurnaroundScaled(i)**3
       bigB                          =+2.0d0                                                           & ! Text after equation of A12 from Shi (2016).
            &                         /    (    +3.0d0+sqrt(1.0d0+4.0d0*bigA)      )
       bigC                          =+                sqrt(1.0d0+4.0d0*bigA)                          & ! Text after equation of A12 from Shi (2016).
            &                         /    (bigA-0.5d0+sqrt(1.0d0+4.0d0*bigA)/2.0d0)
       self%timeTurnaroundScaled  (i)=+    (    +1.0d0+sqrt(1.0d0+4.0d0*bigA)      )                   & ! Equation A12 of Shi (2016).
            &                         /sqrt(bigA-0.5d0+sqrt(1.0d0+4.0d0*bigA)/2.0d0)                   &
            &                         *(                                                               &
            &                           +Elliptic_Integral_Pi(bigC,bigB)                               &
            &                           -Elliptic_Integral_K (bigC     )                               &
            &                          )
       ! Find the epoch corresponding to reaching half of the turnaround radius.
       radiusTurnaroundScaled__=self%radiusTurnaroundScaled(i)
       timeTurnaroundScaled__  =self%timeTurnaroundScaled  (i)
       call finder%rangeExpand(                                                             &
            &                  rangeExpandUpward            =1.1d0                        , &
            &                  rangeExpandDownward          =0.5d0                        , &
            &                  rangeExpandType              =rangeExpandMultiplicative    , &
            &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                  rangeDownwardLimit           =timeTurnaroundScaled__         &
            &                 )
       self%timeHRTAScaled           (i)=finder%find(rootGuess=timeTurnaroundScaled__)
       self%expansionFactorHRTAScaled(i)=expansionFactorFromTimeScaled(self%timeHRTAScaled(i))
    end do
    ! Find enclosed masses. Note that this is found using the self-similarity assumption that the mass of a shell should scale
    ! as the expansion factor at "half-radius-turnaround" to the power of the growth index, i.e. M ∝ uₕᵣₜₐˢ. In regions
    ! where the cosmological constant is negligible (high overdensities, which collapse at early times), this will also give
    ! the expected scaling with overdensity, M ∝ β⁻ˢ. But, at lower overdensities, which collapse when the
    ! cosmological constant is non-negligible, the scaling with overdensity will change. What we do here seems to be consistent
    ! with what Shi (2016) assumes, and has the nice feature that it ensures the mass-epoch relation is the simple scale-free
    ! expectation at all times.
    expansionFactorHRTAScaled         =expansionFactorFromTimeScaled(self%timeHRTAScaled)
    self%massEnclosedInitialOriginal  =(expansionFactorHRTAScaled/self%expansionFactorScaled      )**growthIndex
    self%radiusComovingInitialOriginal=                           self%massEnclosedInitialOriginal **(1.0d0/3.0d0)
    ! Build a variety of interpolators for different radii and masses as functions of scaled time.
    allocate(self%interpolatorRadiusTurnaroundOriginal             )
    allocate(self%interpolatorRadiusHRTAOriginal                   )
    allocate(self%interpolatorRadiusComovingInitialOriginal        )
    allocate(self%interpolatorMassTurnaroundScaled                 )
    allocate(self%interpolatorMassHRTAOriginal                     )
    allocate(self%interpolatorScaleFactorHalfRadiusTurnaroundScaled)
    self%interpolatorRadiusTurnaroundOriginal             =interpolator(              self%timeTurnaroundScaled ,              self%radiusTurnaroundScaled       *self%radiusComovingInitialOriginal/self%cosmologicalConstantScaled**(1.0d0/3.0d0)                             )
    self%interpolatorRadiusComovingInitialOriginal        =interpolator(              self%timeTurnaroundScaled ,              self%radiusComovingInitialOriginal                                                                                                               )
    self%interpolatorMassTurnaroundScaled                 =interpolator(              self%timeTurnaroundScaled ,              self%massEnclosedInitialOriginal                                                                                                                 )
    self%interpolatorRadiusHRTAOriginal                   =interpolator(              self%timeHRTAScaled       ,              self%radiusTurnaroundScaled       *self%radiusComovingInitialOriginal/self%cosmologicalConstantScaled**(1.0d0/3.0d0)/self%ratioRadiusTurnaroundVirial )
    self%interpolatorMassHRTAOriginal                     =interpolator(              self%timeHRTAScaled       ,              self%massEnclosedInitialOriginal                                                                                                                 )
    self%interpolatorScaleFactorHalfRadiusTurnaroundScaled=interpolator(Array_Reverse(self%overdensityScaled   ),Array_Reverse(self%expansionFactorHRTAScaled                                                                                                                  ))
    ! Compute present epoch turnaround and half-turnaround radii.
    self%radiusTurnaroundNowOriginal=self%interpolatorRadiusTurnaroundOriginal%interpolate(self%timeNowScaled)
    self%radiusHRTANowOriginal      =self%interpolatorRadiusHRTAOriginal      %interpolate(self%timeNowScaled)
    ! Make an estimate of the splashback radius (in units of the turnaround radius, using eqn. 2 of Shi 2016 which is for an
    ! Einstein-de Sitter universe). An approximate value is acceptable here as this is used only on the first iteration. It's
    ! not clear if this is precisely what Shi (2016) chose, but it should not matter.
    if (growthIndex <= 1.5d0) then
       self%radiusSplashbackTurnaround=+1.0d0/3.0d0**(2.0d0/3.0d0+2.0d0*growthIndex/9.0d0)
    else
       self%radiusSplashbackTurnaround=+1.0d0/(1.0d0+4.0d0*(4.0d0*growthIndex/9.0d0+1.0d0/3.0d0)/sqrt(Pi))
    end if
    ! Find the splashback radius in units of the present day half-turnaround radius.       
    self%ratioRadiusSplashbackHRTA=self%radiusSplashbackTurnaround*self%radiusTurnaroundNowOriginal/self%radiusHRTANowOriginal
    ! For the first iteration, adopt a mass profile solution in the multi-stream region that has the form f(x)=x (Shi 2016,
    ! section 2.1).
    allocate(self%interpolatorMassMultiStreamScaleHRTA)
    self%radiusMultistreamMinimumScaledHRTA=radiusMultiStreamFractionalSmall*self%ratioRadiusSplashbackHRTA
    self%radiusMultistreamMaximumScaledHRTA=                                 self%ratioRadiusSplashbackHRTA
    self%interpolatorMassMultiStreamScaleHRTA=interpolator([self%radiusMultistreamMinimumScaledHRTA,self%radiusMultistreamMaximumScaledHRTA],[self%radiusMultistreamMinimumScaledHRTA,self%radiusMultistreamMaximumScaledHRTA])
    ! Begin iterating to find a solution.
    iteration           =0
    multistreamConverged=.false.
    do while (iteration < iterationMaximum .and. .not.multistreamConverged)
       iteration=iteration+1
       call displayIndent(var_str('multistream mass profile iteration ')//iteration,verbosity=verbosityLevelWorking)
       ! Iterate over all overdensity shells solving for their radial position and velocity at the present epoch.
       do i=1,countRadii
          call displayCounter(int(100.0d0*dble(i-1)/dble(countRadii)),isNew=i==1,verbosity=verbosityLevelWorking)
          ! Solve the dynamical ODEs to get the scaled radius at the final time.
          radiusComovingInitialOriginal=self%radiusComovingInitialOriginal(i)
          massEnclosedInitialOriginal  =self%massEnclosedInitialOriginal  (i)
          radiusGrowthRateInitialScaled=+sqrt(                                                                                      & ! Equation A8 from Shi (2016).
               &                              +1.0d0/     radiusInitialScaled                                                       &
               &                              +           radiusInitialScaled   **2                                                 &
               &                              -3.0d0*self%overdensityScaled  (i)                                                    &
               &                              /2.0d0                            **(2.0d0/3.0d0)                                     &
               &                             )
          timeInitialScaled            =+2.0d0                                                                                      & ! Equation A6 from Shi (2016).
               &                        /3.0d0                                                                                      &
               &                        *asinh(                                               expansionFactorScaledInitial **1.5d0)
          call radiusScaledSolver(timeInitialScaled,self%timeNowScaled,radiusInitialScaled,radiusGrowthRateInitialScaled,self%radiusScaled(i),self%radiusGrowthRateScaled(i))
       end do
       call displayCounterClear(verbosity=verbosityLevelWorking)
       ! Build an interpolator for the scaled radius as a function of overdensity.
       if (allocated(self%interpolatorRadiusScaled)) deallocate(self%interpolatorRadiusScaled)
       allocate(self%interpolatorRadiusScaled)
       self%interpolatorRadiusScaled=interpolator(Array_Reverse(self%overdensityScaled),Array_Reverse(self%radiusScaled))
       ! Construct the mass and density profile by ordering the shells in radius. Also find the splashback radius.
       radiusOriginal               =abs(self%radiusScaled)*self%radiusComovingInitialOriginal/self%cosmologicalConstantScaled**(1.0d0/3.0d0)
       order                        =sortIndex(radiusOriginal)
       self%radiusSplashBackOriginal=-huge(0.0d0)
       firstZeroFound               =.false.
       do i=1,countRadii
          ! Order radii.
          self%radiusOrderedOriginal(i)=radiusOriginal(order(i))
          ! Compute mass and density in this shell.
          if (order(i) == 1) then
             self%massShellOrderedOriginal(i)=self%massEnclosedInitialOriginal(order(i))
             self%densityOrderedOriginal  (i)=self%massShellOrderedOriginal(i)*3.0d0/4.0d0/Pi/ self%radiusOrderedOriginal(i)**3
          else
             self%massShellOrderedOriginal(i)=self%massEnclosedInitialOriginal(order(i))-self%massEnclosedInitialOriginal(order(i)-1)
             self%densityOrderedOriginal  (i)=self%massShellOrderedOriginal(i)*3.0d0/4.0d0/Pi/(self%radiusOrderedOriginal(i)**3-self%radiusOrderedOriginal(i-1)**3)
          end if
          ! Find the mass enclosed by this shell.
          self%massEnclosedOrderedOriginal(i)=sum(self%massShellOrderedOriginal(1:i))
       end do
       ! Find the splashback radius.
       do i=countRadii,1,-1
          firstZeroFound=firstZeroFound.or.self%radiusScaled(i) <= 0.0d0
          if (firstZeroFound .and. radiusOriginal(i) > self%radiusSplashBackOriginal) then
             self%radiusSplashbackOriginal=         radiusOriginal(i)
             self%radiusSplashbackScaled  =abs(self%radiusScaled  (i))
          end if
       end do
       deallocate(radiusOriginal)
       ! Find the splashback radius in units of the turnaround radius, and the half-turnaround radius.       
       self%radiusSplashbackTurnaround=self%radiusSplashBackOriginal/self%radiusTurnaroundNowOriginal
       self%ratioRadiusSplashbackHRTA =self%radiusSplashBackOriginal/self%radiusHRTANowOriginal
       ! Build a new interpolator for the mass profile in the multistream region, and compare it to the previous one.
       allocate(interpolatorMassMultiStreamNewScaleHRTA)
       interpolatorMassMultiStreamNewScaleHRTA=interpolator(self%radiusOrderedOriginal/self%radiusHRTANowOriginal,self%massEnclosedOrderedOriginal)
       radiusMultistreamMinimumNewScaledHRTA  =self%radiusOrderedOriginal(         1)/self%radiusHRTANowOriginal
       radiusMultistreamMaximumNewScaledHRTA  =self%radiusOrderedOriginal(countRadii)/self%radiusHRTANowOriginal
       radiusCompareMinimumScaledHRTA         =max(radiusMultistreamMinimumNewScaledHRTA,self%radiusMultistreamMinimumScaledHRTA)
       radiusCompareMaximumScaledHRTA         =min(radiusMultistreamMaximumNewScaledHRTA,self%radiusMultistreamMaximumScaledHRTA)
       changeRelativeMaximum                  =-huge(0.0d0)
       do j=1,countCompare
          radiusScaleHRTA      =+                                radiusCompareMinimumScaledHRTA  &
               &                +(radiusCompareMaximumScaledHRTA-radiusCompareMinimumScaledHRTA) &
               &                *dble(           j)                                              &
               &                /dble(countCompare)
          massMultiStream      =self%interpolatorMassMultiStreamScaleHRTA   %interpolate(radiusScaleHRTA)
          massMultiStreamNew   =     interpolatorMassMultiStreamNewScaleHRTA%interpolate(radiusScaleHRTA)
          changeRelative       =+abs(massMultiStream-massMultiStreamNew) &
               &                /   (massMultiStream+massMultiStreamNew) &
               &                /0.5d0
          changeRelativeMaximum=max(changeRelative,changeRelativeMaximum)
       end do
       multistreamConverged=changeRelativeMaximum <= multistreamToleranceRelative
       deallocate(interpolatorMassMultiStreamNewScaleHRTA)
       ! Replace the interpolator for the mass profile in the multistream region with the updated one.
       deallocate(self%interpolatorMassMultiStreamScaleHRTA)
       allocate  (self%interpolatorMassMultiStreamScaleHRTA)
       self%interpolatorMassMultiStreamScaleHRTA=interpolator(self%radiusOrderedOriginal/self%radiusHRTANowOriginal,self%massEnclosedOrderedOriginal)
       self%radiusMultistreamMinimumScaledHRTA         =self%radiusOrderedOriginal(         1)/self%radiusHRTANowOriginal
       self%radiusMultistreamMaximumScaledHRTA         =self%radiusOrderedOriginal(countRadii)/self%radiusHRTANowOriginal
       write (label,'(e8.2)') changeRelativeMaximum
       call displayUnindent(var_str('done [fractional change = ')//trim(adjustl(label))//']',verbosity=verbosityLevelWorking)
    end do
    if (.not.multistreamConverged) call Error_Report('failed to reach convergence in the multistream region'//{introspection:location})
    ! Compute properties along the stream.
    allocate(radiusOriginal  (countRadii))
    allocate(densityOriginal (countRadii))
    allocate(velocityOriginal(countRadii))
    !! Find the radii and velocities in the stream.
    radiusOriginal  =self%radiusScaled          *self%radiusComovingInitialOriginal/self%cosmologicalConstantScaled**(1.0d0/3.0d0)
    velocityOriginal=self%radiusGrowthRateScaled*self%radiusComovingInitialOriginal/self%cosmologicalConstantScaled**(1.0d0/3.0d0)
    !! Compute densities.
    do i=1,countRadii
       if (i == 1) then
          densityOriginal(i)= self%massEnclosedInitialOriginal(i)                                       *3.0d0/4.0d0/Pi/abs(radiusOriginal(i))**3 
       else
          densityOriginal(i)=(self%massEnclosedInitialOriginal(i)-self%massEnclosedInitialOriginal(i-1))*3.0d0/4.0d0/Pi/abs(radiusOriginal(i) **3-radiusOriginal(i-1)**3)
       end if
    end do
    ! Find the shells at their turnaround radius and "half" of their turnaround radius (this defines the virial mass/radius), and
    ! also the shell which is about to make its first passage through zero radius.
    iVirial    =-1
    iTurnaround=-1
    iFirstZero =-1
    do i=countRadii,1,-1
       if (iVirial     < 0 .and. self%radiusGrowthRateScaled(i) <  0.0d0 .and. self%radiusScaled(i) <= self%radiusTurnaroundScaled(i)/self%ratioRadiusTurnaroundVirial) &
            & iVirial    =i
       if (iTurnaround < 0 .and. self%radiusGrowthRateScaled(i) <= 0.0d0                                                                                              ) &
            & iTurnaround=i
       if (iFirstZero  < 0 .and. self%radiusGrowthRateScaled(i) <  0.0d0 .and. self%radiusScaled(i) <  0.0d0                                                          ) &
            & iFirstZero =i+1
    end do
    ! Interpolate to get a more precise virial radius.
    h      =+(1.0d0                       /self%ratioRadiusTurnaroundVirial           -self%radiusScaled(iVirial)/self%radiusTurnaroundScaled(iVirial)) &
         &  /(self%radiusScaled(iVirial+1)/self%radiusTurnaroundScaled(iVirial+1)-self%radiusScaled(iVirial)/self%radiusTurnaroundScaled(iVirial))
    radiusVirialOriginal=+(                                                                                      &
         &                 +self%radiusScaled(iVirial  )*self%radiusComovingInitialOriginal(iVirial  )*(1.0d0-h) &
         &                 +self%radiusScaled(iVirial+1)*self%radiusComovingInitialOriginal(iVirial+1)*       h  &
         &                )                                                                                      &
         &               /self%cosmologicalConstantScaled**(1.0d0/3.0d0)
    massVirialOriginal  =+self%interpolatorMassHRTAOriginal        %interpolate(                     self%timeNowScaled        ) &
         &               *self%interpolatorMassMultiStreamScaleHRTA%interpolate(radiusVirialOriginal/self%radiusHRTANowOriginal)
    ! Scale radii, velocities, and densities to the virial radius/mass.
    allocate(radiusPhysical  (countRadii))
    allocate(velocityPhysical(countRadii))
    allocate(densityPhysical (countRadii))
    radiusPhysical  =+radiusOriginal                                                    &
         &           * self_%radiusVirial/radiusVirialOriginal
    velocityPhysical=+velocityOriginal                                                  &
         &           * self_%radiusVirial/radiusVirialOriginal                          &
         &           * self_%timeNowScaled                    /self%time                &
         &           *megaParsec                                                        &
         &           /gigaYear                                                          &
         &           /kilo
    densityPhysical =+densityOriginal                                                   &
         &           * self%mass                              /massVirialOriginal       &
         &           /(self_%radiusVirial/radiusVirialOriginal)**3
    ! Build interpolators into the infall stream.
    if (allocated(self%interpolatorDensityPhysical )) deallocate(self%interpolatorDensityPhysical )
    if (allocated(self%interpolatorVelocityPhysical)) deallocate(self%interpolatorVelocityPhysical)
    allocate(self%interpolatorDensityPhysical )
    allocate(self%interpolatorVelocityPhysical)
    self%interpolatorDensityPhysical =interpolator(radiusPhysical(iFirstZero:countRadii),densityPhysical (iFirstZero:countRadii),interpolationType=gsl_interp_cspline)
    self%interpolatorVelocityPhysical=interpolator(radiusPhysical(iFirstZero:countRadii),velocityPhysical(iFirstZero:countRadii),interpolationType=gsl_interp_cspline)
    self%radiusMinimumPhysical       =radiusPhysical(iFirstZero)
    self%radiusMaximumPhysical       =radiusPhysical(countRadii)
    ! Compute the analytic solution in the single stream regime.
    allocate(radiusSingleStreamAnalyticPhysical  (countRadii-iVirial-1))
    allocate(densitySingleStreamAnalyticPhysical (countRadii-iVirial-1))
    allocate(velocitySingleStreamAnalyticPhysical(countRadii-iVirial-1))
    do i=iVirial+1,countRadii-1
       radiusSingleStreamAnalyticPhysical  (i-iVirial)=+(                                                                                                                                  &
            &                                            +3.0d0                                                                                                                            &
            &                                            *self%mass                                                                                                                        &
            &                                            /4.0d0                                                                                                                            &
            &                                            /Pi                                                                                                                               &
            &                                            /self%cosmologyFunctions_%matterDensityEpochal(self%time)                                                                         &
            &                                           )**(1.0d0/3.0d0)                                                                                                                   &
            &                                          *(                                                                                                                                  &
            &                                            +self%expansionFactorHRTAScaled(i)                                                                                                &
            &                                            /self%expansionFactorScaled                                                                                                       &
            &                                           )**(growthIndex/3.0d0)                                                                                                             &
            &                                          *self%radiusScaled(i)                                                                                                               &
            &                                          /self%expansionFactorScaled
       densitySingleStreamAnalyticPhysical (i-iVirial)=+self%cosmologyFunctions_%matterDensityEpochal(self%time)                                                                           &
            &                                          *(                                                                                                                                  &
            &                                            +self%expansionFactorScaled                                                                                                       &
            &                                            /self%radiusScaled         (i)                                                                                                    &
            &                                           )**3                                                                                                                               &
            &                                          /(                                                                                                                                  &
            &                                            +1.0d0                                                                                                                            &
            &                                            +3.0d0                                                                                                                            &
            &                                            /growthIndex                                                                                                                      &
            &                                            /(self%interpolatorScaleFactorHalfRadiusTurnaroundScaled%derivative(self%overdensityScaled(i))/self%expansionFactorHRTAScaled(i)) &
            &                                            *(self%interpolatorRadiusScaled                         %derivative(self%overdensityScaled(i))/self%radiusScaled             (i)) &
            &                                          )
       velocitySingleStreamAnalyticPhysical(i-iVirial)=+self%radiusGrowthRateScaled            (i        )                                                                                 &
            &                                          *     radiusSingleStreamAnalyticPhysical(i-iVirial)                                                                                 &
            &                                          /self%radiusScaled                      (i        )                                                                                 &
            &                                          *sqrt(1.0d0-self%cosmologyFunctions_%OmegaMatterEpochal    (self%time))                                                             &
            &                                          *           self%cosmologyFunctions_%HubbleParameterEpochal(self%time)
    end do
    ! Check that the numerical solution matches the analytic solution in the single stream region.
    radiusSplashbackPhysical=+self%radiusSplashBackOriginal &
         &                   *self%radiusVirial             &
         &                   /     radiusVirialOriginal
    do i=iVirial+1,countRadii-1
       if     (                                                                                                                             &
            &   radiusPhysical(i) > radiusSplashBackPhysical                                                                                &
            &  .and.                                                                                                                        &
            &   Values_Differ(                                                                                                              &
            &                        radiusSingleStreamAnalyticPhysical(i-iVirial)/radiusSingleStreamAnalyticPhysical(iTurnaround-iVirial), &
            &                        radiusPhysical                    (i        )/radiusPhysical                    (iTurnaround        ), &
            &                 relTol=1.0d-6                                                                                                 &
            &                )                                                                                                              &
            & )                                                                                                                             &
            & call Error_Report('numerical and analytic solutions disagree in single stream region'//{introspection:location})
    end do
    return
  end subroutine shi2016Solve

  double precision function halfRadiusTurnAroundRoot(timeFinalScaled)
    !!{
    Root function used in finding the epoch at which a shell reaches a radius equal to half of its turnaround radius, $y^*$.
    !!}
    implicit none
    double precision, intent(in   ) :: timeFinalScaled
    double precision                :: radiusScaled   , radiusGrowthRateScaled

    ! Integrate the dynamical equations governing the evolution of the shell radius starting from the turnaround time, I⋆, at
    ! which the radius is y⋆, and the rate of change of radius is zero (by definition).
    noShellCrossing=.true.
    call radiusScaledSolver(timeTurnaroundScaled__,timeFinalScaled,radiusTurnaroundScaled__,0.0d0,radiusScaled,radiusGrowthRateScaled)
    noShellCrossing=.false.
    halfRadiusTurnAroundRoot=+      radiusScaled                &
         &                   -      radiusTurnaroundScaled__    &
         &                   /self_%ratioRadiusTurnaroundVirial
    return
  end function halfRadiusTurnAroundRoot
  
  subroutine radiusScaledSolver(timeInitialScaled,timeFinalScaled,radiusInitialScaled,radiusGrowthRateInitialScaled,radiusScaled,radiusGrowthRateScaled)
    !!{
    Compute the scaled radius (and its growth rate) as a function of the initial state and final time.
    !!}
    use :: Interface_GSL        , only : GSL_Success
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    double precision           , intent(in   ) :: timeInitialScaled         , timeFinalScaled                      , &
         &                                        radiusInitialScaled       , radiusGrowthRateInitialScaled
    double precision           , intent(  out) :: radiusScaled              , radiusGrowthRateScaled
    double precision           , dimension(2)  :: odeVariables
    double precision           , parameter     :: odeToleranceAbsolute=0.0d0, odeToleranceRelative         =1.0d-9
    type            (odeSolver)                :: solver
    double precision                           :: timeScaled

    solver      =odeSolver(2_c_size_t,dynamicalODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)
    odeVariables=[radiusInitialScaled,radiusGrowthRateInitialScaled]
    timeScaled  = timeInitialScaled
    call solver%solve(timeScaled,timeFinalScaled,odeVariables)
    radiusScaled          =odeVariables(1)
    radiusGrowthRateScaled=odeVariables(2)
    return
  end subroutine radiusScaledSolver

  integer function dynamicalODES(timeScaled,odeVariables,odeVariablesGrowthRate)
    !!{
    The dynamical equation describing the motion of a shell of matter in scaled variables \citep[][eqn.~A7]{shi_outer_2016}.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision              , intent(in   ) :: timeScaled
    double precision, dimension(:), intent(in   ) :: odeVariables
    double precision, dimension(:), intent(  out) :: odeVariablesGrowthRate
    double precision                              :: radiusScaled          , radiusGrowthRateScaled, &
         &                                           radiusHRTA            , radiusSplashback      , &
         &                                           radius                , massHRTA              , &
         &                                           massEnclosedRatio
    !$GLC attributes unused :: timeScaled
    
    ! Extract ODE variables into named variables for clarity.
    radiusScaled          =odeVariables(1)
    radiusGrowthRateScaled=odeVariables(2)
    ! Determine the ratio of the mass enclosed within the current radius to the mass enclosed within the shell at the initial
    ! time.
    if (noShellCrossing .or. timeScaled < self_%timeHRTAScaled(1)) then
       ! If shell crossing is being ignored, or if the current time is less than the earliest time for which we have the
       ! half-turnaround radius tabulated, simply assume a mass ratio of unity.
       massEnclosedRatio=1.0d0
    else
       ! Our solution for the mass in the multistream regime is the self-similar solution expressed in units of the radius of the
       ! shell at its half-turnaround radius, and the mass within that shell (assuming no shell crossing). Compute that radius and
       ! mass at the present epoch.
       radiusHRTA=self_%interpolatorRadiusHRTAOriginal%interpolate(timeScaled)
       massHRTA  =self_%interpolatorMassHRTAOriginal  %interpolate(timeScaled)
       ! Find the current splashback radius by scaling the ratio of splashback to half-turnaround radius to the current epoch.
       radiusSplashback=+self_%ratioRadiusSplashbackHRTA &
            &           *      radiusHRTA
       ! Find the radius of the current shell in unscaled units.
       radius=abs(radiusScaled)*radiusComovingInitialOriginal/self_%cosmologicalConstantScaled**(1.0d0/3.0d0)
       ! Determine where in the stream our shell is.
       if (radius < radiusSplashback) then
          ! The shell is within the splashback radius.
          if (radius > self_%radiusMultistreamMinimumScaledHRTA*radiusHRTA) then
             ! Shell is outside the minimum radius that we have tabulated for the multistream region. Therefore, simply
             ! interpolate to get the mass ratio.
             massEnclosedRatio=+massHRTA                                                                  &
                  &            /massEnclosedInitialOriginal                                               &
                  &            *self_%interpolatorMassMultiStreamScaleHRTA%interpolate(radius/radiusHRTA)
          else
             ! Shell is inside the minimum radius that we have tabulated for the multistream region. In this region we assume that
             ! the mass enclosed by the innermost tabulated point is distributed with uniform density. The mass enclosed therefore
             ! grows as the cube of radius.
             massEnclosedRatio=+massHRTA                                                                                         &
                  &            /massEnclosedInitialOriginal                                                                      &
                  &            *self_%interpolatorMassMultiStreamScaleHRTA%interpolate(self_%radiusMultistreamMinimumScaledHRTA) &
                  &            *(radius/radiusHRTA/self_%radiusMultistreamMinimumScaledHRTA)**3
          end if
       else
          ! Mass is outside the splashback radius, so no shell crossing has occurred. The mass ratio is therefore unity.
          massEnclosedRatio=1.0d0
       end if
    end if
    ! Set ODE rates of change.
    !! Radius rate of change is just the velocity.
    odeVariablesGrowthRate(1)=+radiusGrowthRateScaled
    !! Velocity rate of change is given by equation (A7) of Shi (2016).
    if (radiusScaled == 0.0d0) then
       odeVariablesGrowthRate(2)=+0.0d0
    else
       odeVariablesGrowthRate(2)=-0.5d0*massEnclosedRatio*sign(1.0d0,radiusScaled)/radiusScaled**2 &
            &                    +                                                 radiusScaled
    end if
    dynamicalODES=GSL_Success
    return
  end function dynamicalODES

  elemental double precision function expansionFactorFromTimeScaled(timeScaled) result(expansionFactor)
    !!{
    Compute the scaled expansion factor from the scaled time using equation~(A6) of \cite{shi_outer_2016}.
    !!}
    implicit none
    double precision, intent(in   ) :: timeScaled

    expansionFactor=sinh(1.5d0*timeScaled)**(2.0d0/3.0d0)
    return
  end function expansionFactorFromTimeScaled

