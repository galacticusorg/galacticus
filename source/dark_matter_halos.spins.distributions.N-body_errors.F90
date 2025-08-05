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
  An implementation of the dark matter halo spin distribution which modifies another spin distribution to account for the
  effects of particle noise errors in spins measured in N-body simulations.
  !!}

  use :: Cosmology_Functions              , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales          , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales       , only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use :: Halo_Mass_Functions              , only : haloMassFunctionClass
  use :: Root_Finder                      , only : rootFinder
  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <haloSpinDistribution name="haloSpinDistributionNbodyErrors">
   <description>
    A halo spin distribution which modifies another spin distribution to account for the effects of particle noise errors
    in spins measured in N-body simulations.
   </description>
  </haloSpinDistribution>
  !!]
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionNbodyErrors
     !!{
     A dark matter halo spin distribution class which modifies another spin distribution to account for the effects of particle
     noise errors in spins measured in N-body simulations.
     !!}
     private
     class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_           => null()
     class           (haloSpinDistributionClass        ), pointer                     :: distributionIntrinsic         => null()
     class           (nbodyHaloMassErrorClass          ), pointer                     :: nbodyHaloMassError_           => null()
     class           (haloMassFunctionClass            ), pointer                     :: haloMassFunction_             => null()
     class           (darkMatterHaloScaleClass         ), pointer                     :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileScaleRadiusClass), pointer                     :: darkMatterProfileScaleRadius_ => null()
     type            (rootFinder                       )                              :: finder
     double precision                                                                 :: massParticle                           , time       , &
          &                                                                              logNormalRange                         , redshift
     logical                                                                          :: fixedPoint
     integer                                                                          :: particleCountMinimum
     integer                                                                          :: spinCount                              , massCount
     double precision                                                                 :: spinMinimum                            , spinMaximum, &
          &                                                                              massMinimum                            , massMaximum, &
          &                                                                              massDelta                              , spinDelta  , &
          &                                                                              energyEstimateParticleCountMaximum     , spinFixed
     double precision                                   , allocatable, dimension(:  ) :: massWeight
     double precision                                   , allocatable, dimension(:,:) :: distributionTable
   contains
     !![
     <methods>
       <method description="Return the spin distribution function averaged over all halos above the given {\normalfont \ttfamily massLimit}." method="distributionAveraged" />
       <method description="Return the spin distribution function at a fixed point in intrinsic mass and spin." method="distributionFixedPoint" />
       <method description="Tabulate the spin distribution as a fuction of spin and halo mass. Ensure that the table spans the {\normalfont \ttfamily massRequired} and {\normalfont \ttfamily spinRequireed} if provided." method="tabulate" />
     </methods>
     !!]
     final     ::                           nbodyErrorsDestructor
     procedure :: sample                 => nbodyErrorsSample
     procedure :: distribution           => nbodyErrorsDistribution
     procedure :: distributionFixedPoint => nbodyErrorsDistributionFixedPoint
     procedure :: distributionAveraged   => nbodyErrorsDistributionAveraged
     procedure :: tabulate               => nbodyErrorsTabulate
  end type haloSpinDistributionNbodyErrors

  interface haloSpinDistributionNbodyErrors
     !!{
     Constructors for the \refClass{haloSpinDistributionNbodyErrors} dark matter halo spin
     distribution class.
     !!}
     module procedure nbodyErrorsConstructorParameters
     module procedure nbodyErrorsConstructorInternal
  end interface haloSpinDistributionNbodyErrors

  ! Tabulation parameters.
  integer         , parameter :: spinPointsPerDecade=33    , massPointsPerDecade=5
  double precision, parameter :: spinMaximum        =0.5d+0, massMaximum        =1.0d15
  double precision, parameter :: spinMinimum        =3.0d-4, massMinimum        =1.0d11

  ! Module-scope variable used in root finding.
  double precision :: nonCentralChiSquareChi
  !$omp threadprivate(nonCentralChiSquareChi)

contains

  function nbodyErrorsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloSpinDistributionNbodyErrors} dark matter halo spin
    distribution class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloSpinDistributionNbodyErrors  )                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (haloSpinDistributionClass        ), pointer       :: distributionIntrinsic
    class           (nbodyHaloMassErrorClass          ), pointer       :: nbodyHaloMassError_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (haloMassFunctionClass            ), pointer       :: haloMassFunction_
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass), pointer       :: darkMatterProfileScaleRadius_
    double precision                                                   :: massParticle                 , redshift                          , &
         &                                                                time                         , energyEstimateParticleCountMaximum, &
         &                                                                logNormalRange
    integer                                                            :: particleCountMinimum

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>Mass of particle in the simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>particleCountMinimum</name>
      <source>parameters</source>
      <description>Minimum number of particles per halo in the N-body simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>energyEstimateParticleCountMaximum</name>
      <source>parameters</source>
      <description>Maximum number of particles used in estimating the potential energy of halos. Set to a very large number if no such maximum was used.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <description>Redshift at which the spin distribution should be evaluated.</description>
    </inputParameter>
    <inputParameter>
      <name>logNormalRange</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <defaultSource>A large range which will include (almost) the entirety of the distribution.</defaultSource>
      <description>The multiplicative range of the log-normal distribution used to model the distribution of the mass and energy terms in the spin parameter. Specifically, the lognormal distribution is truncated outside the range $(\lambda_\mathrm{m}/R,\lambda_\mathrm{m} R$, where $\lambda_\mathrm{m}$ is the measured spin, and $R=${\normalfont \ttfamily [logNormalRange]}</description>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution"         name="distributionIntrinsic"         source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    ! Find the time corresponding to the given redshift.
    time=  cosmologyFunctions_ %cosmicTime                 (          &
         &  cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                   redshift &
         &                                                  )         &
         &                                                 )
    ! Construct the object.
    self=nbodyErrorsConstructorInternal(                                    &
         &                              distributionIntrinsic             , &
         &                              massParticle                      , &
         &                              particleCountMinimum              , &
         &                              energyEstimateParticleCountMaximum, &
         &                              logNormalRange                    , &
         &                              time                              , &
         &                              nbodyHaloMassError_               , &
         &                              cosmologyFunctions_               , &
         &                              haloMassFunction_                 , &
         &                              darkMatterHaloScale_              , &
         &                              darkMatterProfileScaleRadius_       &
         &                             )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="distributionIntrinsic"        />
    <objectDestructor name="nbodyHaloMassError_"          />
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function nbodyErrorsConstructorParameters

  function nbodyErrorsConstructorInternal(distributionIntrinsic,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,logNormalRange,time,nbodyHaloMassError_,cosmologyFunctions_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileScaleRadius_) result(self)
    !!{
    Internal constructor for the \refClass{haloSpinDistributionNbodyErrors} dark matter halo spin distribution class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (haloSpinDistributionNbodyErrors  )                        :: self
    class           (haloSpinDistributionClass        ), intent(in   ), target :: distributionIntrinsic
    class           (nbodyHaloMassErrorClass          ), intent(in   ), target :: nbodyHaloMassError_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (haloMassFunctionClass            ), intent(in   ), target :: haloMassFunction_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass), intent(in   ), target :: darkMatterProfileScaleRadius_
    double precision                                   , intent(in   )         :: massParticle                      , time          , &
         &                                                                        energyEstimateParticleCountMaximum, logNormalRange
    integer                                            , intent(in   )         :: particleCountMinimum
    !![
    <constructorAssign variables="*distributionIntrinsic, massParticle, particleCountMinimum, energyEstimateParticleCountMaximum, logNormalRange, time, *nbodyHaloMassError_, *cosmologyFunctions_, *haloMassFunction_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_"/>
    !!]

    ! Store the redshift.
    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Set default ranges of spin and mass for tabulation.
    self%fixedPoint =.false.
    self%spinFixed  =-huge(0.0d0)
    self%spinMinimum=spinMinimum
    self%spinMaximum=spinMaximum
    self%massMinimum=massMinimum
    self%massMaximum=massMaximum
    ! Build a root finder.
    self%finder=rootFinder(                                                                      &
         &                 rootFunction                 =nbodyErrorsNonCentralChiSquareModeRoot, &
         &                 toleranceRelative            =1.0d-3                                , &
         &                 rangeExpandUpward            =2.0d0                                 , &
         &                 rangeExpandDownward          =0.5d0                                 , &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive         , &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative         , &
         &                 rangeExpandType              =rangeExpandMultiplicative               &
         &                )
     return
  end function nbodyErrorsConstructorInternal

  subroutine nbodyErrorsTabulate(self,massRequired,spinRequired,massFixed,spinFixed,spinFixedMeasuredMinimum,spinFixedMeasuredMaximum,tree)
    !!{
    Tabulate the halo spin distribution.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Calculations_Resets   , only : Calculations_Reset
    use :: Error                 , only : Error_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic                     , nodeComponentDarkMatterProfile, nodeComponentSpin, treeNode, &
         &                                mergerTree
    use :: Numerical_Integration , only : integrator
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout)                    :: self
    double precision                                 , intent(in   )         , optional :: massRequired                        , spinRequired               , &
         &                                                                                 massFixed                           , spinFixed                  , &
         &                                                                                 spinFixedMeasuredMinimum            , spinFixedMeasuredMaximum
    type            (mergerTree                     ), intent(in   ), target , optional :: tree
    type            (treeNode                       )               , pointer           :: node
    class           (nodeComponentBasic             )               , pointer           :: nodeBasic
    class           (nodeComponentSpin              )               , pointer           :: nodeSpin
    class           (nodeComponentDarkMatterProfile )               , pointer           :: nodeDarkMatterProfile
    double precision                                 , parameter                        :: massIntegrationRange         =10.0d0
    type            (integrator                     )                                   :: integratorMass                      , integratorSpin             , &
         &                                                                                 integratorMassSpin                  , integratorProduct
    logical                                                                             :: retabulate                          , fixedPoint
    integer                                                                             :: iSpin                               , iMass
    double precision                                                                    :: spinMeasured                        , massMeasured               , &
         &                                                                                 densityRatioInternalToSurface       , massError                  , &
         &                                                                                 radiusHalo                          , densityOuterRadius         , &
         &                                                                                 massMinimum                         , massMaximum                , &
         &                                                                                 logNormalWidth                      , logNormalLogMean           , &
         &                                                                                 errorSpinDependent                  , errorSpinIndependent       , &
         &                                                                                 errorSpinIndependent1D              , nonCentrality              , &
         &                                                                                 particleNumber                      , massEnergyEstimateMaximum  , &
         &                                                                                 energyEstimateErrorCorrection       , massAverage                , &
         &                                                                                 massSpinAverage

    ! Validate optional parameter combinations.
    if     (          (present(massRequired).or. present(spinRequired))                                                               &
         &  .and.                                                                                                                     &
         &       .not.(present(massRequired).and.present(spinRequired))                                                               &
         & ) call Error_Report('both massRequired and spinRequired must be provided if either is provided'//{introspection:location})
    if     (          (present(massFixed   ).or. present(spinFixed   ))                                                               &
         &  .and.                                                                                                                     &
         &       .not.(present(massFixed   ).and.present(spinFixed   ))                                                               &
         & ) call Error_Report('both massFixed and spinFixed  must be provided if either is provided'//{introspection:location})
    if     (                                                                                                                          &
         &             present(massRequired).and.present(massFixed   )                                                                &
         & ) call Error_Report('either required or fixed points must be specified - not both'//{introspection:location})
    fixedPoint=present(massFixed)
    ! Determine if the tabulation needs to be rebuilt.
    if (allocated(self%distributionTable)) then
       if (fixedPoint) then
          retabulate=.not.self%fixedPoint
          if (self%massMinimum /= massFixed) then
             retabulate      =.true.
             self%massMinimum=massFixed
          end if
          if (self%spinFixed /= spinFixed) then
             retabulate      =.true.
             self%spinFixed  =spinFixed
          end if
          if (retabulate) then
             self%spinMinimum=min(self%spinMinimum,0.5d0*spinFixed)
             self%spinMaximum=max(self%spinMaximum,2.0d0*spinFixed)
          end if
          if (present(spinFixedMeasuredMinimum)) then
             if (spinFixedMeasuredMinimum < self%spinMinimum) then
                retabulate      =.true.
                self%spinMinimum=spinFixedMeasuredMinimum
             end if
             if (spinFixedMeasuredMaximum > self%spinMaximum) then
                retabulate      =.true.
                self%spinMaximum=spinFixedMeasuredMaximum
             end if
          end if
       else
          retabulate=self%fixedPoint
          if (present(massRequired)) then
             retabulate=massRequired < self%massMinimum .or. massRequired > self%massMaximum
             if (retabulate) then
                self%massMinimum=min(self%massMinimum,0.5d0*massRequired)
                self%massMaximum=max(self%massMaximum,2.0d0*massRequired)
             end if
          end if
          if (present(spinRequired)) then
             retabulate=spinRequired < self%spinMinimum .or. spinRequired > self%spinMaximum
             if (retabulate) then
                self%spinMinimum=min(self%spinMinimum,0.5d0*spinRequired)
                self%spinMaximum=max(self%spinMaximum,2.0d0*spinRequired)
             end if
          end if
       end if
       if (retabulate)  deallocate(self%distributionTable)
    else
       retabulate=.true.
    end if
    if (.not.retabulate) return
    self%fixedPoint=fixedPoint
    ! Allocate array for the distribution.
    if (self%fixedPoint) then
       ! For a fixed point distribution we tabulate at a single mass.
       self%massCount=    1
       self%massDelta=    0.0d0
    else
       self%massCount=int(log10(self%massMaximum/self%massMinimum)*dble(massPointsPerDecade))+1
       self%massDelta=    log10(self%massMaximum/self%massMinimum)/dble(self%massCount-1)
    end if
    self   %spinCount=int(log10(self%spinMaximum/self%spinMinimum)*dble(spinPointsPerDecade))+1
    self   %spinDelta=    log10(self%spinMaximum/self%spinMinimum)/dble(self%spinCount-1)
    allocate(self%distributionTable(self%massCount,self%spinCount))
    ! Build a work node.
    node                  => treeNode                  (                 )
    nodeBasic             => node    %basic            (autoCreate=.true.)
    nodeSpin              => node    %spin             (autoCreate=.true.)
    nodeDarkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
    call nodeBasic%timeSet(self%time)
    if (present(tree)) node%hostTree => tree
    ! Build integrators.
    integratorMass    =integrator(massIntegral               ,toleranceAbsolute=1.0d-30,toleranceRelative=1.0d-3)
    integratorMassSpin=integrator(massSpinIntegral           ,toleranceAbsolute=1.0d-30,toleranceRelative=1.0d-3)
    integratorSpin    =integrator(spinIntegral                                         ,toleranceRelative=1.0d-3)
    integratorProduct =integrator(productDistributionIntegral                          ,toleranceRelative=1.0d-3)
    ! Tabulate the distribution.
    do iMass=1,self%massCount
       massMeasured=10.0d0**(dble(iMass-1)*self%massDelta+log10(self%massMinimum))
       ! Estimate the mass error.
       call nodeBasic%massSet (massMeasured)
       call nodeDarkMatterProfile%scaleSet(self%darkMatterProfileScaleRadius_%radius(node))
       call Calculations_Reset(node)
       massError=  +self%nbodyHaloMassError_%errorFractional(node) &
            &      *massMeasured
       ! Determine the correction to fractional error in potential energy if a maximum number of particles was used for estimation of potential energy.
       massEnergyEstimateMaximum=self%massParticle*self%energyEstimateParticleCountMaximum
       if (massMeasured < massEnergyEstimateMaximum) then
          ! Full sample of halo particles used in estimating potential energy - no correction needed.
          energyEstimateErrorCorrection=1.0d0
       else
          ! Reduced sample of halo particles used in estimating potential energy - compute correction factor.
          call nodeBasic%massSet(massEnergyEstimateMaximum)
          call nodeDarkMatterProfile%scaleSet(self%darkMatterProfileScaleRadius_%radius(node))
          call Calculations_Reset(node)
          energyEstimateErrorCorrection=+self%nbodyHaloMassError_%errorFractional(node) &
               &                        /(                                              &
               &                          +massError                                    &
               &                          /massMeasured                                 &
               &                         )
          call nodeBasic%massSet(massMeasured             )
          call nodeDarkMatterProfile%scaleSet(self%darkMatterProfileScaleRadius_%radius(node))
          call Calculations_Reset(node)
       end if
       ! Integrate over a range of masses corresponding to a fixed number of mass errors around the measured mass.
       massMinimum=max(                                &
            &          dble(self%particleCountMinimum) &
            &          *self%massParticle             ,&
            &          +massMeasured                   &
            &          -massIntegrationRange           &
            &          *massError                      &
            &         )
       massMaximum=    +massMeasured                   &
            &          +massIntegrationRange           &
            &          *massError
       ! Integrate the normalizing factor over halo mass.
       if (self%fixedPoint) then
          massAverage=1.0d0
       else
          massAverage=integratorMass%integrate(massMinimum,massMaximum)
       end if
       do iSpin=1,self%spinCount
          ! Evaluate the spin at this grid point - this corresponds to the measured spin in the N-body simulation.
          spinMeasured=10.0d0**(dble(iSpin-1)*self%spinDelta+log10(self%spinMinimum))
          ! Compute the distribution at the current value of mass and spin.
          if (self%fixedPoint) then
             self%distributionTable(iMass,iSpin)=massSpinIntegral(massMeasured)
          else
             ! Integrate over the intrinsic spin distribution to find the measured distribution at this measured mass and
             ! spin.
             massSpinAverage=integratorMassSpin%integrate(massMinimum,massMaximum)
             self%distributionTable(iMass,iSpin)=+massSpinAverage &
                  &                              /massAverage
          end if
       end do
    end do
    ! Clean up our work node.
    call node%destroy()
    deallocate(node)
    return

  contains

    double precision function massIntegral(massIntrinsic)
      !!{
      Integral over the halo mass function and halo mass error distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: massIntrinsic

      ! Set the mass and compute the mass error.
      call nodeBasic%massSet(massIntrinsic)
      call nodeDarkMatterProfile%scaleSet(self%darkMatterProfileScaleRadius_%radius(node))
      call Calculations_Reset(node)
      massError   =+self%nbodyHaloMassError_%errorFractional(node) &
           &       *massIntrinsic
      ! Evaluate the integrand. Normalization of the distribution function is neglected here since we are average over mass this
      ! will cancel out in numerator and denominator integrals.
      massIntegral=+self%haloMassFunction_%differential(nodeBasic%time(),massIntrinsic) &
           &       *exp(                   &
           &            -0.5d0             &
           &            *(                 &
           &              +(               &
           &                +massMeasured  &
           &                -massIntrinsic &
           &               )               &
           &              /massError       &
           &             )**2              &
           &           )
      return
    end function massIntegral

    double precision function massSpinIntegral(massIntrinsic)
      !!{
      Integral over the halo mass function, spin distribution, halo mass error distribution, and spin error distribution.
      !!}
      use :: Display                   , only : displayMagenta             , displayReset
      use :: Error                     , only : Error_Report               , Warn         , errorStatusFail, errorStatusSuccess
      use :: Input_Parameters          , only : inputParameters
      use :: Numerical_Constants_Math  , only : Pi
      use :: Mass_Distributions        , only : massDistributionClass
      use :: Galactic_Structure_Options, only : componentTypeDarkMatterOnly, massTypeDark
      use :: Coordinates               , only : coordinateSpherical        , assignment(=)
      implicit none
      double precision                       , intent(in   ) :: massIntrinsic
      class           (massDistributionClass), pointer       :: massDistribution_
      double precision                       , parameter     :: rangeIntegration                                           =1.0000d1 ! Range of integration in units of error.
      double precision                       , parameter     :: radiusVelocityDispersionMeanOverSpinSpecificAngularMomentum=0.4175d0 ! Ratio of mean velocity dispersion-radius product to "spin" angular
                                                                                                                                     ! momentum (i.e. the normalizing angular momentum appearing in the definition of
                                                                                                                                     ! halo spin): γ = Mσⱼ/Jₛ, Jₛ = GM²˙⁵/|E⁰˙⁵|. This parameter corresponds to γ.
      integer                                                :: errorStatus
      double precision                                       :: logSpinMinimum , logSpinMaximum      , &
           &                                                    errorMaximum   , scaleAbsolute       , &
           &                                                    tolerance      , massSpinIntegralMass
      character       (len=8                )                :: label          , labelMass           , &
           &                                                    labelSpin
      type            (inputParameters      )                :: descriptor
      type            (coordinateSpherical  )                :: coordinatesHalo

      ! Evaluate the halo mass part of the integrand, unless evaluating the distribution at a fixed point.
      if (self%fixedPoint) then
         massSpinIntegralMass       =  1.0d0
      else
         massSpinIntegralMass       =  massIntegral(massIntrinsic)
      end if
      ! Compute the particle number.
      particleNumber                =  massIntrinsic /self%massParticle
      ! Evaluate the root-variance of the spin-independent error term which arises from the random walk in angular momentum
      ! space. Note that the root-variance that goes into non-central χ-square distribution is the width of the Gaussian for a
      ! single dimension, leading to a factor of √3 in the following.
      errorSpinIndependent          =  +radiusVelocityDispersionMeanOverSpinSpecificAngularMomentum &
           &                           /sqrt(particleNumber)
      errorSpinIndependent1D        =  +errorSpinIndependent                                        &
           &                           /sqrt(3.0d0)
      ! Get the outer radius of the halo.
      radiusHalo                    =  +self%darkMatterHaloScale_ %radiusVirial(node)
      coordinatesHalo               =   [radiusHalo,0.0d0,0.0d0]
      ! Get the density at the edge of the halo.
      massDistribution_             =>  node             %massDistribution(componentTypeDarkMatterOnly,massTypeDark)
      densityOuterRadius            =  +massDistribution_%density         (coordinatesHalo                         )
      !![
      <objectDestructor name="massDistribution_"/>
      !!]
      ! Find the ratio of the mean interior density in the halo to the density at the halo outer radius.
      densityRatioInternalToSurface =  +3.0d0                 &
           &                           *massIntrinsic         &
           &                           /4.0d0                 &
           &                           /Pi                    &
           &                           /radiusHalo        **3 &
           &                           /densityOuterRadius
      ! Evaluate the distribution.
      if (self%fixedPoint) then
         massSpinIntegral=spinErrorIntegral(spinFixed)
      else
         ! Evaluate an estimate of the absolute scale of the spin distribution for use in setting an absolute precision level on the
         ! integration.
         call nodeSpin%angularMomentumSet(spinMeasured*Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_))
         call Calculations_Reset(node)
         scaleAbsolute=self%distributionIntrinsic%distribution(node)
         ! Evaluate the integral over the spin distribution. We integrate from ±ασ around the measured spin, with σ being the larger
         ! of the expected spin-independent and spin-dependent errors, and α a parameter typically set to 10.
         errorMaximum  =max(                                    &
              &             errorSpinIndependent              , &
              &             errorsSpinDependent (spinMeasured)  &
              &            )
         logSpinMinimum=log(                       &
              &             max(                   &
              &                 +self%spinMinimum, &
              &                 +spinMeasured      &
              &                 -rangeIntegration  &
              &                 *errorMaximum      &
              &                )                   &
              &            )
         logSpinMaximum=log(                       &
              &             min(                   &
              &                 +self%spinMaximum, &
              &                 +spinMeasured      &
              &                 +rangeIntegration  &
              &                 *errorMaximum      &
              &                )                   &
              &            )
         errorStatus=errorStatusFail
         tolerance  =1.0d-4
         do while (tolerance < 1.0d0 .and. errorStatus /= errorStatusSuccess)
            tolerance       =+10.0d0    &
                 &           *tolerance
            call integratorSpin%toleranceSet(toleranceAbsolute=1.0d-9*spinMeasured,toleranceRelative=tolerance)
            massSpinIntegral=+massSpinIntegralMass                                                       &
                 &           *integratorSpin%integrate(logSpinMinimum,logSpinMaximum,status=errorStatus) &
                 &           /2.0d0                                                                      & ! <= Partial combined normalization term for the lognormal and non-central chi-square distributions - brought
                 &           /Pi                                                                         & !    outside of integrand since constant. Each contributes √(2π).
                 &           *2.0d0                                                                      & ! <= Factor 2 appears due to change of variables from λ² to λ in the non-central χ² distribution.
                 &           /errorSpinIndependent1D**2                                                    ! <= Partial normalization term for the non-central chi-square distribution - brought outside of integrand since constant.
            if (errorStatus /= errorStatusSuccess) then
               write (label    ,'(e8.1)') tolerance
               write (labelMass,'(i4)'  ) iMass
               write (labelSpin,'(i4)'  ) iSpin
               descriptor=inputParameters()
               call self%distributionIntrinsic%descriptor(descriptor)
               call Warn(                                                                                               &
                    &       displayMagenta()//'WARNING:'//displayReset()//' failed to reach required tolerance ['    // &
                    &       trim(label    )                                                                          // &
                    &       '] in massSpinIntegral [iMass,iSpin='                                                    // &
                    &       trim(labelMass)                                                                          // &
                    &       ','                                                                                      // &
                    &       trim(labelSpin)                                                                          // &
                    &       '] - report on intrinsic spin distribution follows - reattempting with reduced tolerance'// &
                    &       char(10)                                                                                 // &
                    &       descriptor%serializeToString()                                                              &
                 &      )
               if (tolerance >= 1.0d0) call Error_Report('integral failed to reach any tolerance'//{introspection:location})
            end if
         end do
      end if
      return
    end function massSpinIntegral

    double precision function spinIntegral(logSpinIntrinsic)
      !!{
      Integral over the intrinsic spin distribution, and spin error distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: logSpinIntrinsic
      double precision                :: spinIntrinsic

      ! Compute intrinsic spin.
      spinIntrinsic=exp(logSpinIntrinsic)
      ! Set the intrinsic spin.
      call nodeSpin%angularMomentumSet(spinIntrinsic*Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_))
      call Calculations_Reset(node)
      ! Compute the integrand.
      spinIntegral=+self%distributionIntrinsic%distribution(node         ) & ! Weight by the intrinsic spin distribution.
           &       *spinIntrinsic                                          & ! Multiply by spin since our integration variable is log(spin).
           &       *spinErrorIntegral                      (spinIntrinsic)   ! Multiply by the integral over the product distribution of spin-dependent and spin-independent errors.
      return
    end function spinIntegral

    double precision function spinErrorIntegral(spinIntrinsic)
      !!{
      Integral over the spin error distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: spinIntrinsic
      double precision, parameter     :: logNormalMean                              =1.0000d0  ! Mean of the log-normal distribution of mass/energy errors. We assume an unbiased
                                                                                               ! measurement, so the mean is unity.
      double precision, parameter     :: rangeIntegration                           =1.0000d1  ! Integration range (in ~σ - the width of each distribution).
      double precision                :: logSpinMinimum                                      , logSpinMaximum                     , &
           &                             nonCentralChiSquaredMode                            , nonCentralChiSquaredModeLogarithmic, &
           &                             nonCentralChiSquaredRootVarianceLogarithmic

      ! Evaluate the root-variance of the spin-dependent error term which arises from errors in the measurement of mass and
      ! energy.
      errorSpinDependent=errorsSpinDependent(spinIntrinsic)
      ! Evaluate non-centrality parameter of the non-central χ-square distribution used to model the spin-independent
      ! (i.e. random walk in angular momentum space) part of the measured spin distribution.
      nonCentrality=(                        &
           &         +spinIntrinsic          &
           &         /errorSpinIndependent1D &
           &        )**2
      ! Evaluate the (log-)mean and width of the log-normal distribution used to model the multiplicative errors arising from
      ! errors in the measurement of mass and energy. We choose values such that we obtain the desired mean scaling factor for
      ! the spin, and a width which ensures we reproduce the variance predicted by the model.
      logNormalWidth  =+sqrt(                            &
           &                 +log(                       &
           &                      +1.0d0                 &
           &                      +(                     &
           &                        +errorSpinDependent  &
           &                        /spinIntrinsic       &
           &                       )**2                  &
           &                      /logNormalMean         &
           &                     )                       &
           &                )
      logNormalLogMean=+log(logNormalMean)               &
           &           -0.5d0                            &
           &           *logNormalWidth**2
      ! Find the logarithm of the spin corresponding to the mode of the non-central chi-squared distribution.
      nonCentralChiSquaredMode                   =nbodyErrorsNonCentralChiSquareMode(nonCentrality)
      nonCentralChiSquaredModeLogarithmic        =log (sqrt(nonCentralChiSquaredMode)*errorSpinIndependent1D)
      nonCentralChiSquaredRootVarianceLogarithmic=0.5d0*sqrt(2.0d0*(3.0d0+2.0d0*nonCentrality))/nonCentralChiSquaredMode
      ! Determine a suitable integration range. The range is centered on the mode of the non-central chi-squared distribution
      ! offset by the mean of the log-normal distribution to ensure that we always include the peak of the distribution. The
      ! integration range around this mode is then chosen such that the log-normal distribution will have decayed by a large
      ! factor at the edges of the range to ensure that we are not missing any significant contributions from outside of the
      ! integration range.
      logSpinMinimum=min(                                                                                                                   &
           &             nonCentralChiSquaredModeLogarithmic                 -rangeIntegration*nonCentralChiSquaredRootVarianceLogarithmic, &
           &             nonCentralChiSquaredModeLogarithmic-logNormalLogMean-rangeIntegration*logNormalWidth                               &
           &            )
      logSpinMaximum=max(                                                                                                                   &
           &             nonCentralChiSquaredModeLogarithmic                 +rangeIntegration*nonCentralChiSquaredRootVarianceLogarithmic, &
           &             nonCentralChiSquaredModeLogarithmic-logNormalLogMean+rangeIntegration*logNormalWidth                               &
           &            )
      ! Limit the range of the log-normal distribution over which we integrate. This is useful to prevent inclusion of halos which
      ! would be far from equilibrium as judged by an energy condition.
      logSpinMinimum=max(logSpinMinimum,log(spinMeasured/self%logNormalRange))
      logSpinMaximum=min(logSpinMaximum,log(spinMeasured*self%logNormalRange))
      ! Evaluate the integrand.
      call integratorProduct%toleranceSet(toleranceAbsolute=1.0d-9*spinIntrinsic,toleranceRelative=1.0d-3)
      spinErrorIntegral=+integratorProduct%integrate(                & ! Multiply by the integral which gives us the product distribution
           &                                         logSpinMinimum, & ! of log-normal and non-central χ-square distributions.
           &                                         logSpinMaximum  &
           &                                        )                &
           &            /logNormalWidth                              & ! <= Partial normalization term for the lognormal distribution - brought outside of integrand since constant.
           &            /sqrt(nonCentrality)                           ! <= Partial normalization term for the non-central chi-square distribution - brought outside of integrand since constant.      
      return
    end function spinErrorIntegral

    double precision function nbodyErrorsNonCentralChiSquareMode(chi)
      !!{
      Computes the mode of the degree-3 non-central chi-squared distribution function for given $\chi$.
      !!}
      implicit none
      double precision, intent(in   ) :: chi
      
      nonCentralChiSquareChi            =                                     chi
      nbodyErrorsNonCentralChiSquareMode=self%finder%find(rootGuess=max(1.0d0,chi))
      return
    end function nbodyErrorsNonCentralChiSquareMode

    double precision function productDistributionIntegral(logSpinUnscaled)
      !!{
      Product distribution integrand.
      !!}
      implicit none
      double precision, intent(in   ) :: logSpinUnscaled
      double precision                :: x                     , y                    , &
           &                             distributionNonCentral, distributionLogNormal, &
           &                             spinUnscaled          , sinhArgument

      ! Find the unscaled spin. This is the spin resulting from the intrinsic spin plus spin-independent (i.e. random walk in
      ! angular momentum space) error but without any scaling due to errors in mass/energy.
      spinUnscaled=exp(logSpinUnscaled)
      ! Evaluate the non-central χ-square distribution. Note that the non-central distribution is divided by the 1D
      ! spin-independent error squared, since it is the distribution for the scaled parameter "x=λ²/σ₁²" and we need to convert
      ! it back to be a distribution for the unscaled spin. Note also that the distribution function represents a simplified form
      ! of the general non-central χ-square distribution which is specific to the case of k=3 degrees of freedom. The exponential
      ! terms (i.e. the product of the exp() and sinh() functions) are handled separately such that if the argument to sinh()
      ! becomes very large we approximate sinh(z)=exp(z)/2 and combine this into the exp() function to avoid floating point
      ! overflows.
      x                     =(                        &
           &                  +spinUnscaled           &
           &                  /errorSpinIndependent1D &
           &                 )**2
      sinhArgument          =+sqrt(                   &
           &                       +x                 &
           &                       *nonCentrality     &
           &                      )
      if (sinhArgument < 50.0d0) then
         distributionNonCentral=+exp (                     &
              &                       -0.5d0               &
              &                       *    (               &
              &                             +x             &
              &                             +nonCentrality &
              &                            )               &
              &                      )                     &
              &                 *sinh(                     &
              &                       +sqrt(               &
              &                             +x             &
              &                             *nonCentrality &
              &                            )               &
              &                      )
      else
         distributionNonCentral=+exp (                     &
              &                       -0.5d0               &
              &                       *    (               &
              &                             +x             &
              &                             +nonCentrality &
              &                            )               &
              &                       +sinhArgument        &
              &                      )                     &
              &                 /2.0d0
      end if
      ! No need to compute the log-normal distribution if the non-central chi-squared distribution is zero.
      if (distributionNonCentral <= 0.0d0) then
         productDistributionIntegral=0.0d0
      else
         ! Evaluate the log-normal distribution. The parameter "y" here is the ratio of the measured spin to the unscaled spin. Note
         ! that constant normalization terms are excluded here, and are instead included external to the integrand
         y                    =+spinMeasured            &
              &                /spinUnscaled
         distributionLogNormal=+exp(                    &
              &                     -(                  &
              &                       +log(y)           &
              &                       -logNormalLogMean &
              &                      )**2               &
              &                     /2.0d0              &
              &                     /logNormalWidth**2  &
              &                    )                    &
              &                /y
         ! Construct the integrand of the product distribution. There are various factor of "spinUnscaled" that should appear here:
         !  1) In the product distribution we divide through by "spinUnscaled".;
         !  2) Our integration variable is log(spinUnscaled) so we multiply by "spinUnscaled";
         !  3) The non-central χ² distribution is for λ² - changing variables to the distribution for λ requires multiplying by
         !     "spinUnscaled" (and dividing by 2 which we take out of the integration as a constant).
         productDistributionIntegral=+distributionNonCentral &
              &                      *distributionLogNormal  &
              &                      *spinUnscaled
      end if
      return
    end function productDistributionIntegral

    double precision function errorsSpinDependent(spinIntrinsic)
      !!{
      Compute the spin-dependent error term.
      !!}
      implicit none
      double precision, intent(in   ) :: spinIntrinsic
      double precision, parameter     :: spinDependentErrorCoefficient=0.4025d0 ! Coefficient of the spin-dependent error term, β.

      errorsSpinDependent=+(                                   &
           &                +2.5d0                             &
           &                +energyEstimateErrorCorrection     &
           &                *(                                 &
           &                  +1.0d0                           &
           &                  +(                               &
           &                    +densityRatioInternalToSurface &
           &                    /3.0d0                         &
           &                   )                               &
           &                  /2.0d0                           &
           &                 )                                 &
           &               )                                   &
           &              *spinDependentErrorCoefficient       &
           &              *spinIntrinsic                       &
           &              *(                                   &
           &                +massError                         &
           &                /nodeBasic%mass()                  &
           &               )
      return
    end function errorsSpinDependent

  end subroutine nbodyErrorsTabulate

  subroutine nbodyErrorsDestructor(self)
    !!{
    Destructor for the \refClass{haloSpinDistributionNbodyErrors} dark matter halo spin distribution class.
    !!}
    implicit none
    type(haloSpinDistributionNbodyErrors), intent(inout) :: self

    !![
    <objectDestructor name="self%distributionIntrinsic"        />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%nbodyHaloMassError_"          />
    <objectDestructor name="self%haloMassFunction_"            />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine nbodyErrorsDestructor

  double precision function nbodyErrorsSample(self,node)
    !!{
    Sample from the halo spin distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(haloSpinDistributionNbodyErrors), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nbodyErrorsSample=0.0d0
    call Error_Report('sampling from distribution is not implemented'//{introspection:location})
    return
  end function nbodyErrorsSample

  double precision function nbodyErrorsDistribution(self,node)
    !!{
    Compute the spin distribution.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentBasic                     , nodeComponentSpin, treeNode
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    class           (nodeComponentBasic             ), pointer       :: nodeBasic
    class           (nodeComponentSpin              ), pointer       :: nodeSpin
    double precision                                                 :: mass     , spin , &
         &                                                              hMass    , hSpin
    integer                                                             iMass    , iSpin, &
         &                                                              jMass    , jSpin

    ! Extract the mass and spin of the halo.
    nodeBasic =>  node     %basic          ()
    nodeSpin  =>  node     %spin           ()
    mass      =   nodeBasic%mass           ()
    spin      =  +nodeSpin %angularMomentum()                                             &
         &       /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_)
    ! Ensure the table has sufficient extent.
    call self%tabulate(massRequired=mass,spinRequired=spin,tree=node%hostTree)
    ! Find the interpolating factors.
    hMass=log10(mass/self%massMinimum)/self%massDelta
    hSpin=log10(spin/self%spinMinimum)/self%spinDelta
    iMass=int(hMass)+1
    iSpin=int(hSpin)+1
    hMass=hMass-dble(iMass-1)
    hSpin=hSpin-dble(iSpin-1)
    jSpin=min(iSpin+1,self%spinCount)
    jMass=min(iMass+1,self%massCount)
    ! Perform the interpolation.
    nbodyErrorsDistribution=+self%distributionTable(iMass,iSpin)*(1.0d0-hMass)*(1.0d0-hSpin) &
         &                  +self%distributionTable(iMass,jSpin)*(1.0d0-hMass)*       hSpin  &
         &                  +self%distributionTable(jMass,iSpin)*       hMass *(1.0d0-hSpin) &
         &                  +self%distributionTable(jMass,jSpin)*       hMass *       hSpin
    return
  end function nbodyErrorsDistribution

  double precision function nbodyErrorsDistributionFixedPoint(self,node,spinMeasured,spinMeasuredMinimum,spinMeasuredMaximum)
    !!{
    Compute the spin distribution for a fixed point in intrinsic mass and spin.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentBasic                     , nodeComponentSpin, treeNode
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: spinMeasured       , spinMeasuredMinimum, &
         &                                                              spinMeasuredMaximum
    class           (nodeComponentBasic             ), pointer       :: nodeBasic
    class           (nodeComponentSpin              ), pointer       :: nodeSpin
    double precision                                                 :: spin               , hSpin              , &
         &                                                              mass
    integer                                                             iSpin              , jSpin

    ! Extract the mass and spin of the halo.
    nodeBasic =>  node     %basic          ()
    nodeSpin  =>  node     %spin           ()
    mass      =   nodeBasic%mass           ()
    spin      =  +nodeSpin %angularMomentum()                                             &
         &       /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_)
    ! Ensure the table has sufficient extent.
    call self%tabulate(massFixed=mass,spinFixed=spin,spinFixedMeasuredMinimum=spinMeasuredMinimum,spinFixedMeasuredMaximum=spinMeasuredMaximum,tree=node%hostTree)
    ! Find the interpolating factors.
    hSpin=log10(spinMeasured/self%spinMinimum)/self%spinDelta
    iSpin=int(hSpin)+1
    hSpin=hSpin-dble(iSpin-1)
    jSpin=min(iSpin+1,self%spinCount)
    ! Perform the interpolation.
    nbodyErrorsDistributionFixedPoint=+self%distributionTable(1,iSpin)*(1.0d0-hSpin) &
         &                            +self%distributionTable(1,jSpin)*       hSpin
    return
  end function nbodyErrorsDistributionFixedPoint

  double precision function nbodyErrorsDistributionAveraged(self,node,massLimit)
    !!{
    Compute the spin distribution averaged over all halos more massive than the given {\normalfont \ttfamily massLimit}.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Galacticus_Nodes   , only : nodeComponentBasic, treeNode
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: massLimit
    class           (nodeComponentBasic             ), pointer       :: nodeBasic
    double precision                                                 :: mass     , massOriginal, &
         &                                                              massLow  , massHigh    , &
         &                                                              weight   , weightTotal
    integer                                                          :: iMass

    ! Sum the spin distribution over all tabulated points in mass, weighting by the number of halos in that mass range.
    nodeBasic                       => node     %basic()
    massOriginal                    =  nodeBasic%mass ()
    nbodyErrorsDistributionAveraged =  0.0d0
    weightTotal                     =  0.0d0
    do iMass=1,self%massCount
       mass    =    10.0d0**((dble(iMass)-1.0d0)*self%massDelta+log10(self%massMinimum))
       massLow =max(10.0d0**((dble(iMass)-1.5d0)*self%massDelta+log10(self%massMinimum)),massLimit)
       massHigh=    10.0d0**((dble(iMass)-0.5d0)*self%massDelta+log10(self%massMinimum))
       if (massHigh > massLow) then
          call nodeBasic%massSet(mass)
          call Calculations_Reset(node)
          weight                         =self%haloMassFunction_%integrated(nodeBasic%time(),massLow,massHigh)
          nbodyErrorsDistributionAveraged=nbodyErrorsDistributionAveraged+self%distribution(node)*weight
          weightTotal                    =weightTotal                    +                        weight
       end if
    end do
    ! Normalize by the total weight.
    if (weightTotal > 0.0d0) nbodyErrorsDistributionAveraged=nbodyErrorsDistributionAveraged/weightTotal
    call nodeBasic%massSet(massOriginal)
    call Calculations_Reset(node)
    return
  end function nbodyErrorsDistributionAveraged
  
  double precision function nbodyErrorsNonCentralChiSquareModeRoot(x)
    !!{
    Root function used in finding the mode of the degree-3 non-central chi-squared distribution function.
    !!}
    implicit none
    double precision, intent(in   ) :: x
    
    nbodyErrorsNonCentralChiSquareModeRoot=+      sqrt(                        &
         &                                             +x                      &
         &                                             *nonCentralChiSquareChi &
         &                                            )                        &
         &                                 *tanh(                              &
         &                                       +sqrt(                        &
         &                                             +x                      &
         &                                             *nonCentralChiSquareChi &
         &                                            )                        &
         &                                      )                              &
         &                                 -            nonCentralChiSquareChi
    return
  end function nbodyErrorsNonCentralChiSquareModeRoot
