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

  !% An implementation of the dark matter halo spin distribution which modifies another spin distribution to account for the
  !% effects of particle noise errors in spins measured in N-body simulations.
  
  use Statistics_NBody_Halo_Mass_Errors
  use Cosmology_Functions
  use Halo_Mass_Functions
  use Dark_Matter_Profiles
  use Dark_Matter_Halo_Scales

  !# <haloSpinDistribution name="haloSpinDistributionNbodyErrors">
  !#  <description>
  !#   A halo spin distribution which modifies another spin distribution to account for the effects of particle noise errors
  !#   in spins measured in N-body simulations.
  !#  </description>
  !# </haloSpinDistribution>
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionNbodyErrors
     !% A dark matter halo spin distribution class which modifies another spin distribution to account for the effects of particle
     !% noise errors in spins measured in N-body simulations.
     private
     class           (haloSpinDistributionClass), pointer                     :: distributionIntrinsic => null()
     class           (nbodyHaloMassErrorClass  ), pointer                     :: nbodyHaloMassError_   => null()
     class           (haloMassFunctionClass    ), pointer                     :: haloMassFunction_     => null()
     class           (darkMatterHaloScaleClass ), pointer                     :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileClass   ), pointer                     :: darkMatterProfile_    => null()
     double precision                                                         :: massParticle                      , time
     integer                                                                  :: particleCountMinimum
     integer                                                                  :: spinCount                         , massCount
     double precision                                                         :: spinMinimum                       , spinMaximum, &
          &                                                                      massMinimum                       , massMaximum, &
          &                                                                      massDelta                         , spinDelta  , &
          &                                                                      energyEstimateParticleCountMaximum
     double precision                           , allocatable, dimension(:  ) :: massWeight
     double precision                           , allocatable, dimension(:,:) :: distributionTable
   contains
     !@ <objectMethods>
     !@   <object>haloSpinDistributionNbodyErrors</object>
     !@   <objectMethod>
     !@     <method>distributionAveraged</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\arginout, \doublezero\ massLimit\argin</arguments>
     !@     <description>Return the spin distribution function averaged over all halos above the given {\normalfont \ttfamily massLimit}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ [massRequired]\argin, \doublezero\ [spinRequired]\argin</arguments>
     !@     <description>Tabulate the spin distribution as a fuction of spin and halo mass. Ensure that the table spans the {\normalfont \ttfamily massRequired} and {\normalfont \ttfamily spinRequireed} if provided.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                         nbodyErrorsDestructor
     procedure :: sample               => nbodyErrorsSample
     procedure :: distribution         => nbodyErrorsDistribution
     procedure :: distributionAveraged => nbodyErrorsDistributionAveraged
     procedure :: tabulate             => nbodyErrorsTabulate
  end type haloSpinDistributionNbodyErrors

  interface haloSpinDistributionNbodyErrors
     !% Constructors for the {\normalfont \ttfamily nbodyErrors} dark matter halo spin
     !% distribution class.
     module procedure nbodyErrorsConstructorParameters
     module procedure nbodyErrorsConstructorInternal
  end interface haloSpinDistributionNbodyErrors

  ! Tabulation parameters.
  integer         , parameter :: nbodyErrorsSpinPointsPerDecade=33    , nbodyErrorsMassPointsPerDecade=5
  double precision, parameter :: nbodyErrorsSpinMaximum        =0.5d+0, nbodyErrorsMassMaximum        =1.0d15
  double precision, parameter :: nbodyErrorsSpinMinimum        =3.0d-4, nbodyErrorsMassMinimum        =1.0d11

contains

  function nbodyErrorsConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily nbodyErrors} dark matter halo spin
    !% distribution class which takes a parameter list as input.
    use Input_Parameters2
    implicit none
    type            (haloSpinDistributionNbodyErrors)                :: nbodyErrorsConstructorParameters
    type            (inputParameters                ), intent(inout) :: parameters
    class           (haloSpinDistributionClass      ), pointer       :: distributionIntrinsic
    class           (nbodyHaloMassErrorClass        ), pointer       :: nbodyHaloMassError_
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (haloMassFunctionClass          ), pointer       :: haloMassFunction_
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileClass         ), pointer       :: darkMatterProfile_
    double precision                                                 :: massParticle                    , redshift                          , &
         &                                                              time                            , energyEstimateParticleCountMaximum
    integer                                                          :: particleCountMinimum
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>massParticle</variable>
    !#   <description>Mass of particle in the simulation.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>particleCountMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>particleCountMinimum</variable>
    !#   <description>Minimum number of particles per halo in the N-body simulation.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>energyEstimateParticleCountMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>energyEstimateParticleCountMaximum</variable>
    !#   <description>Maximum number of particles used in estimating the potential energy of halos. Set to a very large number if no such maximum was used.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <variable>redshift</variable>
    !#   <description>Redshift at which the spin distribution should be evaluated.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="haloSpinDistribution" name="distributionIntrinsic" source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError"   name="nbodyHaloMassError_"   source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !# <objectBuilder class="haloMassFunction"     name="haloMassFunction_"     source="parameters"/>
    darkMatterHaloScale_ => darkMatterHaloScale()
    darkMatterProfile_   => darkMatterProfile  ()
    ! Find the time corresponding to the given redshift.
    time=  cosmologyFunctions_ %cosmicTime                 (          &
         &  cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                   redshift &
         &                                                  )         &
         &                                                 )
    !# <objectDestructor name="cosmologyFunctions_"/>
    ! Construct the object.
    nbodyErrorsConstructorParameters=nbodyErrorsConstructorInternal(                                    &
         &                                                          distributionIntrinsic             , &
         &                                                          massParticle                      , &
         &                                                          particleCountMinimum              , &
         &                                                          energyEstimateParticleCountMaximum, &
         &                                                          time                              , &
         &                                                          nbodyHaloMassError_               , &
         &                                                          haloMassFunction_                 , &
         &                                                          darkMatterHaloScale_              , &
         &                                                          darkMatterProfile_                  &
         &                                                         )
    return
  end function nbodyErrorsConstructorParameters

  function nbodyErrorsConstructorInternal(distributionIntrinsic,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,time,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfile_)
    !% Internal constructor for the {\normalfont \ttfamily nbodyErrors} dark matter halo spin distribution class.
    implicit none
    type            (haloSpinDistributionNbodyErrors)                        :: nbodyErrorsConstructorInternal
    class           (haloSpinDistributionClass      ), intent(in   ), target :: distributionIntrinsic
    class           (nbodyHaloMassErrorClass        ), intent(in   ), target :: nbodyHaloMassError_
    class           (haloMassFunctionClass          ), intent(in   ), target :: haloMassFunction_
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileClass         ), intent(in   ), target :: darkMatterProfile_
    double precision                                 , intent(in   )         :: massParticle                      , time, &
         &                                                                      energyEstimateParticleCountMaximum
    integer                                          , intent(in   )         :: particleCountMinimum
    ! Store properties.
    nbodyErrorsConstructorInternal%distributionIntrinsic               => distributionIntrinsic
    nbodyErrorsConstructorInternal%nbodyHaloMassError_                 => nbodyHaloMassError_
    nbodyErrorsConstructorInternal%haloMassFunction_                   => haloMassFunction_
    nbodyErrorsConstructorInternal%darkMatterHaloScale_                => darkMatterHaloScale_
    nbodyErrorsConstructorInternal%darkMatterProfile_                  => darkMatterProfile_
    nbodyErrorsConstructorInternal%massParticle                        =  massParticle
    nbodyErrorsConstructorInternal%particleCountMinimum                =  particleCountMinimum
    nbodyErrorsConstructorInternal%energyEstimateParticleCountMaximum  =  energyEstimateParticleCountMaximum
    nbodyErrorsConstructorInternal%time                                =  time
    ! Set default ranges of spin and mass for tabulation.
    nbodyErrorsConstructorInternal%spinMinimum=nbodyErrorsSpinMinimum
    nbodyErrorsConstructorInternal%spinMaximum=nbodyErrorsSpinMaximum
    nbodyErrorsConstructorInternal%massMinimum=nbodyErrorsMassMinimum
    nbodyErrorsConstructorInternal%massMaximum=nbodyErrorsMassMaximum
    ! Tabulate the distribution function.
    call nbodyErrorsConstructorInternal%tabulate()
    return
  end function nbodyErrorsConstructorInternal

  subroutine nbodyErrorsTabulate(self,massRequired,spinRequired)
    !% Tabulate the halo spin distribution.
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Integration
    use               Memory_Management
    use               Galacticus_Nodes
    use               Numerical_Constants_Math
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout)           :: self
    double precision                                 , intent(in   ), optional :: massRequired                        , spinRequired
    type            (treeNode                       ), pointer                 :: node
    class           (nodeComponentBasic             ), pointer                 :: nodeBasic
    class           (nodeComponentSpin              ), pointer                 :: nodeSpin
    class           (nodeComponentDarkMatterProfile ), pointer                 :: nodeDarkMatterProfile
    double precision                                 , parameter               :: massIntegrationRange         =10.0d0
    type            (fgsl_function                  )                          :: integrandFunctionMass               , integrandFunctionSpin      , &
         &                                                                        integrandFunctionMassSpin           , integrandFunctionProduct
    type            (fgsl_integration_workspace     )                          :: integrationWorkspaceMass            , integrationWorkspaceSpin   , &
         &                                                                        integrationWorkspaceMassSpin        , integrationWorkspaceProduct
    logical                                                                    :: retabulate
    integer                                                                    :: iSpin                               , iMass
    double precision                                                           :: spinMeasured                        , massMeasured               , &
         &                                                                        densityRatioInternalToSurface       , massError                  , &
         &                                                                        radiusHalo                          , densityOuterRadius         , &
         &                                                                        massMinimum                         , massMaximum                , &
         &                                                                        logNormalWidth                      , logNormalLogMean           , &
         &                                                                        errorSpinDependent                  , errorSpinIndependent       , &
         &                                                                        errorSpinIndependent1D              , nonCentrality              , &
         &                                                                        particleNumber                      , massEnergyEstimateMaximum  , &
         &                                                                        energyEstimateErrorCorrection

    ! Determine if the tabulation needs to be rebuilt.
    if (allocated(self%distributionTable)) then
       retabulate=.false.
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
       if (retabulate)  call deallocateArray(self%distributionTable)
    else
       retabulate=.true.
    end if
    if (.not.retabulate) return
    ! Allocate array for the distribution.
    self%massCount=int(log10(self%massMaximum/self%massMinimum)*dble(nbodyErrorsMassPointsPerDecade))+1
    self%spinCount=int(log10(self%spinMaximum/self%spinMinimum)*dble(nbodyErrorsSpinPointsPerDecade))+1
    self%massDelta=    log10(self%massMaximum/self%massMinimum)/dble(self%massCount-1)
    self%spinDelta=    log10(self%spinMaximum/self%spinMinimum)/dble(self%spinCount-1)
    call allocateArray(self%distributionTable,[self%massCount,self%spinCount])
    ! Build a work node.
    node                  => treeNode                  (                 )
    nodeBasic             => node    %basic            (autoCreate=.true.)
    nodeSpin              => node    %spin             (autoCreate=.true.)
    nodeDarkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
    call nodeBasic%timeSet(self%time)
    ! Tabulate the distribution.
    do iMass=1,self%massCount
       massMeasured=10.0d0**(dble(iMass-1)*self%massDelta+log10(self%massMinimum))
       ! Estimate the mass error.
       call nodeBasic%massSet(massMeasured)
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
          energyEstimateErrorCorrection=+self%nbodyHaloMassError_%errorFractional(node) &
               &                        /(                                              &
               &                          +massError                                    &
               &                          /massMeasured                                 &
               &                         )
          call nodeBasic%massSet(massMeasured             )
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
       do iSpin=1,self%spinCount
          ! Evaluate the spin at this grid point - this corresponds to the measured spin in the N-body simulation.
          spinMeasured=10.0d0**(dble(iSpin-1)*self%spinDelta+log10(self%spinMinimum))      
          ! Integrate over the intrinsic spin distribution to find the measured distribution at this measureed mass and
          ! spin. 
          self%distributionTable(iMass,iSpin)=+Integrate(                                     &
               &                                         massMinimum                        , &
               &                                         massMaximum                        , &
               &                                         massSpinIntegral                   , &
               &                                         integrandFunctionMassSpin          , &
               &                                         integrationWorkspaceMassSpin       , &
               &                                         toleranceAbsolute           =0.0d+0, &
               &                                         toleranceRelative           =1.0d-3  &
               &                                        )                                     &
               &                              /Integrate(                                     &
               &                                         massMinimum                        , &
               &                                         massMaximum                        , &
               &                                         massIntegral                       , &
               &                                         integrandFunctionMass              , &
               &                                         integrationWorkspaceMass           , &
               &                                         toleranceAbsolute           =0.0d+0, &
               &                                         toleranceRelative           =1.0d-3  &
               &                                        )
          call Integrate_Done(integrandFunctionMassSpin,integrationWorkspaceMassSpin)
          call Integrate_Done(integrandFunctionMass    ,integrationWorkspaceMass    )
       end do
    end do
    ! Clean up our work node.
    call node%destroy()
    deallocate(node)
    return

  contains

    double precision function massIntegral(massIntrinsic)
      !% Integral over the halo mass function and halo mass error distribution.
      implicit none
      double precision, intent(in   ) :: massIntrinsic

      ! Set the mass and compute the mass error.
      call nodeBasic%massSet(massIntrinsic)
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
      !% Integral over the halo mass function, spin distribution, halo mass error distribution, and spin error distribution.
      implicit none
      double precision, intent(in   ) :: massIntrinsic
      double precision, parameter     :: rangeIntegration                                           =1.0000d1 ! Range of integration in units of error.
      double precision, parameter     :: radiusVelocityDispersionMeanOverSpinSpecificAngularMomentum=0.4175d0 ! Ratio of mean velocity dispersion-radius product to "spin" angular
                                                                                                              ! momentum (i.e. the normalizing angular momentum appearing in the definition of
                                                                                                              ! halo spin): γ = Mσⱼ/Jₛ, Jₛ = GM²˙⁵/|E⁰˙⁵|. This parameter corresponds to γ.
      double precision                :: logSpinMinimum, logSpinMaximum, &
           &                             errorMaximum
      
      ! Evaluate the halo mass part of the integrand.
      massSpinIntegral             =massIntegral(massIntrinsic)
      ! Compute the particle number.
      particleNumber               =massIntrinsic/self%massParticle
      ! Evaluate the root-variance of the spin-independent error term which arises from the random walk in angular momentum
      ! space. Note that the root-variance that goes into non-central χ-square distribution is the width of the Gaussian for a
      ! single dimension, leading to a factor of √3 in the following.
      errorSpinIndependent         =+radiusVelocityDispersionMeanOverSpinSpecificAngularMomentum &
           &                        /sqrt(particleNumber)       
      errorSpinIndependent1D       =+errorSpinIndependent                                        &
           &                        /sqrt(3.0d0)
      ! Get the outer radius of the halo.
      radiusHalo                   =+self%darkMatterHaloScale_%virialRadius(node           )
      ! Get the density at the edge of the halo.
      densityOuterRadius           =+self%darkMatterProfile_  %density     (node,radiusHalo)
      ! Find the ratio of the mean interior density in the halo to the density at the halo outer radius.
      densityRatioInternalToSurface=+3.0d0                 &
           &                        *massIntrinsic         &
           &                        /4.0d0                 &
           &                        /Pi                    &
           &                        /radiusHalo        **3 &
           &                        /densityOuterRadius
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
      massSpinIntegral=+massSpinIntegral                           &
           &           *Integrate(                                 &
           &                      logSpinMinimum                 , &
           &                      logSpinMaximum                 , &
           &                      spinIntegral                   , &
           &                      integrandFunctionSpin          , &
           &                      integrationWorkspaceSpin       , &
           &                      toleranceAbsolute       =1.0d-8, &
           &                      toleranceRelative       =1.0d-3  &
           &                    )                                  & 
           &           /2.0d0                                      & ! <= Partial combined normalization term for the lognormal and non-central chi-square distributions - brought
           &           /Pi                                         & !    outside of integrand since constant. Each contributes √(2π).
           &           *2.0d0                                      & ! <= Factor 2 appears due to change of variables in from λ² to λ in non-central χ² distribution.
           &           /errorSpinIndependent1D**2                    ! <= Partial normalization term for the non-central chi-square distribution - brought outside of integrand since constant.
      call Integrate_Done(integrandFunctionSpin,integrationWorkspaceSpin)
      return
    end function massSpinIntegral

    double precision function spinIntegral(logSpinIntrinsic)
      !% Integral over the intrinsic spin distribution, and spin error distribution.
      implicit none
      double precision, intent(in   ) :: logSpinIntrinsic
      double precision, parameter     :: logNormalMean   =1.0000d0                    ! Mean of the log-normal distribution of mass/energy errors. We assume an unbiased
                                                                                      ! measurement, so the mean is unity.
      double precision, parameter     :: rangeIntegration=1.0000d1                    ! Integration range (in ~σ - the width of each distribution).
      double precision                :: spinIntrinsic            , logSpinMinimum, &
           &                             logSpinMaximum

      ! Compute intrinsic spin.
      spinIntrinsic=exp(logSpinIntrinsic)
      ! Set the intrinsic spin.
      call nodeSpin%spinSet(spinIntrinsic)
      ! Evaluate the root-variance of the spin-dependent error term which arises from errors in the measurement of mass and
      ! energy.
      errorSpinDependent=errorsSpinDependent(spinIntrinsic)
      ! Evaluate non-centraility parameter of the non-central χ-square distribution used to model the spin-independent
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
      ! Evaluate the integrand.
      logSpinMinimum=    max(                                                                         &
           &                 log(self%spinMinimum )                                                 , &
           &                 log(     spinMeasured)-logNormalLogMean-rangeIntegration*logNormalWidth  &
           &                )
      logSpinMaximum=min(                                                                             &
           &             min(                                                                         &
           &                 log(self%spinMaximum )                                                 , &
           &                 log(     spinMeasured)-logNormalLogMean+rangeIntegration*logNormalWidth  &
           &                )                                                                       , &
           &                 rangeIntegration*errorSpinIndependent1D*nonCentrality                    &
           &            )
      spinIntegral=+self%distributionIntrinsic%distribution(node)      & ! Weight by the intrinsic spin distribution.
           &       *spinIntrinsic                                      & ! Multiply by spin since our integration variable is log(spin).
           &       *Integrate(                                         & ! Multiply by the integral which gives us the product distribution
           &                  logSpinMinimum                         , & ! of log-normal and non-central χ-square distributions.
           &                  logSpinMaximum                         , &
           &                  productDistributionIntegral            , &
           &                  integrandFunctionProduct               , &
           &                  integrationWorkspaceProduct            , &
           &                  toleranceAbsolute       =+1.0d-6         &
           &                                           /spinIntrinsic, &
           &                  toleranceRelative       =+1.0d-3         &
           &                 )                                         & 
           &       /logNormalWidth                                     & ! <= Partial normalization term for the lognormal distribution - brought outside of integrand since constant.
           &       /sqrt(nonCentrality)                                  ! <= Partial normalization term for the non-central chi-square distribution - brought outside of integrand since constant.
      call Integrate_Done(integrandFunctionProduct,integrationWorkspaceProduct)
      return
    end function spinIntegral

    double precision function productDistributionIntegral(logSpinUnscaled)
      !% Product distribution integrand.
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
      sinhArgument          =+sqrt(                  &
           &                       +x                &
           &                       *nonCentrality    &
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
      return
    end function productDistributionIntegral

    double precision function errorsSpinDependent(spinIntrinsic)
      !% Compute the spin-dependent error term.
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
    !% Destructor for the {\normalfont \ttfamily nbodyErrors} dark matter halo spin distribution class.
    implicit none
    type(haloSpinDistributionNbodyErrors), intent(inout) :: self

    !# <objectDestructor name="self%distributionIntrinsic"/>
    !# <objectDestructor name="self%nbodyHaloMassError_"  />
    !# <objectDestructor name="self%haloMassFunction_"    />
    return
  end subroutine nbodyErrorsDestructor

  double precision function nbodyErrorsSample(self,node)
    !% Sample from the halo spin distribution.
    use Galacticus_Error
    implicit none
    class(haloSpinDistributionNbodyErrors), intent(inout)          :: self
    type (treeNode                       ), intent(inout), pointer :: node
    !GCC$ attributes unused :: self, node

    nbodyErrorsSample=0.0d0
    call Galacticus_Error_Report('nbodyErrorsSample','sampling from distribution is not implemented')
    return
  end function nbodyErrorsSample

  double precision function nbodyErrorsDistribution(self,node)
    !% Compute the spin distribution.
    use Galacticus_Nodes
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout)          :: self
    type            (treeNode                       ), intent(inout), pointer :: node
    class           (nodeComponentBasic             )               , pointer :: nodeBasic
    class           (nodeComponentSpin              )               , pointer :: nodeSpin
    double precision                                                          :: mass     , spin , &
         &                                                                       hMass    , hSpin
    integer                                                                      iMass    , iSpin, &
         &                                                                       jMass    , jSpin

    ! Extract the mass and spin of the halo.
    nodeBasic => node     %basic()
    nodeSpin  => node     %spin ()
    mass      =  nodeBasic%mass ()
    spin      =  nodeSpin %spin ()
    ! Ensure the table has sufficient extent.
    call self%tabulate(mass,spin)
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

  double precision function nbodyErrorsDistributionAveraged(self,node,massLimit)
    !% Compute the spin distribution averaged over all halos more massive than the given {\normalfont \ttfamily massLimit}.
    use Galacticus_Nodes
    implicit none
    class           (haloSpinDistributionNbodyErrors), intent(inout)          :: self
    type            (treeNode                       ), intent(inout), pointer :: node
    double precision                                 , intent(in   )          :: massLimit
    class           (nodeComponentBasic             )               , pointer :: nodeBasic
    double precision                                                          :: mass     , massOriginal, &
         &                                                                       massLow  , massHigh    , &
         &                                                                       weight   , weightTotal
    integer                                                                   :: iMass

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
          weight                         =self%haloMassFunction_%integrated(nodeBasic%time(),massLow,massHigh)
          nbodyErrorsDistributionAveraged=nbodyErrorsDistributionAveraged+self%distribution(node)*weight
          weightTotal                    =weightTotal                    +                        weight
       end if
    end do
    ! Normalize by the total weight.
    if (weightTotal > 0.0d0) nbodyErrorsDistributionAveraged=nbodyErrorsDistributionAveraged/weightTotal
    call nodeBasic%massSet(massOriginal)
    return
  end function nbodyErrorsDistributionAveraged
