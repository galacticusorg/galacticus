!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An implementation of virial orbits using the \cite{jiang_orbital_2014} orbital parameter distribution.

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Tables                 , only : table1DLinearLinear
  use :: Virial_Density_Contrast, only : virialDensityContrastFixed

  !# <virialOrbit name="virialOrbitJiang2014">
  !#  <description>
  !#   A virial orbits class which selects orbital parameters randomly from the distribution given by \cite{jiang_orbital_2014},
  !#   including the mass and mass-ratio dependence of the distributions. If the virial density contrast definition differs from
  !#   that used by \cite{jiang_orbital_2014} then the orbit is assigned based on \cite{jiang_orbital_2014}'s definition and then
  !#   propagated to the virial radius relevant to the current definition of density contrast.
  !#  </description>
  !#  <deepCopy>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </stateStorable>
  !# </virialOrbit>
  type, extends(virialOrbitClass) :: virialOrbitJiang2014
     !% A virial orbit class using the \cite{jiang_orbital_2014} orbital parameter distribution.
     private
     class           (darkMatterHaloScaleClass  ), pointer        :: darkMatterHaloScale_    => null()
     class           (cosmologyParametersClass  ), pointer        :: cosmologyParameters_    => null()
     class           (cosmologyFunctionsClass   ), pointer        :: cosmologyFunctions_     => null()
     type            (virialDensityContrastFixed), pointer        :: virialDensityContrast_  => null()
     double precision                            , dimension(3,3) :: B                                , gamma                        , &
          &                                                          sigma                            , mu                           , &
          &                                                          velocityTangentialMean_          , velocityTotalRootMeanSquared_
     type            (table1DLinearLinear       ), dimension(3,3) :: voightDistributions
   contains
     !@ <objectMethods>
     !@   <object>virialOrbitJiang2014</object>
     !@   <objectMethod>
     !@     <method>parametersSelect</method>
     !@     <arguments>\doublezero\ massHost\argin, \doublezero\ massSatellite\argin, \intzero\ i\argout, \intzero\ j\argout</arguments>
     !@     <type>\void</type>
     !@     <description>Select the parameter set to use for this satellite/host pairing.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                    jiang2014Destructor
     procedure :: orbit                           => jiang2014Orbit
     procedure :: densityContrastDefinition       => jiang2014DensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => jiang2014VelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => jiang2014VelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => jiang2014AngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => jiang2014AngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => jiang2014VelocityTotalRootMeanSquared
     procedure :: energyMean                      => jiang2014EnergyMean
     procedure :: parametersSelect                => jiang2014ParametersSelect
  end type virialOrbitJiang2014

  interface virialOrbitJiang2014
     !% Constructors for the {\normalfont \ttfamily jiang2014} virial orbit class.
     module procedure jiang2014ConstructorParameters
     module procedure jiang2014ConstructorInternal
  end interface virialOrbitJiang2014

  ! Module-scope variables used in root finding.
  class           (virialOrbitJiang2014), pointer :: jiang2014Self
  double precision                                :: jiang2014XTotal                       , jiang2014XRadial              , &
       &                                             jiang204ProbabilityRadialNormalization, jiang2014VelocityTotalInternal
  integer                                         :: jiang2014I                            , jiang2014J
  !$omp threadprivate(jiang2014Self,jiang2014XTotal,jiang2014XRadial,jiang204ProbabilityRadialNormalization,jiang2014VelocityTotalInternal,jiang2014I,jiang2014J)

contains

  function jiang2014ConstructorParameters(parameters) result(self)
    !% Generic constructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitJiang2014    )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                          , dimension(3)  :: bRatioLow           , bRatioIntermediate    , bRatioHigh    , &
         &                                                       gammaRatioLow       , gammaRatioIntermediate, gammaRatioHigh, &
         &                                                       sigmaRatioLow       , sigmaRatioIntermediate, sigmaRatioHigh, &
         &                                                       muRatioLow          , muRatioIntermediate   , muRatioHigh

    !# <inputParameter>
    !#   <name>bRatioLow</name>
    !#   <defaultValue>[+0.049d0,+0.548d0,+1.229d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bRatioIntermediate</name>
    !#   <defaultValue>[+1.044d0,+1.535d0,+3.396d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bRatioHigh</name>
    !#   <defaultValue>[+2.878d0,+3.946d0,+2.982d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioLow</name>
    !#   <defaultValue>[+0.109d0,+0.114d0,+0.110d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioIntermediate</name>
    !#   <defaultValue>[+0.098d0,+0.087d0,+0.050d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioHigh</name>
    !#   <defaultValue>[+0.071d0,+0.030d0,-0.012d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioLow</name>
    !#   <defaultValue>[+0.077d0,+0.094d0,+0.072d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioIntermediate</name>
    !#   <defaultValue>[+0.073d0,+0.083d0,+0.118d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioHigh</name>
    !#   <defaultValue>[+0.091d0,+0.139d0,+0.187d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioLow</name>
    !#   <defaultValue>[+1.220d0,+1.231d0,+1.254d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioIntermediate</name>
    !#   <defaultValue>[+1.181d0,+1.201d0,+1.236d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioHigh</name>
    !#   <defaultValue>[+1.100d0,+1.100d0,+1.084d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=virialOrbitJiang2014(bRatioLow,bRatioIntermediate,bRatioHigh,gammaRatioLow,gammaRatioIntermediate,gammaRatioHigh,sigmaRatioLow,sigmaRatioIntermediate,sigmaRatioHigh,muRatioLow,muRatioIntermediate,muRatioHigh,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    return
  end function jiang2014ConstructorParameters

  function jiang2014ConstructorInternal(bRatioLow,bRatioIntermediate,bRatioHigh,gammaRatioLow,gammaRatioIntermediate,gammaRatioHigh,sigmaRatioLow,sigmaRatioIntermediate,sigmaRatioHigh,muRatioLow,muRatioIntermediate,muRatioHigh,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    use :: Numerical_Integration   , only : integrator                  , GSL_Integ_Gauss61
    use :: Statistics_Distributions, only : distributionFunction1DVoight
    use :: Virial_Density_Contrast , only : fixedDensityTypeCritical
    implicit none
    type            (virialOrbitJiang2014        )                                :: self
    double precision                              , intent(in   ), dimension(3  ) :: bRatioLow                          , bRatioIntermediate            , bRatioHigh                          , &
         &                                                                           gammaRatioLow                      , gammaRatioIntermediate        , gammaRatioHigh                      , &
         &                                                                           sigmaRatioLow                      , sigmaRatioIntermediate        , sigmaRatioHigh                      , &
         &                                                                           muRatioLow                         , muRatioIntermediate           , muRatioHigh
    class           (darkMatterHaloScaleClass    ), intent(in   ), target         :: darkMatterHaloScale_
    class           (cosmologyParametersClass    ), intent(in   ), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass     ), intent(in   ), target         :: cosmologyFunctions_
    integer                                       , parameter                     :: tableCount                 =1000
    integer                                                                       :: i                                  , j                             , k
    type            (distributionFunction1DVoight)                                :: voightDistribution
    logical                                       , save        , dimension(3,3)  :: previousInitialized        =.false.
    double precision                              , save        , dimension(3,3)  :: previousB                          , previousGamma                 , previousSigma                       , &
         &                                                                           previousMu                         , previousVelocityTangentialMean, previousVelocityTotalRootMeanSquared
    type            (table1DLinearLinear         ), save        , dimension(3,3)  :: previousVoightDistributions
    double precision                                                              :: limitLower                         , limitUpper                    , halfWidthHalfMaximum                , &
         &                                                                           fullWidthHalfMaximumLorentzian     , fullWidthHalfMaximumGaussian
    logical                                                                       :: reUse                              , limitFound
    type            (integrator                  )                                :: integratorTangential               , integratorTotal
    !# <constructorAssign variables="*darkMatterHaloScale_, *cosmologyParameters_, *cosmologyFunctions_"/>

    ! Assign parameters of the distribution.
    self%B    (:,1)=    bRatioLow
    self%B    (:,2)=    bRatioIntermediate
    self%B    (:,3)=    bRatioHigh
    self%gamma(:,1)=gammaRatioLow
    self%gamma(:,2)=gammaRatioIntermediate
    self%gamma(:,3)=gammaRatioHigh
    self%sigma(:,1)=sigmaRatioLow
    self%sigma(:,2)=sigmaRatioIntermediate
    self%sigma(:,3)=sigmaRatioHigh
    self%mu   (:,1)=   muRatioLow
    self%mu   (:,2)=   muRatioIntermediate
    self%mu   (:,3)=   muRatioHigh
    ! Tabulate Voight distribution functions for speed.
    integratorTangential=integrator(jiang2014DistributionVelocityTangential  ,toleranceRelative=1.0d-6,integrationRule=GSL_Integ_Gauss61)
    integratorTotal     =integrator(jiang2014DistributionVelocityTotalSquared,toleranceRelative=1.0d-6,integrationRule=GSL_Integ_Gauss61)
    ! Build the distribution function for total velocity.
    do i=1,3
       do j=1,3
          ! Check if this distribution was built previously.
          !$omp critical(virialOrbitJiang2014ReUse)
          reuse  =                                             &
               &                      previousInitialized(i,j) &
               &  .and.                                        &
               &   self%B    (i,j) == previousB          (i,j) &
               &  .and.                                        &
               &   self%gamma(i,j) == previousGamma      (i,j) &
               &  .and.                                        &
               &   self%sigma(i,j) == previousSigma      (i,j) &
               &  .and.                                        &
               &   self%mu   (i,j) == previousMu         (i,j)
          if (reuse) then
             self%voightDistributions          (i,j)=previousVoightDistributions         (i,j)
             self%velocityTangentialMean_      (i,j)=previousVelocityTangentialMean      (i,j)
             self%velocityTotalRootMeanSquared_(i,j)=previousVelocityTotalRootMeanSquared(i,j)
          end if
          !$omp end critical(virialOrbitJiang2014ReUse)
          ! Build the distribution.
          if (.not.reuse) then
             ! Set the lower and upper limit of the distribution to +/-5 times the half-width at half-maximum below/above the mean
             ! (limited also to 0). This avoids attempting to evaluate the distribution far from the mean (where it is small, but
             ! the numerical evaluation of the hypergeometric function used in the CDF is unstable). The half-width at
             ! half-maximum is estimated using the approximation of Olivero (1977; Journal of Quantitative Spectroscopy and
             ! Radiative Transfer; 17; 233; http://adsabs.harvard.edu/abs/1977JQSRT..17..233O).
             fullWidthHalfMaximumLorentzian=+2.0d0*self%gamma(i,j)
             fullWidthHalfMaximumGaussian  =+2.0d0*self%sigma(i,j)*sqrt(2.0d0*log(2.0d0))
             halfWidthHalfMaximum          =+0.5d0                                     &
                  &                         *(                                         &
                  &                           +      0.5346d0                          &
                  &                           *      fullWidthHalfMaximumLorentzian    &
                  &                           +sqrt(                                   &
                  &                                 +0.2166d0                          &
                  &                                 *fullWidthHalfMaximumLorentzian**2 &
                  &                                 +fullWidthHalfMaximumGaussian  **2 &
                  &                                )                                   &
                  &                          )
             limitLower                    =max(self%mu(i,j)-5.0d0*halfWidthHalfMaximum,0.0d0)
             limitUpper                    =    self%mu(i,j)+5.0d0*halfWidthHalfMaximum
             voightDistribution=distributionFunction1DVoight(                        &
                  &                                          self%gamma(i,j)       , &
                  &                                          self%mu   (i,j)       , &
                  &                                          self%sigma(i,j)       , &
                  &                                          limitLower            , &
                  &                                          limitUpper              &
                  &                                         )
             ! Tabulate the cumulative distribution.
             call self%voightDistributions(i,j)%create(limitLower,limitUpper,tableCount)
             !$omp parallel do
             do k=2,tableCount-1
                call self%voightDistributions(i,j)%populate(min(1.0d0,max(0.0d0,voightDistribution%cumulative(self%voightDistributions(i,j)%x(k)))),k)
             end do
             !$omp end parallel do
             ! Ensure tabulation starts at 0 and reaches 1.
             call self%voightDistributions(i,j)%populate(0.0d0,         1)
             call self%voightDistributions(i,j)%populate(1.0d0,tableCount)
             ! Ensure that 0 and 1 are unique within the table, by making any duplicated table entries slightly beyond these
             ! values. This ensures that when seeking values in this table we do not exceed the plausible range.
             limitFound=.false.
             do k=1,tableCount
                if (limitFound)  call self%voightDistributions(i,j)%populate(1.0d0+1.0d-6,k)
                if (self%voightDistributions(i,j)%y(k) >= 1.0d0) limitFound=.true.
             end do
             limitFound=.false.
             do k=tableCount,1,-1
                if (limitFound) call self%voightDistributions(i,j)%populate(0.0d0-1.0d-6,k)
                if (self%voightDistributions(i,j)%y(k) <= 0.0d0) limitFound=.true.
             end do
             ! Compute the mean magnitude of tangential velocity.
             self%velocityTangentialMean_      (i,j)=     integratorTangential%integrate(limitLower,limitUpper)
             ! Compute the root mean squared total velocity.
             self%velocityTotalRootMeanSquared_(i,j)=sqrt(integratorTotal     %integrate(limitLower,limitUpper))
             ! Store this table for later reuse.
             !$omp critical(virialOrbitJiang2014ReUse)
             previousInitialized(i,j)=.true.
             previousB                           (i,j)=self%B                            (i,j)
             previousGamma                       (i,j)=self%gamma                        (i,j)
             previousSigma                       (i,j)=self%sigma                        (i,j)
             previousMu                          (i,j)=self%mu                           (i,j)
             previousVoightDistributions         (i,j)=self%voightDistributions          (i,j)
             previousVelocityTangentialMean      (i,j)=self%velocityTangentialMean_      (i,j)
             previousVelocityTotalRootMeanSquared(i,j)=self%velocityTotalRootMeanSquared_(i,j)
             !$omp end critical(virialOrbitJiang2014ReUse)
          end if
       end do
    end do
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrast_)
    !# <referenceConstruct isResult="yes" owner="self" object="virialDensityContrast_" constructor="virialDensityContrastFixed(200.0d0,fixedDensityTypeCritical,2.0d0,self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    return

  contains

    double precision function jiang2014DistributionVelocityTangential(velocityTotal)
      !% The integrand used to average the tangential velocity over the distribution function for total velocity.
      use :: Bessel_Functions        , only : Bessel_Function_I1
      use :: Numerical_Constants_Math, only : Pi
      use :: Struve_Functions        , only : Struve_Function_L1
      implicit none
      double precision, intent(in   ) :: velocityTotal

      jiang2014DistributionVelocityTangential=+voightDistribution%density(velocityTotal) &
           &                                  *Pi                                        &
           &                                  *                           velocityTotal  &
           &                                  *(                                         &
           &                                    +                         self%B(i,j)    &
           &                                    -2.0d0*Bessel_Function_I1(self%B(i,j))   &
           &                                    -2.0d0*Struve_Function_L1(self%B(i,j))   &
           &                                   )                                         &
           &                                  /4.0d0                                     &
           &                                  /(                                         &
           &                                    +1.0d0                                   &
           &                                    +                         self%B(i,j)    &
           &                                    -      exp               (self%B(i,j))   &
           &                                  )
      return
    end function jiang2014DistributionVelocityTangential

    double precision function jiang2014DistributionVelocityTotalSquared(velocityTotal)
      !% The integrand used to average the squared total velocity over the distribution function for total velocity.
      implicit none
      double precision, intent(in   ) :: velocityTotal

      jiang2014DistributionVelocityTotalSquared=+voightDistribution%density(velocityTotal)    &
           &                                    *                           velocityTotal **2
      return
    end function jiang2014DistributionVelocityTotalSquared

  end function jiang2014ConstructorInternal

  subroutine jiang2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    implicit none
    type(virialOrbitJiang2014), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"  />
    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%virialDensityContrast_"/>
    return
  end subroutine jiang2014Destructor

  function jiang2014Orbit(self,node,host,acceptUnboundOrbits)
    !% Return jiang2014 orbital parameters for a satellite.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Root_Finder                         , only : rangeExpandMultiplicative          , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    type            (keplerOrbit               )                        :: jiang2014Orbit
    class           (virialOrbitJiang2014      ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                   , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: hostBasic              , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrast_
    integer                                     , parameter             :: attemptsMaximum        =10000
    double precision                            , parameter             :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                            , parameter             :: toleranceAbsolute      =0.0d+0
    double precision                            , parameter             :: toleranceRelative      =1.0d-3
    type            (rootFinder                )                        :: totalFinder                   , radialFinder
    double precision                                                    :: velocityHost                  , radiusHost                , &
         &                                                                 massHost                      , massSatellite             , &
         &                                                                 energyInternal                , radiusHostSelf            , &
         &                                                                 velocityRadialInternal        , velocityTangentialInternal
    logical                                                             :: foundOrbit
    integer                                                             :: attempts

    ! Get basic components.
    basic     => node%basic()
    hostBasic => host%basic()
    ! Find virial density contrast under Jiang et al. (2014) definition.
    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>    
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Jiang et al. (2014) virial density contrast
    ! definition. We limit the satellite mass to be less than the host mass here. This is necessary because the correction to the
    ! mass under the definition of Jiang et al. (2014) can sometimes lead to large satellite masses if the satellite is moving to
    ! a new host (and was therefore last isolated at a much earlier time and so can be much denser than the host).
    massHost     =             Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHostSelf,velocityHost)
    massSatellite=min(massHost,Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                            ))
    !# <objectDestructor name="virialDensityContrast_"/>
    ! Select parameters appropriate for this host-satellite pair.
    call self%parametersSelect(massHost,massSatellite,jiang2014I,jiang2014J)
    ! Configure finder objects.
    if (.not.totalFinder %isInitialized()) then
       call totalFinder %rootFunction(jiang2014TotalVelocityCDF          )
       call totalFinder %tolerance   (toleranceAbsolute,toleranceRelative)
       call totalFinder %rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    if (.not.radialFinder%isInitialized()) then
       call radialFinder%rootFunction(jiang2014RadialVelocityCDF         )
       call radialFinder%tolerance   (toleranceAbsolute,toleranceRelative)
       call radialFinder%rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    ! Select an orbit.
    foundOrbit    =  .false.
    attempts      =  0
    jiang2014Self => self
    do while (.not.foundOrbit .and. attempts < attemptsMaximum)
       ! Increment number of attempts.
       attempts=attempts+1
       ! Reset the orbit.
       call jiang2014Orbit%reset()
       ! Set basic properties of the orbit.
       call jiang2014Orbit%massesSet(massSatellite,massHost      )
       call jiang2014Orbit%radiusSet(              radiusHostSelf)
       ! Solve for the total velocity.
       jiang2014XTotal               =node%hostTree%randomNumberGenerator_%uniformSample()
       jiang2014VelocityTotalInternal=totalFinder %find(rootGuess=1.0d0)
       ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
       ! bound that later rounding errors will not make it appear unbound.
       foundOrbit=.true.
       if (.not.acceptUnboundOrbits) then
          energyInternal=-1.0d0+0.5d0*jiang2014VelocityTotalInternal**2*jiang2014Orbit%specificReducedMass()
          foundOrbit=(energyInternal < -boundTolerance)
       end if
       if (.not.foundOrbit) cycle
       ! Solve for the radial velocity.
       jiang2014XRadial                      =+0.0d0
       jiang204ProbabilityRadialNormalization=+1.0d0
       jiang204ProbabilityRadialNormalization=+1.0d0                                                        &
            &                                 /(                                                            &
            &                                   +jiang2014RadialVelocityCDF(jiang2014VelocityTotalInternal) &
            &                                   -jiang2014RadialVelocityCDF(0.0d0                         ) &
            &                                  )
       jiang2014XRadial                      =node%hostTree%randomNumberGenerator_%uniformSample()
       velocityRadialInternal                =radialFinder%find(rootGuess=sqrt(2.0d0)*jiang2014VelocityTotalInternal)
       ! Compute tangential velocity.
       velocityTangentialInternal=sqrt(max(0.0d0,jiang2014VelocityTotalInternal**2-velocityRadialInternal**2))
       call jiang2014Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call jiang2014Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%virialRadius(host)
       foundOrbit=.false.
       if (jiang2014Orbit%radiusApocenter() >= radiusHost .and. jiang2014Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call jiang2014Orbit%propagate(radiusHost  ,infalling=.true.)
          call jiang2014Orbit%massesSet(basic%mass(),hostBasic%mass())
       end if
    end do
    ! If too many iterations were required to find an orbit, abort.
    if (attempts >= attemptsMaximum) call Galacticus_Error_Report('maximum number of attempts exceeded'//{introspection:location})
    return
  end function jiang2014Orbit

  double precision function jiang2014TotalVelocityCDF(velocityTotal)
    !% Cumulative distribution function for the total velocity.
    implicit none
    double precision, intent(in   ) :: velocityTotal

    jiang2014TotalVelocityCDF=jiang2014Self%voightDistributions(jiang2014I,jiang2014J)%interpolate(velocityTotal)-jiang2014XTotal
    return
  end function jiang2014TotalVelocityCDF

  double precision function jiang2014RadialVelocityCDF(velocityRadial)
    !% Cumulative distribution function for the radial velocity.
    implicit none
    double precision, intent(in   ) :: velocityRadial

    jiang2014RadialVelocityCDF=+jiang204ProbabilityRadialNormalization                         &
         &                     *(                                                              &
         &                       +       jiang2014VelocityTotalInternal                        &
         &                       /       jiang2014Self%B               (jiang2014I,jiang2014J) &
         &                       *(                                                            &
         &                         +exp(                                                       &
         &                              +jiang2014Self%B               (jiang2014I,jiang2014J) &
         &                              *         velocityRadial                               &
         &                              /jiang2014VelocityTotalInternal                        &
         &                             )                                                       &
         &                         -1.0d0                                                      &
         &                        )                                                            &
         &                       -velocityRadial                                               &
         &                      )                                                              &
         &                     -jiang2014XRadial
    return
  end function jiang2014RadialVelocityCDF

  function jiang2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{jiang_orbital_2014} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: jiang2014DensityContrastDefinition
    class(virialOrbitJiang2014      ), intent(inout) :: self

    jiang2014DensityContrastDefinition => self%virialDensityContrast_
    return
  end function jiang2014DensityContrastDefinition

  double precision function jiang2014VelocityTangentialMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the tangential velocity.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitJiang2014      ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                  , host
    class           (nodeComponentBasic        ), pointer       :: hostBasic             , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: massHost              , radiusHost   , &
         &                                                         velocityHost          , massSatellite
    integer                                                     :: i                     , j

    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>
    basic         => node%basic()
    hostBasic     => host%basic()
    massHost      =  Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    massSatellite =  Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                        )
    !# <objectDestructor name="virialDensityContrast_"/>
    call self%parametersSelect(massHost,massSatellite,i,j)
    jiang2014VelocityTangentialMagnitudeMean=+self%velocityTangentialMean_(i,j) &
         &                                   *velocityHost
    return
  end function jiang2014VelocityTangentialMagnitudeMean

  function jiang2014VelocityTangentialVectorMean(self,node,host)
    !% Return the mean of the vector tangential velocity.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                      , dimension(3)  :: jiang2014VelocityTangentialVectorMean
    class           (virialOrbitJiang2014), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node                                  , host
    !$GLC attributes unused :: self, node, host

    jiang2014VelocityTangentialVectorMean=0.0d0
    call Galacticus_Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function jiang2014VelocityTangentialVectorMean

  double precision function jiang2014AngularMomentumMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the angular momentum.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitJiang2014), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node        , host
    class           (nodeComponentBasic  ), pointer       :: basic       , hostBasic
    double precision                                      :: massHost    , radiusHost, &
         &                                                   velocityHost

    basic                                 =>  node%basic()
    hostBasic                             =>  host%basic()
    massHost                              =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    jiang2014AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                   *radiusHost                                      &
         &                                   /(                                               & ! Account for reduced mass.
         &                                     +1.0d0                                         &
         &                                     +basic    %mass()                              &
         &                                     /hostBasic%mass()                              &
         &                                    )
    return
  end function jiang2014AngularMomentumMagnitudeMean

  function jiang2014AngularMomentumVectorMean(self,node,host)
    !% Return the mean of the vector angular momentum.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                      , dimension(3)  :: jiang2014AngularMomentumVectorMean
    class           (virialOrbitJiang2014), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node                               , host
    !$GLC attributes unused :: self, node, host

    jiang2014AngularMomentumVectorMean=0.0d0
    call Galacticus_Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function jiang2014AngularMomentumVectorMean

  double precision function jiang2014VelocityTotalRootMeanSquared(self,node,host)
    !% Return the root mean squared total velocity.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitJiang2014      ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                  , host
    class           (nodeComponentBasic        ), pointer       :: hostBasic             , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: massHost              , radiusHost   , &
         &                                                         velocityHost          , massSatellite
    integer                                                     :: i                     , j

    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>
    basic         => node%basic()
    hostBasic     => host%basic()
    massHost      =  Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    massSatellite =  Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                        )
    !# <objectDestructor name="virialDensityContrast_"/>
    call self%parametersSelect(massHost,massSatellite,i,j)
    jiang2014VelocityTotalRootMeanSquared=+self%velocityTotalRootMeanSquared_(i,j) &
         &                                *velocityHost
    return
  end function jiang2014VelocityTotalRootMeanSquared

  double precision function jiang2014EnergyMean(self,node,host)
    !% Return the mean energy of the orbits.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstantGalacticus
    implicit none
    class           (virialOrbitJiang2014), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node        , host
    class           (nodeComponentBasic  ), pointer       :: basic       , hostBasic
    double precision                                      :: massHost    , radiusHost, &
         &                                                   velocityHost

    basic               =>  node%basic()
    hostBasic           =>  host%basic()
    massHost            =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    jiang2014EnergyMean =  +0.5d0                                           &
         &                 *self%velocityTotalRootMeanSquared(node,host)**2 &
         &                 /(                                               & ! Account for reduced mass.
         &                   +1.0d0                                         &
         &                   +basic    %mass()                              &
         &                   /hostBasic%mass()                              &
         &                  )                                               &
         &                 -gravitationalConstantGalacticus                 &
         &                 *massHost                                        &
         &                 /radiusHost
    return
  end function jiang2014EnergyMean

  subroutine jiang2014ParametersSelect(self,massHost,massSatellite,i,j)
    !% Select the parameter set to use for this satellite/host pairing.
    implicit none
    class           (virialOrbitJiang2014), intent(inout) :: self
    double precision                      , intent(in   ) :: massHost, massSatellite
    integer                               , intent(  out) :: i       , j
    !$GLC attributes unused :: self

    if      (              massHost < 10.0d0**12.5d0) then
       i=1
    else if (              massHost < 10.0d0**13.5d0) then
       i=2
    else
       i=3
    end if
    if      (massSatellite/massHost < 0.005d0       ) then
       j=1
    else if (massSatellite/massHost < 0.050d0       ) then
       j=2
    else
       j=3
    end if
    return
  end subroutine jiang2014ParametersSelect
