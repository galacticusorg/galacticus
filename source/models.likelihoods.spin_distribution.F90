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
  Implementation of a posterior sampling likelihood class which implements a likelihood for halo spin distributions.
  !!}

  use :: Cosmology_Functions              , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales          , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales       , only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO         , only : darkMatterProfileDMOClass
  use :: Halo_Mass_Functions              , only : haloMassFunctionClass
  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <enumeration>
   <name>spinDistributionType</name>
   <description>Used to specify the type of intrinsic spin distribution.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="logNormal"/>
   <entry label="bett2007" />
  </enumeration>
  !!]

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodSpinDistribution">
   <description>A posterior sampling likelihood class which implements a likelihood for halo spin distributions.</description>
   <runTimeFileDependencies paths="fileName"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodSpinDistribution
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for fitting dark matter halo spin distributions.
     !!}
     private
     class           (cosmologyFunctionsClass            ), pointer                     :: cosmologyFunctions_           => null()
     class           (haloMassFunctionClass              ), pointer                     :: haloMassFunction_             => null()
     class           (nbodyHaloMassErrorClass            ), pointer                     :: nbodyHaloMassError_           => null()
     class           (darkMatterHaloScaleClass           ), pointer                     :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileScaleRadiusClass  ), pointer                     :: darkMatterProfileScaleRadius_ => null()
     double precision                                     , dimension(:  ), allocatable :: spin                                   , distribution                      , &
          &                                                                                spinMinimum                            , spinMaximum                       , &
          &                                                                                distributionError
     double precision                                                                   :: time                                   , massParticle                      , &
          &                                                                                massHaloMinimum                        , energyEstimateParticleCountMaximum, &
          &                                                                                redshift                               , logNormalRange
     integer                                                                            :: particleCountMinimum
     type            (enumerationSpinDistributionTypeType)                              :: distributionType
     type            (varying_string                     )                              :: fileName
   contains
     final     ::                    spinDistributionDestructor
     procedure :: evaluate        => spinDistributionEvaluate
     procedure :: functionChanged => spinDistributionFunctionChanged
  end type posteriorSampleLikelihoodSpinDistribution

  interface posteriorSampleLikelihoodSpinDistribution
     !!{
     Constructors for the {\normalfont \ttfamily spinDistribution} posterior sampling convergence class.
     !!}
     module procedure spinDistributionConstructorParameters
     module procedure spinDistributionConstructorInternal
  end interface posteriorSampleLikelihoodSpinDistribution

contains

  function spinDistributionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily spinDistribution} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodSpinDistribution)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class           (haloMassFunctionClass                    ), pointer       :: haloMassFunction_
    class           (nbodyHaloMassErrorClass                  ), pointer       :: nbodyHaloMassError_
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass        ), pointer       :: darkMatterProfileScaleRadius_
    type            (varying_string                           )                :: fileName                     , distributionType
    double precision                                                           :: redshift                     , massHaloMinimum                   , &
         &                                                                        massParticle                 , energyEstimateParticleCountMaximum, &
         &                                                                        logNormalRange
    integer                                                                    :: particleCountMinimum

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the target spin distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>distributionType</name>
      <description>The name of the spin distribution to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift at which to compute the spin distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massParticle</name>
      <description>The mass of a particle in the N-body simulation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMinimum</name>
      <description>The minimum halo mass over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>particleCountMinimum</name>
      <description>The minimum particle count used in N-body halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>energyEstimateParticleCountMaximum</name>
      <description>The maximum number of N-body particles used in estimating halo energies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logNormalRange</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <defaultSource>A large range which will include (almost) the entirety of the distribution.</defaultSource>
      <description>The multiplicative range of the log-normal distribution used to model the distribution of the mass and energy terms in the spin parameter. Specifically, the lognormal distribution is truncated outside the range $(\lambda_\mathrm{m}/R,\lambda_\mathrm{m} R$, where $\lambda_\mathrm{m}$ is the measured spin, and $R=${\normalfont \ttfamily [logNormalRange]}</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodSpinDistribution(char(fileName),enumerationSpinDistributionTypeEncode(char(distributionType),includesPrefix=.false.),redshift,logNormalRange,massHaloMinimum,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,cosmologyFunctions_,haloMassFunction_,nbodyHaloMassError_,darkMatterHaloScale_,darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="nbodyHaloMassError_"          />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function spinDistributionConstructorParameters

  function spinDistributionConstructorInternal(fileName,distributionType,redshift,logNormalRange,massHaloMinimum,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,cosmologyFunctions_,haloMassFunction_,nbodyHaloMassError_,darkMatterHaloScale_,darkMatterProfileScaleRadius_) result(self)
    !!{
    Constructor for ``spinDistribution'' posterior sampling likelihood class.
    !!}
    use :: HDF5_Access      , only : hdf5Access
    use :: IO_HDF5          , only : hdf5Object
    implicit none
    type            (posteriorSampleLikelihoodSpinDistribution)                        :: self
    character       (len=*                                    ), intent(in   )         :: fileName
    double precision                                           , intent(in   )         :: redshift                     , massHaloMinimum                   , &
         &                                                                                massParticle                 , energyEstimateParticleCountMaximum, &
         &                                                                                logNormalRange
    integer                                                    , intent(in   )         :: particleCountMinimum
    type            (enumerationSpinDistributionTypeType      ), intent(in   )         :: distributionType
    class           (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class           (haloMassFunctionClass                    ), intent(in   ), target :: haloMassFunction_
    class           (nbodyHaloMassErrorClass                  ), intent(in   ), target :: nbodyHaloMassError_
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass        ), intent(in   ), target :: darkMatterProfileScaleRadius_
    type            (hdf5Object                               )                        :: spinDistributionFile
    double precision                                                                   :: spinIntervalLogarithmic
    integer                                                                            :: i
    !![
    <constructorAssign variables="fileName, distributionType, redshift, logNormalRange, massHaloMinimum, massParticle, particleCountMinimum, energyEstimateParticleCountMaximum, *cosmologyFunctions_, *haloMassFunction_, *nbodyHaloMassError_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_"/>
    !!]

    ! Convert redshift to time.
    self%time=self%cosmologyFunctions_ %cosmicTime                 (          &
         &     self%cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                           redshift &
         &                                                          )         &
         &                                                         )
    ! Read the target spin distribution from file.
    !$ call hdf5Access%set()
    call spinDistributionFile%openFile   (trim(fileName),readOnly=.true.)
    call spinDistributionFile%readDataset("spinParameter"    ,self%spin             )
    call spinDistributionFile%readDataset("distribution"     ,self%distribution     )
    call spinDistributionFile%readDataset("distributionError",self%distributionError)
    call spinDistributionFile%close()
    !$ call hdf5Access%unset()
    ! Compute spin ranges for bins.
    spinIntervalLogarithmic=+log(                            &
         &                       +self%spin(size(self%spin)) &
         &                       /self%spin(              1) &
         &                      )                            &
         &                  /dble(                           &
         &                        +size(self%spin)           &
         &                        -1                         &
         &                       )
    allocate(self%spinMinimum,mold=self%spin)
    allocate(self%spinMaximum,mold=self%spin)
    do i=1,size(self%spin)
       self%spinMinimum(i)=self%spin(i)*exp(-0.5d0*spinIntervalLogarithmic)
       self%spinMaximum(i)=self%spin(i)*exp(+0.5d0*spinIntervalLogarithmic)
    end do
    return
  end function spinDistributionConstructorInternal

  subroutine spinDistributionDestructor(self)
    !!{
    Destructor for ``spinDistribution'' posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodSpinDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%haloMassFunction_"            />
    <objectDestructor name="self%nbodyHaloMassError_"          />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine spinDistributionDestructor

  double precision function spinDistributionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo spin distribution likelihood function.
    !!}
    use :: Error                         , only : Error_Report
    use :: Galacticus_Nodes              , only : nodeComponentBasic             , nodeComponentSpin            , treeNode
    use :: Halo_Spin_Distributions       , only : haloSpinDistributionBett2007   , haloSpinDistributionLogNormal, haloSpinDistributionNbodyErrors
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Numerical_Integration         , only : integrator
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodSpinDistribution), intent(inout), target       :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(inout), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature           , logLikelihoodCurrent    , &
         &                                                                                      logPriorCurrent       , logPriorProposed
    real                                                       , intent(inout)               :: timeEvaluate
    double precision                                           , intent(  out), optional     :: logLikelihoodVariance
    logical                                                    , intent(inout), optional     :: forceAcceptance
    double precision                                           , allocatable  , dimension(:) :: stateVector           , distribution
    type            (treeNode                                 ), pointer                     :: node
    class           (nodeComponentBasic                       ), pointer                     :: nodeBasic
    class           (nodeComponentSpin                        ), pointer                     :: nodeSpin
    class           (haloSpinDistributionLogNormal            ), pointer                     :: distributionLogNormal
    class           (haloSpinDistributionBett2007             ), pointer                     :: distributionBett2007
    type            (haloSpinDistributionNbodyErrors          )                              :: distributionNbody
    type            (integrator                               )                              :: integrator_
    integer                                                                                  :: i
    !$GLC attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       spinDistributionEvaluate=0.0d0
       return
    end if
    ! Extract parameters of the spin distribution.
    stateVector=simulationState%get()
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    ! Build the halo spin distribution object.
    select case (self%distributionType%ID)
    case (spinDistributionTypeLogNormal%ID)
       if (size(stateVector) /= 2)                                                                  &
            & call Error_Report(                                                                    &
            &                   '2 parameters are required for the "lognormal" spin distribution'// &
            &                   {introspection:location}                                            &
            &                  )
       allocate(haloSpinDistributionLogNormal :: distributionLogNormal)
       select type (distributionLogNormal)
       type is (haloSpinDistributionLogNormal)
          distributionLogNormal=haloSpinDistributionLogNormal  (                                         &
               &                                                stateVector(1)                         , &
               &                                                stateVector(2)                         , &
               &                                                self%darkMatterHaloScale_                &
               &                                               )
          distributionNbody    =haloSpinDistributionNbodyErrors(                                         &
               &                                                distributionLogNormal                  , &
               &                                                self%massParticle                      , &
               &                                                self%particleCountMinimum              , &
               &                                                self%energyEstimateParticleCountMaximum, &
               &                                                self%logNormalRange                    , &
               &                                                self%time                              , &
               &                                                self%nbodyHaloMassError_               , &
               &                                                self%cosmologyFunctions_               , &
               &                                                self%haloMassFunction_                 , &
               &                                                self%darkMatterHaloScale_              , &
               &                                                self%darkMatterProfileScaleRadius_       &
               &                                               )
       end select
    case (spinDistributionTypeBett2007%ID)
       if (size(stateVector) /= 2)                                                                 &
            & call Error_Report(                                                                   &
            &                   '2 parameters are required for the "bett2007" spin distribution'// &
            &                   {introspection:location}                                           &
            &                  )
       allocate(haloSpinDistributionBett2007 :: distributionBett2007)
       select type (distributionBett2007)
       type is (haloSpinDistributionBett2007)
          distributionBett2007=haloSpinDistributionBett2007    (                                         &
               &                                                stateVector(1)                         , &
               &                                                stateVector(2)                         , &
               &                                                self%darkMatterHaloScale_                &
               &                                               )
          distributionNbody    =haloSpinDistributionNbodyErrors(                                         &
               &                                                distributionBett2007                   , &
               &                                                self%massParticle                      , &
               &                                                self%particleCountMinimum              , &
               &                                                self%energyEstimateParticleCountMaximum, &
               &                                                self%logNormalRange                    , &
               &                                                self%time                              , &
               &                                                self%nbodyHaloMassError_               , &
               &                                                self%cosmologyFunctions_               , &
               &                                                self%haloMassFunction_                 , &
               &                                                self%darkMatterHaloScale_              , &
               &                                                self%darkMatterProfileScaleRadius_       &
               &                                               )
       end select
    end select
    ! Nullify objects to avoid them being destroyed.
    nullify(distributionLogNormal)
    ! Create a work node and set the appropriate cosmological time.
    node      => treeNode      (                 )
    nodeBasic => node    %basic(autoCreate=.true.)
    nodeSpin  => node    %spin (autoCreate=.true.)
    call nodeBasic%timeSet(self%time)
    ! Compute the spin distribution. The spin distribution is averaged over the width of each bin.
    allocate(distribution(size(self%spin)))
    integrator_=integrator(spinDistributionIntegrate,toleranceRelative=1.0d-6)
    do i=1,size(self%spin)
       distribution(i)=+integrator_%integrate( self%spinMinimum(i),self%spinMaximum(i)) &
            &          /log10                (+self%spinMaximum(i)/self%spinMinimum(i))
    end do
    ! Evaluate the log-likelihood.
    spinDistributionEvaluate=-0.5d0*sum(((distribution-self%distribution)/self%distributionError)**2,mask=self%distributionError > 0.0d0)
    ! Clean up.
    call node%destroy()
    deallocate(node        )
    deallocate(stateVector )
    deallocate(distribution)
    return

  contains

    double precision function spinDistributionIntegrate(spinPrime)
      !!{
      Integrand function used to find cumulative spin distribution over a bin.
      !!}
      use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
      implicit none
      double precision, intent(in   ) :: spinPrime

      call nodeSpin%angularMomentumSet(spinPrime*Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_))
      spinDistributionIntegrate=distributionNbody%distributionAveraged(node,self%massHaloMinimum)
      return
    end function spinDistributionIntegrate

  end function spinDistributionEvaluate

  subroutine spinDistributionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodSpinDistribution), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine spinDistributionFunctionChanged
