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
  Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodHaloMassFunction">
   <description>A posterior sampling likelihood class which implements a likelihood for halo mass functions.</description>
   <runTimeFileDependencies paths="fileNames"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodBaseParameters) :: posteriorSampleLikelihoodHaloMassFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer                       :: cosmologyFunctions_                => null()
     class           (criticalOverdensityClass     ), pointer                       :: criticalOverdensity_               => null()
     class           (cosmologicalMassVarianceClass), pointer                       :: cosmologicalMassVariance_          => null()
     class           (linearGrowthClass            ), pointer                       :: linearGrowth_                      => null()
     class           (randomNumberGeneratorClass   ), pointer                       :: randomNumberGenerator_             => null()
     double precision                               , dimension(:    ), allocatable :: massMinimum                                 , massMaximum              , &
          &                                                                            mass                                        , countConversionFactor
     double precision                               , dimension(:,:,:), allocatable :: covarianceMatrix
     double precision                               , dimension(:,:  ), allocatable :: massFunction
     integer         (c_size_t                     ), dimension(:,:  ), allocatable :: countHalos
     double precision                                                               :: varianceFractionalModelDiscrepancy          , massParticle             , &
          &                                                                            massRangeMinimum                            , massRangeMaximum
     logical                                                                        :: likelihoodPoisson                           , includeDiscrepancyChecked, &
          &                                                                            report                                      , binAverage               , &
          &                                                                            includeCorrelations
     integer                                                                        :: binCountMinimum                             , indexDiscrepancy
     type            (matrix                       ), dimension(:    ), allocatable :: covariance
     type            (varying_string               ), dimension(:    ), allocatable :: fileNames
     double precision                               , dimension(:    ), allocatable :: redshifts                                   , times
   contains
     final     ::                    haloMassFunctionDestructor
     procedure :: evaluate        => haloMassFunctionEvaluate
     procedure :: functionChanged => haloMassFunctionFunctionChanged
  end type posteriorSampleLikelihoodHaloMassFunction

  interface posteriorSampleLikelihoodHaloMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily haloMassFunction} posterior sampling convergence class.
     !!}
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface posteriorSampleLikelihoodHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloMassFunction} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodHaloMassFunction)                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    type            (inputParameters                          ), pointer                     :: parametersModel
    class           (cosmologyFunctionsClass                  ), pointer                     :: cosmologyFunctions_
    class           (criticalOverdensityClass                 ), pointer                     :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass            ), pointer                     :: cosmologicalMassVariance_
    class           (linearGrowthClass                        ), pointer                     :: linearGrowth_
    class           (randomNumberGeneratorClass               ), pointer                     :: randomNumberGenerator_
    type            (varying_string                           ), allocatable  , dimension(:) :: changeParametersFileNames         , fileNames
    double precision                                           , allocatable  , dimension(:) :: redshifts
    type            (varying_string                           )                              :: baseParametersFileName
    double precision                                                                         :: massRangeMinimum                  , massRangeMaximum   , &
         &                                                                                      varianceFractionalModelDiscrepancy
    integer                                                                                  :: binCountMinimum
    logical                                                                                  :: likelihoodPoisson                 , report             , &
         &                                                                                      binAverage                        , includeCorrelations

    if (.not.parameters%isPresent('fileNames')) call Error_Report('`fileNames` parameter is not present'//{introspection:location})
    if (.not.parameters%isPresent('redshifts')) call Error_Report('`redshifts` parameter is not present'//{introspection:location})
    allocate(fileNames(parameters%count('fileNames')))
    allocate(redshifts(parameters%count('redshifts')))
    !![
    <inputParameter>
      <name>baseParametersFileName</name>
      <description>The base set of parameters to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fileNames</name>
      <description>The names of the files containing the halo mass functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshifts</name>
      <description>The redshifts at which to evaluate the halo mass functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massRangeMinimum</name>
      <description>The minimum halo mass to include in the likelihood evaluation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massRangeMaximum</name>
      <description>The maximum halo mass to include in the likelihood evaluation.</description>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>binCountMinimum</name>
      <description>The minimum number of halos per bin required to permit bin to be included in likelihood evaluation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>likelihoodPoisson</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, likelihood is computed assuming a Poisson distribution for the number of halos in each bin (with no covariance between bins). Otherwise a multivariate normal is assumed when computing likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>varianceFractionalModelDiscrepancy</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The fractional variance due to model discrepancy.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>report</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, give detailed reporting on likelihood calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeCorrelations</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, account for correlations between halo mass functions measured at different redshifts.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>binAverage</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, the mass function is averaged over each bin.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    allocate(changeParametersFileNames(parameters%count('changeParametersFileNames',zeroIfNotPresent=.true.)))
    if (size(changeParametersFileNames) > 0) then
       !![
       <inputParameter>
	 <name>changeParametersFileNames</name>
	 <description>The names of files containing parameter changes to be applied.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    allocate(parametersModel)
    parametersModel=inputParameters(baseParametersFileName,noOutput=.true.,changeFiles=changeParametersFileNames)
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parametersModel"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parametersModel"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parametersModel"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parametersModel"/>
    <objectBuilder class="randomNumberGenerator"    name="randomNumberGenerator_"    source="parameters"     />
    !!]
    self=posteriorSampleLikelihoodHaloMassFunction(fileNames,redshifts,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,binAverage,includeCorrelations,report,parametersModel,changeParametersFileNames,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="randomNumberGenerator_"   />
    !!]
    self%baseParametersFileName=baseParametersFileName
    nullify(parametersModel)
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(fileNames,redshifts,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,binAverage,includeCorrelations,report,parametersModel,changeParametersFileNames,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_,randomNumberGenerator_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily haloMassFunction} posterior sampling likelihood class.
    !!}
    use :: Display                 , only : displayMessage  , displayMagenta, displayReset
    use :: Error                   , only : Error_Report
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: ISO_Varying_String      , only : char
    use :: Linear_Algebra          , only : assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    use :: File_Utilities          , only : File_Name_Expand
    implicit none
    type            (posteriorSampleLikelihoodHaloMassFunction)                                :: self
    type            (varying_string                           ), intent(in   ), dimension(:  ) :: fileNames
    double precision                                           , intent(in   ), dimension(:  ) :: redshifts
    double precision                                           , intent(in   )                 :: massRangeMinimum                  , massRangeMaximum   , &
         &                                                                                        varianceFractionalModelDiscrepancy
    integer                                                    , intent(in   )                 :: binCountMinimum
    logical                                                    , intent(in   )                 :: likelihoodPoisson                 , report             , &
         &                                                                                        binAverage                        , includeCorrelations
    type            (inputParameters                          ), intent(inout), target         :: parametersModel
    class           (cosmologyFunctionsClass                  ), intent(inout), target         :: cosmologyFunctions_
    class           (criticalOverdensityClass                 ), intent(inout), target         :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass            ), intent(inout), target         :: cosmologicalMassVariance_
    class           (linearGrowthClass                        ), intent(inout), target         :: linearGrowth_
    class           (randomNumberGeneratorClass               ), intent(in   ), target         :: randomNumberGenerator_
    type            (varying_string                           ), intent(in   ), dimension(:  ) :: changeParametersFileNames
    double precision                                           , allocatable  , dimension(:  ) :: eigenValueArray                   , massOriginal       , &
         &                                                                                        massFunctionOriginal
    integer         (c_size_t                                 ), allocatable  , dimension(:  ) :: massFunctionCountOriginal
    double precision                                           , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
    character       (len=12                                   )                                :: redshiftLabel
    type            (hdf5Object                               )                                :: massFunctionFile                  , simulationGroup
    integer                                                                                    :: i                                 , j               , &
         &                                                                                        ii                                , jj              , &
         &                                                                                        massCountReduced                  , iRedshift
    double precision                                                                           :: massIntervalLogarithmic
    type            (matrix                                   )                                :: eigenVectors
    type            (vector                                   )                                :: eigenValues
    !![
    <constructorAssign variables="fileNames, redshifts, binCountMinimum, massRangeMinimum, massRangeMaximum, likelihoodPoisson, varianceFractionalModelDiscrepancy, binAverage, includeCorrelations, report, changeParametersFileNames, *parametersModel, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_, *randomNumberGenerator_"/>
    !!]

    ! Convert redshifts to times.
    allocate(self%times                ,mold=redshifts)
    allocate(self%countConversionFactor,mold=redshifts)
    do iRedshift=1,size(redshifts)
       self%times(iRedshift)=cosmologyFunctions_%cosmicTime                 (                      &
            &                cosmologyFunctions_%expansionFactorFromRedshift (                     &
            &                                                                 redshifts(iRedshift) &
            &                                                                )                     &
            &                                                               )
    end do
    ! Record that we have not yet checked if model discrepancy terms should be included.
    self%includeDiscrepancyChecked=.false.
    ! Read the halo mass function files.
    do iRedshift=1,size(redshifts)
       write (redshiftLabel,'(f6.3)') redshifts(iRedshift)
       !$ call hdf5Access%set()
       call massFunctionFile%openFile(char(File_Name_Expand(char(fileNames(iRedshift)))),readOnly=.true.)
       simulationGroup=massFunctionFile%openGroup('simulation0001')
       call simulationGroup %readDataset("mass"        ,massOriginal             )
       call simulationGroup %readDataset("massFunction",massFunctionOriginal     )
       call simulationGroup %readDataset("count"       ,massFunctionCountOriginal)
       call simulationGroup %close      (                                        )
       call massFunctionFile%close      (                                        )
       !$ call hdf5Access%unset()    
       ! Compute quantities needed for likelihood calculations.
       if (self%likelihoodPoisson) then
          ! Find a reduced mass function excluding bins below the mass threshold.
          massCountReduced=0
          do i=1,size(massOriginal)       
             if     (                                    &
                  &   massOriginal(i) < massRangeMinimum &
                  &  .or.                                &
                  &   massOriginal(i) > massRangeMaximum &
                  & ) cycle
             massCountReduced=massCountReduced+1
          end do
          if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileNames(iRedshift))//"'"//{introspection:location})
          if (iRedshift == 1) then
             allocate(self%mass        (massCountReduced                ))
             allocate(self%massFunction(massCountReduced,size(redshifts)))
             allocate(self%countHalos  (massCountReduced,size(redshifts)))
          end if
          ii=0
          do i=1,size(massOriginal)
             if     (                                    &
                  &   massOriginal(i) < massRangeMinimum &
                  &  .or.                                &
                  &   massOriginal(i) > massRangeMaximum &
                  & ) cycle
             ii=ii+1
             self%mass        (ii          )=massOriginal             (i)
             self%massFunction(ii,iRedshift)=massFunctionOriginal     (i)
             self%countHalos  (ii,iRedshift)=massFunctionCountOriginal(i)
          end do
          ! Compute the conversion factor between halo count per bin and the mass function.
          self%countConversionFactor(iRedshift)=+     sum  (dble(self%countHalos(:,iRedshift))/self%massFunction(:,iRedshift),mask=self%massFunction(:,iRedshift) > 0.0d0)  &
               &                                /dble(count(                                                                  mask=self%massFunction(:,iRedshift) > 0.0d0))
       else
          ! Construct the covariance matrix.
          allocate(massFunctionCovarianceOriginal(size(massOriginal),size(massOriginal)))
          massFunctionCovarianceOriginal=0.0d0
          do i=1,size(massOriginal)
             do j=1,size(massOriginal)
                if   (                                           &
                     &   massFunctionCountOriginal(i) > 0_c_size_t &
                     &  .and.                                      &
                     &   massFunctionCountOriginal(j) > 0_c_size_t &
                     & ) then
                   ! Compute the Poisson contribution.
                   if (i == j) massFunctionCovarianceOriginal(i,j)=+     massFunctionCovarianceOriginal(i,j)     &
                        &                                          +     massFunctionOriginal          (i  ) **2 &
                        &                                          /dble(massFunctionCountOriginal     (i  ))
                end if
             end do
          end do
          ! Find a reduced mass function excluding any empty bins.
          massCountReduced=0
          do i=1,size(massOriginal)       
             if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
             if (massOriginal             (i) <  massRangeMinimum) cycle
             if (massOriginal             (i) >  massRangeMaximum) cycle
             if (massFunctionCountOriginal(i) <  binCountMinimum ) cycle
             massCountReduced=massCountReduced+1
          end do
          if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileNames(iRedshift))//"'"//{introspection:location})
          if (iRedshift == 1) then
             allocate(self%mass            (massCountReduced                                 ))
             allocate(self%massFunction    (massCountReduced                 ,size(redshifts)))
             allocate(self%covarianceMatrix(massCountReduced,massCountReduced,size(redshifts)))
          end if
          ii=0
          do i=1,size(massOriginal)
             if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
             if (massOriginal             (i) <  massRangeMinimum) cycle
             if (massOriginal             (i) >  massRangeMaximum) cycle
             if (massFunctionCountOriginal(i) <  binCountMinimum ) cycle
             ii=ii+1
             self%mass        (ii          )=massOriginal        (i)
             self%massFunction(ii,iRedshift)=massFunctionOriginal(i)
             jj=0
             do j=1,size(massOriginal)
                if (massFunctionOriginal     (j) <= 0.0d0           ) cycle
                if (massOriginal             (j) <  massRangeMinimum) cycle
                if (massOriginal             (j) >  massRangeMaximum) cycle
                if (massFunctionCountOriginal(j) <  binCountMinimum ) cycle
                jj=jj+1
                self%covarianceMatrix(ii,jj,iRedshift)=massFunctionCovarianceOriginal(i,j)
             end do
          end do
          ! Find the covariance matrices.
          self%covariance(iRedshift)=self%covarianceMatrix(:,:,iRedshift)
          ! Get eigenvalues and vectors of the covariance matrix.
          if (iRedshift == 1) allocate(eigenValueArray(size(self%mass)))
          call self%covariance(iRedshift)%eigenSystem(eigenVectors,eigenValues)
          eigenValueArray=eigenValues
          if (any(eigenValueArray < 0.0d0)) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' inverse covariance matrix is not semi-positive definite')
          deallocate(eigenValueArray               )
       end if
       ! Compute mass ranges for bins.
       if (iRedshift == 1) then
          massIntervalLogarithmic=+log(                                  &
               &                       +massOriginal(size(massOriginal)) &
               &                       /massOriginal(                 1) &
               &                      )                                  &
               &                  /dble(                                 &
               &                        +size(massOriginal)              &
               &                        -1                               &
               &                       )
          allocate(self%massMinimum,mold=self%mass)
          allocate(self%massMaximum,mold=self%mass)
          do i=1,size(self%mass)
             self%massMinimum(i)=self%mass(i)*exp(-0.5d0*massIntervalLogarithmic)
             self%massMaximum(i)=self%mass(i)*exp(+0.5d0*massIntervalLogarithmic)
          end do
       end if
    end do
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily haloMassFunction} posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/> 
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%randomNumberGenerator_"   />
   !!]
    if (associated(self%parametersModel)) then
       call self%parametersModel%destroy()
       deallocate(self%parametersModel)
    end if
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance) result(logLikelihood)
    !!{
    Return the log-likelihood for the halo mass function likelihood function. If {\normalfont \ttfamily [likelihoodPoisson]=false}
    then Gaussian statistics are assumed, otherwise, Poisson statistics are assumed. No covariance between bins is assumed.
    !!}

use :: MPI_Utilities, only : mpiSelf

    
    use :: Error                            , only : Error_Report
    use :: Display                          , only : displayMessage                                , displayIndent                         , displayUnindent
    use :: Halo_Mass_Functions              , only : haloMassFunctionClass
    use :: Interface_GSL                    , only : GSL_Success                                   , GSL_ETol                              , GSL_EMaxIter
    use :: Linear_Algebra                   , only : assignment(=)                                 , operator(*)
    use :: Models_Likelihoods_Constants     , only : logImpossible                                 , logImprobable
    use :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass                       , nbodyHaloMassErrorPowerLaw            , nbodyHaloMassErrorSOHaloFinder     , nbodyHaloMassErrorTrenti2010
    use :: Statistics_Distributions         , only : distributionFunction1DNormal                  , distributionFunctionMultivariateNormal
    use :: Statistics_Distributions_Discrete, only : distributionFunctionDiscrete1DNegativeBinomial, distributionFunctionDiscrete1DPoisson , distributionFunctionDiscrete1DClass
    use :: Factorials                       , only : Logarithmic_Factorial
    use :: Gamma_Functions                  , only : Gamma_Function_Logarithmic
    use :: Galacticus_Nodes                 , only : treeNode                       , nodeComponentBasic
    implicit none
    class           (posteriorSampleLikelihoodHaloMassFunction), intent(inout), target         :: self
    class           (posteriorSampleStateClass                ), intent(inout)                 :: simulationState
    type            (modelParameterList                       ), intent(inout), dimension(:)   :: modelParametersActive_                        , modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)                 :: simulationConvergence
    double precision                                           , intent(in   )                 :: temperature                                   , logLikelihoodCurrent     , &
         &                                                                                        logPriorCurrent                               , logPriorProposed
    real                                                       , intent(inout)                 :: timeEvaluate
    double precision                                           , intent(  out), optional       :: logLikelihoodVariance
    logical                                                    , intent(inout), optional       :: forceAcceptance
    double precision                                           , allocatable  , dimension(:  ) :: stateVector                                   , meanCopula               , &
         &                                                                                        xCopulaLow                                    , xCopulaHigh
    double precision                                           , allocatable  , dimension(:,:) :: countHalosMeanPerBin                          , massFunction             , &
         &                                                                                        likelihoodPerBin                              , correlationCopula
    double precision                                           , parameter                     :: errorFractionalMaximum                 =1.0d+1
    double precision                                           , parameter                     :: countStandardDeviationsMaximum         =1.0d+1
    class           (haloMassFunctionClass                    ), pointer                       :: haloMassFunction_
    class           (distributionFunctionDiscrete1DClass      ), allocatable                   :: distributionFunctionDiscrete1D_
    type            (distributionFunction1DNormal             )                                :: distributionFunctionNormal_
    type            (distributionFunctionMultivariateNormal   )                                :: distributionFunctionMultivariateNormal_
    type            (vector                                   )                                :: difference
    type            (treeNode                                 ), pointer                       :: node
    class           (nodeComponentBasic                       ), pointer                       :: basic
    logical                                                                                    :: evaluationFailed
    integer                                                                                    :: i                                             , status                              , &
         &                                                                                        iTime                                         , jTime                               , &
         &                                                                                        iCopula                                       , jCopula                             , &
         &                                                                                        dimensionCopula                               , j
    double precision                                                                           :: countHalosMean                                , stoppingTimeParameter               , &
         &                                                                                        varianceFractionalModelDiscrepancy            , logLikelihood_                      , &
         &                                                                                        deviationMaximum                              , logLikelihoodUncorrelated           , &
         &                                                                                        linearGrowthFactorLate                        , linearGrowthFactorEarly             , &
         &                                                                                        rootVarianceLate                              , rootVarianceEarly                   , &
         &                                                                                        rootVarianceLogarithmicGradientLate           , rootVarianceLogarithmicGradientEarly, &
         &                                                                                        widthBinMassLogarithmic                       , peakHeightEffective                 , &
         &                                                                                        fractionMassProgenitors                       , correlation                         , &
         &                                                                                        yCopula
    type            (varying_string                           )                                :: message
    character       (len=17                                   )                                :: label
    real                                                                                       :: timeBegin                                     , timeEnd
    !$GLC attributes unused :: simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Initialize log-likelihood.
    logLikelihood=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) return
    ! Ensure pointers into the base parameters are initialized.
    call self%initialize(modelParametersActive_,modelParametersInactive_)
    ! Record start time.
    call CPU_Time(timeBegin)
    ! Get states for all chains.
    allocate(stateVector(simulationState%dimension()))
    stateVector=simulationState%get()
    ! Update parameter values.
    call self%update(simulationState,modelParametersActive_,modelParametersInactive_,stateVector)
    if (.not.self%includeDiscrepancyChecked) then
       self%includeDiscrepancyChecked=.true.
       self%indexDiscrepancy         =-1
       do i=1,size(self%modelParametersActive_)
          if (self%modelParametersActive_(i)%definition == "varianceFractionalModelDiscrepancy") then
             if (self%indexDiscrepancy > 0) call Error_Report('multiple instances of parameter [varianceFractionalModelDiscrepancy] found'//{introspection:location})
             self%indexDiscrepancy=i
          end if
       end do
    end if
    ! Get the halo mass function object.
    !![
    <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="self%parametersModel"/>
    !!]
    ! Determine the model discrepancy variance term.
    if (self%indexDiscrepancy > 0) then
       varianceFractionalModelDiscrepancy=modelParametersActive_(self%indexDiscrepancy)%modelParameter_%unmap(stateVector(self%indexDiscrepancy))
    else
       varianceFractionalModelDiscrepancy=self%varianceFractionalModelDiscrepancy
    end if
    if (varianceFractionalModelDiscrepancy > 0.0d0) then
       stoppingTimeParameter=1.0d0/varianceFractionalModelDiscrepancy
    else
       stoppingTimeParameter=0.0d0
    end if
    ! Create a node.
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    ! If correlations are to be included, allocate space for the relevant calculations.
    if (self%includeCorrelations) then
       dimensionCopula=+size(self%mass ) &
            &          *size(self%times)
       if (self%likelihoodPoisson) then
          allocate(xCopulaLow (dimensionCopula))
          allocate(xCopulaHigh(dimensionCopula))
          distributionFunctionNormal_=distributionFunction1DNormal(0.0d0,1.0d0)
       end if
    end if
    ! Iterate over times.
    allocate(massFunction(size(self%mass),size(self%times)))
    do iTime=1,size(self%times)
       call basic%timeSet(self%times(iTime))
       ! Compute the mass function.
       evaluationFailed=.false.
       do i=1,size(self%mass)
          if (self%binAverage) then
             massFunction(i,iTime)=+haloMassFunction_%integrated  (                                &
                  &                                                       self%times      (iTime), &
                  &                                                       self%massMinimum(i    ), &
                  &                                                       self%massMaximum(i    ), &
                  &                                                node  =     node              , &
                  &                                                status=     status              &
                  &                                               )                                &
                  &                /log(                                                           &
                  &                                             +         self%massMaximum(i    )  &
                  &                                             /         self%massMinimum(i    )  &
                  &                    )             
             if (status /= errorStatusSuccess) then
                logLikelihood   =logImprobable
                evaluationFailed=.true.
                exit
             end if
          else
            massFunction(i,iTime)=+haloMassFunction_%differential(                                &
                  &                                                       self%times      (iTime), &
                  &                                                       self%mass       (i    ), &
                  &                                                node  =     node                &
                  &                                               )                                &
                  &                *                                      self%mass       (i    )
          end if
       end do
       if (evaluationFailed) exit
    end do
    call node%destroy()
    deallocate(node)
    ! Allocate array for per bin likelihood if needed.
    if (self%report) then
       allocate(likelihoodPerBin    (size(self%mass),size(self%times)))
       allocate(countHalosMeanPerBin(size(self%mass),size(self%times)))
       likelihoodPerBin    =+0.0d0
       countHalosMeanPerBin=-1.0d0
    else
       allocate(likelihoodPerBin    (             0 ,               0))
       allocate(countHalosMeanPerBin(             0 ,               0))
    end if
    ! Evaluate the log-likelihood.
    if (.not.evaluationFailed) then
       do iTime=1,size(self%times)
          logLikelihood_=0.0d0
          if (self%likelihoodPoisson) then
             ! Assume Poisson/negative binomial statistics. We treat each bin as independent with a pure Poisson/negative binomial distribution.
             do i=1,size(self%mass)
                ! Find the mean number of halos expected in this bin based on our model mass function.
                countHalosMean=+self%countConversionFactor(  iTime) &
                     &         *     massFunction         (i,iTime)
                if (self%report) countHalosMeanPerBin(i,iTime)=countHalosMean
                ! If the expected mean is zero, and the measured number is non-zero, this is impossible.
                if (countHalosMean <= 0.0d0) then
                   if (self%countHalos(i,iTime) > 0) then
                      logLikelihood_=logImprobable
                      if (self%report) likelihoodPerBin=logImprobable
                      exit
                   else if (self%includeCorrelations.and.self%likelihoodPoisson) then
                      iCopula             =(iTime-1)*size(self%mass)+i
                      xCopulaLow (iCopula)=-huge(0.0d0)
                      xCopulaHigh(iCopula)=+huge(0.0d0)
                   end if
                else
                   if (varianceFractionalModelDiscrepancy <= 0.0d0) then
                      ! Evaluate the Poisson likelihood (zero model discrepancy term).
                      logLikelihood_=+dble                 (    self%countHalos    (i,iTime))  &
                           &         *log                  (         countHalosMean         )  &
                           &         -                               countHalosMean            &
                           &         -Logarithmic_Factorial(int(self%countHalos    (i,iTime)))
                      if (self%includeCorrelations) then
                         allocate(distributionFunctionDiscrete1DPoisson :: distributionFunctionDiscrete1D_)
                         distributionFunctionDiscrete1D_=distributionFunctionDiscrete1DPoisson(countHalosMean)
                      end if
                   else
                      ! Evaluate the negative binomial likelihood (non-zero model discrepancy term). Here the negative binomial
                      ! distribution (which is used in the alternative parameterization as given by, e.g.,
                      ! https://en.wikipedia.org/wiki/Negative_binomial_distribution#Poisson_distribution) represents an
                      ! over-dispersed Poisson distribution.
                      logLikelihood_=+dble(self%countHalos           (i,iTime))*log                       (                                               countHalosMean          ) &
                           &         -                                          Logarithmic_Factorial     (                                    +int (self%countHalos    (i,iTime))) &
                           &         +                                          Gamma_Function_Logarithmic(               stoppingTimeParameter+dble(self%countHalos    (i,iTime))) &
                           &         -                                          Gamma_Function_Logarithmic(               stoppingTimeParameter                                   ) &
                           &         -dble(self%countHalos           (i,iTime))*log                       (               stoppingTimeParameter+          countHalosMean          ) &
                           &         -          stoppingTimeParameter          *log                       (countHalosMean/stoppingTimeParameter+1.0d0                             )
                      if (self%includeCorrelations) then
                         allocate(distributionFunctionDiscrete1DNegativeBinomial :: distributionFunctionDiscrete1D_)
                         distributionFunctionDiscrete1D_=distributionFunctionDiscrete1DNegativeBinomial(stoppingTimeParameter/(stoppingTimeParameter+countHalosMean),stoppingTimeParameter)
                     end if
                   end if
                   if (self%includeCorrelations) then
                      ! Evaluate the cumulative probability bounds for this marginal distribution.
                      iCopula         =(iTime-1)*size(self%mass)+i
                      deviationMaximum=countStandardDeviationsMaximum*max(sqrt(countHalosMean),sqrt(max(0.0d0,varianceFractionalModelDiscrepancy))*countHalosMean)
                      !! Evaluate the upper bound.
                      !!! First find the cumulative probability for our marginal distribution.
                      yCopula                   =+distributionFunctionDiscrete1D_%cumulative             (int(self%countHalos(i,iTime)  ),status)
                      if (status /= GSL_Success) then
                         logLikelihood_=logImprobable
                         exit
                      end if
                      if    (yCopula < 1.0d0) then
                         ! The marginal probability is less than 1 - find the argument in a standard normal with the same
                         ! cumulative probability.
                         xCopulaHigh   (iCopula)=+distributionFunctionNormal_    %inverse                (         yCopula                      )
                      else
                         ! The marginal probability is 1 to numerical precision. In this case, we can exploit the fact that the
                         ! standard normal distribution is symmetric. We find the complementary cumulative probabilty (i.e. 1
                         ! minus the usual cumulative probability) - which can be determined more accurately in this regime - then
                         ! find the argument of the standard normal that gives this same probability, and then negate it so that
                         ! we find the argument that we would have found if we were able to evaluate the regular cumulative
                         ! probability to better precision.
                         yCopula                =+distributionFunctionDiscrete1D_%cumulativeComplementary(int(self%countHalos(i,iTime)  ),status)
                         if (status /= GSL_Success) then
                            logLikelihood_=logImprobable
                            exit
                         end if
                         xCopulaHigh   (iCopula)=-distributionFunctionNormal_    %inverse                (         yCopula                      )
                      end if
                      !! Evaluate the lower bound.
                      if (self%countHalos(i,iTime) > 0.0d0) then
                         ! The count of halos in the target is greater than zero, so evaluate the lower bound as one less than this.
                         !! First find the cumulative probability for our marginal distribution.
                         yCopula                =+distributionFunctionDiscrete1D_%cumulative             (int(self%countHalos(i,iTime)-1),status)
                         if (status /= GSL_Success) then
                            logLikelihood_=logImprobable
                            exit
                         end if
                         if (yCopula < 1.0d0) then
                            xCopulaLow (iCopula)=+distributionFunctionNormal_    %inverse                (         yCopula                      )
                         else
                            ! The marginal probabiltiy is 1 to numerical precision. In this case, we can exploit the fact that the
                            ! standard normal distribution is symmetric. We find the complementary cumulative probabilty (i.e. 1
                            ! minus the usual cumulative probability) - which can be determined more accurately in this regime - then
                            ! find the argument of the standard normal that gives this same probability, and then negate it so that
                            ! we find the argument that we would have found if we were able to evaluate the regular cumulative
                            ! probability to better precision.
                            yCopula             =+distributionFunctionDiscrete1D_%cumulativeComplementary(int(self%countHalos(i,iTime)-1),status)
                            if (status /= GSL_Success) then
                               logLikelihood_=logImprobable
                               exit
                            end if
                            xCopulaLow (iCopula)=-distributionFunctionNormal_    %inverse                (         yCopula                      )
                         end if
                      else
                         ! The count of halos in the target is zero, so the lower bound of the copula must correspond to negative infinity.
                         xCopulaLow    (iCopula)=-huge(0.0d0)
                      end if
                      ! Ensure that the copula bounds are correctly ordered. (Mis-orderings can happen due to numerical precision.)
                      xCopulaHigh(iCopula)=max(xCopulaLow(iCopula),xCopulaHigh(iCopula))
                      deallocate(distributionFunctionDiscrete1D_)                         
                   end if
                   ! Accumulate the likelihood.
                   logLikelihood=+logLikelihood  &
                        &        +logLikelihood_
                   if (self%report) likelihoodPerBin(i,iTime)=logLikelihood_
                end if
             end do
          else
             ! Assume Gaussian statistics.
             difference   =+     massFunction(:,iTime)                               &
                  &        -self%massFunction(:,iTime)
             logLikelihood=+logLikelihood                                            &
                    &      -0.5d0                                                    &
                    &      *self%covariance  (  iTime)%covarianceProduct(difference)
          end if
          if (logLikelihood_ <= logImprobable) then
             logLikelihood=logImprobable
             exit
          end if
       end do
       if (logLikelihood > logImprobable) then
          ! If correlations are to be included, evaluate the likelihood for the multivariate distribution.
          if (self%includeCorrelations) then
             logLikelihoodUncorrelated=logLikelihood
             if (self%likelihoodPoisson) then
                ! Construct the copula and evaluate the probability.
                dimensionCopula=+size(self%mass ) &
                     &          *size(self%times)
                allocate(meanCopula       (dimensionCopula                ))
                allocate(correlationCopula(dimensionCopula,dimensionCopula))
                meanCopula       =0.0d0
                correlationCopula=0.0d0
                widthBinMassLogarithmic=log(self%mass(2)/self%mass(1))
                do i=1,dimensionCopula
                   correlationCopula(i,i)=1.0d0
                end do
                do iTime=1,size(self%times)
                   linearGrowthFactorLate       =self%linearGrowth_%value(self%times(iTime))
                   do i=2,size(self%mass)
                      iCopula      =(iTime-1)*size(self%mass)+i
                      call       self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(self%mass(i),self%times(iTime),rootVarianceLate ,rootVarianceLogarithmicGradientLate )
                      rootVarianceLate=rootVarianceLate/linearGrowthFactorLate
                      do jTime=1,size(self%times)
                         if (self%times(jTime) >= self%times(iTime)) cycle
                         linearGrowthFactorEarly=self%linearGrowth_%value(self%times(jTime))
                         do j=1,i-1
                            jCopula=(jTime-1)*size(self%mass)+j
                            call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(self%mass(j),self%times(jTime),rootVarianceEarly,rootVarianceLogarithmicGradientEarly)
                            rootVarianceEarly=rootVarianceEarly/linearGrowthFactorEarly
                            ! Evaluate the model for the correlation between these two bins.
                            if (massFunction(j,jTime) > 0.0d0) then
                               peakHeightEffective    =+    (                                                                            &
                                    &                        +self%criticalOverdensity_%value(self%times(jTime))/linearGrowthFactorEarly &
                                    &                        -self%criticalOverdensity_%value(self%times(iTime))/linearGrowthFactorLate  &
                                    &                       )                                                                            &
                                    &                  /sqrt(                                                                            &
                                    &                        +rootVarianceEarly**2                                                       &
                                    &                        -rootVarianceLate **2                                                       &
                                    &                       )
                               fractionMassProgenitors=+0.4d0                                     & ! "Global fit" from Cole et al. (2008; https://ui.adsabs.harvard.edu/abs/2008MNRAS.383..546C).
                                    &                  *    ( peakHeightEffective**0.75d0       ) &
                                    &                  *exp (-peakHeightEffective**3     /10.0d0)
                               correlation            =+min(                                                              & ! Limit correlation to maximum of 1.0 - this can be exceeded
                                    &                       +1.0d0                                                      , & ! if the current model mass function has unrealistic redshift
                                    &                       +sqrt(                                                        & ! evolution.
                                    &                             +self%mass        (i      )/self%mass        (j      )  &
                                    &                             *     massFunction(i,iTime)/     massFunction(j,jTime)  &
                                    &                             *    fractionMassProgenitors                            &
                                    &                             *abs(rootVarianceLogarithmicGradientLate)               &
                                    &                             *    widthBinMassLogarithmic                            &
                                    &                            )                                                        &
                                    &                      )
                            else
                               correlation            =+0.0d0
                            end if
                            correlationCopula(iCopula,jCopula)=correlation
                            correlationCopula(jCopula,iCopula)=correlation
                         end do
                      end do
                   end do
                end do
                distributionFunctionMultivariateNormal_=distributionFunctionMultivariateNormal            (meanCopula,correlationCopula,self%randomNumberGenerator_              )
                logLikelihood                          =distributionFunctionMultivariateNormal_%cumulative(xCopulaLow,xCopulaHigh      ,logarithmic=.true.         ,status=status)
                select case (status)
                case (GSL_Success             )
                   ! Successful evaluation - proceed.
                case (GSL_ETol   ,GSL_EMaxIter)
                   ! Tolerable errors - proceed.
                case default
                   ! Intolerable error - assume an improbable likelihood.
                   logLikelihood=logImprobable
                   return
                end select
             else
                ! Correlations not supported for Gaussian statistics.
                call Error_Report('inclusion of correlations is not supported for Gaussian statistics'//{introspection:location})
             end if
          end if
       end if
    end if
    ! Record timing information.
    call CPU_Time(timeEnd)
    timeEvaluate=timeEnd-timeBegin
    if (self%report) then
       write (label,'(e17.10)') logLikelihood
       if (evaluationFailed) then
          call displayMessage("model evaluation failed - no likelihood computed")
       else
          call displayMessage("logℒ = "//trim(label))
          if (self%includeCorrelations) then
             write (label,'(e17.10)') logLikelihoodUncorrelated
             call displayMessage("logℒ = "//trim(label)//" (assuming uncorrelated data)")
          end if
          if (self%likelihoodPoisson) then
             if (varianceFractionalModelDiscrepancy <= 0.0d0) then
                call displayMessage("Likelihood model: Poisson"          )
             else
                call displayMessage("Likelihood model: Negative binomial")
             end if
          else
             call    displayMessage("Likelihood model: normal"           )
          end if
          do iTime=1,size(self%times)
             call displayIndent("Likelihood report for: "//char(self%fileNames(iTime)))      
             if (self%likelihoodPoisson) then
                call displayIndent("per bin likelihoods")
                do i=1,size(self%mass)
                   write (label,'(i4)') i
                   message=trim(label)//": "
                   write (label,'(e17.10)') log10(self%mass                (i      ))
                   message=message//"log₁₀(M/M☉) = " //trim(label)//"; "
                   write (label,'(e17.10)')            likelihoodPerBin    (i,iTime)
                   message=message//"logℒ = "        //trim(label)//"; "
                   write (label,'(i07)'   )       self%countHalos          (i,iTime)
                   message=message//"Nhalo(target) = "//trim(label)//"; "
                   write (label,'(e17.10)')            countHalosMeanPerBin(i,iTime)
                   message=message//"Nhalo(model) = " //trim(label)
                   call displayMessage(message)
                end do
                call displayUnindent("done")
             end if
             call displayUnindent("done")
          end do
       end if
    end if
    ! Clean up.
    !![
    <objectDestructor name="haloMassFunction_"/>
    !!]
    deallocate(stateVector )
    deallocate(massFunction)
    return
  end function haloMassFunctionEvaluate
  
  subroutine haloMassFunctionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine haloMassFunctionFunctionChanged
