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

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  
  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodHaloMassFunction">
   <description>A posterior sampling likelihood class which implements a likelihood for halo mass functions.</description>
   <runTimeFileDependencies paths="fileName"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodBaseParameters) :: posteriorSampleLikelihoodHaloMassFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_                => null()
     double precision                               , dimension(:  ), allocatable :: mass                                        , massFunction             , &
          &                                                                          massMinimum                                 , massMaximum
     double precision                               , dimension(:,:), allocatable :: covarianceMatrix
     integer         (c_size_t                     ), dimension(:  ), allocatable :: countHalos
     double precision                                                             :: time                                        , massParticle             , &
          &                                                                          massRangeMinimum                            , massRangeMaximum         , &
          &                                                                          countConversionFactor                       , redshift                 , &
          &                                                                          varianceFractionalModelDiscrepancy
     logical                                                                      :: likelihoodPoisson                           , includeDiscrepancyChecked, &
          &                                                                          report                                      , binAverage
     integer                                                                      :: binCountMinimum                             , indexDiscrepancy
     type            (vector                       )                              :: means
     type            (matrix                       )                              :: covariance
     type            (varying_string               )                              :: fileName
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
    type            (varying_string                           ), allocatable  , dimension(:) :: changeParametersFileNames
    type            (varying_string                           )                              :: fileName                 , baseParametersFileName
    double precision                                                                         :: redshift                 , massRangeMinimum                  , &
         &                                                                                      massRangeMaximum         , varianceFractionalModelDiscrepancy
    integer                                                                                  :: binCountMinimum
    logical                                                                                  :: likelihoodPoisson        , report                            , &
         &                         binAverage

    !![
    <inputParameter>
      <name>baseParametersFileName</name>
      <description>The base set of parameters to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift at which to evaluate the halo mass function.</description>
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
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parametersModel"/>
    !!]
    self=posteriorSampleLikelihoodHaloMassFunction(char(fileName),redshift,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,binAverage,report,parametersModel,changeParametersFileNames,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    self%baseParametersFileName=baseParametersFileName
    nullify(parametersModel)
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(fileName,redshift,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,binAverage,report,parametersModel,changeParametersFileNames,cosmologyFunctions_) result(self)
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
    character       (len=*                                    ), intent(in   )                 :: fileName
    double precision                                           , intent(in   )                 :: redshift                      , massRangeMinimum                  , &
         &                                                                                        massRangeMaximum              , varianceFractionalModelDiscrepancy
    integer                                                    , intent(in   )                 :: binCountMinimum
    logical                                                    , intent(in   )                 :: likelihoodPoisson             , report                             , &
         &                           binAverage
    type            (inputParameters                          ), intent(inout), target         :: parametersModel
    class           (cosmologyFunctionsClass                  ), intent(inout), target         :: cosmologyFunctions_
    type            (varying_string                           ), intent(in   ), dimension(:  ) :: changeParametersFileNames
    double precision                                           , allocatable  , dimension(:  ) :: eigenValueArray               , massOriginal                      , &
         &                                                                                        massFunctionOriginal
    integer         (c_size_t                                 ), allocatable  , dimension(:  ) :: massFunctionCountOriginal
    double precision                                           , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
    character       (len=12                                   )                                :: redshiftLabel
    type            (hdf5Object                               )                                :: massFunctionFile              , simulationGroup
    integer                                                                                    :: i                             , j                                 , &
         &                                                                                        ii                            , jj                                , &
         &                                                                                        massCountReduced
    double precision                                                                           :: massIntervalLogarithmic
    type            (matrix                                   )                                :: eigenVectors
    type            (vector                                   )                                :: eigenValues
    !![
    <constructorAssign variables="fileName, redshift, binCountMinimum, massRangeMinimum, massRangeMaximum, likelihoodPoisson, varianceFractionalModelDiscrepancy, binAverage, report, changeParametersFileNames, *parametersModel, *cosmologyFunctions_"/>
    !!]

    ! Convert redshift to time.
    self%time=cosmologyFunctions_%cosmicTime                 (          &
         &    cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                     redshift &
         &                                                    )         &
         &                                                   )
    ! Record that we have not yet checked if model discrepancy terms should be included.
    self%includeDiscrepancyChecked=.false.
    ! Read the halo mass function file.
    write (redshiftLabel,'(f6.3)') redshift
    !$ call hdf5Access%set()
    call massFunctionFile%openFile(char(File_Name_Expand(trim(fileName))),readOnly=.true.)
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
       if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileName)//"'"//{introspection:location})
       allocate(self%mass        (massCountReduced))
       allocate(self%massFunction(massCountReduced))
       allocate(self%countHalos  (massCountReduced))
       ii=0
       do i=1,size(massOriginal)
          if     (                                    &
               &   massOriginal(i) < massRangeMinimum &
               &  .or.                                &
               &   massOriginal(i) > massRangeMaximum &
               & ) cycle
          ii=ii+1
          self%mass        (ii)=massOriginal             (i)
          self%massFunction(ii)=massFunctionOriginal     (i)
          self%countHalos  (ii)=massFunctionCountOriginal(i)
       end do
       ! Compute the conversion factor between halo count per bin and the mass function.
       self%countConversionFactor=+     sum  (dble(self%countHalos)/self%massFunction,mask=self%massFunction > 0.0d0)  &
            &                     /dble(count(                                        mask=self%massFunction > 0.0d0))
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
       if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileName)//"'"//{introspection:location})
       allocate(self%mass            (massCountReduced                 ))
       allocate(self%massFunction    (massCountReduced                 ))
       allocate(self%covarianceMatrix(massCountReduced,massCountReduced))
       ii=0
       do i=1,size(massOriginal)
          if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
          if (massOriginal             (i) <  massRangeMinimum) cycle
          if (massOriginal             (i) >  massRangeMaximum) cycle
          if (massFunctionCountOriginal(i) <  binCountMinimum ) cycle
          ii=ii+1
          self%mass        (ii)=massOriginal        (i)
          self%massFunction(ii)=massFunctionOriginal(i)
          jj=0
          do j=1,size(massOriginal)
             if (massFunctionOriginal     (j) <= 0.0d0           ) cycle
             if (massOriginal             (j) <  massRangeMinimum) cycle
             if (massOriginal             (j) >  massRangeMaximum) cycle
             if (massFunctionCountOriginal(j) <  binCountMinimum ) cycle
             jj=jj+1
             self%covarianceMatrix(ii,jj)=massFunctionCovarianceOriginal(i,j)
          end do
       end do
       ! Find the covariance matrix.
       self%covariance=self%covarianceMatrix
       ! Get eigenvalues and vectors of the covariance matrix.
       allocate(eigenValueArray(size(self%mass)))
       call self%covariance%eigenSystem(eigenVectors,eigenValues)
       eigenValueArray=eigenValues
       if (any(eigenValueArray < 0.0d0)) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' inverse covariance matrix is not semi-positive definite')
       deallocate(eigenValueArray               )
    end if
    ! Compute mass ranges for bins.
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
    if (allocated(massOriginal                  )) deallocate(massOriginal                  )
    if (allocated(massFunctionOriginal          )) deallocate(massFunctionOriginal          )
    if (allocated(massFunctionCovarianceOriginal)) deallocate(massFunctionCovarianceOriginal)
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily haloMassFunction} posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    if (associated(self%parametersModel)) then
       call self%parametersModel%destroy()
       deallocate(self%parametersModel)
    end if
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo mass function likelihood function. If {\normalfont \ttfamily [likelihoodPoisson]=false}
    then Gaussian statistics are assumed, otherwise, Poisson statistics are assumed. No covariance between bins is assumed.
    !!}
    use :: Error                            , only : Error_Report
    use :: Display                          , only : displayMessage                 , displayIndent             , displayUnindent
    use :: Halo_Mass_Functions              , only : haloMassFunctionClass
    use :: Linear_Algebra                   , only : assignment(=)                  , operator(*)
    use :: Models_Likelihoods_Constants     , only : logImpossible                  , logImprobable
    use :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass        , nbodyHaloMassErrorPowerLaw, nbodyHaloMassErrorSOHaloFinder, nbodyHaloMassErrorTrenti2010
    use :: Factorials                       , only : Logarithmic_Factorial
    use :: Gamma_Functions                  , only : Gamma_Function_Logarithmic
    use :: Galacticus_Nodes                 , only : treeNode                       , nodeComponentBasic
    implicit none
    class           (posteriorSampleLikelihoodHaloMassFunction), intent(inout), target       :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(inout), dimension(:) :: modelParametersActive_                   , modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature                              , logLikelihoodCurrent    , &
         &                                                                                      logPriorCurrent                          , logPriorProposed
    real                                                       , intent(inout)               :: timeEvaluate
    double precision                                           , intent(  out), optional     :: logLikelihoodVariance
    logical                                                    , intent(inout), optional     :: forceAcceptance
    double precision                                           , allocatable  , dimension(:) :: stateVector                              , massFunction            , &
         &                                                                                      likelihoodPerBin                         , countHalosMeanPerBin
    double precision                                           , parameter                   :: errorFractionalMaximum            =1.0d+1
    class           (haloMassFunctionClass                    ), pointer                     :: haloMassFunction_
    type            (vector                                   )                              :: difference
    type            (treeNode                                 ), pointer                     :: node
    class           (nodeComponentBasic                       ), pointer                     :: basic
    logical                                                                                  :: evaluationFailed
    integer                                                                                  :: i                                        , status
    double precision                                                                         :: countHalosMean                           , stoppingTimeParameter   , &
         &                                                                                      varianceFractionalModelDiscrepancy       , logLikelihood
    type            (varying_string                           )                              :: message
    character       (len=17                                   )                              :: label
    real                                                                                     :: timeBegin                                , timeEnd
    !$GLC attributes unused :: simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       haloMassFunctionEvaluate=0.0d0
       return
    end if
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
    call basic%timeSet(self%time)
    ! Compute the mass function.
    allocate(massFunction(size(self%mass)))
    evaluationFailed=.false.
    do i=1,size(self%mass)
       if (self%binAverage) then
          massFunction(i)=+haloMassFunction_%integrated  (                            &
               &                                                 self%time          , &
               &                                                 self%massMinimum(i), &
               &                                                 self%massMaximum(i), &
               &                                          node  =     node          , &
               &                                          status=     status          &
               &                                         )                            &
               &          /log(                                                       &
               &                                       +         self%massMaximum(i)  &
               &                                       /         self%massMinimum(i)  &
               &              )
          if (status /= errorStatusSuccess) then
             haloMassFunctionEvaluate=logImprobable
             evaluationFailed        =.true.
             exit
          end if
       else
          massFunction(i)=+haloMassFunction_%differential(                            &
               &                                                 self%time          , &
               &                                                 self%mass       (i), &
               &                                          node  =     node            &
               &                                         )                            &
               &          *                                      self%mass       (i)
       end if
    end do
    call node%destroy()
    deallocate(node)
    ! Evaluate the log-likelihood.
    if (.not.evaluationFailed) then
       if (self%likelihoodPoisson) then
          ! Allocate array for per bin likelihood if needed.
          if (self%report) then
             allocate(likelihoodPerBin    (size(self%mass)))
             allocate(countHalosMeanPerBin(size(self%mass)))
             likelihoodPerBin    =+0.0d0
             countHalosMeanPerBin=-1.0d0
          else
             allocate(likelihoodPerBin    (             0 ))
             allocate(countHalosMeanPerBin(             0 ))
          end if
          ! Assume Poisson/negative binomial statistics. We treat each bin as independent with a pure Poisson/negative binomial distribution.
          haloMassFunctionEvaluate=0.0d0
          do i=1,size(self%mass)
             ! Find the mean number of halos expected in this bin based on our model mass function.
             countHalosMean          =+                          self%countConversionFactor         &
                  &                   *                               massFunction            (i)
             if (self%report) countHalosMeanPerBin(i) =countHalosMean
             ! If the expected mean is zero, and the measured number is non-zero, this is impossible.
             if (countHalosMean <= 0.0d0) then
                if (self%countHalos(i) > 0) then
                   haloMassFunctionEvaluate=logImprobable
                   if (self%report) likelihoodPerBin=logImprobable
                   exit
                end if
             else
                if (varianceFractionalModelDiscrepancy <= 0.0d0) then
                   ! Evaluate the Poisson likelihood (zero model discrepancy term).
                   logLikelihood=+dble                 (    self%countHalos    (i))  &
                        &        *log                  (         countHalosMean   )  &
                        &        -                               countHalosMean      &
                        &        -Logarithmic_Factorial(int(self%countHalos    (i)))
                else
                   ! Evaluate the negative binomial likelihood (non-zero model discrepancy term). Here the negative binomial
                   ! distribution (which is used in the alternative parameterization as given by, e.g.,
                   ! https://en.wikipedia.org/wiki/Negative_binomial_distribution#Poisson_distribution) represents an
                   ! over-dispersed Poisson distribution.
                   logLikelihood=+dble(self%countHalos           (i))*log                       (                                               countHalosMean    ) &
                        &        -                                    Logarithmic_Factorial     (                                    +int (self%countHalos    (i))) &
                        &        +                                    Gamma_Function_Logarithmic(               stoppingTimeParameter+dble(self%countHalos    (i))) &
                        &        -                                    Gamma_Function_Logarithmic(               stoppingTimeParameter                             ) &
                        &        -dble(self%countHalos           (i))*log                       (               stoppingTimeParameter+          countHalosMean    ) &
                        &        -          stoppingTimeParameter    *log                       (countHalosMean/stoppingTimeParameter+1.0d0                       )
                   if (self%report) likelihoodPerBin(i)=logLikelihood
                end if
                haloMassFunctionEvaluate=+haloMassFunctionEvaluate &
                     &                   +logLikelihood
                if (self%report) likelihoodPerBin(i)=logLikelihood
             end if
          end do
       else
          ! Assume Gaussian statistics.
          difference              =+massFunction                                  &
               &                   -self%massFunction
          haloMassFunctionEvaluate=-0.5d0                                         &
               &                   *self%covariance%covarianceProduct(difference)
       end if
    end if
    ! Record timing information.
    call CPU_Time(timeEnd)
    timeEvaluate=timeEnd-timeBegin
    if (self%report) then
       call displayIndent("Likelihood report for: "//char(self%fileName))
       write (label,'(e17.10)') haloMassFunctionEvaluate
       if (self%likelihoodPoisson) then
          if (varianceFractionalModelDiscrepancy <= 0.0d0) then
             call displayMessage("Likelihood model: Poisson"          )
          else
             call displayMessage("Likelihood model: Negative binomial")
          end if
       else
          call    displayMessage("Likelihood model: normal"           )
       end if
       if (evaluationFailed) then
          call displayMessage("model evaluation failed - no likelihood computed")
       else
          call displayMessage("logℒ = "//trim(label))
          if (self%likelihoodPoisson) then
             call displayIndent("per bin likelihoods")
             do i=1,size(self%mass)
                write (label,'(i4)') i
                message=trim(label)//": "
                write (label,'(e17.10)') log10(self%mass                (i))
                message=message//"log₁₀(M/M☉) = " //trim(label)//"; "
                write (label,'(e17.10)')            likelihoodPerBin    (i)
                message=message//"logℒ = "        //trim(label)//"; "
                write (label,'(i07)'   )       self%countHalos          (i)
                message=message//"Nhalo(target) = "//trim(label)//"; "
                write (label,'(e17.10)')            countHalosMeanPerBin(i)
                message=message//"Nhalo(model) = " //trim(label)
                call displayMessage(message)
             end do
             call displayUnindent("done")
          end if
       end if
       call displayUnindent("done")
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
