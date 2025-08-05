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
Implements a generic 1D volume function (i.e. number density of objects binned by some property, e.g. a
mass function) output analysis class.
!!}

  use               :: Galactic_Filters                        , only : galacticFilterClass
  use   , intrinsic :: ISO_C_Binding                           , only : c_size_t
  use               :: ISO_Varying_String                      , only : varying_string
  use               :: Node_Property_Extractors                , only : nodePropertyExtractorClass
  !$ use            :: OMP_Lib                                 , only : omp_lock_kind
  use               :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerClass
  use               :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorClass
  use               :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorClass
  use               :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorClass
  use               :: Output_Times                            , only : outputTimesClass
  use               :: Output_Analyses_Options                 , only : enumerationOutputAnalysisCovarianceModelType
  !![
  <outputAnalysis name="outputAnalysisVolumeFunction1D">
   <description>
     A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.
  
     In addition to the volume function itself, the covariance matrix, $\mathbf{C}_\mathrm{model}$, of the mass function is also
     computed. The assumptions used when constructing the covariance matrix are controlled by the parameter {\normalfont
     \ttfamily [covarianceModel]}. If set to {\normalfont \ttfamily binomial}, then to construct $\mathbf{C}_\mathrm{model}$ we
     make use of the fact that \glc\ works by sampling a set of tree ``root masses'' from the $z=0$ dark matter halo mass
     function. From each root, a tree is grown, within which the physics of galaxy formation is then solved. Root masses are
     sampled uniformly from the halo mass function. That is, the cumulative halo mass function, $N(M)$, is constructed between
     the maximum and minimum halo masses to be simulated. The number of root masses, $N_\mathrm{r}$, to be used in a model
     evaluation is then determined. Root masses are then chosen such that
     \begin{equation}
      N(M_i) = N(M_\mathrm{min}) {i-1 \over N_\mathrm{r}-1}
     \end{equation}
     for $i=1\ldots N_\mathrm{r}$ (noting that $N(M_\mathrm{max})=0$ by construction). 
  
     Consider first those galaxies which form in the main branch of each tree (i.e. those galaxies which are destined to become
     the central galaxy of the $z=0$ halo). Suppose that we simulate $N_k$ halos of root mass $M_k$ at $z=0$. In such halos the
     main branch galaxies will, at any time, have property values drawn from some distribution $p_k(M_\star|t)$. The number of
     such galaxies contributing to bin $i$ of the mass function is therefore binomially distributed with success probability
     $p_{ik} = \int_{M_{i,\mathrm min}}^{M_{i,\mathrm max}} p_k(M_\star|t) \d M_\star$ and a sample size of $N_k$.

     Generalizing to consider all bins in our volume function, the number of galaxies in each bin will jointly follow a
     \href{https://en.wikipedia.org/wiki/Multinomial_distribution}{multinomial distribution}. The contribution to the covariance
     matrix from these main branch galaxies is therefore:     
     \begin{equation}
      \mathcal{C}_{ij} = \left\{ \begin{array}{ll} p_{ik}(1-p_{ik}) N_k w_k^2 &amp; \hbox{ if } i = j \\ -p_{ik} p_{jk} N_k w_k^2 &amp; \hbox{ otherwise,} \end{array} \right.
     \end{equation}     
     where $w_k$ is the weight to be assigned to each tree. To compute this covariance requires knowledge of the probabilities,
     $p_{ik}$. We estimate these directly from the model. To do this, we bin trees into narrow bins of root mass and assume that
     $p_{ik}$ does not vary significantly across the mass range of each bin. Using all realizations of trees that fall within a
     given bin, $k$, we can directly estimate $p_{ik}$. Similarly, $N_k w_k^2$ is found by accumulating squared weights in bins of
     root mass. In computing $p_{ik}$ and $N_k$, the range of halo masses considered and the fineness of binning in halo mass are
     determined by the parameters {\normalfont \ttfamily [covarianceBinomialMassHaloMinimum]}, {\normalfont \ttfamily
     [covarianceBinomialMassHaloMaximum]}, and {\normalfont \ttfamily [covarianceBinomialBinsPerDecade]}.
  
     If instead, {\normalfont \ttfamily [covarianceModel]}$=${\normalfont \ttfamily Poisson}, the main branch galaxies are
     modeled as being sampled from a Poisson distribution (and so off-diagonal terms in the covariance matrix will be zero).
  
     In addition to the main branch galaxies, each tree will contain a number of other galaxies (these will be ``satellite''
     galaxies at $z=0$, but at higher redshifts may still be central galaxies in their own halos). Tests have established that
     the number of satellites in halos is well described by a Poisson process. Note that, as described above, each galaxy
     contributes a Gaussian distribution to the mass function due to modeling of random errors in property value
     determinations. For main branch galaxies this is simply accounted for when accumulating the probabilities, $p_{ik}$. For
     satellite galaxies, off-diagonal contributions to the covariance matrix arise as a result, $C_{ij} = w_k f_i f_j$, where
     $f_i$ is the fraction of the galaxy contributing to bin $i$ of the mass function.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisVolumeFunction1D
     !!{
     A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output
     analysis class.
     !!}
     private
     type            (varying_string                              )                              :: label                                          , comment                                          , &
          &                                                                                         propertyLabel                                  , propertyComment                                  , &
          &                                                                                         distributionLabel                              , distributionComment                              , &
          &                                                                                         propertyUnits                                  , distributionUnits                                , &
          &                                                                                         xAxisLabel                                     , yAxisLabel                                       , &
          &                                                                                         targetLabel
     double precision                                                                            :: propertyUnitsInSI                              , distributionUnitsInSI
     class           (nodePropertyExtractorClass                  ), pointer                     :: nodePropertyExtractor_                => null()
     class           (outputAnalysisPropertyOperatorClass         ), pointer                     :: outputAnalysisPropertyOperator_       => null()                                                   , &
          &                                                                                         outputAnalysisPropertyUnoperator_     => null()
     class           (outputAnalysisWeightOperatorClass           ), pointer                     :: outputAnalysisWeightOperator_         => null()
     class           (outputAnalysisDistributionOperatorClass     ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (outputAnalysisDistributionNormalizerClass   ), pointer                     :: outputAnalysisDistributionNormalizer_ => null()
     class           (galacticFilterClass                         ), pointer                     :: galacticFilter_                       => null()
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                          => null()
     double precision                                              , dimension(:,:), allocatable :: outputWeight                                   , functionCovariance                               , &
          &                                                                                         weightMainBranch                               , functionCovarianceTarget
     double precision                                              , dimension(:  ), allocatable :: binCenter                                      , functionValue                                    , &
          &                                                                                         functionValueTarget                            , weightMainBranchSquared                          , &
          &                                                                                         functionCovarianceTarget1D
     double precision                                              , dimension(:  ), allocatable :: binMinimum                                     , binMaximum
     integer         (c_size_t                                    )                              :: binCount                                       , bufferCount                                      , &
          &                                                                                         binCountTotal                                  , covarianceModelBinomialBinCount
     type            (enumerationOutputAnalysisCovarianceModelType)                              :: covarianceModel
     integer                                                                                     :: covarianceBinomialBinsPerDecade
     double precision                                                                            :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum                , &
          &                                                                                         covarianceModelHaloMassMinimumLogarithmic      , covarianceModelHaloMassIntervalLogarithmicInverse, &
          &                                                                                         binWidth
     logical                                                                                     :: finalized                                      , xAxisIsLog                                       , &
          &                                                                                         yAxisIsLog                                     , likelihoodNormalize
     !$ integer      (omp_lock_kind                               )                              :: accumulateLock
   contains
     !![
     <methods>
       <method description="Return the results of the volume function operator." method="results"         />
       <method description="Finalize the analysis of this function."             method="finalizeAnalysis"/>
     </methods>
     !!]
     final     ::                     volumeFunction1DDestructor
     procedure :: analyze          => volumeFunction1DAnalyze
     procedure :: finalize         => volumeFunction1DFinalize
     procedure :: results          => volumeFunction1DResults
     procedure :: reduce           => volumeFunction1DReduce
     procedure :: logLikelihood    => volumeFunction1DLogLikelihood
     procedure :: finalizeAnalysis => volumeFunction1DFinalizeAnalysis
  end type outputAnalysisVolumeFunction1D

  interface outputAnalysisVolumeFunction1D
     !!{
     Constructors for the \refClass{outputAnalysisVolumeFunction1D} output analysis class.
     !!}
     module procedure volumeFunction1DConstructorParameters
     module procedure volumeFunction1DConstructorInternal
  end interface outputAnalysisVolumeFunction1D

contains

  function volumeFunction1DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisVolumeFunction1D} output analysis class which takes a parameter set as input.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter                                , inputParameters
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisCovarianceModelEncode
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass               ), pointer                     :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass        ), pointer                     :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    double precision                                           , dimension(:  ), allocatable :: binCenter                            , outputWeight                     , &
         &                                                                                      functionValueTarget                  , functionCovarianceTarget1D
    double precision                                           , dimension(:,:), allocatable :: functionCovarianceTarget
    integer         (c_size_t                                 )                              :: bufferCount
    type            (varying_string                           )                              :: label                                , comment                          , &
         &                                                                                      propertyLabel                        , propertyComment                  , &
         &                                                                                      distributionLabel                    , distributionComment              , &
         &                                                                                      propertyUnits                        , distributionUnits                , &
         &                                                                                      covarianceModel                      , targetLabel                      , &
         &                                                                                      xAxisLabel                           , yAxisLabel
    integer                                                                                  :: covarianceBinomialBinsPerDecade
    type            (inputParameters                          )                              :: unoperatorParameters
    double precision                                                                         :: propertyUnitsInSI                    , distributionUnitsInSI            , &
         &                                                                                      covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum, &
         &                                                                                      binWidth
    logical                                                                                  :: xAxisIsLog                           , yAxisIsLog                       , &
         &                                                                                      likelihoodNormalize

    ! Check and read parameters.
    unoperatorParameters=parameters%subParameters('unoperatorParameters',requireValue=.false.)
    !![
    <objectBuilder class="nodePropertyExtractor"                name="nodePropertyExtractor_"                source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"     source="unoperatorParameters"/>
    <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"         source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionNormalizer" name="outputAnalysisDistributionNormalizer_" source="parameters"          />
    <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"          />
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"          />
    !!]
    allocate(binCenter   (int(parameters%count('binCenter'))                          ))
    allocate(outputWeight(int(parameters%count('binCenter'))*self%outputTimes_%count()))
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*self%outputTimes_%count()) &
         & call Error_Report('incorrect number of output weights provided'//{introspection:location})
    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <variable>label</variable>
      <description>A label for the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>xAxisLabel</name>
      <source>parameters</source>
      <description>A label for the $x$-axis in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>yAxisLabel</name>
      <source>parameters</source>
      <description>A label for the $y$-axis in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>xAxisIsLog</name>
      <source>parameters</source>
      <description>If true, indicates that the $x$-axis should be logarithmic in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>yAxisIsLog</name>
      <source>parameters</source>
      <description>If true, indicates that the $y$-axis should be logarithmic in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>comment</name>
      <source>parameters</source>
      <variable>comment</variable>
      <description>A descriptive comment for the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyLabel</name>
      <source>parameters</source>
      <variable>propertyLabel</variable>
      <description>A label for the property variable.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyComment</name>
      <source>parameters</source>
      <variable>propertyComment</variable>
      <description>A descriptive comment for the property variable.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyUnits</name>
      <source>parameters</source>
      <variable>propertyUnits</variable>
      <description>A human-readable description of the units for the property.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyUnitsInSI</name>
      <source>parameters</source>
      <variable>propertyUnitsInSI</variable>
      <description>A units for the property in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>distributionLabel</name>
      <source>parameters</source>
      <variable>distributionLabel</variable>
      <description>A label for the distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>distributionComment</name>
      <source>parameters</source>
      <variable>distributionComment</variable>
      <description>A descriptive comment for the distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>distributionUnits</name>
      <source>parameters</source>
      <variable>distributionUnits</variable>
      <description>A human-readable description of the units for the distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>distributionUnitsInSI</name>
      <source>parameters</source>
      <variable>distributionUnitsInSI</variable>
      <description>A units for the distribution in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>The value of the property at the center of each bin.</description>
    </inputParameter>
    !!]
    if (size(binCenter) == 1) then
       !![
       <inputParameter>
	 <name>binWidth</name>
	 <source>parameters</source>
	 <variable>binWidth</variable>
	 <description>The width of the bins.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>bufferCount</name>
      <source>parameters</source>
      <variable>bufferCount</variable>
      <description>The number of buffer bins to include below and above the range of actual bins.</description>
    </inputParameter>
    <inputParameter>
      <name>outputWeight</name>
      <source>parameters</source>
      <variable>outputWeight</variable>
      <description>The weight to assign to each bin at each output.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceModel</name>
      <source>parameters</source>
      <variable>covarianceModel</variable>
      <description>The model to use for computing covariances.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true then normalize the likelihood to make it a probability density.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('functionValueTarget')) then
       if (parameters%isPresent('functionCovarianceTarget')) then
          !![
          <inputParameter>
            <name>functionValueTarget</name>
            <source>parameters</source>
            <description>The target function for likelihood calculations.</description>
          </inputParameter>
          <inputParameter>
            <name>functionCovarianceTarget</name>
            <source>parameters</source>
            <variable>functionCovarianceTarget1D</variable>
            <description>The target function covariance for likelihood calculations.</description>
          </inputParameter>
          !!]
          if (size(functionCovarianceTarget1D) == size(functionValueTarget)**2) then
             allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
             functionCovarianceTarget=reshape(functionCovarianceTarget1D,shape(functionCovarianceTarget))
          else
             call Error_Report('functionCovariance has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('functionCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('functionCovariance')) call Error_Report('functionTarget must be specified if functionCovariance is present'//{introspection:location})
    end if
    !![
    <inputParameter>
      <name>targetLabel</name>
      <source>parameters</source>
      <description>A label for the target dataset in a plot of this analysis.</description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    !!]
    ! Build the object.
    !![
    <conditionalCall>
     <call>
      self=outputAnalysisVolumeFunction1D(                                                                                                    &amp;
           &amp;                          label                                                                                             , &amp;
           &amp;                          comment                                                                                           , &amp;
           &amp;                          propertyLabel                                                                                     , &amp;
           &amp;                          propertyComment                                                                                   , &amp;
           &amp;                          propertyUnits                                                                                     , &amp;
           &amp;                          propertyUnitsInSI                                                                                 , &amp;
           &amp;                          distributionLabel                                                                                 , &amp;
           &amp;                          distributionComment                                                                               , &amp;
           &amp;                          distributionUnits                                                                                 , &amp;
           &amp;                          distributionUnitsInSI                                                                             , &amp;
           &amp;                          binCenter                                                                                         , &amp;
           &amp;                          bufferCount                                                                                       , &amp;
           &amp;                          reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),self%outputTimes_%count()]), &amp;
           &amp;                          nodePropertyExtractor_                                                                  , &amp;
           &amp;                          outputAnalysisPropertyOperator_                                                                   , &amp;
           &amp;                          outputAnalysisPropertyUnoperator_                                                                 , &amp;
           &amp;                          outputAnalysisWeightOperator_                                                                     , &amp;
           &amp;                          outputAnalysisDistributionOperator_                                                               , &amp;
           &amp;                          outputAnalysisDistributionNormalizer_                                                             , &amp;
           &amp;                          galacticFilter_                                                                                   , &amp;
           &amp;                          outputTimes_                                                                                      , &amp;
           &amp;                          enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.)      , &amp;
           &amp;                          covarianceBinomialBinsPerDecade                                                                   , &amp;
           &amp;                          covarianceBinomialMassHaloMinimum                                                                 , &amp;
           &amp;                          covarianceBinomialMassHaloMaximum                                                                 , &amp;
           &amp;                          likelihoodNormalize                                                                               , &amp;
           &amp;                          xAxisLabel                                                                                        , &amp;
           &amp;                          yAxisLabel                                                                                        , &amp;
           &amp;                          xAxisIsLog                                                                                        , &amp;
           &amp;                          yAxisIsLog                                                                                        , &amp;
           &amp;                          targetLabel                                                                                         &amp;
           &amp;                          {conditions}                                                                                        &amp;
           &amp;                         )
     </call>
     <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"          />
     <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"          />
     <argument name="binWidth"                 value="binWidth"                        condition="size(binCenter) == 1"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"               />
    <objectDestructor name="outputAnalysisPropertyOperator_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="outputAnalysisWeightOperator_"        />
    <objectDestructor name="outputAnalysisDistributionOperator_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="galacticFilter_"                      />
    <objectDestructor name="outputTimes_"                         />
    !!]
    return
  end function volumeFunction1DConstructorParameters

  function volumeFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,distributionLabel,distributionComment,distributionUnits,distributionUnitsInSI,binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,outputAnalysisDistributionNormalizer_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,xAxisLabel,yAxisLabel,xAxisIsLog,yAxisIsLog,targetLabel,functionValueTarget,functionCovarianceTarget,binWidth) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisVolumeFunction1D} output analysis class for internal use.
    !!}
    use    :: Error                   , only : Error_Report
    use    :: Node_Property_Extractors, only : nodePropertyExtractorClass           , nodePropertyExtractorScalar
    use    :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial
    !$ use :: OMP_Lib                 , only : OMP_Init_Lock
    implicit none
    type            (outputAnalysisVolumeFunction1D              )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                , comment                          , &
         &                                                                                                     propertyLabel                        , propertyComment                  , &
         &                                                                                                     distributionLabel                    , distributionComment              , &
         &                                                                                                     propertyUnits                        , distributionUnits
    type            (varying_string                              ), intent(in   ), optional                 :: xAxisLabel                           , yAxisLabel                       , &
         &                                                                                                     targetLabel
    logical                                                       , intent(in   ), optional                 :: xAxisIsLog                           , yAxisIsLog                       , &
         &                                                                                                     likelihoodNormalize
    double precision                                              , intent(in   )                           :: propertyUnitsInSI                    , distributionUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    class           (nodePropertyExtractorClass                  ), intent(in   ), target                   :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(in   ), target                   :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(in   ), target                   :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass   ), intent(in   ), target                   :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                         ), intent(in   ), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(in   ), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum, &
         &                                                                                                     binWidth
    double precision                                              , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                              , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    integer         (c_size_t                                    )                                          :: i
    !![
    <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, distributionLabel, distributionComment, distributionUnits, distributionUnitsInSI, binCenter, bufferCount, outputWeight, *nodePropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *outputAnalysisDistributionNormalizer_, *galacticFilter_, *outputTimes_, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, xAxisLabel='x', yAxisLabel='y', xAxisIsLog=.false., yAxisIsLog=.false., targetLabel, functionValueTarget, functionCovarianceTarget, binWidth"/>
    !!]

    ! Assign 1D version of the target covariance for use in the descriptor.
    if (present(functionCovarianceTarget)) self%functionCovarianceTarget1D=reshape(functionCovarianceTarget,[size(functionCovarianceTarget)])
    ! Validate.
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    ! Set normalization state for likelihood.
    self%likelihoodNormalize=.true.
    if (present(likelihoodNormalize)) self%likelihoodNormalize=likelihoodNormalize
    ! Count bins.
    self%binCount     =size(binCenter,kind=c_size_t)
    self%binCountTotal=self%binCount+2*bufferCount
    ! Determine bin minima and maxima. Allocate arrays for bins that include buffer regions.
    allocate(self%binMinimum(-bufferCount+1:-bufferCount+self%binCountTotal))
    allocate(self%binMaximum(-bufferCount+1:-bufferCount+self%binCountTotal))
    if (present(binWidth)) then
       self%binMinimum(1:self%binCount)=+binCenter-0.5d0*binWidth
       self%binMaximum(1:self%binCount)=+binCenter+0.5d0*binWidth
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*binWidth
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*binWidth
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*binWidth
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*binWidth
       end do
    else
       do i=1,self%binCount
          if (i == 1) then
             self%binMinimum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i+1))
          else
             self%binMinimum(i)=             +0.5d0*(self%binCenter(i-1)+self%binCenter(i  ))
          end if
          if (i == self%binCount) then
             self%binMaximum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i-1))
          else
             self%binMaximum(i)=             +0.5d0*(self%binCenter(i  )+self%binCenter(i+1))
          end if
       end do
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
       end do
    end if
    ! Allocate and initialize function values.
    allocate(self%functionValue     (self%binCount              ))
    allocate(self%functionCovariance(self%binCount,self%binCount))
    self%functionValue     =0.0d0
    self%functionCovariance=0.0d0
    ! Allocate and initialize binomial covariance model halo arrays if necessary.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       self%covarianceModelBinomialBinCount                  =int(log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)*dble(self%covarianceBinomialBinsPerDecade)+0.5d0)
       self%covarianceModelHaloMassMinimumLogarithmic        =log10(self%covarianceBinomialMassHaloMinimum)
       self%covarianceModelHaloMassIntervalLogarithmicInverse=dble(self%covarianceModelBinomialBinCount)/log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)
       allocate(self%weightMainBranch       (self%binCount,self%covarianceModelBinomialBinCount))
       allocate(self%weightMainBranchSquared(self%covarianceModelBinomialBinCount))
       self%weightMainBranch       =0.0d0
       self%weightMainBranchSquared=0.0d0
    end if
    ! Initialize finalization status.
    self%finalized=.false.
    ! Initialize OpenMP accumulation lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
   return
  end function volumeFunction1DConstructorInternal

  subroutine volumeFunction1DDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisVolumeFunction1D} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisVolumeFunction1D), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"               />
    <objectDestructor name="self%outputAnalysisPropertyOperator_"      />
    <objectDestructor name="self%outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="self%outputAnalysisWeightOperator_"        />
    <objectDestructor name="self%outputAnalysisDistributionOperator_"  />
    <objectDestructor name="self%outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="self%galacticFilter_"                      />
    <objectDestructor name="self%outputTimes_"                         />
    !!]
    ! Destroy OpenMP lock.
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine volumeFunction1DDestructor

  subroutine volumeFunction1DAnalyze(self,node,iOutput)
    !!{
    Implement a volumeFunction1D output analysis.
    !!}
    use    :: Galacticus_Nodes        , only : nodeComponentBasic                   , treeNode
    use    :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    use    :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial, enumerationOutputAnalysisPropertyTypeType, enumerationOutputAnalysisPropertyQuantityType
    !$ use :: OMP_Lib                 , only : OMP_Set_Lock                         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisVolumeFunction1D               ), intent(inout)                 :: self
    type            (treeNode                                     ), intent(inout)                 :: node
    integer         (c_size_t                                     ), intent(in   )                 :: iOutput
    double precision                                               , allocatable  , dimension(:  ) :: distribution
    double precision                                               , allocatable  , dimension(:,:) :: covariance
    class           (nodeComponentBasic                           ), pointer                       :: basic
    double precision                                                                               :: propertyValue         , weightValue  , &
         &                                                                                            propertyValueIntrinsic
    type            (enumerationOutputAnalysisPropertyTypeType    )                                :: propertyType 
    type            (enumerationOutputAnalysisPropertyQuantityType)                                :: propertyQuantity
    integer         (c_size_t                                     )                                :: j                     , k            , &
         &                                                                                            indexHaloMass

    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    ! Filter this node.
    if (.not.self%galacticFilter_%passes(node)) return
    ! Allocate work arrays.
    allocate(distribution(-self%bufferCount+1:self%binCount+self%bufferCount))
    ! Extract the property from the node.
    propertyType             =self%nodePropertyExtractor_%type    (    )
    propertyQuantity         =self%nodePropertyExtractor_%quantity(    )
    select type (extractor_ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       propertyValue         =                           extractor_%extract (node)
    class default
       propertyValue         =0.0d0
    end select
    propertyValueIntrinsic   =propertyValue
    ! Apply property operators.
    propertyValue=self%outputAnalysisPropertyOperator_%operate(propertyValue,node,propertyType,iOutput)
    ! Apply distribution operators.
    distribution=self%outputAnalysisDistributionOperator_%operateScalar(propertyValue,propertyType,self%binMinimum,self%binMaximum,iOutput,node)
    ! Compute the weight.
    weightValue=node%hostTree%volumeWeight
    ! Apply weight operators.
    weightValue=self%outputAnalysisWeightOperator_%operate(weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,iOutput)
    ! Apply weights.
    distribution(1:self%binCount)=+distribution              (1:self%binCount        ) &
         &                        *self        %outputWeight ( :             ,iOutput) &
         &                        *weightValue
    ! Accumulate the property, including weights from both the host tree and the output. Note that we accumulate only the
    ! non-buffer bins of the distribution.
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%functionValue=+self%functionValue+distribution(1:self%binCount)
    !$ call OMP_Unset_Lock(self%accumulateLock)
    ! Accumulate covariance. If using the binomial model for main branch galaxies, handle them separately.
    if (node%isOnMainBranch() .and. self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       ! Find the bin to which this halo mass belongs.
       basic         => node%basic()
       indexHaloMass =  floor((log10(basic%mass())-self%covarianceModelHaloMassMinimumLogarithmic)*self%covarianceModelHaloMassIntervalLogarithmicInverse)+1
       ! Accumulate weights to halo mass arrays.
       if (indexHaloMass >= 1 .and. indexHaloMass <= self%covarianceModelBinomialBinCount) then
          self                     %weightMainBranch       ( :             ,indexHaloMass)=    &
               &  +    self        %weightMainBranch       ( :             ,indexHaloMass)     &
               &  +    distribution                        (1:self%binCount)
          self                     %weightMainBranchSquared(                indexHaloMass)=    &
               &  +    self        %weightMainBranchSquared(                indexHaloMass)     &
               &  +sum(distribution                        (1:self%binCount              ))**2
       end if
    else
       ! Construct contribution to the covariance matrix assuming Poisson statistics.
       allocate(covariance(self%binCount,self%binCount))
       forall(j=1:self%binCount)
          forall(k=j:self%binCount)
             covariance(j,k)=+distribution(j  ) &
                  &          *distribution(  k)
             covariance(k,j)=+covariance  (j,k)
          end forall
       end forall
       ! Accumulate covariance.
       !$ call OMP_Set_Lock(self%accumulateLock)
       self        %functionCovariance= &
            & +self%functionCovariance  &
            & +             covariance
       !$ call OMP_Unset_Lock(self%accumulateLock)
       deallocate(covariance)
    end if
    ! Deallocate workspace.
    deallocate(distribution)
    return
  end subroutine volumeFunction1DAnalyze

  subroutine volumeFunction1DReduce(self,reduced)
    !!{
    Implement a volumeFunction1D output analysis reduction.
    !!}
    use    :: Error                  , only : Error_Report
    use    :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    !$ use :: OMP_Lib                , only : OMP_Set_Lock                         , OMP_Unset_Lock
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout) :: self
    class(outputAnalysisClass           ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisVolumeFunction1D)
       !$ call OMP_Set_Lock(reduced%accumulateLock)
       if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
          reduced%weightMainBranch       =reduced%weightMainBranch       +self%weightMainBranch
          reduced%weightMainBranchSquared=reduced%weightMainBranchSquared+self%weightMainBranchSquared
       end if
       reduced%functionCovariance        =reduced%functionCovariance     +self%functionCovariance
       reduced%functionValue             =reduced%functionValue          +self%functionValue
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine volumeFunction1DReduce

  subroutine volumeFunction1DFinalize(self,groupName)
    !!{
    Implement a volumeFunction1D output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout)           :: self
    type (varying_string                ), intent(in   ), optional :: groupName
    type (hdf5Object                    )               , target   :: analysesGroup, subGroup
    type (hdf5Object                    )               , pointer  :: inGroup
    type (hdf5Object                    )                          :: analysisGroup, dataset

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Output.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup =  inGroup      %openGroup(char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call    analysisGroup%writeAttribute(     char(self%       comment   )                       ,'description'                                                                                                          )
    call    analysisGroup%writeAttribute("function1D"                                            ,'type'                                                                                                                 )
    call    analysisGroup%writeAttribute(     char(self%       xAxisLabel)                       ,'xAxisLabel'                                                                                                           )
    call    analysisGroup%writeAttribute(     char(self%       yAxisLabel)                       ,'yAxisLabel'                                                                                                           )
    call    analysisGroup%writeAttribute(          self%       xAxisIsLog                        ,'xAxisIsLog'                                                                                                           )
    call    analysisGroup%writeAttribute(          self%       yAxisIsLog                        ,'yAxisIsLog'                                                                                                           )
    call    analysisGroup%writeAttribute(     char(self%    propertyLabel)                       ,'xDataset'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%distributionLabel)                       ,'yDataset'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%distributionLabel)//"Covariance"         ,'yCovariance'                                                                                                          )
    if (allocated(self%functionValueTarget)) then
       call analysisGroup%writeAttribute(     char(self%distributionLabel)//"Target"             ,'yDatasetTarget'                                                                                                       )
       call analysisGroup%writeAttribute(     char(self%distributionLabel)//"CovarianceTarget"   ,'yCovarianceTarget'                                                                                                    )
    end if
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (self%binCenter    (1:self%binCount                     ),char(self%    propertyLabel)                   ,char(self%   propertyComment)                  ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%    propertyUnits    )                   ,'units'                                                                                                                )
    call    dataset      %writeAttribute(          self%    propertyUnitsInSI                    ,'unitsInSI'                                                                                                            )
    call    analysisGroup%writeDataset  (self%functionValue(1:self%binCount                     ),char(self%distributionLabel)                   ,char(self%distributionComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%distributionUnits    )                   ,'units'                                                                                                                )
    call    dataset      %writeAttribute(          self%distributionUnitsInSI                    ,'unitsInSI'                                                                                                            )
    call    analysisGroup%writeDataset  (self%functionCovariance(1:self%binCount,1:self%binCount),char(self%distributionLabel)//"Covariance"     ,char(self%distributionComment)//" [covariance]",datasetReturned=dataset)
    call    dataset      %writeAttribute("["//char(self%distributionUnits    )//"]²"             ,'units'                                                                                                                )
    call    dataset      %writeAttribute(          self%distributionUnitsInSI   **2              ,'unitsInSI'                                                                                                            )
    ! If available, include the log-likelihood and target dataset.
    if (allocated(self%functionValueTarget)) then
       call analysisGroup%writeAttribute(          self%logLikelihood()                         ,'logLikelihood'                                                                                                         )
       call analysisGroup%writeAttribute(     char(self%targetLabel          )                  ,'targetLabel'                                                                                                           )
       call analysisGroup%writeDataset  (          self%functionValueTarget                     ,char(self%distributionLabel)//"Target"          ,char(self%distributionComment)                 ,datasetReturned=dataset)
       call dataset      %writeAttribute(     char(self%distributionUnits    )                  ,'units'                                                                                                                 )
       call dataset      %writeAttribute(          self%distributionUnitsInSI                   ,'unitsInSI'                                                                                                             )
       call analysisGroup%writeDataset  (          self%functionCovarianceTarget                ,char(self%distributionLabel)//"CovarianceTarget",char(self%distributionComment)//" [covariance]",datasetReturned=dataset)
       call dataset      %writeAttribute("["//char(self%distributionUnits    )//"]²"            ,'units'                                                                                                                 )
       call dataset      %writeAttribute(          self%distributionUnitsInSI   **2             ,'unitsInSI'                                                                                                             )
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine volumeFunction1DFinalize

  subroutine volumeFunction1DFinalizeAnalysis(self)
    !!{
    Compute final covariances and normalize.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities          , only : mpiSelf
#endif
    use :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout) :: self
    integer         (c_size_t                      )                :: i                    , j, &
         &                                                             m
    double precision                                                :: weightMainBranchTotal

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! If using the binomial covariance model, add on the contribution to covariance from main branch galaxies.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       do m=1,self%covarianceModelBinomialBinCount
          weightMainBranchTotal=sum(self%weightMainBranch(:,m))
          if (weightMainBranchTotal > 0.0d0) then
             do i=1,self%binCount
                ! For on-diagonal terms the variance is just n pᵢ (1-pᵢ). Here, pᵢ is the probability of contributing to bin i in
                ! our function. We find this by dividing the weight of main branch galaxies in bin i by the sum of weights over
                ! all bins. Finally, n is just the sum of squared weights for this halo mass bin.
                self               %functionCovariance     (i,i)=                       &
                     &         self%functionCovariance     (i,i)                        &
                     & +       self%weightMainBranchSquared(  m)                        & ! +   n
                     & *       self%weightMainBranch       (i,m)/weightMainBranchTotal  & !     pᵢ
                     & *(1.0d0-self%weightMainBranch       (i,m)/weightMainBranchTotal)   !  (1-pᵢ)
                do j=1,self%binCount
                   if (i == j) cycle
                   ! For off-diagonal terms the covariance is - n pᵢ pⱼ.
                   self            %functionCovariance     (i,j)=                       &
                        &  +   self%functionCovariance     (i,j)                        &
                        &  -   self%weightMainBranchSquared(  m)                        & ! -   n
                        &  *   self%weightMainBranch       (i,m)/weightMainBranchTotal  & !     pᵢ
                        &  *   self%weightMainBranch       (j,m)/weightMainBranchTotal    !     pⱼ
                end do
             end do
          end if
       end do
    end if
#ifdef USEMPI
    ! If running under MPI, perform a summation reduction across all processes.
    self%functionValue     =mpiSelf%sum(self%functionValue     )
    self%functionCovariance=mpiSelf%sum(self%functionCovariance)
#endif
    ! Apply final distribution operators - pass only the non-buffer bin values here.
    call self%outputAnalysisDistributionNormalizer_%normalize(self%functionValue,self%functionCovariance,self%binMinimum(1:self%binCount),self%binMaximum(1:self%binCount))
    ! Apply any "unoperator" to output bin values. This can be used to reverse transformations (e.g. if masses were converted to
    ! log10 for analysis, that can be undone here).
    do i=1,self%binCount
       self%binCenter(i)=self%outputAnalysisPropertyUnoperator_%operate(self%binCenter(i))
    end do
    return
  end subroutine volumeFunction1DFinalizeAnalysis

  subroutine volumeFunction1DResults(self,binCenter,functionValue,functionCovariance)
    !!{
    Implement a volumeFunction1D output analysis finalization.
    !!}
    implicit none
    class           (outputAnalysisVolumeFunction1D)                             , intent(inout)           :: self
    double precision                                , allocatable, dimension(:  ), intent(inout), optional :: binCenter         , functionValue
    double precision                                , allocatable, dimension(:,:), intent(inout), optional :: functionCovariance

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Return results.
    if (present(binCenter         )) then
       if (allocated(binCenter         )) deallocate(binCenter         )
       allocate(binCenter(size(self%binCenter)))
       binCenter         =self%binCenter
    end if
    if (present(functionValue     )) then
       if (allocated(functionValue     )) deallocate(functionValue     )
       allocate(functionValue(size(self%functionValue)))
       functionValue     =self%functionValue
    end if
    if (present(functionCovariance)) then
       if (allocated(functionCovariance)) deallocate(functionCovariance)
        allocate(functionCovariance(size(self%functionCovariance,dim=1),size(self%functionCovariance,dim=2)))
       !!]
       functionCovariance=self%functionCovariance
    end if
    return
  end subroutine volumeFunction1DResults

  double precision function volumeFunction1DLogLikelihood(self)
    !!{
    Return the log-likelihood of a volumeFunction1D output analysis.
    !!}
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Models_Likelihoods_Constants, only : logImprobable
    use :: Error                       , only : Error_Report
    use :: Interface_GSL               , only : GSL_Success
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout)                 :: self
    double precision                                , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                , allocatable  , dimension(:  ) :: functionValueDifference
    type            (vector                        )                                :: residual
    type            (matrix                        )                                :: covariance
    integer                                                                         :: status

    ! Check for existence of a target distribution.
    if (allocated(self%functionValueTarget)) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! If model has everywhere zero return an improbable likelihood.
       if (all(self%functionValue == 0.0d0) .and. all(self%functionCovariance == 0.0d0)) then
          volumeFunction1DLogLikelihood=logImprobable
       else
          ! Allocate workspaces.
          allocate(functionCovarianceCombined(self%binCount,self%binCount))
          allocate(functionValueDifference   (self%binCount              ))
          ! Find combined covariance and difference between model and target.
          functionValueDifference   =+self%functionValue            &
               &                     -self%functionValueTarget
          functionCovarianceCombined=+self%functionCovariance       &
               &                     +self%functionCovarianceTarget
          residual                  = vector(functionValueDifference   )
          covariance                = matrix(functionCovarianceCombined)
          ! Compute the log-likelihood.
          volumeFunction1DLogLikelihood          =-0.5d0*covariance%covarianceProduct(residual,status)
          if (status == GSL_Success) then
             if (self%likelihoodNormalize)                                                   &
                  & volumeFunction1DLogLikelihood=+volumeFunction1DLogLikelihood             &
                  &                               -0.5d0*covariance%logarithmicDeterminant() &
                  &                               -0.5d0*dble(self%binCount)                 &
                  &                               *log(2.0d0*Pi)
          else
             volumeFunction1DLogLikelihood       =+logImprobable
          end if
       end if
    else
       volumeFunction1DLogLikelihood=0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function volumeFunction1DLogLikelihood
