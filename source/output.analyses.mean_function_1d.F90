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
  Implements a generic 1D mean function (i.e. mean value of some property weighted by number density of
  objects binned by some property) output analysis class.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <outputAnalysis name="outputAnalysisMeanFunction1D">
   <description>A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.</description>
   <deepCopy>
    <functionClass variables="volumeFunctionUnweighted, volumeFunctionWeighted, crossCovariance"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="volumeFunctionUnweighted, volumeFunctionWeighted, crossCovariance"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisMeanFunction1D
     !!{
     A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.
     !!}
     private
     type            (varying_string                              )                              :: label                                          , comment                                         , &
          &                                                                                         propertyLabel                                  , propertyComment                                 , &
          &                                                                                         meanLabel                                      , meanComment                                     , &
          &                                                                                         propertyUnits                                  , meanUnits                                       , &
          &                                                                                         xAxisLabel                                     , yAxisLabel                                      , &
          &                                                                                         targetLabel
     type            (enumerationOutputAnalysisCovarianceModelType)                              :: covarianceModel
     double precision                                                                            :: propertyUnitsInSI                              , meanUnitsInSI
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                          => null()
     class           (nodePropertyExtractorClass                  ), pointer                     :: nodePropertyExtractor_                => null(), outputAnalysisWeightPropertyExtractor_ => null()
     class           (outputAnalysisPropertyOperatorClass         ), pointer                     :: outputAnalysisPropertyOperator_       => null(), outputAnalysisPropertyUnoperator_      => null(), &
          &                                                                                         outputAnalysisWeightPropertyOperator_ => null()
     class           (outputAnalysisWeightOperatorClass           ), pointer                     :: outputAnalysisWeightOperator_         => null()
     class           (outputAnalysisDistributionOperatorClass     ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (galacticFilterClass                         ), pointer                     :: galacticFilter_                       => null()
     type            (outputAnalysisVolumeFunction1D              ), pointer                     :: volumeFunctionUnweighted              => null(), volumeFunctionWeighted                  => null()
     type            (outputAnalysisCrossCorrelator1D             ), pointer                     :: crossCovariance                       => null()
     double precision                                              , allocatable, dimension(:  ) :: binCenter                                      , meanValue                                        , &
          &                                                                                         meanValueTarget                                , meanCovarianceTarget1D                           , &
          &                                                                                         outputWeight
     double precision                                              , allocatable, dimension(:,:) :: meanCovariance                                 , meanCovarianceTarget
     double precision                                                                            :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum                , &
          &                                                                                         binWidth
     logical                                                                                     :: finalized                                      , likelihoodNormalize                              , &
          &                                                                                         xAxisIsLog                                     , yAxisIsLog
     integer                                                                                     :: covarianceBinomialBinsPerDecade
     integer         (c_size_t                                    )                              :: bufferCount
   contains
     !![
     <methods>
       <method description="Return the results of the mean function operator." method="results"         />
       <method description="Finalize analysis of the mean function operator."  method="finalizeAnalysis"/>
     </methods>
     !!]
     final     ::                     meanFunction1DDestructor
     procedure :: analyze          => meanFunction1DAnalyze
     procedure :: finalize         => meanFunction1DFinalize
     procedure :: finalizeAnalysis => meanFunction1DFinalizeAnalysis
     procedure :: reduce           => meanFunction1DReduce
     procedure :: results          => meanFunction1DResults
     procedure :: logLikelihood    => meanFunction1DLogLikelihood
  end type outputAnalysisMeanFunction1D

  interface outputAnalysisMeanFunction1D
     !!{
     Constructors for the \refClass{outputAnalysisMeanFunction1D} output analysis class.
     !!}
     module procedure meanFunction1DConstructorParameters
     module procedure meanFunction1DConstructorInternal
  end interface outputAnalysisMeanFunction1D

contains

  function meanFunction1DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisMeanFunction1D} output analysis class which takes a parameter set as input.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter                                , inputParameters
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisCovarianceModelEncode
    implicit none
    type            (outputAnalysisMeanFunction1D           )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass             ), pointer                     :: nodePropertyExtractor_               , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_     , &
         &                                                                                    outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass      ), pointer                     :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                    ), pointer                     :: galacticFilter_
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    double precision                                         , dimension(:  ), allocatable :: binCenter                            , outputWeight                          , &
         &                                                                                    meanValueTarget                      , meanCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: meanCovarianceTarget
    integer         (c_size_t                               )                              :: bufferCount
    type            (varying_string                         )                              :: label                                , comment                               , &
         &                                                                                    propertyLabel                        , propertyComment                       , &
         &                                                                                    meanLabel                            , meanComment                           , &
         &                                                                                    propertyUnits                        , meanUnits                             , &
         &                                                                                    covarianceModel                      , xAxisLabel                            , &
         &                                                                                    yAxisLabel                           , targetLabel
    integer                                                                                :: covarianceBinomialBinsPerDecade
    type            (inputParameters                        )                              :: unoperatorParameters
    type            (inputParameters                        )                              :: weightParameters
    double precision                                                                       :: propertyUnitsInSI                    , meanUnitsInSI                         , &
         &                                                                                    covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum     , &
         &                                                                                    binWidth
    logical                                                                                :: likelihoodNormalize                  , xAxisIsLog                            , &
         &                                                                                    yAxisIsLog

    !![
    <objectBuilder class="nodePropertyExtractor"                name="nodePropertyExtractor_"                 source="parameters"          />
    <objectBuilder class="nodePropertyExtractor"                name="outputAnalysisWeightPropertyExtractor_" source="weightParameters"    />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"        source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisWeightPropertyOperator_"  source="weightParameters"    />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"      source="unoperatorParameters"/>
    <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"          source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"    source="parameters"          />
    <objectBuilder class="galacticFilter"                       name="galacticFilter_"                        source="parameters"          />
    <objectBuilder class="outputTimes"                          name="outputTimes_"                           source="parameters"          />
    !!]
    unoperatorParameters=parameters%subParameters('unoperator',requireValue=.false.)
    weightParameters    =parameters%subParameters('weight'    ,requireValue=.false.)
    allocate(binCenter   (int(parameters%count('binCenter'))                     ))
    allocate(outputWeight(int(parameters%count('binCenter'))*outputTimes_%count()))
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*outputTimes_%count()) &
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
      <name>meanLabel</name>
      <source>parameters</source>
      <variable>meanLabel</variable>
      <description>A label for the mean.</description>
    </inputParameter>
    <inputParameter>
      <name>meanComment</name>
      <source>parameters</source>
      <variable>meanComment</variable>
      <description>A descriptive comment for the mean.</description>
    </inputParameter>
    <inputParameter>
      <name>meanUnits</name>
      <source>parameters</source>
      <variable>meanUnits</variable>
      <description>A human-readable description of the units for the mean.</description>
    </inputParameter>
    <inputParameter>
      <name>meanUnitsInSI</name>
      <source>parameters</source>
      <variable>meanUnitsInSI</variable>
      <description>A units for the mean in the SI system.</description>
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
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>The value of the property at the center of each bin.</description>
    </inputParameter>
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
   if (parameters%isPresent('meanValueTarget')) then
       if (parameters%isPresent('meanCovarianceTarget')) then
          !![
          <inputParameter>
            <name>meanValueTarget</name>
            <source>parameters</source>
            <description>The target function for likelihood calculations.</description>
          </inputParameter>
          <inputParameter>
            <name>meanCovarianceTarget</name>
            <source>parameters</source>
            <variable>meanCovarianceTarget1D</variable>
            <description>The target function covariance for likelihood calculations.</description>
          </inputParameter>
          !!]
          if (size(meanCovarianceTarget1D) == size(meanValueTarget)**2) then
             allocate(meanCovarianceTarget(size(meanValueTarget),size(meanValueTarget)))
             meanCovarianceTarget=reshape(meanCovarianceTarget1D,shape(meanCovarianceTarget))
          else
             call Error_Report('meanCovariance has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('meanCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('meanCovariance')) call Error_Report('functionTarget must be specified if meanCovariance is present'//{introspection:location})
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
      self=outputAnalysisMeanFunction1D(                                                                                               &amp;
           &amp;                        label                                                                                        , &amp;
           &amp;                        comment                                                                                      , &amp;
           &amp;                        propertyLabel                                                                                , &amp;
           &amp;                        propertyComment                                                                              , &amp;
           &amp;                        propertyUnits                                                                                , &amp;
           &amp;                        propertyUnitsInSI                                                                            , &amp;
           &amp;                        meanLabel                                                                                    , &amp;
           &amp;                        meanComment                                                                                  , &amp;
           &amp;                        meanUnits                                                                                    , &amp;
           &amp;                        meanUnitsInSI                                                                                , &amp;
           &amp;                        binCenter                                                                                    , &amp;
           &amp;                        bufferCount                                                                                  , &amp;
           &amp;                        reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),outputTimes_%count()]), &amp;
           &amp;                        nodePropertyExtractor_                                                                       , &amp;
           &amp;                        outputAnalysisWeightPropertyExtractor_                                                       , &amp;
           &amp;                        outputAnalysisPropertyOperator_                                                              , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                                                        , &amp;
           &amp;                        outputAnalysisPropertyUnoperator_                                                            , &amp;
           &amp;                        outputAnalysisWeightOperator_                                                                , &amp;
           &amp;                        outputAnalysisDistributionOperator_                                                          , &amp;
           &amp;                        galacticFilter_                                                                              , &amp;
           &amp;                        outputTimes_                                                                                 , &amp;
           &amp;                        enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.) , &amp;
           &amp;                        covarianceBinomialBinsPerDecade                                                              , &amp;
           &amp;                        covarianceBinomialMassHaloMinimum                                                            , &amp;
           &amp;                        covarianceBinomialMassHaloMaximum                                                            , &amp;
           &amp;                        likelihoodNormalize                                                                          , &amp;
           &amp;                        xAxisLabel                                                                                   , &amp;
           &amp;                        yAxisLabel                                                                                   , &amp;
           &amp;                        xAxisIsLog                                                                                   , &amp;
           &amp;                        yAxisIsLog                                                                                   , &amp;
           &amp;                        targetLabel                                                                                    &amp;
           &amp;                        {conditions}                                                                                   &amp;
           &amp;                       )
     </call>
     <argument name="meanValueTarget"      value="meanValueTarget"      parameterPresent="parameters"          />
     <argument name="meanCovarianceTarget" value="meanCovarianceTarget" parameterPresent="parameters"          />
     <argument name="binWidth"             value="binWidth"                    condition="size(binCenter) == 1"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"      />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"/>
    <objectDestructor name="outputAnalysisPropertyOperator_"       />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_" />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="outputAnalysisWeightOperator_"         />
    <objectDestructor name="outputAnalysisDistributionOperator_"   />
    <objectDestructor name="galacticFilter_"                       />
    <objectDestructor name="outputTimes_"                          />
    !!]
    return
  end function meanFunction1DConstructorParameters

  function meanFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,meanLabel,meanComment,meanUnits,meanUnitsInSI,binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisWeightPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisWeightPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,xAxisLabel,yAxisLabel,xAxisIsLog,yAxisIsLog,targetLabel,meanValueTarget,meanCovarianceTarget,binWidth) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMeanFunction1D} output analysis class for internal use.
    !!}
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorBoolean       , outputAnalysisPropertyOperatorClass , outputAnalysisPropertyOperatorSequence, propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorClass           , outputAnalysisWeightOperatorProperty, outputAnalysisWeightOperatorSequence  , weightOperatorList
    implicit none
    type            (outputAnalysisMeanFunction1D                )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                    , comment                                , &
         &                                                                                                     propertyLabel                            , propertyComment                        , &
         &                                                                                                     meanLabel                                , meanComment                            , &
         &                                                                                                     propertyUnits                            , meanUnits
    type            (varying_string                              ), intent(in   ), optional                 :: xAxisLabel                               , yAxisLabel                             , &
         &                                                                                                     targetLabel
    double precision                                              , intent(in   )                           :: propertyUnitsInSI                        , meanUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    logical                                                       , intent(in   ), optional                 :: xAxisIsLog                               , yAxisIsLog                             , &
         &                                                                                                     likelihoodNormalize
    double precision                                              , intent(in   ), optional                 :: binWidth
    class           (nodePropertyExtractorClass                  ), intent(inout), target                   :: nodePropertyExtractor_                   , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(inout), target                   :: outputAnalysisPropertyOperator_          , outputAnalysisPropertyUnoperator_      , &
         &                                                                                                     outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(inout), target                   :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(inout), target                   :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                         ), intent(inout), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(inout), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum        , covarianceBinomialMassHaloMaximum
    double precision                                              , intent(in   ), optional, dimension(:  ) :: meanValueTarget
    double precision                                              , intent(in   ), optional, dimension(:,:) :: meanCovarianceTarget
    type            (outputAnalysisDistributionNormalizerIdentity), pointer                                 :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisWeightOperatorSequence        ), pointer                                 :: outputAnalysisWeightOperatorWeighted_    , outputAnalysisWeightOperatorUnweighted_
    type            (weightOperatorList                          ), pointer                                 :: weightOperatorWeight_                    , weightOperatorUnweighted_
    type            (outputAnalysisPropertyOperatorSequence      ), pointer                                 :: weightOperatorUnweightedPropertyOperator_
    type            (propertyOperatorList                        ), pointer                                 :: propertyOperator_                        , propertyOperatorUnweighted_
    type            (outputAnalysisWeightOperatorProperty        ), pointer                                 :: weightOperatorUnweightedProperty_        , weightOperatorWeightProperty_
    type            (outputAnalysisPropertyOperatorBoolean       ), pointer                                 :: propertyOperatorUnweightedBoolean_
    !![
    <constructorAssign variables="binWidth, bufferCount, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, meanLabel, meanComment, meanUnits, meanUnitsInSI, xAxisLabel, yAxisLabel, xAxisIsLog, yAxisIsLog, targetLabel, meanValueTarget, meanCovarianceTarget, *nodePropertyExtractor_, *outputAnalysisWeightPropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisWeightPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *galacticFilter_, *outputTimes_"/>
    !!]

    ! Set 1D covariance for descriptor.
    if (present(meanCovarianceTarget)) &
         & self%meanCovarianceTarget1D=reshape(meanCovarianceTarget,[size(meanCovarianceTarget)])
    self       %outputWeight          =reshape(outputWeight        ,[size(outputWeight        )])
    ! Mark as unfinalized.
    self%finalized=.false.
    ! Set normalization state for likelihood.
    self%likelihoodNormalize=.true.
    if (present(likelihoodNormalize)) self%likelihoodNormalize=likelihoodNormalize
    ! Build an identity distribution normalizer.
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_" constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    !!]
    ! Build weight operator that includes the weight property.
    allocate(outputAnalysisWeightOperatorWeighted_     )
    allocate(weightOperatorWeight_                     )
    allocate(weightOperatorWeight_                %next)
    allocate(weightOperatorWeightProperty_             )
    !![
    <referenceConstruct object="weightOperatorWeightProperty_" constructor="outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractor_,outputAnalysisWeightPropertyOperator_)"/>
    !!]
    weightOperatorWeight_     %operator_ => outputAnalysisWeightOperator_
    weightOperatorWeight_%next%operator_ => weightOperatorWeightProperty_
    !![
    <referenceConstruct object="outputAnalysisWeightOperatorWeighted_" constructor="outputAnalysisWeightOperatorSequence(weightOperatorWeight_)"/>
    !!]
    ! Build weight operator that includes a final boolean operator. We use this in the unweighted case to allow filters to be
    ! applied to the weight property. If any filter-type property operator sets the property value to zero, this boolean operator
    ! will return zero, otherwise it will return unity.
    allocate(outputAnalysisWeightOperatorUnweighted_       )
    allocate(weightOperatorUnweighted_                     )
    allocate(weightOperatorUnweighted_                %next)
    allocate(weightOperatorUnweightedProperty_             )
    allocate(weightOperatorUnweightedPropertyOperator_     )
    allocate(propertyOperatorUnweighted_                   )
    allocate(propertyOperatorUnweighted_              %next)
    allocate(propertyOperatorUnweightedBoolean_            )
    ! Build a sequence property operator for the unweight property. This includes any property operator passed to this
    ! constructor, plus the boolean operator.
    !![
    <referenceConstruct object="propertyOperatorUnweightedBoolean_" constructor="outputAnalysisPropertyOperatorBoolean(preciseZero=.true.)"/>
    !!]
    propertyOperatorUnweighted_     %operator_ => outputAnalysisWeightPropertyOperator_ ! First operator in the sequence is whatever operator was passed to this constructor.
    propertyOperatorUnweighted_%next%operator_ => propertyOperatorUnweightedBoolean_
    !![
    <referenceConstruct object="weightOperatorUnweightedPropertyOperator_" constructor="outputAnalysisPropertyOperatorSequence(propertyOperatorUnweighted_)"/>
    !!]
    ! Build a weight operator sequence - the first operator is whatever weight operator was passed to this constructor, the second
    ! is the operator which weights by the property value.
    !![
    <referenceConstruct object="weightOperatorUnweightedProperty_" constructor="outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractor_,weightOperatorUnweightedPropertyOperator_)"/>
    !!]
    weightOperatorUnweighted_     %operator_ => outputAnalysisWeightOperator_
    weightOperatorUnweighted_%next%operator_ => weightOperatorUnweightedProperty_
    !![
    <referenceConstruct object="outputAnalysisWeightOperatorUnweighted_" constructor="outputAnalysisWeightOperatorSequence(weightOperatorUnweighted_)"/>
    !!]
    ! Build weighted, unweighted, and cross volume function objects.
    allocate(self%volumeFunctionUnweighted)
    allocate(self%volumeFunctionWeighted  )
    allocate(self%crossCovariance         )
    !![
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionUnweighted">
     <constructor>
      outputAnalysisVolumeFunction1D (                                                  &amp;
       &amp;                                   label                                  , &amp;
       &amp;                                   comment                                , &amp;
       &amp;                                   propertyLabel                          , &amp;
       &amp;                                   propertyComment                        , &amp;
       &amp;                                   propertyUnits                          , &amp;
       &amp;                                   propertyUnitsInSI                      , &amp;
       &amp;                                   meanLabel                              , &amp;
       &amp;                                   meanComment                            , &amp;
       &amp;                                   meanUnits                              , &amp;
       &amp;                                   meanUnitsInSI                          , &amp;
       &amp;                                   binCenter                              , &amp;
       &amp;                                   bufferCount                            , &amp;
       &amp;                                   outputWeight                           , &amp;
       &amp;                                   nodePropertyExtractor_                 , &amp;
       &amp;                                   outputAnalysisPropertyOperator_        , &amp;
       &amp;                                   outputAnalysisPropertyUnoperator_      , &amp;
       &amp;                                   outputAnalysisWeightOperatorUnweighted_, &amp;
       &amp;                                   outputAnalysisDistributionOperator_    , &amp;
       &amp;                                   outputAnalysisDistributionNormalizer_  , &amp;
       &amp;                                   galacticFilter_                        , &amp;
       &amp;                                   outputTimes_                           , &amp;
       &amp;                                   covarianceModel                        , &amp;
       &amp;                                   covarianceBinomialBinsPerDecade        , &amp;
       &amp;                                   covarianceBinomialMassHaloMinimum      , &amp;
       &amp;                                   covarianceBinomialMassHaloMaximum      , &amp;
       &amp;                          binWidth=binWidth                                 &amp;
       &amp;                         )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionWeighted">
     <constructor>
      outputAnalysisVolumeFunction1D (                                                  &amp;
       &amp;                                   label                                  , &amp;
       &amp;                                   comment                                , &amp;
       &amp;                                   propertyLabel                          , &amp;
       &amp;                                   propertyComment                        , &amp;
       &amp;                                   propertyUnits                          , &amp;
       &amp;                                   propertyUnitsInSI                      , &amp;
       &amp;                                   meanLabel                              , &amp;
       &amp;                                   meanComment                            , &amp;
       &amp;                                   meanUnits                              , &amp;
       &amp;                                   meanUnitsInSI                          , &amp;
       &amp;                                   binCenter                              , &amp;
       &amp;                                   bufferCount                            , &amp;
       &amp;                                   outputWeight                           , &amp;
       &amp;                                   nodePropertyExtractor_                 , &amp;
       &amp;                                   outputAnalysisPropertyOperator_        , &amp;
       &amp;                                   outputAnalysisPropertyUnoperator_      , &amp;
       &amp;                                   outputAnalysisWeightOperatorWeighted_  , &amp;
       &amp;                                   outputAnalysisDistributionOperator_    , &amp;
       &amp;                                   outputAnalysisDistributionNormalizer_  , &amp;
       &amp;                                   galacticFilter_                        , &amp;
       &amp;                                   outputTimes_                           , &amp;
       &amp;                                   covarianceModel                        , &amp;
       &amp;                                   covarianceBinomialBinsPerDecade        , &amp;
       &amp;                                   covarianceBinomialMassHaloMinimum      , &amp;
       &amp;                                   covarianceBinomialMassHaloMaximum      , &amp;
       &amp;                          binWidth=binWidth                                 &amp;
       &amp;                         )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="crossCovariance">
     <constructor>
      outputAnalysisCrossCorrelator1D(                                                  &amp;
       &amp;                                   binCenter                              , &amp;
       &amp;                                   bufferCount                            , &amp;
       &amp;                                   outputWeight                           , &amp;
       &amp;                                   nodePropertyExtractor_                 , &amp;
       &amp;                                   outputAnalysisPropertyOperator_        , &amp;
       &amp;                                   outputAnalysisPropertyUnoperator_      , &amp;
       &amp;                                   outputAnalysisWeightOperatorUnweighted_, &amp;
       &amp;                                   outputAnalysisWeightOperatorWeighted_  , &amp;
       &amp;                                   outputAnalysisDistributionOperator_    , &amp;
       &amp;                                   outputAnalysisDistributionNormalizer_  , &amp;
       &amp;                                   galacticFilter_                        , &amp;
       &amp;                                   outputTimes_                           , &amp;
       &amp;                                   covarianceModel                        , &amp;
       &amp;                                   covarianceBinomialBinsPerDecade        , &amp;
       &amp;                                   covarianceBinomialMassHaloMinimum      , &amp;
       &amp;                                   covarianceBinomialMassHaloMaximum      , &amp;
       &amp;                                   binWidth                                 &amp;
       &amp;                         )
     </constructor>
    </referenceConstruct>
    <objectDestructor name="outputAnalysisDistributionNormalizer_"    />
    <objectDestructor name="outputAnalysisWeightOperatorWeighted_"    />
    <objectDestructor name="outputAnalysisWeightOperatorUnweighted_"  />
    <objectDestructor name="weightOperatorUnweightedPropertyOperator_"/>
    <objectDestructor name="propertyOperatorUnweightedBoolean_"       />
    <objectDestructor name="weightOperatorWeightProperty_"            />
    <objectDestructor name="weightOperatorUnweightedProperty_"        />
    !!]
    nullify(weightOperatorWeight_      )
    nullify(weightOperatorUnweighted_  )
    nullify(propertyOperator_          )
    nullify(propertyOperatorUnweighted_)
    return
  end function meanFunction1DConstructorInternal

  subroutine meanFunction1DDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMeanFunction1D} output analysis class.
    !!}
    implicit none
    type(outputAnalysisMeanFunction1D), intent(inout) :: self

    !![
    <objectDestructor name="self%volumeFunctionUnweighted"              />
    <objectDestructor name="self%volumeFunctionWeighted"                />
    <objectDestructor name="self%crossCovariance"                       />
    <objectDestructor name="self%nodePropertyExtractor_"                />
    <objectDestructor name="self%outputAnalysisWeightPropertyExtractor_"/>
    <objectDestructor name="self%outputAnalysisPropertyOperator_"       />
    <objectDestructor name="self%outputAnalysisWeightPropertyOperator_" />
    <objectDestructor name="self%outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="self%outputAnalysisWeightOperator_"         />
    <objectDestructor name="self%outputAnalysisDistributionOperator_"   />
    <objectDestructor name="self%galacticFilter_"                       />
    <objectDestructor name="self%outputTimes_"                          />
    !!]
    return
  end subroutine meanFunction1DDestructor

  subroutine meanFunction1DAnalyze(self,node,iOutput)
    !!{
    Implement a meanFunction1D output analysis.
    !!}
    implicit none
    class  (outputAnalysisMeanFunction1D), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(c_size_t                    ), intent(in   ) :: iOutput

    ! Analyze for all three volume functions.
    call self%volumeFunctionUnweighted%analyze(node,iOutput)
    call self%volumeFunctionWeighted  %analyze(node,iOutput)
    call self%crossCovariance         %analyze(node,iOutput)
    return
  end subroutine meanFunction1DAnalyze

  subroutine meanFunction1DReduce(self,reduced)
    !!{
    Implement a volumeFunction1D output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisMeanFunction1D), intent(inout) :: self
    class(outputAnalysisClass         ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisMeanFunction1D)
       call self%volumeFunctionUnweighted%reduce(reduced%volumeFunctionUnweighted)
       call self%volumeFunctionWeighted  %reduce(reduced%volumeFunctionWeighted  )
       call self%crossCovariance         %reduce(reduced%crossCovariance         )
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine meanFunction1DReduce

  subroutine meanFunction1DFinalizeAnalysis(self)
    !!{
    Finalize analysis of a {\normalfont \ttfamily meanFunction1D} output analysis.
    !!}
    implicit none
    class           (outputAnalysisMeanFunction1D), intent(inout)                 :: self
    double precision                              , allocatable  , dimension(:  ) :: unweightedValue     , weightedValue
    double precision                              , allocatable  , dimension(:,:) :: unweightedCovariance, weightedCovariance, &
         &                                                                           crossCovariance
    integer         (c_size_t                    )                                :: i                   , j

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! Retrieve results from our 1-D volume functions.
    call self%volumeFunctionUnweighted%results(binCenter=self%binCenter,functionValue=unweightedValue,functionCovariance=unweightedCovariance)
    call self%volumeFunctionWeighted  %results(                         functionValue=weightedValue  ,functionCovariance=weightedCovariance  )
    call self%crossCovariance         %results(                                                       functionCovariance=crossCovariance     )
    ! Estimate covariance of ratio using Taylor series expansion approach
    ! (e.g. http://math.stackexchange.com/questions/40713/calculating-the-variance-of-the-ratio-of-random-variables).
    allocate(self%meanValue     (size(self%binCenter)                     ))
    allocate(self%meanCovariance(size(self%binCenter),size(self%binCenter)))
    self%meanValue     =0.0d0
    self%meanCovariance=0.0d0
    do i=1,size(self%binCenter)
       if (unweightedValue(i) > 0.0d0) then
          self%meanValue(i)=weightedValue(i)/unweightedValue(i)
          do j=i,size(self%binCenter)
             if (unweightedValue(j) > 0.0d0) then
                self%meanCovariance(i,j)=+(                                                                                                      &
                     &                     +  weightedCovariance(i,j)* unweightedValue(i)                 *unweightedValue(j)                    &
                     &                     +unweightedCovariance(i,j)*                    weightedValue(i)                   *weightedValue(j)   &
                     &                     -     crossCovariance(i,j)*(unweightedValue(i)*weightedValue(j)+unweightedValue(j)*weightedValue(i))  &
                     &                    )                                                                                                      &
                     &                   /unweightedValue(i)**2                                                                                  &
                     &                   /unweightedValue(j)**2
                if (i == j) self%meanCovariance(i,j)=max(self%meanCovariance(i,j),0.0d0)
                self%meanCovariance(j,i)=self%meanCovariance(i,j)
            end if
          end do
       end if
    end do
    return
  end subroutine meanFunction1DFinalizeAnalysis

  subroutine meanFunction1DFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily meanFunction1D} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisMeanFunction1D), intent(inout)           :: self
    type (varying_string              ), intent(in   ), optional :: groupName
    type (hdf5Object                  )               , target   :: analysesGroup, subGroup
    type (hdf5Object                  )               , pointer  :: inGroup
    type (hdf5Object                  )                          :: analysisGroup, dataset

    ! Finalize the analysis.
    call self%finalizeAnalysis()
    ! Output the resulting mean function.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup(char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call    analysisGroup%writeAttribute(     char(self%   comment   )                    ,'description'                                                                                                   )
    call    analysisGroup%writeAttribute("function1D"                                     ,'type'                                                                                                          )
    call    analysisGroup%writeAttribute(     char(self%   xAxisLabel)                    ,'xAxisLabel'                                                                                                    )
    call    analysisGroup%writeAttribute(     char(self%   yAxisLabel)                    ,'yAxisLabel'                                                                                                    )
    call    analysisGroup%writeAttribute(          self%   xAxisIsLog                     ,'xAxisIsLog'                                                                                                    )
    call    analysisGroup%writeAttribute(          self%   yAxisIsLog                     ,'yAxisIsLog'                                                                                                    )
    call    analysisGroup%writeAttribute(     char(self%propertyLabel)                    ,'xDataset'                                                                                                      )
    call    analysisGroup%writeAttribute(     char(self%    meanLabel)                    ,'yDataset'                                                                                                      )
    call    analysisGroup%writeAttribute(     char(self%    meanLabel)//"Covariance"      ,'yCovariance'                                                                                                   )
    if (allocated(self%meanValueTarget)) then
       call analysisGroup%writeAttribute(     char(self%    meanLabel)//"Target"          ,'yDatasetTarget'                                                                                                )
       call analysisGroup%writeAttribute(     char(self%    meanLabel)//"CovarianceTarget",'yCovarianceTarget'                                                                                             )
    end if
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (          self%binCenter                         ,char(self%propertyLabel)                    ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%propertyUnits    )                ,'units'                                                                                                         )
    call    dataset      %writeAttribute(          self%propertyUnitsInSI                 ,'unitsInSI'                                                                                                     )
    call    dataset      %close         (                                                                                                                                                                  )
    call    analysisGroup%writeDataset  (          self%meanValue                         ,char(self%    meanLabel)                    ,char(self%    meanComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
    call    dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
    call    dataset      %close         (                                                                                                                                                                  )
    call    analysisGroup%writeDataset  (          self%meanCovariance                    ,char(self%    meanLabel)//"Covariance"      ,char(self%    meanComment)//" [covariance]",datasetReturned=dataset)
    call    dataset      %writeAttribute("["//char(self%    meanUnits    )//"]²"          ,'units'                                                                                                         )
    call    dataset      %writeAttribute(          self%    meanUnitsInSI   **2           ,'unitsInSI'                                                                                                     )
    call    dataset      %close         (                                                                                                                                                                  )
    ! If available, include the log-likelihood and target dataset.
    if (allocated(self%meanValueTarget)) then
       call analysisGroup%writeAttribute(          self%logLikelihood()                   ,'logLikelihood'                                                                                                 )
       call analysisGroup%writeAttribute(     char(self%targetLabel         )             ,'targetLabel'                                                                                                   )
       call analysisGroup%writeDataset  (          self%meanValueTarget                   ,char(self%    meanLabel)//"Target"          ,char(self%    meanComment)                 ,datasetReturned=dataset)
       call dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
       call dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
       call dataset      %close         (                                                                                                                                                                  )
       call analysisGroup%writeDataset  (          self%meanCovarianceTarget              ,char(self%    meanLabel)//"CovarianceTarget",char(self%    meanComment)//" [covariance]",datasetReturned=dataset)
       call dataset      %writeAttribute("["//char(self%    meanUnits    )//"]²"          ,'units'                                                                                                         )
       call dataset      %writeAttribute(          self%    meanUnitsInSI   **2           ,'unitsInSI'                                                                                                     )
       call dataset      %close         (                                                                                                                                                                  )
    end if
    call    analysisGroup%close         (                                                                                                                                                                  )
    if (present(groupName)) &
         & call subGroup %close         (                                                                                                                                                                  )
    call    analysesGroup%close         (                                                                                                                                                                  )
    !$ call hdf5Access%unset()
    return
  end subroutine meanFunction1DFinalize

  subroutine meanFunction1DResults(self,binCenter,meanValue,meanCovariance)
    !!{
    Implement a meanFunction1D output analysis finalization.
    !!}
    implicit none
    class           (outputAnalysisMeanFunction1D)                             , intent(inout)           :: self
    double precision                              , allocatable, dimension(:  ), intent(inout), optional :: binCenter     , meanValue
    double precision                              , allocatable, dimension(:,:), intent(inout), optional :: meanCovariance

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Return results.
    if (present(binCenter     )) then
       if (allocated(binCenter     )) deallocate(binCenter     )
       allocate(binCenter(size(self%binCenter)))
       binCenter         =self%binCenter
    end if
    if (present(meanValue     )) then
       if (allocated(meanValue     )) deallocate(meanValue     )
       allocate(meanValue(size(self%meanValue)))
       meanValue     =self%meanValue
    end if
    if (present(meanCovariance)) then
       if (allocated(meanCovariance)) deallocate(meanCovariance)
       allocate(meanCovariance(size(self%meanCovariance,dim=1),size(self%meanCovariance,dim=2)))
       meanCovariance=self%meanCovariance
    end if
    return
  end subroutine meanFunction1DResults

  double precision function meanFunction1DLogLikelihood(self)
    !!{
    Return the log-likelihood of a meanFunction1D output analysis.
    !!}
    use :: Error                       , only : Error_Report
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisMeanFunction1D), intent(inout)                 :: self
    double precision                              , allocatable  , dimension(:,:) :: meanCovarianceCombined
    double precision                              , allocatable  , dimension(:  ) :: meanValueDifference
    type            (vector                      )                                :: residual
    type            (matrix                      )                                :: covariance
    integer                                                                       :: status

    ! Check for existence of a target distribution.
    if (allocated(self%meanValueTarget)) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! Allocate workspaces.
       allocate(meanCovarianceCombined(size(self%binCenter),size(self%binCenter)))
       allocate(meanValueDifference   (size(self%binCenter)                     ))
       ! Find combined covariance and difference between model and target.
       meanValueDifference   =+self%meanValue            &
            &                 -self%meanValueTarget
       meanCovarianceCombined=+self%meanCovariance       &
            &                 +self%meanCovarianceTarget
       residual              = vector(meanValueDifference   )
       covariance            = matrix(meanCovarianceCombined)
       ! Compute the log-likelihood.
       meanFunction1DLogLikelihood           =-0.5d0*covariance%covarianceProduct(residual,status)
       if (status == GSL_Success) then
          if (self%likelihoodNormalize)                                                      &
               & meanFunction1DLogLikelihood=+meanFunction1DLogLikelihood                    &
               &                             -0.5d0*covariance%logarithmicDeterminant()      &
               &                             -0.5d0*dble(size(self%binCenter))*log(2.0d0*Pi)
       else
          meanFunction1DLogLikelihood       =+logImprobable
       end if
    else
       meanFunction1DLogLikelihood=0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function meanFunction1DLogLikelihood
