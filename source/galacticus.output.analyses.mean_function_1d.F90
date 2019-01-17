!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Contains a module which implements a generic 1D mean function (i.e. mean value of some property weighted by number density of
  !% objects binned by some property) output analysis class.
  
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Output_Analysis_Property_Extractions
  use               Output_Analysis_Property_Operators
  use               Output_Analysis_Weight_Operators
  use               Output_Analysis_Distribution_Operators
  use               Output_Analysis_Distribution_Normalizers
  use               Output_Analyses_Options
  use               Output_Times
  use               Galactic_Filters

  !# <outputAnalysis name="outputAnalysisMeanFunction1D" defaultThreadPrivate="yes">
  !#  <description>A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisClass) :: outputAnalysisMeanFunction1D
     !% A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.
     type            (varying_string                )                              :: label                   , comment               , &
          &                                                                           propertyLabel           , propertyComment       , &
          &                                                                           meanLabel               , meanComment           , &
          &                                                                           propertyUnits           , meanUnits
     double precision                                                              :: propertyUnitsInSI       , meanUnitsInSI
     type            (outputAnalysisVolumeFunction1D)                              :: volumeFunctionUnweighted, volumeFunctionWeighted, &
          &                                                                           volumeFunctionCross     
     double precision                                , allocatable, dimension(:  ) :: binCenter               , meanValue             , &
          &                                                                           meanValueTarget
     double precision                                , allocatable, dimension(:,:) :: meanCovariance          , meanCovarianceTarget
     logical                                                                       :: finalized               , likelihoodNormalize
   contains
     !@ <objectMethods>
     !@   <object>outputAnalysisMeanFunction1D</object>
     !@   <objectMethod>
     !@     <method>results</method>
     !@     <arguments>\doubleone\ [binCenter]\arginout, \doubletwo\ [functionValue]\arginout, \doubletwo\ [functionCovariance]\arginout</arguments>
     !@     <type>\void</type>
     !@     <description>Return the results of the mean function operator.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: deepCopy      => meanFunction1DDeepCopy
     procedure :: analyze       => meanFunction1DAnalyze
     procedure :: finalize      => meanFunction1DFinalize
     procedure :: reduce        => meanFunction1DReduce
     procedure :: results       => meanFunction1DResults
     procedure :: logLikelihood => meanFunction1DLogLikelihood
  end type outputAnalysisMeanFunction1D

  interface outputAnalysisMeanFunction1D
     !% Constructors for the ``meanFunction1D'' output analysis class.
     module procedure meanFunction1DConstructorParameters
     module procedure meanFunction1DConstructorInternal
  end interface outputAnalysisMeanFunction1D

contains

  function meanFunction1DConstructorParameters(parameters) result(self)
    !% Constructor for the ``meanFunction1D'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (outputAnalysisMeanFunction1D           )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (outputAnalysisPropertyExtractorClass   ), pointer                     :: outputAnalysisPropertyExtractor_     , outputAnalysisWeightPropertyExtractor_
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
         &                                                                                    covarianceModel
    integer                                                                                :: covarianceBinomialBinsPerDecade
    type            (inputParameters                        )                              :: unoperatorParameters
    type            (inputParameters                        )                              :: weightParameters
    double precision                                                                       :: propertyUnitsInSI                    , meanUnitsInSI                         , &
         &                                                                                    covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum
    logical                                                                                :: likelihoodNormalize

    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisPropertyExtractor_"       source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisWeightPropertyExtractor_" source="weightParameters"    />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"        source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisWeightPropertyOperator_"  source="weightParameters"    />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"      source="unoperatorParameters"/>
    !# <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"          source="parameters"          />
    !# <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"    source="parameters"          />
    !# <objectBuilder class="galacticFilter"                       name="galacticFilter_"                        source="parameters"          />
    !# <objectBuilder class="outputTimes"                          name="outputTimes_"                           source="parameters"          />
    unoperatorParameters=parameters%subParameters('unoperator',requireValue=.false.)
    weightParameters    =parameters%subParameters('weight'    ,requireValue=.false.)
    call allocateArray(binCenter   ,[int(parameters%count('binCenter'),kind=c_size_t)                     ])
    call allocateArray(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t)*outputTimes_%count()])
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*outputTimes_%count()) &
         & call Galacticus_Error_Report('incorrect number of output weights provided'//{introspection:location})
    !# <inputParameter>
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <variable>label</variable>
    !#   <description>A label for the analysis.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>comment</name>
    !#   <source>parameters</source>
    !#   <variable>comment</variable>
    !#   <description>A descriptive comment for the analysis.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyLabel</name>
    !#   <source>parameters</source>
    !#   <variable>propertyLabel</variable>
    !#   <description>A label for the property variable.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyComment</name>
    !#   <source>parameters</source>
    !#   <variable>propertyComment</variable>
    !#   <description>A descriptive comment for the property variable.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyUnits</name>
    !#   <source>parameters</source>
    !#   <variable>propertyUnits</variable>
    !#   <description>A human-readable description of the units for the property.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyUnitsInSI</name>
    !#   <source>parameters</source>
    !#   <variable>propertyUnitsInSI</variable>
    !#   <description>A units for the property in the SI system.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLabel</name>
    !#   <source>parameters</source>
    !#   <variable>meanLabel</variable>
    !#   <description>A label for the mean.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanComment</name>
    !#   <source>parameters</source>
    !#   <variable>meanComment</variable>
    !#   <description>A descriptive comment for the mean.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanUnits</name>
    !#   <source>parameters</source>
    !#   <variable>meanUnits</variable>
    !#   <description>A human-readable description of the units for the mean.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanUnitsInSI</name>
    !#   <source>parameters</source>
    !#   <variable>meanUnitsInSI</variable>
    !#   <description>A units for the mean in the SI system.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>    
    !# <inputParameter>
    !#   <name>binCenter</name>
    !#   <source>parameters</source>
    !#   <variable>binCenter</variable>
    !#   <description>The value of the property at the center of each bin.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bufferCount</name>
    !#   <source>parameters</source>
    !#   <variable>bufferCount</variable>
    !#   <description>The number of buffer bins to include below and above the range of actual bins.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputWeight</name>
    !#   <source>parameters</source>
    !#   <variable>outputWeight</variable>
    !#   <description>The weight to assign to each bin at each output.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceModel</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceModel</variable>
    !#   <description>The model to use for computing covariances.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>likelihoodNormalize</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true then normalize the likelihood to make it a probability density.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
   if (parameters%isPresent('meanValueTarget')) then
       if (parameters%isPresent('meanCovarianceTarget')) then
          !# <inputParameter>
          !#   <name>meanValueTarget</name>
          !#   <source>parameters</source>
          !#   <description>The target function for likelihood calculations.</description>
          !#   <type>real</type>
          !#   <cardinality>0..1</cardinality>
          !# </inputParameter> 
          !# <inputParameter>
          !#   <name>meanCovarianceTarget</name>
          !#   <source>parameters</source>
          !#   <variable>meanCovarianceTarget1D</variable>
          !#   <description>The target function covariance for likelihood calculations.</description>
          !#   <type>real</type>
          !#   <cardinality>0..1</cardinality>
          !# </inputParameter>
          if (size(meanCovarianceTarget1D) == size(meanValueTarget)**2) then
             allocate(meanCovarianceTarget(size(meanValueTarget),size(meanValueTarget)))
             meanCovarianceTarget=reshape(meanCovarianceTarget1D,shape(meanCovarianceTarget))
          else
             call Galacticus_Error_Report('meanCovariance has wrong size'//{introspection:location})
          end if
       else
          call Galacticus_Error_Report('meanCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('meanCovariance')) call Galacticus_Error_Report('functionTarget must be specified if meanCovariance is present'//{introspection:location})
    end if
    ! Build the object.
    !# <conditionalCall>
    !#  <call>
    !#   self=outputAnalysisMeanFunction1D(                                                                                               &amp;
    !#        &amp;                        label                                                                                        , &amp;
    !#        &amp;                        comment                                                                                      , &amp;
    !#        &amp;                        propertyLabel                                                                                , &amp;
    !#        &amp;                        propertyComment                                                                              , &amp;
    !#        &amp;                        propertyUnits                                                                                , &amp;
    !#        &amp;                        propertyUnitsInSI                                                                            , &amp;
    !#        &amp;                        meanLabel                                                                                    , &amp;
    !#        &amp;                        meanComment                                                                                  , &amp;
    !#        &amp;                        meanUnits                                                                                    , &amp;
    !#        &amp;                        meanUnitsInSI                                                                                , &amp;
    !#        &amp;                        binCenter                                                                                    , &amp;
    !#        &amp;                        bufferCount                                                                                  , &amp;
    !#        &amp;                        reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),outputTimes_%count()]), &amp;
    !#        &amp;                        outputAnalysisPropertyExtractor_                                                             , &amp;
    !#        &amp;                        outputAnalysisWeightPropertyExtractor_                                                       , &amp;
    !#        &amp;                        outputAnalysisPropertyOperator_                                                              , &amp;
    !#        &amp;                        outputAnalysisWeightPropertyOperator_                                                        , &amp;
    !#        &amp;                        outputAnalysisPropertyUnoperator_                                                            , &amp;
    !#        &amp;                        outputAnalysisWeightOperator_                                                                , &amp;
    !#        &amp;                        outputAnalysisDistributionOperator_                                                          , &amp;
    !#        &amp;                        galacticFilter_                                                                              , &amp;
    !#        &amp;                        outputTimes_                                                                                 , &amp;
    !#        &amp;                        enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.) , &amp;
    !#        &amp;                        covarianceBinomialBinsPerDecade                                                              , &amp;
    !#        &amp;                        covarianceBinomialMassHaloMinimum                                                            , &amp;
    !#        &amp;                        covarianceBinomialMassHaloMaximum                                                            , &amp;                      
    !#        &amp;                        likelihoodNormalize                                                                            &amp;
    !#        &amp;                        {conditions}                                                                                   &amp;
    !#        &amp;                       )
    !#  </call>
    !#  <argument name="meanValueTarget"      value="meanValueTarget"      parameterPresent="parameters"/>
    !#  <argument name="meanCovarianceTarget" value="meanCovarianceTarget" parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    return
  end function meanFunction1DConstructorParameters

  function meanFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,meanLabel,meanComment,meanUnits,meanUnitsInSI,binCenter,bufferCount,outputWeight,outputAnalysisPropertyExtractor_,outputAnalysisWeightPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisWeightPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperatorIn_,outputAnalysisDistributionOperator_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,meanValueTarget,meanCovarianceTarget) result (self)
    !% Constructor for the ``meanFunction1D'' output analysis class for internal use.
    use Memory_Management
    implicit none
    type            (outputAnalysisMeanFunction1D                )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                           , comment                                        , &
         &                                                                                                     propertyLabel                                   , propertyComment                                , &
         &                                                                                                     meanLabel                                       , meanComment                                    , &
         &                                                                                                     propertyUnits                                   , meanUnits
    double precision                                                                                        :: propertyUnitsInSI                               , meanUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    logical                                                       , intent(in   ), optional                 :: likelihoodNormalize
    class           (outputAnalysisPropertyExtractorClass        ), intent(inout), target                   :: outputAnalysisPropertyExtractor_                , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(inout), target                   :: outputAnalysisPropertyOperator_                 , outputAnalysisPropertyUnoperator_              , &
         &                                                                                                     outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(inout), target                   :: outputAnalysisWeightOperatorIn_
    class           (outputAnalysisDistributionOperatorClass     ), intent(inout), target                   :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                         ), intent(inout), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(inout), target                   :: outputTimes_
    integer                                                       , intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum               , covarianceBinomialMassHaloMaximum
    double precision                                              , intent(in   ), optional, dimension(:  ) :: meanValueTarget
    double precision                                              , intent(in   ), optional, dimension(:,:) :: meanCovarianceTarget
    type            (outputAnalysisDistributionNormalizerIdentity), pointer                                 :: outputAnalysisDistributionNormalizerWeighted_   , outputAnalysisDistributionNormalizerUnweighted_, &
         &                                                                                                     outputAnalysisDistributionNormalizerCross_
    type            (outputAnalysisWeightOperatorSequence        ), pointer                                 :: outputAnalysisWeightOperatorWeighted_           , outputAnalysisWeightOperatorCross_             , &
         &                                                                                                     outputAnalysisWeightOperatorUnweighted_
    type            (weightOperatorList                          ), pointer                                 :: weightOperatorWeight_                           , weightOperatorCross_                           , &
         &                                                                                                     weightOperatorUnweighted_
    type            (outputAnalysisPropertyOperatorSequence      ), pointer                                 :: weightOperatorCrossPropertyOperator_            , weightOperatorUnweightedPropertyOperator_
    type            (propertyOperatorList                        ), pointer                                 :: propertyOperator_                               , propertyOperatorUnweighted_
    class           (outputAnalysisPropertyExtractorClass        ), pointer                                 :: outputAnalysisPropertyExtractorUnweighted_      , outputAnalysisPropertyExtractorCross_          , &
         &                                                                                                     outputAnalysisWeightPropertyExtractorUnweighted_, outputAnalysisWeightPropertyExtractorCross_
    class           (outputAnalysisPropertyOperatorClass         ), pointer                                 :: outputAnalysisPropertyOperatorUnweighted_       , outputAnalysisPropertyOperatorCross_           , &
         &                                                                                                     outputAnalysisPropertyUnoperatorUnweighted_     , outputAnalysisPropertyUnoperatorCross_         , &
         &                                                                                                     outputAnalysisWeightPropertyOperatorUnweighted_ , outputAnalysisWeightPropertyOperatorCross_
    class           (outputAnalysisWeightOperatorClass           ), pointer                                 :: outputAnalysisWeightOperatorInUnweighted_       , outputAnalysisWeightOperatorInCross_
    class           (outputAnalysisDistributionOperatorClass     ), pointer                                 :: outputAnalysisDistributionOperatorUnweighted_   , outputAnalysisDistributionOperatorCross_
    class           (galacticFilterClass                         ), pointer                                 :: galacticFilterUnweighted_                       , galacticFilterCross_
    class           (outputTimesClass                            ), pointer                                 :: outputTimesUnweighted_                          , outputTimesCross_
    !# <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, meanLabel, meanComment, meanUnits, meanUnitsInSI, meanValueTarget, meanCovarianceTarget"/>

    ! Mark as unfinalized.
    self%finalized=.false.
    ! Set normalization state for likelihood.
    self%likelihoodNormalize=.true.
    if (present(likelihoodNormalize)) self%likelihoodNormalize=likelihoodNormalize
    ! Build copies of all objects for the weighted and cross volume function 1D objects.
    allocate(outputAnalysisWeightOperatorInUnweighted_       ,mold=outputAnalysisWeightOperatorIn_       )
    allocate(outputAnalysisWeightOperatorInCross_            ,mold=outputAnalysisWeightOperatorIn_       )
    allocate(outputAnalysisPropertyExtractorUnweighted_      ,mold=outputAnalysisPropertyExtractor_      )
    allocate(outputAnalysisPropertyExtractorCross_           ,mold=outputAnalysisPropertyExtractor_      )
    allocate(outputAnalysisPropertyOperatorUnweighted_       ,mold=outputAnalysisPropertyOperator_       )
    allocate(outputAnalysisPropertyOperatorCross_            ,mold=outputAnalysisPropertyOperator_       )
    allocate(outputAnalysisPropertyUnoperatorUnweighted_     ,mold=outputAnalysisPropertyUnoperator_     )
    allocate(outputAnalysisPropertyUnoperatorCross_          ,mold=outputAnalysisPropertyUnoperator_     )
    allocate(outputAnalysisWeightPropertyExtractorUnweighted_,mold=outputAnalysisWeightPropertyExtractor_)
    allocate(outputAnalysisWeightPropertyExtractorCross_     ,mold=outputAnalysisWeightPropertyExtractor_)
    allocate(outputAnalysisWeightPropertyOperatorUnweighted_ ,mold=outputAnalysisWeightPropertyOperator_ )
    allocate(outputAnalysisWeightPropertyOperatorCross_      ,mold=outputAnalysisWeightPropertyOperator_ )
    allocate(outputAnalysisDistributionOperatorUnweighted_   ,mold=outputAnalysisDistributionOperator_   )
    allocate(outputAnalysisDistributionOperatorCross_        ,mold=outputAnalysisDistributionOperator_   )
    allocate(galacticFilterUnweighted_                       ,mold=galacticFilter_                       )
    allocate(galacticFilterCross_                            ,mold=galacticFilter_                       )
    allocate(outputTimesUnweighted_                          ,mold=outputTimes_                          )
    allocate(outputTimesCross_                               ,mold=outputTimes_                          )
    call outputAnalysisWeightOperatorIn_       %deepCopy(outputAnalysisWeightOperatorInUnweighted_       )
    call outputAnalysisWeightOperatorIn_       %deepCopy(outputAnalysisWeightOperatorInCross_            )
    call outputAnalysisPropertyExtractor_      %deepCopy(outputAnalysisPropertyExtractorUnweighted_      )
    call outputAnalysisPropertyExtractor_      %deepCopy(outputAnalysisPropertyExtractorCross_           )
    call outputAnalysisPropertyOperator_       %deepCopy(outputAnalysisPropertyOperatorUnweighted_       )
    call outputAnalysisPropertyOperator_       %deepCopy(outputAnalysisPropertyOperatorCross_            )
    call outputAnalysisPropertyUnoperator_     %deepCopy(outputAnalysisPropertyUnoperatorUnweighted_     )
    call outputAnalysisPropertyUnoperator_     %deepCopy(outputAnalysisPropertyUnoperatorCross_          )
    call outputAnalysisWeightPropertyExtractor_%deepCopy(outputAnalysisWeightPropertyExtractorUnweighted_)
    call outputAnalysisWeightPropertyExtractor_%deepCopy(outputAnalysisWeightPropertyExtractorCross_     )
    call outputAnalysisWeightPropertyOperator_ %deepCopy(outputAnalysisWeightPropertyOperatorUnweighted_ )
    call outputAnalysisWeightPropertyOperator_ %deepCopy(outputAnalysisWeightPropertyOperatorCross_      )
    call outputAnalysisDistributionOperator_   %deepCopy(outputAnalysisDistributionOperatorUnweighted_   )
    call outputAnalysisDistributionOperator_   %deepCopy(outputAnalysisDistributionOperatorCross_        )
    call galacticFilter_                       %deepCopy(galacticFilterUnweighted_                       )
    call galacticFilter_                       %deepCopy(galacticFilterCross_                            )
    call outputTimes_                          %deepCopy(outputTimesUnweighted_                          )
    call outputTimes_                          %deepCopy(outputTimesCross_                               )
    ! Build an identity distribution normalizer.
    allocate(outputAnalysisDistributionNormalizerWeighted_  )
    allocate(outputAnalysisDistributionNormalizerUnweighted_)
    allocate(outputAnalysisDistributionNormalizerCross_     )
    outputAnalysisDistributionNormalizerWeighted_  =outputAnalysisDistributionNormalizerIdentity()
    outputAnalysisDistributionNormalizerUnweighted_=outputAnalysisDistributionNormalizerIdentity()
    outputAnalysisDistributionNormalizerCross_     =outputAnalysisDistributionNormalizerIdentity()
    ! Build weight operator that includes the weight property.
    allocate(                                        outputAnalysisWeightOperatorWeighted_               )
    allocate(                                        weightOperatorWeight_                               )
    allocate(                                        weightOperatorWeight_                %next          )
    allocate(outputAnalysisWeightOperatorProperty :: weightOperatorWeight_                %next%operator_)
    weightOperatorWeight_%operator_       => outputAnalysisWeightOperatorIn_
    select type (operator_ => weightOperatorWeight_%next%operator_)
    type is (outputAnalysisWeightOperatorProperty)
       operator_=outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractor_,outputAnalysisWeightPropertyOperator_)
    end select
    outputAnalysisWeightOperatorWeighted_=outputAnalysisWeightOperatorSequence(weightOperatorWeight_)
    ! Build weight operator that includes sqaure-root of the weight property. We use this to effectively compute the covariance
    ! between numerator and denominator in our ratio of volume functions. If we define:
    !
    ! 〈W〉  = ∑ wᵢ  ; Var(W ) = ∑  wᵢ²
    ! 〈WX〉 = ∑ wᵢ xᵢ; Var(WX) = ∑ (wᵢ xᵢ)²
    !
    ! where wᵢ is the weight, and xᵢ is the property for which we wish to compute the mean, then the above quantities are those
    ! compute (mean and variance) for the unweighted and weighted volume functions (assuming Poisson statistics, but the approach
    ! applies generally). The covariance of the weighted vs. unweighted functions is then:
    !
    ! Cov(WX,W) = ∑ wᵢ wᵢ xᵢ
    !
    ! which is the covariance computed by the volume function class for a property √x, providing x∈[0,∞).
    allocate(                                            outputAnalysisWeightOperatorCross_                      )
    allocate(                                            weightOperatorCross_                                    )
    allocate(                                            weightOperatorCross_                     %next          )
    allocate(outputAnalysisWeightOperatorProperty     :: weightOperatorCross_                     %next%operator_)
    allocate(                                            weightOperatorCrossPropertyOperator_                    )
    allocate(                                            propertyOperator_                                       )
    allocate(                                            propertyOperator_                        %next          )
    allocate(outputAnalysisPropertyOperatorSquareRoot :: propertyOperator_                        %next%operator_)
    ! Build a sequence property operator for the cross-correlation weight property. This includes any property operator passed to
    ! this constructor, plus the square-root operator.
    propertyOperator_%operator_       => outputAnalysisWeightPropertyOperatorCross_ ! First operator in the sequence is whatever operator was passed to this constructor.
    select type (operator_ => propertyOperator_%next%operator_)                ! Second operator in the sequence is the square-root operator.
    type is (outputAnalysisPropertyOperatorSquareRoot)
       operator_=outputAnalysisPropertyOperatorSquareRoot()
    end select
    weightOperatorCrossPropertyOperator_=outputAnalysisPropertyOperatorSequence(propertyOperator_)
    ! Build a weight operator sequence - the first operator is whatever weight operator was passed to this constructor, the second
    ! is the operator which weights by the property value.
    weightOperatorCross_%operator_ => outputAnalysisWeightOperatorInCross_
    select type (operator_ => weightOperatorCross_%next%operator_)
    type is (outputAnalysisWeightOperatorProperty)
       operator_=outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractorCross_,weightOperatorCrossPropertyOperator_)
    end select
    outputAnalysisWeightOperatorCross_=outputAnalysisWeightOperatorSequence(weightOperatorCross_)
    ! Build weight operator that includes a final boolean operator. We use this in the unweighted case to allow filters to be
    ! applied to the weight property. If any filter sets the property value to zero, this boolean operator will return zero,
    ! otherwise it will return unity.
    allocate(                                         outputAnalysisWeightOperatorUnweighted_                 )
    allocate(                                         weightOperatorUnweighted_                               )
    allocate(                                         weightOperatorUnweighted_                %next          )
    allocate(outputAnalysisWeightOperatorProperty  :: weightOperatorUnweighted_                %next%operator_)
    allocate(                                         weightOperatorUnweightedPropertyOperator_               )
    allocate(                                         propertyOperatorUnweighted_                             )
    allocate(                                         propertyOperatorUnweighted_              %next          )
    allocate(outputAnalysisPropertyOperatorBoolean :: propertyOperatorUnweighted_              %next%operator_)
    ! Build a sequence property operator for the unweight property. This includes any property operator passed to this
    ! constructor, plus the boolean operator.
    propertyOperatorUnweighted_%operator_       => outputAnalysisWeightPropertyOperatorUnweighted_ ! First operator in the sequence is whatever operator was passed to this constructor.
    select type (operator_ => propertyOperatorUnweighted_%next%operator_)                ! Second operator in the sequence is the boolean operator.
    type is (outputAnalysisPropertyOperatorBoolean)
       operator_=outputAnalysisPropertyOperatorBoolean()
    end select
    weightOperatorUnweightedPropertyOperator_=outputAnalysisPropertyOperatorSequence(propertyOperatorUnweighted_)
    ! Build a weight operator sequence - the first operator is whatever weight operator was passed to this constructor, the second
    ! is the operator which weights by the property value.
    weightOperatorUnweighted_%operator_ => outputAnalysisWeightOperatorInUnweighted_
    select type (operator_ => weightOperatorUnweighted_%next%operator_)
    type is (outputAnalysisWeightOperatorProperty)
       operator_=outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractorUnweighted_,weightOperatorUnweightedPropertyOperator_)
    end select
    outputAnalysisWeightOperatorUnweighted_=outputAnalysisWeightOperatorSequence(weightOperatorUnweighted_)
    ! Build weighted, unweighted, and cross volume function objects.
    self%volumeFunctionUnweighted=outputAnalysisVolumeFunction1D(                                                 &
         &                                                       label                                          , &
         &                                                       comment                                        , &
         &                                                       propertyLabel                                  , &
         &                                                       propertyComment                                , &
         &                                                       propertyUnits                                  , &
         &                                                       propertyUnitsInSI                              , &
         &                                                       meanLabel                                      , &
         &                                                       meanComment                                    , &
         &                                                       meanUnits                                      , &
         &                                                       meanUnitsInSI                                  , &
         &                                                       binCenter                                      , &
         &                                                       bufferCount                                    , &
         &                                                       outputWeight                                   , &
         &                                                       outputAnalysisPropertyExtractorUnweighted_     , &
         &                                                       outputAnalysisPropertyOperatorUnweighted_      , &
         &                                                       outputAnalysisPropertyUnoperatorUnweighted_    , &
         &                                                       outputAnalysisWeightOperatorUnweighted_        , &
         &                                                       outputAnalysisDistributionOperatorUnweighted_  , &
         &                                                       outputAnalysisDistributionNormalizerUnweighted_, &
         &                                                       galacticFilterUnweighted_                      , &
         &                                                       outputTimesUnweighted_                         , &
         &                                                       covarianceModel                                , &
         &                                                       covarianceBinomialBinsPerDecade                , &
         &                                                       covarianceBinomialMassHaloMinimum              , &
         &                                                       covarianceBinomialMassHaloMaximum                &
         &                                                      )
    self%volumeFunctionWeighted  =outputAnalysisVolumeFunction1D(                                                 &
         &                                                       label                                          , &
         &                                                       comment                                        , &
         &                                                       propertyLabel                                  , &
         &                                                       propertyComment                                , &
         &                                                       propertyUnits                                  , &
         &                                                       propertyUnitsInSI                              , &
         &                                                       meanLabel                                      , &
         &                                                       meanComment                                    , &
         &                                                       meanUnits                                      , &
         &                                                       meanUnitsInSI                                  , &
         &                                                       binCenter                                      , &
         &                                                       bufferCount                                    , &
         &                                                       outputWeight                                   , &
         &                                                       outputAnalysisPropertyExtractor_               , &
         &                                                       outputAnalysisPropertyOperator_                , &
         &                                                       outputAnalysisPropertyUnoperator_              , &
         &                                                       outputAnalysisWeightOperatorWeighted_          , &
         &                                                       outputAnalysisDistributionOperator_            , &
         &                                                       outputAnalysisDistributionNormalizerWeighted_  , &
         &                                                       galacticFilter_                                , &
         &                                                       outputTimes_                                   , &
         &                                                       covarianceModel                                , &
         &                                                       covarianceBinomialBinsPerDecade                , &
         &                                                       covarianceBinomialMassHaloMinimum              , &
         &                                                       covarianceBinomialMassHaloMaximum                &
         &                                                      )
    self%volumeFunctionCross     =outputAnalysisVolumeFunction1D(                                                 &
         &                                                       label                                          , &
         &                                                       comment                                        , &
         &                                                       propertyLabel                                  , &
         &                                                       propertyComment                                , &
         &                                                       propertyUnits                                  , &
         &                                                       propertyUnitsInSI                              , &
         &                                                       meanLabel                                      , &
         &                                                       meanComment                                    , &
         &                                                       meanUnits                                      , &
         &                                                       meanUnitsInSI                                  , &
         &                                                       binCenter                                      , &
         &                                                       bufferCount                                    , &
         &                                                       outputWeight                                   , &
         &                                                       outputAnalysisPropertyExtractorCross_          , &
         &                                                       outputAnalysisPropertyOperatorCross_           , &
         &                                                       outputAnalysisPropertyUnoperatorCross_         , &
         &                                                       outputAnalysisWeightOperatorCross_             , &
         &                                                       outputAnalysisDistributionOperatorCross_       , &
         &                                                       outputAnalysisDistributionNormalizerCross_     , &
         &                                                       galacticFilterCross_                           , &
         &                                                       outputTimesCross_                              , &
         &                                                       covarianceModel                                , &
         &                                                       covarianceBinomialBinsPerDecade                , &
         &                                                       covarianceBinomialMassHaloMinimum              , &
         &                                                       covarianceBinomialMassHaloMaximum                &
         &                                                      )
    ! Nullify objects to avoid destruction.
    nullify(outputAnalysisWeightOperatorWeighted_           )
    nullify(outputAnalysisWeightOperatorCross_              )
    nullify(outputAnalysisWeightOperatorUnweighted_         )
    nullify(weightOperatorWeight_                           )
    nullify(weightOperatorCross_                            )
    nullify(weightOperatorUnweighted_                       )
    nullify(outputAnalysisWeightOperatorUnweighted_         )
    nullify(outputAnalysisWeightOperatorCross_              )
    nullify(outputAnalysisPropertyExtractorUnweighted_      )
    nullify(outputAnalysisPropertyExtractorCross_           )
    nullify(outputAnalysisPropertyOperatorUnweighted_       )
    nullify(outputAnalysisPropertyOperatorCross_            )
    nullify(outputAnalysisPropertyUnoperatorUnweighted_     )
    nullify(outputAnalysisPropertyUnoperatorCross_          )
    nullify(outputAnalysisWeightPropertyExtractorUnweighted_)
    nullify(outputAnalysisWeightPropertyExtractorCross_     )
    nullify(outputAnalysisWeightPropertyOperatorUnweighted_ )
    nullify(outputAnalysisWeightPropertyOperatorCross_      )
    nullify(outputAnalysisDistributionOperatorUnweighted_   )
    nullify(outputAnalysisDistributionOperatorCross_        )
    nullify(galacticFilterUnweighted_                       )
    nullify(galacticFilterCross_                            )
    nullify(outputTimesUnweighted_                          )
    nullify(outputTimesCross_                               )
    return
  end function meanFunction1DConstructorInternal
  
  subroutine meanFunction1DAnalyze(self,node,iOutput)
    !% Implement a meanFunction1D output analysis.
    implicit none
    class  (outputAnalysisMeanFunction1D), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(c_size_t                    ), intent(in   ) :: iOutput

    ! Analyze for all three volume functions.
    call self%volumeFunctionUnweighted%analyze(node,iOutput)
    call self%volumeFunctionWeighted  %analyze(node,iOutput)
    call self%volumeFunctionCross     %analyze(node,iOutput)
    return
  end subroutine meanFunction1DAnalyze

  subroutine meanFunction1DReduce(self,reduced)
    !% Implement a volumeFunction1D output analysis reduction.
    use Galacticus_Error
    implicit none
    class(outputAnalysisMeanFunction1D), intent(inout) :: self
    class(outputAnalysisClass         ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisMeanFunction1D)
       call self%volumeFunctionUnweighted%reduce(reduced%volumeFunctionUnweighted)
       call self%volumeFunctionWeighted  %reduce(reduced%volumeFunctionWeighted  )
       call self%volumeFunctionCross     %reduce(reduced%volumeFunctionCross     )
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine meanFunction1DReduce

  subroutine meanFunction1DFinalizeAnalysis(self)
    !% Finalize analysis of a {\normalfont \ttfamily meanFunction1D} output analysis.
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
    call self%volumeFunctionCross     %results(                                                       functionCovariance=crossCovariance     )
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

  subroutine meanFunction1DFinalize(self)
    !% Implement a {\normalfont \ttfamily meanFunction1D} output analysis finalization.
    use IO_HDF5
    use Galacticus_HDF5
    implicit none
    class(outputAnalysisMeanFunction1D), intent(inout) :: self
    type (hdf5Object                  )                :: analysesGroup, analysisGroup, &
         &                                                dataset

    ! Finalize the analysis.
    call meanFunction1DFinalizeAnalysis(self)
    ! Output the resulting mean function.
    !$ call hdf5Access%set()
    analysesGroup=galacticusOutputFile%openGroup('analyses'                         )
    analysisGroup=analysesGroup       %openGroup(char(self%label),char(self%comment))
    call analysisGroup%writeDataset  (          self%binCenter               ,char(self%propertyLabel)              ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%propertyUnits    )      ,'units'                                                                                                   )
    call dataset      %writeAttribute(          self%propertyUnitsInSI       ,'unitsInSI'                                                                                               )
    call dataset      %close()
    call analysisGroup%writeDataset  (          self%meanValue               ,char(self%    meanLabel)              ,char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )      ,'units'                                                                                                   )
    call dataset      %writeAttribute(          self%meanUnitsInSI           ,'unitsInSI'                                                                                               )
    call dataset      %close()
    call analysisGroup%writeDataset  (          self%meanCovariance          ,char(self%    meanLabel)//"Covariance",char(self%    meanComment)//" [covariance]",datasetReturned=dataset)
    call dataset      %writeAttribute("["//char(self%    meanUnits    )//"]²",'units'                                                                                                   )
    call dataset      %writeAttribute(          self%    meanUnitsInSI   **2 ,'unitsInSI'                                                                                               )
    call dataset      %close()
    call analysisGroup%close()
    call analysesGroup%close()
    !$ call hdf5Access%unset()    
    return
  end subroutine meanFunction1DFinalize

  subroutine meanFunction1DResults(self,binCenter,meanValue,meanCovariance)
    !% Implement a meanFunction1D output analysis finalization.
    use Memory_Management
    implicit none
    class           (outputAnalysisMeanFunction1D)                             , intent(inout)           :: self
    double precision                              , allocatable, dimension(:  ), intent(inout), optional :: binCenter     , meanValue
    double precision                              , allocatable, dimension(:,:), intent(inout), optional :: meanCovariance

    ! Finalize analysis.
    call meanFunction1DFinalizeAnalysis(self)
    ! Return results.
    if (present(binCenter         )) then
       if (allocated(binCenter         )) call deallocateArray(binCenter         )
       call allocateArray(binCenter         ,shape(self%binCenter         ))
       binCenter         =self%binCenter
    end if
    if (present(meanValue     )) then
       if (allocated(meanValue     )) call deallocateArray(meanValue     )
       call allocateArray(meanValue     ,shape(self%meanValue     ))
       meanValue     =self%meanValue
    end if
    if (present(meanCovariance)) then
       if (allocated(meanCovariance)) call deallocateArray(meanCovariance)
       call allocateArray(meanCovariance,shape(self%meanCovariance))
       meanCovariance=self%meanCovariance
    end if
    return
  end subroutine meanFunction1DResults

  double precision function meanFunction1DLogLikelihood(self)
    !% Return the log-likelihood of a meanFunction1D output analysis.
    use Linear_Algebra          , only : vector, matrix, assignment(=), operator(*)
    use Numerical_Constants_Math, only : Pi
    use Galacticus_Error        , only : Galacticus_Error_Report
    implicit none
    class           (outputAnalysisMeanFunction1D), intent(inout)                 :: self
    double precision                              , allocatable  , dimension(:,:) :: meanCovarianceCombined
    double precision                              , allocatable  , dimension(:  ) :: meanValueDifference
    type            (vector                      )                                :: residual
    type            (matrix                      )                                :: covariance            , covarianceInverse
    
    ! Check for existance of a target distribution.
    if (allocated(self%meanValueTarget)) then
       ! Finalize analysis.
       call meanFunction1DFinalizeAnalysis(self)
       ! Allocate workspaces.
       allocate(meanCovarianceCombined(size(self%binCenter),size(self%binCenter)))
       allocate(meanValueDifference   (size(self%binCenter)                     ))
       ! Find combined covariance and difference between model and target.
       meanValueDifference   =+self%meanValue            &
            &                 -self%meanValueTarget
       meanCovarianceCombined=+self%meanCovariance       &
            &                 +self%meanCovarianceTarget
       residual              = meanValueDifference
       covariance            = meanCovarianceCombined       
       ! Find the inverse covariance matrix of the combined model and target covariances.
       covarianceInverse=covariance%invert()
       ! Compute the log-likelihood.
       meanFunction1DLogLikelihood       =-0.5d0*(residual*(covarianceInverse*residual))
       if (self%likelihoodNormalize)                                                      &
            & meanFunction1DLogLikelihood=+meanFunction1DLogLikelihood                    &
            &                             -0.5d0*covariance%determinant()                 &
            &                             -0.5d0*dble(size(self%binCenter))*log(2.0d0*Pi)
    else
       meanFunction1DLogLikelihood=0.0d0
       call Galacticus_Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function meanFunction1DLogLikelihood

  subroutine meanFunction1DDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily meanFunction1D} analysis class.
    use Galacticus_Error
    implicit none
    class(outputAnalysisMeanFunction1D), intent(inout) :: self
    class(outputAnalysisClass         ), intent(  out) :: destination

    call self%outputAnalysisClass%deepCopy(destination)
    select type (destination)
    type is (outputAnalysisMeanFunction1D)
       destination%label              =self%label
       destination%comment            =self%comment
       destination%propertyLabel      =self%propertyLabel
       destination%propertyComment    =self%propertyComment
       destination%meanLabel          =self%meanLabel
       destination%meanComment        =self%meanComment
       destination%propertyUnits      =self%propertyUnits
       destination%meanUnits          =self%meanUnits
       destination%propertyUnitsInSI  =self%propertyUnitsInSI
       destination%meanUnitsInSI      =self%meanUnitsInSI
       destination%finalized          =self%finalized
       destination%likelihoodNormalize=self%likelihoodNormalize
       if (allocated(self%binCenter)) then
          allocate(destination%binCenter           (size(self%binCenter           ,dim=1)                                      ))
          destination%binCenter=self%binCenter
       end if
       if (allocated(self%meanValue)) then
          allocate(destination%meanValue           (size(self%meanValue           ,dim=1)                                      ))
          destination%meanValue=self%meanValue
       end if
       if (allocated(self%meanValueTarget)) then
          allocate(destination%meanValueTarget     (size(self%meanValueTarget     ,dim=1)                                      ))
          destination%meanValueTarget=self%meanValueTarget
       end if
       if (allocated(self%meanCovariance)) then
          allocate(destination%meanCovariance      (size(self%meanCovariance      ,dim=1),size(self%meanCovariance      ,dim=2)))
          destination%meanCovariance=self%meanCovariance
       end if
       if (allocated(self%meanCovarianceTarget)) then
          allocate(destination%meanCovarianceTarget(size(self%meanCovarianceTarget,dim=1),size(self%meanCovarianceTarget,dim=2)))
          destination%meanCovarianceTarget=self%meanCovarianceTarget
       end if
       call self%volumeFunctionUnweighted%deepCopy(destination%volumeFunctionUnweighted)
       call self%volumeFunctionWeighted  %deepCopy(destination%volumeFunctionWeighted  )
       call self%volumeFunctionCross     %deepCopy(destination%volumeFunctionCross     )
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine meanFunction1DDeepCopy
