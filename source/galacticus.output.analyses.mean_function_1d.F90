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
  use               Galactic_Filters

  !# <outputAnalysis name="outputAnalysisMeanFunction1D" defaultThreadPrivate="yes">
  !#  <description>A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisClass) :: outputAnalysisMeanFunction1D
     !% A generic 1D mean function (i.e. mean value of some property weighted by number density of objects binned by some property) output analysis class.
     type            (varying_string                )          :: label                             , comment                         , &
          &                                                       propertyLabel                     , propertyComment                 , &
          &                                                       meanLabel                         , meanComment                     , &
          &                                                       propertyUnits                     , meanUnits
     double precision                                          :: propertyUnitsInSI                 , meanUnitsInSI
     type            (outputAnalysisVolumeFunction1D), pointer :: volumeFunctionUnweighted => null(), volumeFunctionWeighted => null(), &
          &                                                       volumeFunctionCross      => null()
   contains
     final     ::             meanFunction1DDestructor
     procedure :: analyze  => meanFunction1DAnalyze
     procedure :: finalize => meanFunction1DFinalize
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
    use Galacticus_Output_Times
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
    double precision                                         , dimension(:  ), allocatable :: binCenter                            , outputWeight
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

    ! Check and read parameters.
    unoperatorParameters=parameters%subParameters('unoperator',requireValue=.false.)
    weightParameters    =parameters%subParameters('weight'    ,requireValue=.false.)
    call allocateArray(binCenter   ,[int(parameters%count('binCenter'),kind=c_size_t)                               ])
    call allocateArray(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t)*Galacticus_Output_Time_Count()])
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*Galacticus_Output_Time_Count()) &
         & call Galacticus_Error_Report('meanFunction1DConstructorParameters','incorrect number of output weights provided')
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
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisPropertyExtractor_"       source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisWeightPropertyExtractor_" source="weightParameters"    />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"        source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisWeightPropertyOperator_"  source="weightParameters"    />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"      source="unoperatorParameters"/>
    !# <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"          source="parameters"          />
    !# <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"    source="parameters"          />
    !# <objectBuilder class="galacticFilter"                       name="galacticFilter_"                        source="parameters"          />
    ! Build the object.
    self=outputAnalysisMeanFunction1D(                                                                                                         &
         &                            label                                                                                                  , &
         &                            comment                                                                                                , &
         &                            propertyLabel                                                                                          , &
         &                            propertyComment                                                                                        , &
         &                            propertyUnits                                                                                          , &
         &                            propertyUnitsInSI                                                                                      , &
         &                            meanLabel                                                                                              , &
         &                            meanComment                                                                                            , &
         &                            meanUnits                                                                                              , &
         &                            meanUnitsInSI                                                                                          , &
         &                            binCenter                                                                                              , &
         &                            bufferCount                                                                                            , &
         &                            reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),Galacticus_Output_Time_Count()]), &
         &                            outputAnalysisPropertyExtractor_                                                                       , &
         &                            outputAnalysisWeightPropertyExtractor_                                                                 , &
         &                            outputAnalysisPropertyOperator_                                                                        , &
         &                            outputAnalysisWeightPropertyOperator_                                                                  , &
         &                            outputAnalysisPropertyUnoperator_                                                                      , &
         &                            outputAnalysisWeightOperator_                                                                          , &
         &                            outputAnalysisDistributionOperator_                                                                    , &
         &                            galacticFilter_                                                                                        , &
         &                            enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.)           , &
         &                            covarianceBinomialBinsPerDecade                                                                        , &
         &                            covarianceBinomialMassHaloMinimum                                                                      , &
         &                            covarianceBinomialMassHaloMaximum                                                                        &                      
         &                           )
    !# <inputParametersValidate source="parameters"/>
    return
  end function meanFunction1DConstructorParameters

  function meanFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,meanLabel,meanComment,meanUnits,meanUnitsInSI,binCenter,bufferCount,outputWeight,outputAnalysisPropertyExtractor_,outputAnalysisWeightPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisWeightPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,galacticFilter_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !% Constructor for the ``meanFunction1D'' output analysis class for internal use.
    use Memory_Management
    implicit none
    type            (outputAnalysisMeanFunction1D                )                                :: self
    type            (varying_string                              ), intent(in   )                 :: label                                  , comment                                  , &
         &                                                                                           propertyLabel                          , propertyComment                          , &
         &                                                                                           meanLabel                              , meanComment                              , &
         &                                                                                           propertyUnits                          , meanUnits
    double precision                                                                              :: propertyUnitsInSI                      , meanUnitsInSI
    double precision                                              , intent(in   ), dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                 :: bufferCount
    double precision                                              , intent(in   ), dimension(:,:) :: outputWeight
    class           (outputAnalysisPropertyExtractorClass        ), intent(in   ), target         :: outputAnalysisPropertyExtractor_       , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(in   ), target         :: outputAnalysisPropertyOperator_        , outputAnalysisPropertyUnoperator_        , &
         &                                                                                           outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(in   ), target         :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                         ), intent(in   ), target         :: galacticFilter_
    integer                                                       , intent(in   )                 :: covarianceModel
    integer                                                       , intent(in   ), optional       :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional       :: covarianceBinomialMassHaloMinimum      , covarianceBinomialMassHaloMaximum
    type            (outputAnalysisDistributionNormalizerIdentity), pointer                       :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisWeightOperatorSequence        ), pointer                       :: outputAnalysisWeightOperatorWeighted_  , outputAnalysisWeightOperatorCross_       , &
         &                                                                                           outputAnalysisWeightOperatorUnweighted_
    type            (weightOperatorList                          ), pointer                       :: weightOperatorWeight_                  , weightOperatorCross_                     , &
         &                                                                                           weightOperatorUnweighted_
    type            (outputAnalysisPropertyOperatorSequence      ), pointer                       :: weightOperatorCrossPropertyOperator_   , weightOperatorUnweightedPropertyOperator_
    type            (propertyOperatorList                        ), pointer                       :: propertyOperator_                      , propertyOperatorUnweighted_
    !# <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, meanLabel, meanComment, meanUnits, meanUnitsInSI"/>

    ! Build an identity distribution normalizer.
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerIdentity()
    ! Build weight operator that includes the weight property.
    allocate(                                        outputAnalysisWeightOperatorWeighted_               )
    allocate(                                        weightOperatorWeight_                               )
    allocate(                                        weightOperatorWeight_                %next          )
    allocate(outputAnalysisWeightOperatorProperty :: weightOperatorWeight_                %next%operator_)
    weightOperatorWeight_%operator_       => outputAnalysisWeightOperator_
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
    propertyOperator_%operator_       => outputAnalysisWeightPropertyOperator_ ! First operator in the sequence is whatever operator was passed to this constructor.
    select type (operator_ => propertyOperator_%next%operator_)                ! Second operator in the sequence is the square-root operator.
    type is (outputAnalysisPropertyOperatorSquareRoot)
       operator_=outputAnalysisPropertyOperatorSquareRoot()
    end select
    weightOperatorCrossPropertyOperator_=outputAnalysisPropertyOperatorSequence(propertyOperator_)
    ! Build a weight operator sequence - the first operator is whatever weight operator was passed to this constructor, the second
    ! is the operator which weights by the property value.
    weightOperatorCross_%operator_ => outputAnalysisWeightOperator_
    select type (operator_ => weightOperatorCross_%next%operator_)
    type is (outputAnalysisWeightOperatorProperty)
       operator_=outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractor_,weightOperatorCrossPropertyOperator_)
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
    propertyOperatorUnweighted_%operator_       => outputAnalysisWeightPropertyOperator_ ! First operator in the sequence is whatever operator was passed to this constructor.
    select type (operator_ => propertyOperatorUnweighted_%next%operator_)                ! Second operator in the sequence is the boolean operator.
    type is (outputAnalysisPropertyOperatorBoolean)
       operator_=outputAnalysisPropertyOperatorBoolean()
    end select
    weightOperatorUnweightedPropertyOperator_=outputAnalysisPropertyOperatorSequence(propertyOperatorUnweighted_)
    ! Build a weight operator sequence - the first operator is whatever weight operator was passed to this constructor, the second
    ! is the operator which weights by the property value.
    weightOperatorUnweighted_%operator_ => outputAnalysisWeightOperator_
    select type (operator_ => weightOperatorUnweighted_%next%operator_)
    type is (outputAnalysisWeightOperatorProperty)
       operator_=outputAnalysisWeightOperatorProperty(outputAnalysisWeightPropertyExtractor_,weightOperatorUnweightedPropertyOperator_)
    end select
    outputAnalysisWeightOperatorUnweighted_=outputAnalysisWeightOperatorSequence(weightOperatorUnweighted_)
    ! Build weighted, unweighted, and cross volume function objects.
    allocate(self%volumeFunctionWeighted  )
    allocate(self%volumeFunctionUnweighted)
    allocate(self%volumeFunctionCross     )
    self%volumeFunctionUnweighted=outputAnalysisVolumeFunction1D(                                         &
         &                                                       label                                  , &
         &                                                       comment                                , &
         &                                                       propertyLabel                          , &
         &                                                       propertyComment                        , &
         &                                                       propertyUnits                          , &
         &                                                       propertyUnitsInSI                      , &
         &                                                       meanLabel                              , &
         &                                                       meanComment                            , &
         &                                                       meanUnits                              , &
         &                                                       meanUnitsInSI                          , &
         &                                                       binCenter                              , &
         &                                                       bufferCount                            , &
         &                                                       outputWeight                           , &
         &                                                       outputAnalysisPropertyExtractor_       , &
         &                                                       outputAnalysisPropertyOperator_        , &
         &                                                       outputAnalysisPropertyUnoperator_      , &
         &                                                       outputAnalysisWeightOperatorUnweighted_, &
         &                                                       outputAnalysisDistributionOperator_    , &
         &                                                       outputAnalysisDistributionNormalizer_  , &
         &                                                       galacticFilter_                        , &
         &                                                       covarianceModel                        , &
         &                                                       covarianceBinomialBinsPerDecade        , &
         &                                                       covarianceBinomialMassHaloMinimum      , &
         &                                                       covarianceBinomialMassHaloMaximum        &
         &                                                      )
    self%volumeFunctionWeighted  =outputAnalysisVolumeFunction1D(                                         &
         &                                                       label                                  , &
         &                                                       comment                                , &
         &                                                       propertyLabel                          , &
         &                                                       propertyComment                        , &
         &                                                       propertyUnits                          , &
         &                                                       propertyUnitsInSI                      , &
         &                                                       meanLabel                              , &
         &                                                       meanComment                            , &
         &                                                       meanUnits                              , &
         &                                                       meanUnitsInSI                          , &
         &                                                       binCenter                              , &
         &                                                       bufferCount                            , &
         &                                                       outputWeight                           , &
         &                                                       outputAnalysisPropertyExtractor_       , &
         &                                                       outputAnalysisPropertyOperator_        , &
         &                                                       outputAnalysisPropertyUnoperator_      , &
         &                                                       outputAnalysisWeightOperatorWeighted_  , &
         &                                                       outputAnalysisDistributionOperator_    , &
         &                                                       outputAnalysisDistributionNormalizer_  , &
         &                                                       galacticFilter_                        , &
         &                                                       covarianceModel                        , &
         &                                                       covarianceBinomialBinsPerDecade        , &
         &                                                       covarianceBinomialMassHaloMinimum      , &
         &                                                       covarianceBinomialMassHaloMaximum        &
         &                                                      )
    self%volumeFunctionCross     =outputAnalysisVolumeFunction1D(                                         &
         &                                                       label                                  , &
         &                                                       comment                                , &
         &                                                       propertyLabel                          , &
         &                                                       propertyComment                        , &
         &                                                       propertyUnits                          , &
         &                                                       propertyUnitsInSI                      , &
         &                                                       meanLabel                              , &
         &                                                       meanComment                            , &
         &                                                       meanUnits                              , &
         &                                                       meanUnitsInSI                          , &
         &                                                       binCenter                              , &
         &                                                       bufferCount                            , &
         &                                                       outputWeight                           , &
         &                                                       outputAnalysisPropertyExtractor_       , &
         &                                                       outputAnalysisPropertyOperator_        , &
         &                                                       outputAnalysisPropertyUnoperator_      , &
         &                                                       outputAnalysisWeightOperatorCross_     , &
         &                                                       outputAnalysisDistributionOperator_    , &
         &                                                       outputAnalysisDistributionNormalizer_  , &
         &                                                       galacticFilter_                        , &
         &                                                       covarianceModel                        , &
         &                                                       covarianceBinomialBinsPerDecade        , &
         &                                                       covarianceBinomialMassHaloMinimum      , &
         &                                                       covarianceBinomialMassHaloMaximum        &
         &                                                      )
    ! Nullify objects to avoid destruction.
    nullify(outputAnalysisWeightOperatorWeighted_  )
    nullify(outputAnalysisWeightOperatorCross_     )
    nullify(outputAnalysisWeightOperatorUnweighted_)
    nullify(weightOperatorWeight_                  )
    nullify(weightOperatorCross_                   )
    nullify(weightOperatorUnweighted_              )
    return
  end function meanFunction1DConstructorInternal

  subroutine meanFunction1DDestructor(self)
    !% Destructor for  the ``meanFunction1D'' output analysis class.
    type(outputAnalysisMeanFunction1D), intent(inout) :: self
    
    !# <objectDestructor name="self%volumeFunctionUnweighted"/>
    !# <objectDestructor name="self%volumeFunctionWeighted"  />
    return
  end subroutine meanFunction1DDestructor
  
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

  subroutine meanFunction1DFinalize(self)
    !% Implement a meanFunction1D output analysis finalization.
    use IO_HDF5
    use Galacticus_HDF5
    implicit none
    class           (outputAnalysisMeanFunction1D), intent(inout)                 :: self
    double precision                              , allocatable  , dimension(:  ) :: binCenter           , unweightedValue   , &
         &                                                                           weightedValue       , meanValue
    double precision                              , allocatable  , dimension(:,:) :: unweightedCovariance, weightedCovariance, &
         &                                                                           crossCovariance     , meanCovariance
    integer         (c_size_t                    )                                :: i                   , j
    type            (hdf5Object                  )                                :: analysesGroup       , analysisGroup     , &
         &                                                                           dataset

    ! Retrieve results from our 1-D volume functions.
    call self%volumeFunctionUnweighted%results(binCenter=binCenter,functionValue=unweightedValue,functionCovariance=unweightedCovariance)
    call self%volumeFunctionWeighted  %results(                    functionValue=weightedValue  ,functionCovariance=weightedCovariance  )
    call self%volumeFunctionCross     %results(                                                  functionCovariance=crossCovariance     )
    ! Estimate covariance of ratio using Taylor series expansion approach
    ! (e.g. http://math.stackexchange.com/questions/40713/calculating-the-variance-of-the-ratio-of-random-variables).
    allocate(meanValue     (size(binCenter)                ))
    allocate(meanCovariance(size(binCenter),size(binCenter)))
    meanValue     =0.0d0
    meanCovariance=0.0d0
    do i=1,size(binCenter)
       if (unweightedValue(i) > 0.0d0) then
          meanValue(i)=weightedValue(i)/unweightedValue(i)
          do j=i,size(binCenter)
             if (unweightedValue(j) > 0.0d0) then
                meanCovariance(i,j)=+(                                                                                                      &
                     &                +  weightedCovariance(i,j)* unweightedValue(i)                 *unweightedValue(j)                    &
                     &                +unweightedCovariance(i,j)*                    weightedValue(i)                   *weightedValue(j)   &
                     &                -     crossCovariance(i,j)*(unweightedValue(i)*weightedValue(j)+unweightedValue(j)*weightedValue(i))  &
                     &               )                                                                                                      &
                     &              /unweightedValue(i)**2                                                                                  &
                     &              /unweightedValue(j)**2
                if (i == j) meanCovariance(i,j)=max(meanCovariance(i,j),0.0d0)
                meanCovariance(j,i)=meanCovariance(i,j)
            end if
          end do
       end if
    end do
    ! Output the resulting mean function.
    !$omp critical(HDF5_Access)
    analysesGroup=galacticusOutputFile%openGroup('analyses'                         )
    analysisGroup=analysesGroup       %openGroup(char(self%label),char(self%comment))
    call analysisGroup%writeDataset  (binCenter     ,char(self%propertyLabel)              ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%propertyUnits    )      ,'units'                                                                          )
    call dataset      %writeAttribute(          self%propertyUnitsInSI       ,'unitsInSI'                                                                      )
    call dataset      %close()
    call analysisGroup%writeDataset  (meanValue     ,char(self%    meanLabel)              ,char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )      ,'units'                                                                          )
    call dataset      %writeAttribute(          self%    meanUnitsInSI       ,'unitsInSI'                                                                      )
    call dataset      %close()
    call analysisGroup%writeDataset  (meanCovariance,char(self%    meanLabel)//"Covariance",char(self%    meanComment)//" [covariance]",datasetReturned=dataset)
    call dataset      %writeAttribute("["//char(self%    meanUnits    )//"]²",'units'                                                                          )
    call dataset      %writeAttribute(          self%    meanUnitsInSI   **2 ,'unitsInSI'                                                                      )
    call dataset      %close()
    call analysisGroup%close()
    call analysesGroup%close()
    !$omp end critical(HDF5_Access)    
    return
  end subroutine meanFunction1DFinalize
