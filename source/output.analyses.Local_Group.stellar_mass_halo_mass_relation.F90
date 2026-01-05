!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements an output analysis class that computes the stellar mass-halo mass relation in the Local
  Group.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupStellarMassHaloMassRelation">
   <description>An output analysis class for the Local Group stellar mass-halo mass relation.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupStellarMassHaloMassRelation
     !!{
     An output analysis class for the Local Group stellar mass-halo mass relation.
     !!}
     private
     class           (outputAnalysisClass        ), pointer                     :: outputAnalysis_                                 => null()
     class           (outputTimesClass           ), pointer                     :: outputTimes_                                    => null()
     double precision                             , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient                         , systematicErrorPolynomialCoefficient, &
          &                                                                        massStellarSystematicErrorPolynomialCoefficient
     integer                                                                    :: covarianceBinomialBinsPerDecade
     double precision                                                           :: covarianceBinomialMassHaloMinimum                        , covarianceBinomialMassHaloMaximum   , &
          &                                                                        randomErrorMinimum                                       , randomErrorMaximum
     type            (enumerationPositionTypeType)                              :: positionType
   contains
     final     ::                  localGroupStellarMassHaloMassRelationDestructor
     procedure :: analyze       => localGroupStellarMassHaloMassRelationAnalyze
     procedure :: finalize      => localGroupStellarMassHaloMassRelationFinalize
     procedure :: reduce        => localGroupStellarMassHaloMassRelationReduce
     procedure :: logLikelihood => localGroupStellarMassHaloMassRelationLogLikelihood
  end type outputAnalysisLocalGroupStellarMassHaloMassRelation

  interface outputAnalysisLocalGroupStellarMassHaloMassRelation
     !!{
     Constructors for the \refClass{outputAnalysisLocalGroupStellarMassHaloMassRelation} output analysis class.
     !!}
     module procedure localGroupStellarMassHaloMassRelationConstructorParameters
     module procedure localGroupStellarMassHaloMassRelationConstructorInternal
  end interface outputAnalysisLocalGroupStellarMassHaloMassRelation

contains

  function localGroupStellarMassHaloMassRelationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupStellarMassHaloMassRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupStellarMassHaloMassRelation)                              :: self
    type            (inputParameters                                    ), intent(inout)               :: parameters
    class           (outputTimesClass                                   ), pointer                     :: outputTimes_
    double precision                                                     , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient               , systematicErrorPolynomialCoefficient, &
         &                                                                                                massStellarSystematicErrorPolynomialCoefficient
    integer                                                                                            :: covarianceBinomialBinsPerDecade
    double precision                                                                                   :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum   , &
         &                                                                                                randomErrorMinimum                             , randomErrorMaximum
    type            (varying_string                                     )                              :: positionType

    ! Check and read parameters.
    allocate(massStellarSystematicErrorPolynomialCoefficient(max(1,parameters%count('massStellarSystematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(           systematicErrorPolynomialCoefficient(max(1,parameters%count(           'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(               randomErrorPolynomialCoefficient(max(1,parameters%count(               'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for halo masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for halo masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.1d0]</defaultValue>
      <description>The coefficients of the random error polynomial for halo masses.</description>
    </inputParameter>
    <inputParameter>
      <name>massStellarSystematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>massStellarSystematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the stellar mass systematic error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing Local Group stellar mass-halo mass relation covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing Local Group stellar mass-halo mass relation covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing Local Group stellar mass-halo mass relation covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('orbital')</defaultValue>
      <description>The type of position to use in survey geometry filters.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=outputAnalysisLocalGroupStellarMassHaloMassRelation(outputTimes_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,massStellarSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function localGroupStellarMassHaloMassRelationConstructorParameters

  function localGroupStellarMassHaloMassRelationConstructorInternal(outputTimes_,positionType,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,massStellarSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupStellarMassHaloMassRelation} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                        , only : filterList                                          , galacticFilterAll                         , galacticFilterHaloNotIsolated         , galacticFilterHostMassRange                    , &
          &                                                 galacticFilterSurveyGeometry                        , galacticFilterStellarMass                 , enumerationPositionTypeType
    use :: Geometry_Surveys                        , only : surveyGeometryFullSky
    use :: HDF5_Access                             , only : hdf5Access
    use :: IO_HDF5                                 , only : hdf5Object
    use :: Input_Paths                             , only : inputPath                                           , pathTypeDataStatic
    use :: Node_Property_Extractors                , only : nodePropertyExtractorMassStellar                    , nodePropertyExtractorMassBasic
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Constants_Astronomical        , only : massSolar
    use :: Numerical_Ranges                        , only : Make_Range                                          , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10             , outputAnalysisPropertyOperatorLog10       , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupStellarMassHaloMassRelation   )                                :: self
    integer                                                                 , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                        , intent(in   )                 :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum                                 , &
         &                                                                                                     randomErrorMinimum                                         , randomErrorMaximum
    double precision                                                        , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient                              , &
         &                                                                                                     massStellarSystematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                           ), intent(in   )                 :: positionType
    class           (outputTimesClass                                      ), intent(inout), target         :: outputTimes_
    type            (nodePropertyExtractorMassBasic                        )               , pointer        :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellar                      )               , pointer        :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial       )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_           , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorLog10                   )               , pointer        :: outputAnalysisPropertyOperatorLog10_                       , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorSequence                )               , pointer        :: outputAnalysisPropertyOperator_                            , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10               )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling               )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml   )               , pointer        :: outputAnalysisDistributionOperator_
    type            (surveyGeometryFullSky                                 )               , pointer        :: surveyGeometry_
    type            (galacticFilterHaloNotIsolated                         )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                           )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterSurveyGeometry                          )               , pointer        :: galacticFilterSurveyGeometry_
    type            (galacticFilterStellarMass                             )               , pointer        :: galacticFilterStellarMass_
    type            (galacticFilterAll                                     )               , pointer        :: galacticFilter_
    type            (filterList                                            )               , pointer        :: filters_
    type            (propertyOperatorList                                  )               , pointer        :: operators_                                                  , weightPropertyOperators_
    double precision                                                        , allocatable  , dimension(:  ) :: masses                                                      , functionValueTarget                                              , &
         &                                                                                                     massHaloData                                                , massStellarData                                                  , &
         &                                                                                                     massStellarScatterData
    double precision                                                        , allocatable  , dimension(:,:) :: outputWeight                                                , functionCovarianceTarget
    double precision                                                        , parameter                     :: bufferWidthLogarithmic                          =+3.0d+0    , errorZeroPoint                                       =10.00d00   , &
         &                                                                                                     massStellarErrorPolynomialZeroPoint             =-4.0d+0
    integer         (c_size_t                                              ), parameter                     :: bufferCountMinimum                              = 5_c_size_t
    double precision                                                        , parameter                     :: radiusOuter                                     =+3.0d-1
    logical                                                                 , parameter                     :: likelihoodNormalize                             =.false.
    double precision                                                        , parameter                     :: massStellarThreshold                            =+1.0d-3
    integer         (c_size_t                                              )                                :: i                                                           , bufferCount
    type            (hdf5Object                                            )                                :: fileData
    !![
    <constructorAssign variables="*outputTimes_, positionType, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, massStellarSystematicErrorPolynomialCoefficient, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum"/>
    !!]
    
    ! Construct the target distribution.
    !$ call hdf5Access%set  ()
    fileData=hdf5Object(char(inputPath(pathTypeDataStatic))//"observations/stellarHaloMassRelation/stellarHaloMassRelation_Local_Group_Nadler2020.hdf5",readOnly=.true.)
    call fileData%readDataset('massHalo'          ,massHaloData          )
    call fileData%readDataset('massStellar'       ,massStellarData       )
    call fileData%readDataset('massStellarScatter',massStellarScatterData)
    !$ call hdf5Access%unset()
    ! Construct mass bins.
    allocate(masses                  (size(massHaloData)                   ))
    allocate(functionValueTarget     (size(massHaloData)                   ))
    allocate(functionCovarianceTarget(size(massHaloData),size(massHaloData)))
    functionCovarianceTarget=0.0d0
    do i=1,size(massHaloData)
       functionCovarianceTarget(i,i)=(massStellarScatterData(i)/massStellarData(i)/log(10.0))**2
    end do
    masses             =log10(massHaloData   )
    functionValueTarget=log10(massStellarData)
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                 constructor="nodePropertyExtractorMassBasic                (                                                                                    )"/>
    !!]
    ! Create a stellar mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                 constructor="nodePropertyExtractorMassStellar               (                                                                                   )"/>
    !!]
    ! Build a size weight property operator.
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(massStellarErrorPolynomialZeroPoint,massStellarSystematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                                    )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    weightPropertyOperators_     %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    weightPropertyOperators_%next%operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                                                            )"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"                   constructor="outputAnalysisPropertyOperatorLog10            (                                                                                    )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_"       constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorZeroPoint              ,systematicErrorPolynomialCoefficient                   )"/>
    !!]
    allocate(operators_     )
    allocate(operators_%next)
    operators_     %operator_ => outputAnalysisPropertyOperatorLog10_
    operators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                        constructor="outputAnalysisPropertyOperatorSequence         (operators_                                                                          )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                      constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                    )"/>
    !!]
    ! Create a subsampling weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                          constructor="outputAnalysisWeightOperatorSubsampling        (                                                                                    )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
    <constructor>
    outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
         &amp;                                           randomErrorMinimum              , &amp;
         &amp;                                           randomErrorMaximum              , &amp;
         &amp;                                           errorZeroPoint                  , &amp;
         &amp;                                           randomErrorPolynomialCoefficient  &amp;
         &amp;                                          )
     </constructor>
    </referenceConstruct>
    !!]
    !!]
    ! Build survey geometry.
    !! Use a full sky geometry out to the outer radius considered (300 kpc).
    allocate(surveyGeometry_               )
    !![
    <referenceConstruct object="surveyGeometry_"                constructor="surveyGeometryFullSky            (distanceMinimum=0.0d0,distanceMaximum=radiusOuter)"/>
    !!]
    ! Build filters which select satellites in a specified range of host halo mass.
    allocate(galacticFilterHaloNotIsolated_)
    !![
    <referenceConstruct object="galacticFilterHaloNotIsolated_" constructor="galacticFilterHaloNotIsolated    (                                                 )"/>
    !!]
    allocate(galacticFilterSurveyGeometry_ )
    !![
    <referenceConstruct object="galacticFilterSurveyGeometry_"  constructor="galacticFilterSurveyGeometry     (positionType         ,surveyGeometry_            )"/>
    !!]
    allocate(galacticFilterStellarMass_    )
    !![
    <referenceConstruct object="galacticFilterStellarMass_"     constructor="galacticFilterStellarMass        (massStellarThreshold                             )"/>
    !!]
    allocate(galacticFilterHostMassRange_  )
    !![
    <referenceConstruct object="galacticFilterHostMassRange_">
     <constructor>
      galacticFilterHostMassRange(                      &amp;
        &amp;                     massMinimum =1.00d12, &amp;
        &amp;                     massMaximum =2.00d12, &amp;
        &amp;                     useFinalHost=.true.   &amp;
        &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(filters_               )
    allocate(filters_%next          )
    allocate(filters_%next%next     )
    allocate(filters_%next%next%next)
    filters_               %filter_ => galacticFilterHaloNotIsolated_ 
    filters_%next          %filter_ => galacticFilterHostMassRange_
    filters_%next%next     %filter_ => galacticFilterSurveyGeometry_
    filters_%next%next%next%filter_ => galacticFilterStellarMass_
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterAll(filters_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(size(massHaloData),outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(i),0.0d0,absTol=1.0d-6)) then
          outputWeight(:,i)=1.0d0
       else
          outputWeight(:,i)=0.0d0
       end if
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity'//{introspection:location})
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidthLogarithmic/(masses(2)-masses(1)),kind=c_size_t)+1_c_size_t,bufferCountMinimum)
    ! Construct the analysis objects.
    allocate(outputAnalysisMeanFunction1D :: self%outputAnalysis_)
    select type (outputAnalysis_ => self%outputAnalysis_)
    type is (outputAnalysisMeanFunction1D)
       !![
       <referenceConstruct isResult="yes" object="outputAnalysis_">
	 <constructor>
	   outputAnalysisMeanFunction1D(                                                                      &amp;
	   &amp;                        var_str('localGroupStellarMassHaloMassRelation'                    ), &amp;
	   &amp;                        var_str('Stellar mass-halo mass relation of Local Group satellites'), &amp;
	   &amp;                        var_str('massHalo'                                                 ), &amp;
	   &amp;                        var_str('Halo mass at the bin center'                              ), &amp;
	   &amp;                        var_str('M☉'                                                       ), &amp;
	   &amp;                        massSolar                                                           , &amp;
	   &amp;                        var_str('massStellarMean'                                          ), &amp;
           &amp;                        var_str('Mean logarithmic stellar mass; ⟨log₁₀(M★/M☉)⟩'           ), &amp;
           &amp;                        var_str('dimensionless'                                            ), &amp;
           &amp;                        0.0d0                                                               , &amp;
	   &amp;                        masses                                                              , &amp;
	   &amp;                        bufferCount                                                         , &amp;
	   &amp;                        outputWeight                                                        , &amp;
	   &amp;                        nodePropertyExtractor_                                              , &amp;
           &amp;                        outputAnalysisWeightPropertyExtractor_                              , &amp;
	   &amp;                        outputAnalysisPropertyOperator_                                     , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                               , &amp;
	   &amp;                        outputAnalysisPropertyUnoperator_                                   , &amp;
	   &amp;                        outputAnalysisWeightOperator_                                       , &amp;
	   &amp;                        outputAnalysisDistributionOperator_                                 , &amp;
	   &amp;                        galacticFilter_                                                     , &amp;
	   &amp;                        outputTimes_                                                        , &amp;
	   &amp;                        outputAnalysisCovarianceModelBinomial                               , &amp;
	   &amp;                        covarianceBinomialBinsPerDecade                                     , &amp;
	   &amp;                        covarianceBinomialMassHaloMinimum                                   , &amp;
	   &amp;                        covarianceBinomialMassHaloMaximum                                   , &amp;
           &amp;                        likelihoodNormalize                                                 , &amp;
           &amp;                        var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'                       ), &amp;
           &amp;                        var_str('$\langle\log_{10}(M_\star/\mathrm{M}_\odot)\rangle$'      ), &amp;
           &amp;                        .true.                                                              , &amp;
           &amp;                        .false.                                                             , &amp;
           &amp;                        var_str('Nadler et al. (2020)'                                     ), &amp;
           &amp;                        functionValueTarget                                                 , &amp;
           &amp;                        functionCovarianceTarget                                              &amp;
	   &amp;                       )
	 </constructor>
       </referenceConstruct>
       !!]         
    end select
    !![
    <objectDestructor name="nodePropertyExtractor_"                                />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                        />
    <objectDestructor name="galacticFilterHostMassRange_"                          />
    <objectDestructor name="galacticFilterStellarMass_"                            />
    <objectDestructor name="galacticFilterSurveyGeometry_"                         />
    <objectDestructor name="galacticFilter_"                                       />
    <objectDestructor name="surveyGeometry_"                                       />
    !!]
    nullify(filters_                )
    nullify(operators_              )
    nullify(weightPropertyOperators_)
    return
  end function localGroupStellarMassHaloMassRelationConstructorInternal

  subroutine localGroupStellarMassHaloMassRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLocalGroupStellarMassHaloMassRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupStellarMassHaloMassRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"/>
    <objectDestructor name="self%outputTimes_"   />
    !!]
    return
  end subroutine localGroupStellarMassHaloMassRelationDestructor

  subroutine localGroupStellarMassHaloMassRelationAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassHaloMassRelation} output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupStellarMassHaloMassRelation), intent(inout) :: self
    type   (treeNode                                           ), intent(inout) :: node
    integer(c_size_t                                           ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine localGroupStellarMassHaloMassRelationAnalyze

  subroutine localGroupStellarMassHaloMassRelationReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassHaloMassRelation} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisLocalGroupStellarMassHaloMassRelation), intent(inout) :: self
    class(outputAnalysisClass                                ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupStellarMassHaloMassRelation)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine localGroupStellarMassHaloMassRelationReduce

  subroutine localGroupStellarMassHaloMassRelationFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassHaloMassRelation} output analysis finalization.
    !!}
    implicit none
    class(outputAnalysisLocalGroupStellarMassHaloMassRelation), intent(inout)           :: self
    type (varying_string                                     ), intent(in   ), optional :: groupName

    call self%outputAnalysis_%finalize(groupName)
    return
  end subroutine localGroupStellarMassHaloMassRelationFinalize

  double precision function localGroupStellarMassHaloMassRelationLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily localGroupStellarMassHaloMassRelation} output analysis.
    !!}
    implicit none
    class(outputAnalysisLocalGroupStellarMassHaloMassRelation), intent(inout) :: self

    localGroupStellarMassHaloMassRelationLogLikelihood=self%outputAnalysis_%logLikelihood()
    return
  end function localGroupStellarMassHaloMassRelationLogLikelihood
