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
  Implements an output analysis class that computes the stellar mass-halo mass relation in the Local
  Group.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupOccupationFraction">
   <description>An output analysis class for the Local Group stellar mass-halo mass relation.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupOccupationFraction
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
     final     ::                  localGroupOccupationFractionDestructor
     procedure :: analyze       => localGroupOccupationFractionAnalyze
     procedure :: finalize      => localGroupOccupationFractionFinalize
     procedure :: reduce        => localGroupOccupationFractionReduce
     procedure :: logLikelihood => localGroupOccupationFractionLogLikelihood
  end type outputAnalysisLocalGroupOccupationFraction

  interface outputAnalysisLocalGroupOccupationFraction
     !!{
     Constructors for the ``localGroupOccupationFraction'' output analysis class.
     !!}
     module procedure localGroupOccupationFractionConstructorParameters
     module procedure localGroupOccupationFractionConstructorInternal
  end interface outputAnalysisLocalGroupOccupationFraction

contains

  function localGroupOccupationFractionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``localGroupOccupationFraction'' output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupOccupationFraction)                              :: self
    type            (inputParameters                           ), intent(inout)               :: parameters
    class           (outputTimesClass                          ), pointer                     :: outputTimes_
    double precision                                            , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient               , systematicErrorPolynomialCoefficient, &
         &                                                                                       massStellarSystematicErrorPolynomialCoefficient
    integer                                                                                   :: covarianceBinomialBinsPerDecade
    double precision                                                                          :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum   , &
         &                                                                                       randomErrorMinimum                             , randomErrorMaximum
    type            (varying_string                            )                              :: positionType

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
    self=outputAnalysisLocalGroupOccupationFraction(outputTimes_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,massStellarSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function localGroupOccupationFractionConstructorParameters

  function localGroupOccupationFractionConstructorInternal(outputTimes_,positionType,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,massStellarSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !!{
    Constructor for the ``localGroupOccupationFraction'' output analysis class for internal use.
    !!}
    use :: Galactic_Filters                        , only : filterList                                          , galacticFilterAll                           , galacticFilterHaloNotIsolated         , galacticFilterHostMassRange                    , &
          &                                                 galacticFilterSurveyGeometry                        , enumerationPositionTypeType
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
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10             , outputAnalysisPropertyOperatorLog10         , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 outputAnalysisPropertyOperatorBoolean               , outputAnalysisPropertyOperatorFilterHighPass, propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupOccupationFraction         )                                :: self
    integer                                                              , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )                 :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum                                 , &
         &                                                                                                  randomErrorMinimum                                         , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient                              , &
         &                                                                                                  massStellarSystematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                        ), intent(in   )                 :: positionType
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    type            (nodePropertyExtractorMassBasic                     )               , pointer        :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellar                   )               , pointer        :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_           , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorLog10                )               , pointer        :: outputAnalysisPropertyOperatorLog10_                       , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorFilterHighPass       )               , pointer        :: outputAnalysisWeightPropertyOperatorFilterHighPass_
    type            (outputAnalysisPropertyOperatorBoolean              )               , pointer        :: outputAnalysisWeightPropertyOperatorBoolean_
    type            (outputAnalysisPropertyOperatorSequence             )               , pointer        :: outputAnalysisPropertyOperator_                            , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10            )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling            )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer        :: outputAnalysisDistributionOperator_
    type            (surveyGeometryFullSky                              )               , pointer        :: surveyGeometry_
    type            (galacticFilterHaloNotIsolated                      )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                        )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterSurveyGeometry                       )               , pointer        :: galacticFilterSurveyGeometry_
    type            (galacticFilterAll                                  )               , pointer        :: galacticFilter_
    type            (filterList                                         )               , pointer        :: filters_
    type            (propertyOperatorList                               )               , pointer        :: operators_                                                  , weightPropertyOperators_
    double precision                                                     , allocatable  , dimension(:  ) :: masses                                                      , functionValueTarget                                              , &
         &                                                                                                  massHaloData                                                , fractionOccupationData
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                                , functionCovarianceTarget
    double precision                                                     , parameter                     :: bufferWidthLogarithmic                          =+3.0d+0    , errorZeroPoint                                       =10.00d00   , &
         &                                                                                                  massStellarErrorPolynomialZeroPoint             =-4.0d+0
    integer         (c_size_t                                           ), parameter                     :: bufferCountMinimum                              = 5_c_size_t
    double precision                                                     , parameter                     :: radiusOuter                                     =+3.0d-1
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    double precision                                                     , parameter                     :: massStellarThreshold                            =+1.0d+0
    integer         (c_size_t                                           )                                :: i                                                           , bufferCount
    type            (hdf5Object                                         )                                :: fileData
    !![
    <constructorAssign variables="*outputTimes_, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, massStellarSystematicErrorPolynomialCoefficient, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, randomErrorMinimum, randomErrorMaximum, positionType"/>
    !!]
    
    ! Construct the target distribution.
    !$ call hdf5Access%set  ()
    call fileData%openFile(char(inputPath(pathTypeDataStatic))//"observations/stellarHaloMassRelation/fractionOccupation_Local_Group_Nadler2020.hdf5",readOnly=.true.)
    call fileData%readDataset('massHalo'          ,massHaloData          )
    call fileData%readDataset('fractionOccupation',fractionOccupationData)
    call fileData%close      (                                           )
    !$ call hdf5Access%unset()
    ! Construct mass bins.
    allocate(masses                  (size(massHaloData)                   ))
    allocate(functionValueTarget     (size(massHaloData)                   ))
    allocate(functionCovarianceTarget(size(massHaloData),size(massHaloData)))
    functionCovarianceTarget=0.0d0
    masses             =log10(massHaloData          )
    functionValueTarget=      fractionOccupationData
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
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10              (                                                                                 )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorFilterHighPass_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorFilterHighPass_"    constructor="outputAnalysisPropertyOperatorFilterHighPass    (log10(massStellarThreshold)                                                       )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorBoolean_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorBoolean_"           constructor="outputAnalysisPropertyOperatorBoolean           (.false.                                                                           )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    allocate(weightPropertyOperators_%next%next                    )
    allocate(weightPropertyOperators_%next%next%next               )
    weightPropertyOperators_               %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    weightPropertyOperators_%next          %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    weightPropertyOperators_%next%next     %operator_ => outputAnalysisWeightPropertyOperatorFilterHighPass_
    weightPropertyOperators_%next%next%next%operator_ => outputAnalysisWeightPropertyOperatorBoolean_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                                                           )"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"                   constructor="outputAnalysisPropertyOperatorLog10            (                                                                                   )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_"       constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorZeroPoint              ,systematicErrorPolynomialCoefficient                  )"/>
    !!]
    allocate(operators_     )
    allocate(operators_%next)
    operators_     %operator_ => outputAnalysisPropertyOperatorLog10_
    operators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                        constructor="outputAnalysisPropertyOperatorSequence         (operators_                                                                         )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                      constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                   )"/>
    !!]
    ! Create a subsampling weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                          constructor="outputAnalysisWeightOperatorSubsampling        (                                                                                   )"/>
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
    allocate(filters_          )
    allocate(filters_%next     )
    allocate(filters_%next%next)
    filters_          %filter_ => galacticFilterHaloNotIsolated_ 
    filters_%next     %filter_ => galacticFilterHostMassRange_
    filters_%next%next%filter_ => galacticFilterSurveyGeometry_
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
	   outputAnalysisMeanFunction1D(                                                                  &amp;
	   &amp;                        var_str('localGroupOccupationFraction'                         ), &amp;
	   &amp;                        var_str('Subhalo occupation fraction of Local Group satellites'), &amp;
	   &amp;                        var_str('massHalo'                                             ), &amp;
	   &amp;                        var_str('Halo mass at the bin center'                          ), &amp;
	   &amp;                        var_str('M☉'                                                   ), &amp;
	   &amp;                        massSolar                                                       , &amp;
	   &amp;                        var_str('fractionOccupation'                                   ), &amp;
           &amp;                        var_str('Occupation fraction'                                  ), &amp;
           &amp;                        var_str('dimensionless'                                        ), &amp;
           &amp;                        0.0d0                                                           , &amp;
	   &amp;                        masses                                                          , &amp;
	   &amp;                        bufferCount                                                     , &amp;
	   &amp;                        outputWeight                                                    , &amp;
	   &amp;                        nodePropertyExtractor_                                          , &amp;
           &amp;                        outputAnalysisWeightPropertyExtractor_                          , &amp;
	   &amp;                        outputAnalysisPropertyOperator_                                 , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                           , &amp;
	   &amp;                        outputAnalysisPropertyUnoperator_                               , &amp;
	   &amp;                        outputAnalysisWeightOperator_                                   , &amp;
	   &amp;                        outputAnalysisDistributionOperator_                             , &amp;
	   &amp;                        galacticFilter_                                                 , &amp;
	   &amp;                        outputTimes_                                                    , &amp;
	   &amp;                        outputAnalysisCovarianceModelBinomial                           , &amp;
	   &amp;                        covarianceBinomialBinsPerDecade                                 , &amp;
	   &amp;                        covarianceBinomialMassHaloMinimum                               , &amp;
	   &amp;                        covarianceBinomialMassHaloMaximum                               , &amp;
           &amp;                        likelihoodNormalize                                             , &amp;
           &amp;                        var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'                   ), &amp;
           &amp;                        var_str('$f_\mathrm{occupied}$'                                ), &amp;
           &amp;                        .true.                                                          , &amp;
           &amp;                        .false.                                                         , &amp;
           &amp;                        var_str('Nadler et al. (2020)'                                 ), &amp;
           &amp;                        functionValueTarget                                             , &amp;
           &amp;                        functionCovarianceTarget                                          &amp;
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
    <objectDestructor name="outputAnalysisWeightPropertyOperatorFilterHighPass_"   />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorBoolean_"          />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                        />
    <objectDestructor name="galacticFilterHostMassRange_"                          />
    <objectDestructor name="galacticFilterSurveyGeometry_"                         />
    <objectDestructor name="galacticFilter_"                                       />
    <objectDestructor name="surveyGeometry_"                                       />
    !!]
    nullify(filters_                )
    nullify(operators_              )
    nullify(weightPropertyOperators_)
    return
  end function localGroupOccupationFractionConstructorInternal

  subroutine localGroupOccupationFractionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily localGroupOccupationFraction} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupOccupationFraction), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"/>
    <objectDestructor name="self%outputTimes_"   />
    !!]
    return
  end subroutine localGroupOccupationFractionDestructor

  subroutine localGroupOccupationFractionAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily localGroupOccupationFraction} output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupOccupationFraction), intent(inout) :: self
    type   (treeNode                                  ), intent(inout) :: node
    integer(c_size_t                                  ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine localGroupOccupationFractionAnalyze

  subroutine localGroupOccupationFractionReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily localGroupOccupationFraction} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisLocalGroupOccupationFraction), intent(inout) :: self
    class(outputAnalysisClass                       ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupOccupationFraction)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine localGroupOccupationFractionReduce

  subroutine localGroupOccupationFractionFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily localGroupOccupationFraction} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisLocalGroupOccupationFraction), intent(inout)           :: self
    type (varying_string                            ), intent(in   ), optional :: groupName
    type (hdf5Object                                )               , target   :: analysesGroup, subGroup
    type (hdf5Object                                )               , pointer  :: inGroup
    type (hdf5Object                                )                          :: analysisGroup

    call self%outputAnalysis_%finalize(groupName)
    ! Overwrite log-likelihood with our own calculation.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('localGroupOccupationFraction','Subhalo occupation fraction of Local Group satellites')
    call analysisGroup%writeAttribute(self%logLikelihood(),'logLikelihood')
    call analysisGroup%close         (                                    )
    if (present(groupName)) &
         & call subGroup%close       (                                    )
    call analysesGroup%close         (                                    )
    !$ call hdf5Access%unset()
    return
  end subroutine localGroupOccupationFractionFinalize

  double precision function localGroupOccupationFractionLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily localGroupOccupationFraction} output analysis.
    !!}
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisLocalGroupOccupationFraction), intent(inout)              :: self
    integer                                                                                  :: i
    double precision                                            , dimension(:) , allocatable :: massHalo                                , fractionOccupation                           , &
         &                                                                                      massHaloLogarithmic
    double precision                                            , parameter                  :: massHaloHalfLogarithmicConstraint=7.51d0, massHaloHalfLogarithmicConstraintWidth=0.21d0
    double precision                                                                         :: massHaloHalfLogarithmic

    localGroupOccupationFractionLogLikelihood=logImprobable
    select type (outputAnalysis_ => self%outputAnalysis_)
    type is (outputAnalysisMeanFunction1D)
       ! Get the results.
       call outputAnalysis_%results(massHalo,fractionOccupation)
       massHaloLogarithmic=log10(massHalo)
       ! Find the 50% occupation mass.
       do i=size(massHalo),1,-1
          if (fractionOccupation(i) < 0.5d0) then
             if (i == size(massHalo)) then
                ! 50% occupation reached at the highest halo mass - return improbable.
                localGroupOccupationFractionLogLikelihood=logImprobable
             else
                ! Interpolate to get the precise halo mass for 50% occupation.
                massHaloHalfLogarithmic=+                           massHaloLogarithmic(i)  &
                     &                  +(+massHaloLogarithmic(i+1)-massHaloLogarithmic(i)) &
                     &                  /(+fractionOccupation (i+1)-fractionOccupation (i)) &
                     &                  *(+0.5d0                   -fractionOccupation (i))
                ! Nadler et al. (2020) report that the 50% mass is log₁₀(M₅₀/M☉)=7.51⁺⁰˙²¹₋₀.₀₀ - that is it is unconstrained below
                ! 7.51. Therefore, we compute our likelihood as a half-normal distribution.
                if (massHaloHalfLogarithmic > massHaloHalfLogarithmicConstraint) then
                   localGroupOccupationFractionLogLikelihood=-0.5d0                                                               &
                        &                                    *(massHaloHalfLogarithmic-massHaloHalfLogarithmicConstraint     )**2 &
                        &                                    /                         massHaloHalfLogarithmicConstraintWidth **2
                else
                   localGroupOccupationFractionLogLikelihood=+0.0d0
                end if
             end if
             exit
          end if
       end do
    end select
    return
  end function localGroupOccupationFractionLogLikelihood
