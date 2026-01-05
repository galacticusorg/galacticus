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
  Implements a stellar vs halo mass relation analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisMorphologicalFractionGAMAMoffett2016">
   <description>A morphological fraction output analysis class for the analysis of \cite{moffett_galaxy_2016}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisMorphologicalFractionGAMAMoffett2016
     !!{
     A morphological fraction output analysis class.
     !!}
     private
     double precision                         , allocatable, dimension(:) :: countAllTarget                                , countEarlyTarget                , &
          &                                                                  functionErrorLowerTarget                      , functionErrorUpperTarget        , &
          &                                                                  systematicErrorPolynomialCoefficient          , randomErrorPolynomialCoefficient
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_                  => null()
     double precision                                                     :: ratioEarlyType                                , ratioEarlyTypeError             , &
          &                                                                  randomErrorMinimum                            , randomErrorMaximum
   contains
     final     ::                  morphologicalFractionGAMAMoffett2016Destructor
     procedure :: finalize      => morphologicalFractionGAMAMoffett2016Finalize
     procedure :: logLikelihood => morphologicalFractionGAMAMoffett2016LogLikelihood
  end type outputAnalysisMorphologicalFractionGAMAMoffett2016

  interface outputAnalysisMorphologicalFractionGAMAMoffett2016
     !!{
     Constructors for the \refClass{outputAnalysisMorphologicalFractionGAMAMoffett2016} output analysis class.
     !!}
     module procedure morphologicalFractionGAMAMoffett2016ConstructorParameters
     module procedure morphologicalFractionGAMAMoffett2016ConstructorInternal
  end interface outputAnalysisMorphologicalFractionGAMAMoffett2016

contains

  function morphologicalFractionGAMAMoffett2016ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMorphologicalFractionGAMAMoffett2016} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter    , inputParameters
    implicit none
    type            (outputAnalysisMorphologicalFractionGAMAMoffett2016)                              :: self
    type            (inputParameters                                   ), intent(inout)               :: parameters
    double precision                                                    , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                           ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                  ), pointer                     :: outputTimes_
    double precision                                                                                  :: ratioEarlyType                      , ratioEarlyTypeError             , &
         &                                                                                               randomErrorMinimum                  , randomErrorMaximum


    ! Check and read parameters.
    allocate(systematicErrorPolynomialCoefficient(max(1,parameters%count('systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorPolynomialCoefficient(max(1,parameters%count(    'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>ratioEarlyType</name>
      <defaultValue>0.5d0</defaultValue>
      <description>The minimum spheroid-to-total ratio for a galaxy to be classified as ``early-type'' when constructing the \gls{gama} early-type fraction function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ratioEarlyTypeError</name>
      <defaultValue>0.3d0</defaultValue>
      <description>The error in spheroid fraction to be used when constructing the \gls{gama} early-type fraction function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the random error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="outputTimes"        name="outputTimes_"        source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisMorphologicalFractionGAMAMoffett2016(ratioEarlyType,ratioEarlyTypeError,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="outputTimes_"       />
    !!]
    return
  end function morphologicalFractionGAMAMoffett2016ConstructorParameters

  function morphologicalFractionGAMAMoffett2016ConstructorInternal(ratioEarlyType,ratioEarlyTypeError,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMorphologicalFractionGAMAMoffett2016} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                            , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass
    use :: Input_Paths                           , only : inputPath                                          , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryBaldry2012GAMA
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassStellar                   , nodePropertyExtractorMassStellarMorphology
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorNormal, &
          &                                               outputAnalysisPropertyOperatorSequence             , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    use :: Statistics_Distributions              , only : distributionFunction1DBeta
    implicit none
    type            (outputAnalysisMorphologicalFractionGAMAMoffett2016   )                                :: self
    double precision                                                       , intent(in   )                 :: ratioEarlyType                                                         , ratioEarlyTypeError                     , &
         &                                                                                                    randomErrorMinimum                                                     , randomErrorMaximum
    double precision                                                       , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                                   , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                              ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                     ), intent(inout), target         :: outputTimes_
    integer                                                                , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                       , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.000d08             , covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                       , allocatable  , dimension(:  ) :: masses                                                                 , functionValueTarget
    double precision                                                       , allocatable  , dimension(:,:) :: outputWeight                                                           , functionCovarianceTarget
    type            (galacticFilterStellarMass                            ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml  ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity                 ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence               ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                  ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10              ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorNormal                 ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (nodePropertyExtractorMassStellar                     ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellarMorphology           ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc      ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial      ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                            ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                       ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                                 ), pointer                       :: propertyOperators_
    type            (surveyGeometryBaldry2012GAMA                         ), pointer                       :: surveyGeometry_
    logical                                                                , parameter                     :: likelihoodNormalize                             =.false.
    double precision                                                       , parameter                     :: errorPolynomialZeroPoint                        =11.300d00
    double precision                                                       , parameter                     :: confidenceLevel                                 = 0.683d00 ! 1-sigma confidence level
    double precision                                                       , parameter                     :: alpha                                           = 1.0d0-confidenceLevel
    integer         (c_size_t                                             ), parameter                     :: bufferCount                                     =10
    type            (distributionFunction1DBeta                           )                                :: betaDistributionLower                                                  , betaDistributionUpper
    integer         (c_size_t                                             )                                :: iBin                                                                   , binCount
    type            (hdf5Object                                           )                                :: dataFile
    double precision                                                                                       :: probit                                                                 , sqrtArg
    !![
    <constructorAssign variables="ratioEarlyType, ratioEarlyTypeError, systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, *cosmologyFunctions_"/>
    !!]
    
    ! Read masses at which fraction was measured.
    !$ call hdf5Access%set()
    dataFile=hdf5Object(char(inputPath(pathTypeDataStatic))//"observations/morphology/earlyTypeFractionGAMA.hdf5",readOnly=.true.)
    call dataFile%readDataset("mass"      ,masses               )
    call dataFile%readDataset("countEarly",self%countEarlyTarget)
    !$ call hdf5Access%unset()
    binCount=size(masses,kind=c_size_t)
    ! Compute confidence intervals on data. In each mass bin the quantity of interest is the probability, p, of a galaxy being
    ! early type. Since the observations look at N galaxies in the bin, and find k early types, we therefore expect k to be
    ! binomially distributed. Our estimate of p is just the fraction of galaxies observed that are early types. To get confidence
    ! intervals on this we need to find confidence intervals of the binomial distribution probability, p. Here we follow the
    ! Clopper-Pearson interval method
    ! (https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper–Pearson_interval) which determines the
    ! confidence interval using the β distribution.
    allocate(     functionValueTarget     (binCount         ))
    allocate(self%functionErrorLowerTarget(binCount         ))
    allocate(self%functionErrorUpperTarget(binCount         ))
    allocate(     functionCovarianceTarget(binCount,binCount))
    functionCovarianceTarget=0.0d0
    do iBin=1,binCount
       functionValueTarget          (iBin)=+self%countEarlyTarget(iBin) &
            &                              /self%countAllTarget  (iBin)
       if (self%countAllTarget(iBin) < 100.0d0) then
          ! Use the Clopper-Pearson interval: Clopper, C.; Pearson, E. S. (1934), Biometrika. 26: 04–413. doi:10.1093/biomet/26.4.404 ;
          ! https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval
          betaDistributionLower              = distributionFunction1DBeta(self%countEarlyTarget(iBin)      ,self%countAllTarget(iBin)-self%countEarlyTarget(iBin)+1.0d0)
          betaDistributionUpper              = distributionFunction1DBeta(self%countEarlyTarget(iBin)+1.0d0,self%countAllTarget(iBin)-self%countEarlyTarget(iBin)      )
          if (self%countEarlyTarget(iBin) <= 0.0d0                    ) then
             self%functionErrorLowerTarget(iBin)= 0.0d0
          else
             self%functionErrorLowerTarget(iBin)= betaDistributionLower%inverse(      0.5d0*alpha)
          end if
          if (self%countEarlyTarget(iBin) >= self%countAllTarget(iBin)) then
             self%functionErrorUpperTarget(iBin)= 1.0d0
          else
             self%functionErrorUpperTarget(iBin)= betaDistributionUpper%inverse(1.0d0-0.5d0*alpha)
          end if
       else
          ! Use Wilson score interval with continuity correction: Newcombe, R. G. (1998), Statistics in Medicine. 17 (8):
          ! 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. ;
          ! https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval_with_continuity_correction)
          probit=1.0d0-alpha/2.0d0
          sqrtArg=probit**2-1.0d0/self%countAllTarget(iBin)+4.0d0*self%countAllTarget(iBin)*functionValueTarget(iBin)*(1.0d0-functionValueTarget(iBin))+(4.0d0*functionValueTarget(iBin)-2.0d0)
          if (sqrtArg >= 0.0d0) then
             self%functionErrorLowerTarget(iBin)=max(0.0d0,(2*self%countAllTarget(iBin)*functionValueTarget(iBin)+probit**2-(probit*sqrt(sqrtArg)+1.0d0))/2.0d0/(self%countAllTarget(iBin)+probit**2))
          else
             self%functionErrorLowerTarget(iBin)=0.0d0
          end if
          sqrtArg=probit**2-1.0d0/self%countAllTarget(iBin)+4.0d0*self%countAllTarget(iBin)*functionValueTarget(iBin)*(1.0d0-functionValueTarget(iBin))-(4.0d0*functionValueTarget(iBin)-2.0d0)
          if (sqrtArg >= 0.0d0) then
             self%functionErrorUpperTarget(iBin)=min(1.0d0,(2*self%countAllTarget(iBin)*functionValueTarget(iBin)+probit**2+(probit*sqrt(sqrtArg)+1.0d0))/2.0d0/(self%countAllTarget(iBin)+probit**2))
          else
             self%functionErrorUpperTarget(iBin)=1.0d0
          end if
       end if
       functionCovarianceTarget(iBin,iBin)=(0.5d0*(self%functionErrorUpperTarget(iBin)-self%functionErrorLowerTarget(iBin)))**2
    end do
    self%functionErrorLowerTarget=-self%functionErrorLowerTarget+functionValueTarget
    self%functionErrorUpperTarget=+self%functionErrorUpperTarget-functionValueTarget
    ! Construct survey geometry.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryBaldry2012GAMA(cosmologyFunctions_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(binCount,outputTimes_%count()))
    do iBin=1,binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,masses(iBin))
    end do
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple     (                            &amp;
        &amp;                        OmegaMatter    = 0.30000d0, &amp;
        &amp;                        OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                        HubbleConstant =70.00000d0, &amp;
        &amp;                        temperatureCMB = 2.72548d0, &amp;
        &amp;                        OmegaBaryon    = 0.00000d0  &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
      cosmologyFunctionsMatterLambda(                            &amp;
        &amp;                        cosmologyParametersData     &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilter_                                       )
    !![
    <referenceConstruct object="galacticFilter_"                                  constructor="galacticFilterStellarMass                           (massThreshold=1.0d8                                                       )"/>
    !!]
     ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorIdentity                (                                                                          )"/>
    !!]
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc     (cosmologyFunctions_     ,cosmologyFunctionsData              ,outputTimes_)"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial     (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient             )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10                 (                                                                          )"/>
    !!]
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_          %operator_  => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next     %operator_  => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next%operator_  => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence              (propertyOperators_                                                        )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
     outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a weight property operator.
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_">
     <constructor>
     outputAnalysisPropertyOperatorNormal                (                                  &amp;
        &amp;                                             rangeLower  =ratioEarlyType     , &amp;
        &amp;                                             rangeUpper  =1.0d0              , &amp;
        &amp;                                             extentLower =0.0d0              , &amp;
        &amp;                                             extentUpper =1.0d0              , &amp;
        &amp;                                             rootVariance=ratioEarlyTypeError  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorAntiLog10   (                                                                          )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                      )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorMassStellar          (                                                                          )"/>
    !!]
    ! Create a morphology weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"           constructor="nodePropertyExtractorMassStellarMorphology(                                                                          )"/>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                 &
         &                                                         var_str('morphologicalFractionGAMAMoffett2016'), &
         &                                                         var_str('Early-type fraction'                 ), &
         &                                                         var_str('massStellar'                         ), &
         &                                                         var_str('Stellar mass'                        ), &
         &                                                         var_str('M☉'                                  ), &
         &                                                         massSolar                                      , &
         &                                                         var_str('earlyTypeFraction'                   ), &
         &                                                         var_str('Early-type fraction'                 ), &
         &                                                         var_str(' '                                   ), &
         &                                                         0.0d0                                          , &
         &                                                         log10(masses)                                  , &
         &                                                         bufferCount                                    , &
         &                                                         outputWeight                                   , &
         &                                                         nodePropertyExtractor_                         , &
         &                                                         outputAnalysisWeightPropertyExtractor_         , &
         &                                                         outputAnalysisPropertyOperator_                , &
         &                                                         outputAnalysisWeightPropertyOperator_          , &
         &                                                         outputAnalysisPropertyUnoperator_              , &
         &                                                         outputAnalysisWeightOperator_                  , &
         &                                                         outputAnalysisDistributionOperator_            , &
         &                                                         galacticFilter_                                , &
         &                                                         outputTimes_                                   , &
         &                                                         outputAnalysisCovarianceModelBinomial          , &
         &                                                         covarianceBinomialBinsPerDecade                , &
         &                                                         covarianceBinomialMassHaloMinimum              , &
         &                                                         covarianceBinomialMassHaloMaximum              , &
         &                                                         likelihoodNormalize                            , &
         &                                                         var_str('$M_\star/\mathrm{M}_\odot$')          , &
         &                                                         var_str('$f_\mathrm{early}$'        )          , &
         &                                                         .true.                                         , &
         &                                                         .false.                                        , &
         &                                                         var_str('Moffett et al. (2016)')               , &
         &                                                         functionValueTarget                            , &
         &                                                         functionCovarianceTarget                         &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                 />
    <objectDestructor name="galacticFilter_"                                 />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"           />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"          />
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="cosmologyParametersData"                         />
    <objectDestructor name="cosmologyFunctionsData"                          />
    !!]
    nullify(propertyOperators_)
    return
  end function morphologicalFractionGAMAMoffett2016ConstructorInternal

  subroutine morphologicalFractionGAMAMoffett2016Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMorphologicalFractionGAMAMoffett2016} output analysis class.
    !!}
    implicit none
    type(outputAnalysisMorphologicalFractionGAMAMoffett2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%outputTimes_"       />
    !!]
    return
  end subroutine morphologicalFractionGAMAMoffett2016Destructor

  double precision function morphologicalFractionGAMAMoffett2016LogLikelihood(self)
    !!{
    Return the log-likelihood of a morphologicalFractionGAMAMoffett2016 output analysis.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class  (outputAnalysisMorphologicalFractionGAMAMoffett2016), intent(inout) :: self
    integer                                                                    :: i

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Compute the log-likelihood. This assumes that the number of early types in each bin follows a binomial distribution. We do
    ! not account for the variance in the model expectation here - doing so would require evaluating the likelihood of the
    ! binomial distribution compounded with the distribution of the model mean. The distribution of the model mean could be
    ! described by a β distribution (since it too should follow a binomial [or maybe Poisson-binomial] distribution), and then
    ! the compound distribution is the β-binomial distribution. However, providing the model mean is determined with high
    ! precision, neglecting the distribution of the model mean here should not matter.
    morphologicalFractionGAMAMoffett2016LogLikelihood=0.0d0
    do i=1,size(self%countAllTarget)
       if (self%countEarlyTarget(i) > 0.0d0) then
          if (      +self%meanValue(i) > 0.0d0) then
             morphologicalFractionGAMAMoffett2016LogLikelihood=+morphologicalFractionGAMAMoffett2016LogLikelihood                               &
                  &                                            +                        self%countEarlyTarget(i) *log(      +self%meanValue(i))
          else
             morphologicalFractionGAMAMoffett2016LogLikelihood=+logImpossible
             exit
          end if
       end if
       if (self%countAllTarget(i)-self%countEarlyTarget(i) > 0.0d0) then
          if (+1.0d0-self%meanValue(i) > 0.0d0) then
             morphologicalFractionGAMAMoffett2016LogLikelihood=+morphologicalFractionGAMAMoffett2016LogLikelihood                               &
                  &                                            +(self%countAllTarget(i)-self%countEarlyTarget(i))*log(+1.0d0-self%meanValue(i))
          else
             morphologicalFractionGAMAMoffett2016LogLikelihood=+logImpossible
             exit
          end if
       end if
    end do
    return
  end function morphologicalFractionGAMAMoffett2016LogLikelihood

  subroutine morphologicalFractionGAMAMoffett2016Finalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily morphologicalFractionGAMAMoffett2016} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisMorphologicalFractionGAMAMoffett2016), intent(inout)           :: self
    type (varying_string                                    ), intent(in   ), optional :: groupName
    type (hdf5Object                                        )                          :: analysesGroup, analysisGroup, &
         &                                                                                dataset

    ! Finalize the analysis.
    call self%finalizeAnalysis()
    ! Output the resulting mean function.
    !$ call hdf5Access%set()
    analysesGroup=outputFile   %openGroup('analyses'                         )
    analysisGroup=analysesGroup%openGroup(char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute(     char(self%   comment   )                    ,'description'                                                                                                   )
    call analysisGroup%writeAttribute("function1D"                                     ,'type'                                                                                                          )
    call analysisGroup%writeAttribute(     char(self%   xAxisLabel)                    ,'xAxisLabel'                                                                                                    )
    call analysisGroup%writeAttribute(     char(self%   yAxisLabel)                    ,'yAxisLabel'                                                                                                    )
    call analysisGroup%writeAttribute(          self%   xAxisIsLog                     ,'xAxisIsLog'                                                                                                    )
    call analysisGroup%writeAttribute(          self%   yAxisIsLog                     ,'yAxisIsLog'                                                                                                    )
    call analysisGroup%writeAttribute(     char(self%propertyLabel)                    ,'xDataset'                                                                                                      )
    call analysisGroup%writeAttribute(     char(self%    meanLabel)                    ,'yDataset'                                                                                                      )
    call analysisGroup%writeAttribute(     char(self%    meanLabel)//"Target"          ,'yDatasetTarget'                                                                                                )
    call analysisGroup%writeAttribute(     char(self%    meanLabel)//"Covariance"      ,'yCovariance'                                                                                                   )
    call analysisGroup%writeAttribute(     char(self%    meanLabel)//"ErrorLowerTarget",'yErrorLowerTarget'                                                                                             )
    call analysisGroup%writeAttribute(     char(self%    meanLabel)//"ErrorUpperTarget",'yErrorUpperTarget'                                                                                             )
    ! Write computed datasets.
    call analysisGroup%writeDataset  (          self%binCenter                         ,char(self%propertyLabel)                    ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%propertyUnits    )                ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%propertyUnitsInSI                 ,'unitsInSI'                                                                                                     )
    call analysisGroup%writeDataset  (          self%meanValue                         ,char(self%    meanLabel)                    ,char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
    call analysisGroup%writeDataset  (          self%meanCovariance                    ,char(self%    meanLabel)//"Covariance"      ,char(self%    meanComment)//" [covariance]",datasetReturned=dataset)
    call dataset      %writeAttribute("["//char(self%    meanUnits    )//"]²"          ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%    meanUnitsInSI   **2           ,'unitsInSI'                                                                                                     )
    ! Include the log-likelihood and target dataset.
    call analysisGroup%writeAttribute(          self%logLikelihood()                   ,'logLikelihood'                                                                                                 )
    call analysisGroup%writeAttribute(     char(self%targetLabel      )                ,'targetLabel'                                                                                                   )
    call analysisGroup%writeDataset  (          self%meanValueTarget                   ,char(self%    meanLabel)//"Target"          ,char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
    call analysisGroup%writeDataset  (          self%functionErrorLowerTarget          ,char(self%    meanLabel)//"ErrorLowerTarget",char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
    call analysisGroup%writeDataset  (          self%functionErrorUpperTarget          ,char(self%    meanLabel)//"ErrorUpperTarget",char(self%    meanComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    meanUnits    )                ,'units'                                                                                                         )
    call dataset      %writeAttribute(          self%meanUnitsInSI                     ,'unitsInSI'                                                                                                     )
    !$ call hdf5Access%unset()
    return
  end subroutine morphologicalFractionGAMAMoffett2016Finalize
