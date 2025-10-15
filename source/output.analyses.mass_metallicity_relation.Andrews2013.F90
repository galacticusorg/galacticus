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
  Implements a mass-metallicity relation analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisMassMetallicityAndrews2013">
   <description>A mass-metallicity relation output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisMassMetallicityAndrews2013
     !!{
     A mass-metallicity relation output analysis class.
     !!}
     private
     double precision                                           , allocatable, dimension(:) :: systematicErrorPolynomialCoefficient                     , randomErrorPolynomialCoefficient, &
          &                                                                                    metallicitySystematicErrorPolynomialCoefficient
     class           (cosmologyFunctionsClass                  ), pointer                   :: cosmologyFunctions_                             => null()
     class           (starFormationRateDisksClass              ), pointer                   :: starFormationRateDisks_                         => null()
     class           (starFormationRateSpheroidsClass          ), pointer                   :: starFormationRateSpheroids_                     => null()
     class           (starFormationRateNuclearStarClustersClass), pointer                   :: starFormationRateNuclearStarClusters_           => null()
     double precision                                                                       :: randomErrorMinimum                                       , randomErrorMaximum              , &
          &                                                                                    fractionGasThreshold
   contains
     final :: massMetallicityAndrews2013Destructor
  end type outputAnalysisMassMetallicityAndrews2013
  
  interface outputAnalysisMassMetallicityAndrews2013
     !!{
     Constructors for the \refClass{outputAnalysisMassMetallicityAndrews2013} output analysis class.
     !!}
     module procedure massMetallicityAndrews2013ConstructorParameters
     module procedure massMetallicityAndrews2013ConstructorInternal
  end interface outputAnalysisMassMetallicityAndrews2013

contains

  function massMetallicityAndrews2013ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassMetallicityAndrews2013} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions                       , only : cosmologyFunctions             , cosmologyFunctionsClass
    use :: Input_Parameters                          , only : inputParameter                 , inputParameters
    use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
    use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
    use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
    implicit none
    type            (outputAnalysisMassMetallicityAndrews2013 )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    double precision                                           , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient           , randomErrorPolynomialCoefficient, &
         &                                                                                      metallicitySystematicErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                  ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass              ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), pointer                     :: starFormationRateNuclearStarClusters_
    double precision                                                                         :: randomErrorMinimum                             , randomErrorMaximum              , &
         &                                                                                      fractionGasThreshold

    ! Check and read parameters.
    allocate(metallicitySystematicErrorPolynomialCoefficient(max(1,parameters%count('metallicitySystematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(           systematicErrorPolynomialCoefficient(max(1,parameters%count(           'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(               randomErrorPolynomialCoefficient(max(1,parameters%count(               'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>metallicitySystematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>metallicitySystematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the metallicity systematic error polynomial.</description>
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
      <description>The maximum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionGasThreshold</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The minimum gas fraction to include in the sample.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisMassMetallicityAndrews2013(metallicitySystematicErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,fractionGasThreshold,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function massMetallicityAndrews2013ConstructorParameters

  function massMetallicityAndrews2013ConstructorInternal(metallicitySystematicErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,fractionGasThreshold,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassMetallicityAndrews2013} output analysis class for internal use.
    !!}
    use :: Abundances_Structure                      , only : Abundances_Index_From_Name                         , abundances
    use :: Atomic_Data                               , only : Atomic_Mass
    use :: Cosmology_Functions                       , only : cosmologyFunctionsClass                            , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                      , only : cosmologyParametersSimple
    use :: Display                                   , only : displayGreen                                       , displayReset
    use :: Galactic_Filters                          , only : filterList                                         , galacticFilterAll                              , galacticFilterGasFractionISM                   , galacticFilterStarFormationRate    , &
          &                                                   galacticFilterStellarMass
    use :: Error                                     , only : Error_Report
    use :: Input_Paths                               , only : inputPath                                          , pathTypeDataStatic
    use :: Geometry_Surveys                          , only : surveyGeometryLiWhite2009SDSS
    use :: HDF5_Access                               , only : hdf5Access
    use :: IO_HDF5                                   , only : hdf5Object
    use :: Node_Property_Extractors                  , only : nodePropertyExtractorMassStellar                   , nodePropertyExtractorMetallicityISM
    use :: Numerical_Constants_Astronomical          , only : massSolar
    use :: Output_Analyses_Options                   , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators    , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators        , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorFilterHighPass   , outputAnalysisPropertyOperatorLog10, &
          &                                                   outputAnalysisPropertyOperatorMetallicity12LogNH   , outputAnalysisPropertyOperatorSequence         , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Utilities                 , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators          , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                              , only : outputTimesClass
    use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
    use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
    use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
    use :: String_Handling                           , only : stringXMLFormat
    implicit none
    type            (outputAnalysisMassMetallicityAndrews2013           )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum                                            , &
         &                                                                                                  fractionGasThreshold
    double precision                                                     , intent(in   ), dimension(:  ) :: metallicitySystematicErrorPolynomialCoefficient         , systematicErrorPolynomialCoefficient                          , &
         &                                                                                                  randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    class           (starFormationRateDisksClass                        ), intent(in   ), target         :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target         :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target         :: starFormationRateNuclearStarClusters_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum                      =1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: masses                                                  , functionValueTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                            , functionCovarianceTarget
    type            (galacticFilterStellarMass                          ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterStarFormationRate                    ), pointer                       :: galacticFilterStarFormationRate_
    type            (galacticFilterGasFractionISM                       ), pointer                       :: galacticFilterGasFractionISM_
    type            (galacticFilterAll                                  ), pointer                       :: galacticFilter_
    type            (filterList                                         ), pointer                       :: filters_                                       , filter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_                         , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorFilterHighPass       ), pointer                       :: outputAnalysisPropertyOperatorFilterHighPass_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorMetallicity12LogNH   ), pointer                       :: outputAnalysisPropertyOperatorMetallicity12LogNH_
    type            (nodePropertyExtractorMassStellar                   ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMetallicityISM                ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_        , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_                                      , weightPropertyOperators_
    type            (surveyGeometryLiWhite2009SDSS                      ), pointer                       :: surveyGeometry_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.3d00
    double precision                                                     , parameter                     :: metallicityErrorPolynomialZeroPoint             = 8.8d00
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    integer         (c_size_t                                           )                                :: iBin                                                    , binCount
    type            (hdf5Object                                         )                                :: dataFile
    integer                                                                                              :: indexOxygen
    !![
    <constructorAssign variables="metallicitySystematicErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, fractionGasThreshold, *cosmologyFunctions_, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]
    
    ! Read masses at which fraction was measured.
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(inputPath(pathTypeDataStatic))//"observations/abundances/gasPhaseMetallicityAndrews2013.hdf5",readOnly=.true.                  )
    call dataFile%readDataset("mass"                                                                                            ,         masses                  )
    call dataFile%readDataset("metallicity"                                                                                     ,         functionValueTarget     )
    call dataFile%readDataset("metallicityCovariance"                                                                           ,         functionCovarianceTarget)
    call dataFile%close      (                                                                                                                                    )
    !$ call hdf5Access%unset()
    ! Convert masses fro logarithmic.
    masses=10.0d0**masses
    ! Construct survey geometry. Use a lower redshift limit than actually used by Andrews et al. to ensure that low mass bins have non-zero weight.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryLiWhite2009SDSS(redshiftMinimum=0.0d0,redshiftMaximum=0.25d0,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    binCount=size(masses,kind=c_size_t)
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
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample and with a
    ! selection for the star forming main sequence (based on that proposed by Renzini & Peng; 2015;
    ! http://adsabs.harvard.edu/abs/2015ApJ...801L..29R with a downward shift of 0.6dex to allow for the width of the
    ! distribution).
    allocate(galacticFilterStellarMass_                            )
    !![
    <referenceConstruct object="galacticFilterStellarMass_" constructor="galacticFilterStellarMass(massThreshold=1.00d7)"/>
    !!]
    allocate(galacticFilterStarFormationRate_                      )
    !![
    <referenceConstruct object="galacticFilterStarFormationRate_">
     <constructor>
      galacticFilterStarFormationRate(                                                                             &amp;
        &amp;                         logM0                                =0.00d0                               , &amp;
        &amp;                         logSFR0                              =0.76d0                               , &amp;
        &amp;                         logSFR1                              =0.76d0                               , &amp;
        &amp;                         starFormationRateDisks_              =starFormationRateDisks_              , &amp;
        &amp;                         starFormationRateSpheroids_          =starFormationRateSpheroids_          , &amp;
        &amp;                         starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_  &amp;        
        &amp;                        )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(galacticFilterGasFractionISM_                         )
    !![
    <referenceConstruct object="galacticFilterGasFractionISM_" constructor="galacticFilterGasFractionISM(fractionGasThreshold=fractionGasThreshold)"/>
    !!]
    allocate(filters_                                              )
    filter_ => filters_
    filter_%filter_ =>galacticFilterStellarMass_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterStarFormationRate_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterGasFractionISM_
    allocate(galacticFilter_                                       )
    !![
    <referenceConstruct object="galacticFilter_"                                   constructor="galacticFilterAll                                      (filters_                                                                  )"/>
    !!]
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                     constructor="outputAnalysisWeightOperatorIdentity                   (                                                                          )"/>
    !!]
    ! Build luminosity distance, systematic, and log10() property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"  constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc        (cosmologyFunctions_     ,cosmologyFunctionsData              ,outputTimes_)"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_"  constructor="outputAnalysisPropertyOperatorSystmtcPolynomial        (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient             )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"              constructor="outputAnalysisPropertyOperatorLog10                    (                                                                          )"/>
    !!]
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_          %operator_           => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next     %operator_           => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next%operator_           => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                   constructor="outputAnalysisPropertyOperatorSequence                 (propertyOperators_                                                        )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorRandomErrorPlynml(                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a metallicity weight property operator.
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial (metallicityErrorPolynomialZeroPoint,metallicitySystematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(outputAnalysisPropertyOperatorMetallicity12LogNH_     )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorMetallicity12LogNH_">
     <constructor>
      outputAnalysisPropertyOperatorMetallicity12LogNH(                            &amp;
        &amp;                                          Atomic_Mass(shortLabel="O") &amp;
        &amp;                                         )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(outputAnalysisPropertyOperatorFilterHighPass_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorFilterHighPass_"          constructor="outputAnalysisPropertyOperatorFilterHighPass    (0.0d0                                                        )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    allocate(weightPropertyOperators_%next%next                    )
    weightPropertyOperators_          %operator_ => outputAnalysisPropertyOperatorMetallicity12LogNH_
    weightPropertyOperators_%next     %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    weightPropertyOperators_%next%next%operator_ => outputAnalysisPropertyOperatorFilterHighPass_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence          (weightPropertyOperators_                                     )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                      constructor="outputAnalysisPropertyOperatorAntiLog10         (                                                             )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                      )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                 constructor="nodePropertyExtractorMassStellar                (                                                             )"/>
    !!]
    ! Find the index for the oxygen abundance.
    indexOxygen=Abundances_Index_From_Name("O")
    if (indexOxygen < 0)                                                                                           &
         & call Error_Report(                                                                                      &
         &                   'oxygen abundance is required for this analysis'                         //char(10)// &
         &                   displayGreen()//'HELP:'//displayReset()                                            // &
         &                   ' you can track oxygen abundance by including:'                //char(10)//char(10)// &
         &                   stringXMLFormat('<elementsToTrack value="O"/>',indentInitial=6)//char(10)//char(10)// &
         &                   ' in your parameter file'                                                          // &
         &                   {introspection:location}                                                              &
         &                  )
    ! Create an ISM metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                 constructor="nodePropertyExtractorMetallicityISM             (Abundances_Index_From_Name('O')                              )"/>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                         &
         &                                                         var_str('massMetallicityAndrews2013'                  ), &
         &                                                         var_str('Mass-metallicity relation'                   ), &
         &                                                         var_str('massStellar'                                 ), &
         &                                                         var_str('Stellar mass'                                ), &
         &                                                         var_str('M☉'                                          ), &
         &                                                         massSolar                                              , &
         &                                                         var_str('metallicityMean'                             ), &
         &                                                         var_str('Mean metallicity'                            ), &
         &                                                         var_str('dimensionless'                               ), &
         &                                                         0.0d0                                                  , &
         &                                                         log10(masses)                                          , &
         &                                                         bufferCount                                            , &
         &                                                         outputWeight                                           , &
         &                                                         nodePropertyExtractor_                                 , &
         &                                                         outputAnalysisWeightPropertyExtractor_                 , &
         &                                                         outputAnalysisPropertyOperator_                        , &
         &                                                         outputAnalysisWeightPropertyOperator_                  , &
         &                                                         outputAnalysisPropertyUnoperator_                      , &
         &                                                         outputAnalysisWeightOperator_                          , &
         &                                                         outputAnalysisDistributionOperator_                    , &
         &                                                         galacticFilter_                                        , &
         &                                                         outputTimes_                                           , &
         &                                                         outputAnalysisCovarianceModelBinomial                  , &
         &                                                         covarianceBinomialBinsPerDecade                        , &
         &                                                         covarianceBinomialMassHaloMinimum                      , &
         &                                                         covarianceBinomialMassHaloMaximum                      , &
         &                                                         likelihoodNormalize                                    , &
         &                                                         var_str('$M_\star/\mathrm{M}_\odot$')                  , &
         &                                                         var_str('$\langle 12+[\mathrm{O}/\mathrm{H}] \rangle$'), &
         &                                                         .true.                                                 , &
         &                                                         .false.                                                , &
         &                                                         var_str('Andrews \& Martini (2013)')                   , &
         &                                                         functionValueTarget                                    , &
         &                                                         functionCovarianceTarget                                 &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="galacticFilter_"                                       />
    <objectDestructor name="galacticFilterStellarMass_"                            />
    <objectDestructor name="galacticFilterStarFormationRate_"                      />
    <objectDestructor name="galacticFilterGasFractionISM_"                         />
    <objectDestructor name="surveyGeometry_"                                       />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"      />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyOperatorMetallicity12LogNH_"     />
    <objectDestructor name="outputAnalysisPropertyOperatorFilterHighPass_"         />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="nodePropertyExtractor_"                                />
    <objectDestructor name="cosmologyParametersData"                               />
    <objectDestructor name="cosmologyFunctionsData"                                />
    !!]
    nullify(propertyOperators_      )
    nullify(weightPropertyOperators_)
    nullify(filter_                 )
    nullify(filters_                )
    return
  end function massMetallicityAndrews2013ConstructorInternal

  subroutine massMetallicityAndrews2013Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMassMetallicityAndrews2013} output analysis class.
    !!}
    implicit none
    type(outputAnalysisMassMetallicityAndrews2013), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                  />
    <objectDestructor name="self%outputTimes_"                         />
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine massMetallicityAndrews2013Destructor
