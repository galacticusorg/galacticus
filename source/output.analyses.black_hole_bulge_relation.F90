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
  Implements a black hole-bulge mass relation analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisBlackHoleBulgeRelation">
   <description>A black hole-bulge mass relation output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisBlackHoleBulgeRelation
     !!{
     A black hole-bulge mass relation output analysis class.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_                  => null()
     double precision                         , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient          , randomErrorPolynomialCoefficient
     double precision                                                       :: randomErrorMinimum                            , randomErrorMaximum
   contains
     final :: blackHoleBulgeRelationDestructor
  end type outputAnalysisBlackHoleBulgeRelation

  interface outputAnalysisBlackHoleBulgeRelation
     !!{
     Constructors for the \refClass{outputAnalysisBlackHoleBulgeRelation} output analysis class.
     !!}
     module procedure blackHoleBulgeRelationConstructorParameters
     module procedure blackHoleBulgeRelationConstructorInternal
  end interface outputAnalysisBlackHoleBulgeRelation

contains

  function blackHoleBulgeRelationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisBlackHoleBulgeRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter    , inputParameters
    implicit none
    type            (outputAnalysisBlackHoleBulgeRelation)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    double precision                                      , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass             ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                    ), pointer                     :: outputTimes_
    double precision                                                                    :: randomErrorMinimum                  , randomErrorMaximum


    ! Check and read parameters.
    allocate(systematicErrorPolynomialCoefficient(max(1,parameters%count('systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorPolynomialCoefficient(max(1,parameters%count(    'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
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
      <defaultValue>[0.09d0]</defaultValue>
      <description>The coefficients of the random error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.09d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.09d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="outputTimes"        name="outputTimes_"        source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisBlackHoleBulgeRelation(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="outputTimes_"       />
    !!]
    return
  end function blackHoleBulgeRelationConstructorParameters

  function blackHoleBulgeRelationConstructorInternal(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisBlackHoleBulgeRelation} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                            , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterSpheroidStellarMass
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                          , pathTypeDataStatic
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassBlackHole                 , nodePropertyExtractorMassStellarSpheroid
    use :: Numerical_Comparison                  , only : Values_Agree
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorMinMax, &
          &                                               outputAnalysisPropertyOperatorSequence             , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    implicit none
    type            (outputAnalysisBlackHoleBulgeRelation               )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                    , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: masses                                                  , functionValueTarget                     , &
         &                                                                                                  functionErrorTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                            , functionCovarianceTarget
    type            (galacticFilterSpheroidStellarMass                  ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_                         , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_                    , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorMinMax               ), pointer                       :: outputAnalysisWeightPropertyOperatorMinMax_
    type            (nodePropertyExtractorMassStellarSpheroid           ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassBlackHole                 ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_                                      , weightPropertyOperators_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.3d00
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    integer         (c_size_t                                           )                                :: iOutput                                                 , i
    type            (hdf5Object                                         )                                :: dataFile
    type            (varying_string                                     )                                :: targetLabel
    !![
    <constructorAssign variables="systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, *cosmologyFunctions_, *outputTimes_"/>
    !!]
    
    !$ call hdf5Access%set()
    call dataFile%openFile     (char(inputPath(pathTypeDataStatic)//'/observations/blackHoles/blackHoleMassVsBulgeMass_KormendyHo2013.hdf5'),readOnly=.true.             )
    call dataFile%readDataset  ('massBulgeBinned'                                                                                           ,         masses             )
    call dataFile%readAttribute('label'                                                                                                     ,         targetLabel        )
    call dataFile%readDataset  ('massBlackHoleMean'                                                                                         ,         functionValueTarget)
    call dataFile%readDataset  ('massBlackHoleMeanError'                                                                                    ,         functionErrorTarget)
    call dataFile%close        (                                                                                                                                         )
    !$ call hdf5Access%unset()
    allocate(functionCovarianceTarget(size(functionErrorTarget),size(functionErrorTarget)))
    functionCovarianceTarget=0.0d0
    do i=1,size(functionErrorTarget)
       functionCovarianceTarget(i,i)=functionErrorTarget(i)**2
    end do
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(size(masses),outputTimes_%count()))
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('zero redshift output is required'//{introspection:location})
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple     (                            &amp;
        &amp;                        OmegaMatter    = 0.30000d0, &amp;
        &amp;                        OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                        HubbleConstant =70.50000d0, &amp;
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
    <referenceConstruct object="galacticFilter_"                                  constructor="galacticFilterSpheroidStellarMass              (massThreshold=1.0d8                                                       )"/>
    !!]
     ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorIdentity           (                                                                          )"/>
    !!]
    ! Build luminosity distance, systematic, and log10() property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_     ,cosmologyFunctionsData              ,outputTimes_)"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient             )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                          )"/>
    !!]
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_               %operator_      => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next          %operator_      => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next     %operator_      => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (propertyOperators_                                                        )"/>
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
    ! Build weight property operators.
    allocate(outputAnalysisWeightPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"       constructor="outputAnalysisPropertyOperatorLog10            (                                                                          )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorMinMax_           )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorMinMax_"      constructor="outputAnalysisPropertyOperatorMinMax           (thresholdMinimum=1.0d1,thresholdMaximum=huge(0.0d0)                       )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    weightPropertyOperators_         %operator_      => outputAnalysisWeightPropertyOperatorMinMax_
    weightPropertyOperators_%next    %operator_      => outputAnalysisWeightPropertyOperatorLog10_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"            constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                                                  )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                          )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                      )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorMassStellarSpheroid       (                                                                          )"/>
    !!]
    ! Create an ISM metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"           constructor="nodePropertyExtractorMassBlackHole             (                                                                          )"/>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                   &
         &                                                         var_str('blackHoleBulgeRelation'                                ), &
         &                                                         var_str('Black hole mass-bulge mass relation'                   ), &
         &                                                         var_str('massStellarSpheroid'                                   ), &
         &                                                         var_str('Stellar mass of spheroid'                              ), &
         &                                                         var_str('M☉'                                                    ), &
         &                                                         massSolar                                                        , &
         &                                                         var_str('massBlackHole'                                         ), &
         &                                                         var_str('Mean logarithmic (base-10) mass of central black hole' ), &
         &                                                         var_str('M☉'                                                    ), &
         &                                                         massSolar                                                        , &
         &                                                         masses                                                           , &
         &                                                         bufferCount                                                      , &
         &                                                         outputWeight                                                     , &
         &                                                         nodePropertyExtractor_                                           , &
         &                                                         outputAnalysisWeightPropertyExtractor_                           , &
         &                                                         outputAnalysisPropertyOperator_                                  , &
         &                                                         outputAnalysisWeightPropertyOperator_                            , &
         &                                                         outputAnalysisPropertyUnoperator_                                , &
         &                                                         outputAnalysisWeightOperator_                                    , &
         &                                                         outputAnalysisDistributionOperator_                              , &
         &                                                         galacticFilter_                                                  , &
         &                                                         outputTimes_                                                     , &
         &                                                         outputAnalysisCovarianceModelBinomial                            , &
         &                                                         covarianceBinomialBinsPerDecade                                  , &
         &                                                         covarianceBinomialMassHaloMinimum                                , &
         &                                                         covarianceBinomialMassHaloMaximum                                , &
         &                                                         likelihoodNormalize                                              , &
         &                                                         var_str('$M_{\star,\mathrm{bulge}}$ [M$_\odot$]'                ), &
         &                                                         var_str('$\langle \log_{10} M_\bullet/\mathrm{M}_\odot \rangle$'), &
         &                                                         .true.                                                           , &
         &                                                         .false.                                                          , &
         &                                                         targetLabel                                                      , &
         &                                                         functionValueTarget                                              , &
         &                                                         functionCovarianceTarget                                           &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="galacticFilter_"                                 />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorMinMax_"     />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"      />
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
    nullify(propertyOperators_      )
    nullify(weightPropertyOperators_)
    return
  end function blackHoleBulgeRelationConstructorInternal

  subroutine blackHoleBulgeRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisBlackHoleBulgeRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisBlackHoleBulgeRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine blackHoleBulgeRelationDestructor
