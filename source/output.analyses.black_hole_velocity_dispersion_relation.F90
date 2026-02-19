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

!+    Contributions to this file made by Sachi Weerasooriya
  
  !!{
  Implements a black hole-velocity dispersion mass relation analysis class using the data from \cite{mcconnell_revisiting_2013}.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  !![
  <outputAnalysis name="outputAnalysisBlackHoleVelocityDispersionRelation">
   <description>A black hole-velocity dispersion mass relation output analysis class using the data from \cite{mcconnell_revisiting_2013}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisBlackHoleVelocityDispersionRelation
     !!{
     A black hole-velocity dispersion mass relation output analysis class using the data from \cite{mcconnell_revisiting_2013}.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer                     :: darkMatterHaloScale_                        => null()
     class           (cosmologyFunctionsClass ), pointer                     :: cosmologyFunctions_                         => null()
     double precision                          , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient                 , randomErrorPolynomialCoefficient
     double precision                                                        :: randomErrorMinimum                                   , randomErrorMaximum
     double precision                                                        :: toleranceRelative                   =1.0d-3
     logical                                                                 :: includeRadii                                         , integrationFailureIsFatal
  contains
     final :: blackHoleVelocityDispersionRelationDestructor
  end type outputAnalysisBlackHoleVelocityDispersionRelation

  interface outputAnalysisBlackHoleVelocityDispersionRelation
     !!{
     Constructors for the \refClass{outputAnalysisBlackHoleVelocityDispersionRelation} output analysis class.
     !!}
     module procedure blackHoleVelocityDispersionRelationConstructorParameters
     module procedure blackHoleVelocityDispersionRelationConstructorInternal
  end interface outputAnalysisBlackHoleVelocityDispersionRelation

contains

  function blackHoleVelocityDispersionRelationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisBlackHoleVelocityDispersionRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (outputAnalysisBlackHoleVelocityDispersionRelation)                              :: self
    type            (inputParameters                                  ), intent(inout)               :: parameters
    double precision                                                   , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient       , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                          ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                 ), pointer                     :: outputTimes_
    class           (darkMatterHaloScaleClass                         ), pointer                     :: darkMatterHaloScale_
    double precision                                                                                 :: randomErrorMinimum                         , randomErrorMaximum
    double precision                                                   , parameter                   :: toleranceRelative                   =1.0d-3

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
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the random error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.01d0</defaultValue>
      <description>The minimum random error for velocity dispersions.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.01d0</defaultValue>
      <description>The minimum random error for velocity dispersions.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisBlackHoleVelocityDispersionRelation(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_,toleranceRelative=1.0d-3,darkMatterHaloScale_=darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function blackHoleVelocityDispersionRelationConstructorParameters

  function blackHoleVelocityDispersionRelationConstructorInternal(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_,toleranceRelative,darkMatterHaloScale_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisBlackHoleVelocityDispersionRelation} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterSpheroidStellarMass
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                          , pathTypeDataStatic
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassBlackHole                 , nodePropertyExtractorVelocityDispersion        , nodePropertyExtractorScalarizer
    use :: Numerical_Comparison                  , only : Values_Agree
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Numerical_Constants_Prefixes          , only : kilo
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorMinMax, &
          &                                               outputAnalysisPropertyOperatorSequence             , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    implicit none
    type            (outputAnalysisBlackHoleVelocityDispersionRelation  )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum                                , &
         &                                                                                                  toleranceRelative
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                    , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    class           (darkMatterHaloScaleClass                           ), intent(in   ), target         :: darkMatterHaloScale_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum          =1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: velocities                                              , functionValueTarget                               , &
         &                                                                                                  functionErrorTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                            , functionCovarianceTarget
    type            (galacticFilterSpheroidStellarMass                  ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_                         , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_                    , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorMinMax               ), pointer                       :: outputAnalysisWeightPropertyOperatorMinMax_
    type            (nodePropertyExtractorVelocityDispersion            ), pointer                       :: nodePropertyExtractorVelocityDispersion_
    type            (nodePropertyExtractorMassBlackHole                 ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_                                      , weightPropertyOperators_
    type            (nodePropertyExtractorScalarizer     )               , pointer                       :: nodePropertyExtractor_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.3d0
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    integer         (c_size_t                                           )                                :: iOutput                                                 , i
    logical                                                                                              :: includeRadii                                            , integrationFailureisFatal
    type            (hdf5Object                                         )                                :: dataFile
    type            (varying_string                                     )                                :: targetLabel
    type            (varying_string                                     )               , dimension(1  ) :: radiusSpecifiers
    !![
    <constructorAssign variables="systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, *cosmologyFunctions_, *outputTimes_, toleranceRelative, *darkMatterHaloScale_"/>
    !!]
    
    !$ call hdf5Access%set()
    call dataFile%openFile     (char(inputPath(pathTypeDataStatic)//'/observations/blackHoles/blackHoleMassVsVelocityDispersion_McConnellMa2013.hdf5'),readOnly=.true.             )
    call dataFile%readDataset  ('velocityDispersionBinned'                                                                                           ,          velocities         )
    call dataFile%readAttribute('label'                                                                                                              ,          targetLabel        )
    call dataFile%readDataset  ('massBlackHoleMean'                                                                                                  ,          functionValueTarget)
    call dataFile%readDataset  ('massBlackHoleMeanError'                                                                                             ,          functionErrorTarget)
    call dataFile%close        (                                                                                                                                                   )
    !$ call hdf5Access%unset()
    allocate(functionCovarianceTarget(size(functionErrorTarget),size(functionErrorTarget)))
    functionCovarianceTarget=0.0d0
    do i=1,size(functionErrorTarget)
       functionCovarianceTarget(i,i)=functionErrorTarget(i)**2
    end do
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(size(velocities),outputTimes_%count()))
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
        &amp;                        OmegaMatter    = 0.27400d0, &amp;
        &amp;                        OmegaDarkEnergy= 0.72600d0, &amp;
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
    <referenceConstruct object="galacticFilter_"                                  constructor="galacticFilterSpheroidStellarMass              (massThreshold=1.0d4                                                       )"/>
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
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"       constructor="outputAnalysisPropertyOperatorLog10            (                                                                                             )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorMinMax_           )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorMinMax_"      constructor="outputAnalysisPropertyOperatorMinMax           (thresholdMinimum=1.0d1,thresholdMaximum=huge(0.0d0)                                          )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    weightPropertyOperators_         %operator_ => outputAnalysisWeightPropertyOperatorMinMax_
    weightPropertyOperators_%next    %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"            constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                                                                     )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                             )"/>
    !!]
    ! Create a velocity dispersion property extractor.
    allocate(nodePropertyExtractorVelocityDispersion_              )
    radiusSpecifiers         (1)='spheroidHalfMassRadius:all:stellar:lineOfSight:1.0'
    includeRadii                =.false.
    integrationFailureIsFatal   =.false.
    !![
    <referenceConstruct object="nodePropertyExtractorVelocityDispersion_"         constructor="nodePropertyExtractorVelocityDispersion       (radiusSpecifiers,includeRadii,integrationFailureIsFatal,toleranceRelative,darkMatterHaloScale_)"/>
    !!]
    allocate(nodePropertyExtractor_                                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorScalarizer               (1,1,nodePropertyExtractorVelocitydispersion_                                                  )"/>
    !!]
    ! Create an ISM metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"           constructor="nodePropertyExtractorMassBlackHole             (                                                                                              )"/>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                   &
         &                                                         var_str('blackHoleVelocityDispersionRelation'                   ), &
         &                                                         var_str('Black hole mass-VelocityDispersion mass relation'      ), &
         &                                                         var_str('velocityDispersion'                                    ), &
         &                                                         var_str('Velocity dispersion of spheroid'                       ), &
         &                                                         var_str('km/s'                                                  ), &
         &                                                         kilo                                                             , &
         &                                                         var_str('massBlackHole'                                         ), &
         &                                                         var_str('Mean logarithmic (base-10) mass of central black hole' ), &
         &                                                         var_str('M☉'                                                    ), &
         &                                                         massSolar                                                        , &
         &                                                         velocities                                                       , &
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
         &                                                         var_str('$\sigma_{\star,\mathrm{spheroid}}$ [km/s]'             ), &
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
    <objectDestructor name="nodePropertyExtractorVelocityDispersion_"        />
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="cosmologyParametersData"                         />
    <objectDestructor name="cosmologyFunctionsData"                          />
    !!]
    nullify(propertyOperators_      )
    nullify(weightPropertyOperators_)
    return
  end function blackHoleVelocityDispersionRelationConstructorInternal

  subroutine blackHoleVelocityDispersionRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisBlackHoleVelocityDispersionRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisBlackHoleVelocityDispersionRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine blackHoleVelocityDispersionRelationDestructor

