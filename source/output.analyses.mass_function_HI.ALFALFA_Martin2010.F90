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
Implements an ALFALFA HI mass function output analysis class.
!!}


  !![
  <outputAnalysis name="outputAnalysisMassFunctionHIALFALFAMartin2010">
   <description>
    An ALFALFA HI $z\approx 0.0$ mass function output analysis class measured by \cite{martin_arecibo_2010}. HI mass estimates
    can be affected by HI self-absorption for highly inclined galaxies. \cite[][see also
    \protect\citealt{zwaan_hipass_2005}]{zwaan_h_1997} estimate that this effect would lead to a mean underestimation of HI
    masses by a factor $1.1$ for a randomly oriented galaxy sample. Therefore, a value of $-0.0414$ for the systematic parameter
    {\normalfont \ttfamily [systematicErrorPolynomialCoefficient]} is recommended.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMassFunctionHI) :: outputAnalysisMassFunctionHIALFALFAMartin2010
     !!{
     An ALFALFA HI mass function output analysis class.
     !!}
     private
     class           (cosmologyParametersClass               ), pointer                     :: cosmologyParameters_                           => null()
     class           (gravitationalLensingClass              ), pointer                     :: gravitationalLensing_                          => null()
     class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperatorRandomError_ => null()
     double precision                                         , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient
     double precision                                                                       :: sizeSourceLensing
   contains
     final :: massFunctionHIALFALFAMartin2010Destructor
  end type outputAnalysisMassFunctionHIALFALFAMartin2010

  interface outputAnalysisMassFunctionHIALFALFAMartin2010
     !!{
     Constructors for the \refClass{outputAnalysisMassFunctionHIALFALFAMartin2010} output analysis class.
     !!}
     module procedure massFunctionHIALFALFAMartin2010ConstructorParameters
     module procedure massFunctionHIALFALFAMartin2010ConstructorInternal
  end interface outputAnalysisMassFunctionHIALFALFAMartin2010

contains

  function massFunctionHIALFALFAMartin2010ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassFunctionHIALFALFAMartin2010} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters            , only : cosmologyParameters         , cosmologyParametersClass
    use :: Input_Parameters                , only : inputParameter              , inputParameters
    use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatio, outputAnalysisMolecularRatioClass
    implicit none
    type            (outputAnalysisMassFunctionHIALFALFAMartin2010)                              :: self
    type            (inputParameters                              ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                      ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                             ), pointer                     :: outputTimes_
    class           (cosmologyParametersClass                     ), pointer                     :: cosmologyParameters_
    class           (gravitationalLensingClass                    ), pointer                     :: gravitationalLensing_
    class           (outputAnalysisMolecularRatioClass            ), pointer                     :: outputAnalysisMolecularRatio_
    class           (outputAnalysisDistributionOperatorClass      ), pointer                     :: outputAnalysisDistributionOperatorRandomError_
    double precision                                               , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient
    integer                                                                                      :: covarianceBinomialBinsPerDecade
    double precision                                                                             :: covarianceBinomialMassHaloMinimum             , covarianceBinomialMassHaloMaximum, &
         &                                                                                          sizeSourceLensing

    ! Check and read parameters.
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                   ))
    end if
    !![
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for ALFALFA HI masses.</description>
    </inputParameter>
    <inputParameter>
      <name>sizeSourceLensing</name>
      <source>parameters</source>
      <variable>sizeSourceLensing</variable>
      <defaultValue>2.0d-3</defaultValue>
      <description>The characteristic source size for gravitational lensing calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                            source="parameters"/>
    <objectBuilder class="outputTimes"                        name="outputTimes_"                                   source="parameters"/>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                           source="parameters"/>
    <objectBuilder class="gravitationalLensing"               name="gravitationalLensing_"                          source="parameters"/>
    <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperatorRandomError_" source="parameters"/>
    <objectBuilder class="outputAnalysisMolecularRatio"       name="outputAnalysisMolecularRatio_"                  source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisMassFunctionHIALFALFAMartin2010(cosmologyFunctions_,cosmologyParameters_,outputAnalysisDistributionOperatorRandomError_,outputAnalysisMolecularRatio_,gravitationalLensing_,outputTimes_,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                           />
    <objectDestructor name="outputTimes_"                                  />
    <objectDestructor name="cosmologyParameters_"                          />
    <objectDestructor name="gravitationalLensing_"                         />
    <objectDestructor name="outputAnalysisDistributionOperatorRandomError_"/>
    <objectDestructor name="outputAnalysisMolecularRatio_"                 />
    !!]
    return
  end function massFunctionHIALFALFAMartin2010ConstructorParameters

  function massFunctionHIALFALFAMartin2010ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,outputAnalysisDistributionOperatorRandomError_,outputAnalysisMolecularRatio_,gravitationalLensing_,outputTimes_,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassFunctionHIALFALFAMartin2010} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                        , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass                       , cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterISMMass
    use :: Input_Paths                           , only : inputPath                                      , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryMartin2010ALFALFA
    use :: Gravitational_Lensing                 , only : gravitationalLensingClass
    use :: Output_Analysis_Distribution_Operators, only : distributionOperatorList                       , outputAnalysisDistributionOperatorClass, outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Molecular_Ratios      , only : outputAnalysisMolecularRatioClass
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    implicit none
    type            (outputAnalysisMassFunctionHIALFALFAMartin2010  )                              :: self
    class           (cosmologyFunctionsClass                        ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                               ), intent(inout), target       :: outputTimes_
    class           (cosmologyParametersClass                       ), intent(in   ), target       :: cosmologyParameters_
    class           (gravitationalLensingClass                      ), intent(in   ), target       :: gravitationalLensing_
    class           (outputAnalysisMolecularRatioClass              ), intent(in   ), target       :: outputAnalysisMolecularRatio_
    class           (outputAnalysisDistributionOperatorClass        ), intent(in   ), target       :: outputAnalysisDistributionOperatorRandomError_
    double precision                                                 , intent(in   )               :: sizeSourceLensing
    double precision                                                 , intent(in   ), dimension(:) :: systematicErrorPolynomialCoefficient
    integer                                                          , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                 , intent(in   )               :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum
    type            (galacticFilterISMMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryMartin2010ALFALFA                )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial)               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng  )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence     )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                      )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                 )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                       )               , pointer      :: distributionOperatorSequence
    double precision                                                 , parameter                   :: errorPolynomialZeroPoint                            =11.3d+0
    !![
    <constructorAssign variables="sizeSourceLensing, systematicErrorPolynomialCoefficient, *cosmologyParameters_, *gravitationalLensing_, *outputAnalysisDistributionOperatorRandomError_"/>
    !!]
    
    ! Build a filter which select galaxies with ISM mass 10⁴M☉ or greater.
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterISMMass(massThreshold=1.0d4)"/>
    !!]
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
    <constructor>
    cosmologyParametersSimple     (                            &amp;
        &amp;                      OmegaMatter    = 0.30000d0, &amp;
        &amp;                      OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                      HubbleConstant =70.00000d0, &amp;
        &amp;                      temperatureCMB = 2.72548d0, &amp;
        &amp;                      OmegaBaryon    = 0.04550d0  &amp;
        &amp;                     )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
    <constructor>
    cosmologyFunctionsMatterLambda(                            &amp;
         &amp;                     cosmologyParametersData     &amp;
         &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the ALFALFA survey geometry of Martin et al. (2010).
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryMartin2010ALFALFA(cosmologyParameters_)"/>
    !!]
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorGrvtnlLnsng_">
    <constructor>
      outputAnalysisDistributionOperatorGrvtnlLnsng(                                    &amp;
         &amp;                                      gravitationalLensing_             , &amp;
         &amp;                                      outputTimes_                      , &amp;
         &amp;                                      sizeSourceLensing                   &amp;
         &amp;                                     )
     </constructor>
    </referenceConstruct>
    !!]
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence                )
    allocate(distributionOperatorSequence           %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomError_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
    <constructor>
      outputAnalysisDistributionOperatorSequence   (                                     &amp;
         &amp;                                      distributionOperatorSequence         &amp;
         &amp;                                     )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the object.
    self%outputAnalysisMassFunctionHI=                                                                                                           &
         & outputAnalysisMassFunctionHI(                                                                                                         &
         &                              var_str('Martin2010ALFALFA'                                             )                              , &
         &                              var_str('HI mass function for the Martin et al. (2010) ALFALFA analysis')                              , &
         &                              char(inputPath(pathTypeDataStatic)//'/observations/massFunctionsHI/HI_Mass_Function_ALFALFA_2010.hdf5'), &
         &                              galacticFilter_                                                                                        , &
         &                              surveyGeometry_                                                                                        , &
         &                              cosmologyFunctions_                                                                                    , &
         &                              cosmologyFunctionsData                                                                                 , &
         &                              outputAnalysisPropertyOperator_                                                                        , &
         &                              outputAnalysisDistributionOperator_                                                                    , &
         &                              outputAnalysisMolecularRatio_                                                                          , &
         &                              outputTimes_                                                                                           , &
         &                              covarianceBinomialBinsPerDecade                                                                        , &
         &                              covarianceBinomialMassHaloMinimum                                                                      , &
         &                              covarianceBinomialMassHaloMaximum                                                                        &
         &                             )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                               />
    <objectDestructor name="outputAnalysisPropertyOperator_"               />
    <objectDestructor name="galacticFilter_"                               />
    <objectDestructor name="cosmologyParametersData"                       />
    <objectDestructor name="cosmologyFunctionsData"                        />
    <objectDestructor name="outputAnalysisDistributionOperator_"           />
    <objectDestructor name="outputAnalysisDistributionOperatorGrvtnlLnsng_"/>
    !!]
    nullify(distributionOperatorSequence)
    return
  end function massFunctionHIALFALFAMartin2010ConstructorInternal

  subroutine massFunctionHIALFALFAMartin2010Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMassFunctionHIALFALFAMartin2010} output analysis class.
    !!}
    implicit none
    type(outputAnalysisMassFunctionHIALFALFAMartin2010), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"                          />
    <objectDestructor name="self%gravitationalLensing_"                         />
    <objectDestructor name="self%outputAnalysisDistributionOperatorRandomError_"/>
    !!]
    return
  end subroutine massFunctionHIALFALFAMartin2010Destructor
  
