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
Implements a luminosity function output analysis class for the \cite{stefanon_evolution_2013} analysis.
!!}

  !![
  <outputAnalysis name="outputAnalysisLuminosityFunctionStefanonMarchesini2013">
   <description>A  luminosity function output analysis class for the \cite{stefanon_evolution_2013} analysis.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisLuminosityFunction) :: outputAnalysisLuminosityFunctionStefanonMarchesini2013
     !!{
     A luminosity function output analysis class for the \cite{stefanon_evolution_2013} analysis.
     !!}
     private
     class           (gravitationalLensingClass), pointer                     :: gravitationalLensing_            => null()
     double precision                           , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                         :: randomErrorMinimum                        , randomErrorMaximum                  , &
          &                                                                      sizeSourceLensing
     character       (len=1                    )                              :: band
     integer                                                                  :: redshiftInterval
   contains
     final :: luminosityFunctionStefanonMarchesini2013Destructor
  end type outputAnalysisLuminosityFunctionStefanonMarchesini2013

  interface outputAnalysisLuminosityFunctionStefanonMarchesini2013
     !!{
     Constructors for the \refClass{outputAnalysisLuminosityFunctionStefanonMarchesini2013} output analysis class.
     !!}
     module procedure luminosityFunctionStefanonMarchesini2013ConstructorParameters
     module procedure luminosityFunctionStefanonMarchesini2013ConstructorInternal
  end interface outputAnalysisLuminosityFunctionStefanonMarchesini2013

contains

  function luminosityFunctionStefanonMarchesini2013ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunctionStefanonMarchesini2013} output analysis class which takes a parameter set as input.
    !!}
    use :: Gravitational_Lensing, only : gravitationalLensing, gravitationalLensingClass
    use :: Input_Parameters     , only : inputParameter      , inputParameters
    implicit none
    type            (outputAnalysisLuminosityFunctionStefanonMarchesini2013)                              :: self
    type            (inputParameters                                       ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                               ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                      ), pointer                     :: outputTimes_
    class           (gravitationalLensingClass                             ), pointer                     :: gravitationalLensing_
    double precision                                                        , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                                      :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                                   randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                                   sizeSourceLensing
    character       (len=1                                                 )                              :: band
    integer                                                                                               :: redshiftInterval

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                       ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    !![
    <inputParameter>
      <name>band</name>
      <source>parameters</source>
      <description>The band (J or H) for which the luminosity function should be computed.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftInterval</name>
      <source>parameters</source>
      <variable>redshiftInterval</variable>
      <description>The redshift interval (0-3) to use.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>No estimate of photometric uncertainty is provided by the reference paper---a tiny value is adopted.</defaultSource>
      <description>The minimum random error for absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>No estimate of photometric uncertainty is provided by the reference paper---a tiny value is adopted.</defaultSource>
      <description>The minimum random error for absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.01d0]</defaultValue>
      <defaultSource>No estimate of photometric uncertainty is provided by the reference paper---a tiny value is adopted.</defaultSource>
      <description>The coefficients of the random error polynomial for absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for absolute magnitudes.</description>
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
      <description>The number of bins per decade of halo mass to use when constructing luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing SDSS luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing SDSS luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisLuminosityFunctionStefanonMarchesini2013(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,band,redshiftInterval)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="outputTimes_"         />
    <objectDestructor name="gravitationalLensing_"/>
    !!]
    return
  end function luminosityFunctionStefanonMarchesini2013ConstructorParameters

  function luminosityFunctionStefanonMarchesini2013ConstructorInternal(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,band,redshiftInterval) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunctionStefanonMarchesini2013} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                        , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                      , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: Gravitational_Lensing                 , only : gravitationalLensingClass
    use :: Output_Analysis_Distribution_Operators, only : distributionOperatorList                       , outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorRandomErrorPlynml, outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    use :: String_Handling                       , only : String_Lower_Case                              , operator(//)
    implicit none
    type            (outputAnalysisLuminosityFunctionStefanonMarchesini2013)                              :: self
    class           (cosmologyFunctionsClass                               ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                                      ), intent(inout), target       :: outputTimes_
    class           (gravitationalLensingClass                             ), intent(in   ), target       :: gravitationalLensing_
    double precision                                                        , intent(in   )               :: randomErrorMinimum                                  , randomErrorMaximum                  , &
         &                                                                                                   sizeSourceLensing
    double precision                                                        , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient
    integer                                                                 , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                        , intent(in   )               :: covarianceBinomialMassHaloMinimum                   , covarianceBinomialMassHaloMaximum
    character       (len=1                                                 ), intent(in   )               :: band
    integer                                                                 , intent(in   )               :: redshiftInterval
    type            (galacticFilterStellarMass                             )               , pointer      :: galacticFilter_
    type            (surveyGeometryFullSky                                 )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial       )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml   )               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng         )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence            )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                             )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                        )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                              )               , pointer      :: distributionOperatorSequence
    double precision                                                                                      :: errorPolynomialZeroPoint
    type            (varying_string                                        )                              :: fileName
    double precision                                                                                      :: redshiftMinimum                                     , redshiftMaximum
    character       (len=7                                                 )                              :: redshiftLabel
    !![
    <constructorAssign variables="randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, sizeSourceLensing, redshiftInterval, *gravitationalLensing_"/>
    !!]

    ! Validate redshift interval.
    if (redshiftInterval < 0    .or.  redshiftInterval >  3  ) call Error_Report('0 ≤ redshiftInterval ≤ 3 is required'//{introspection:location})
    ! Validate band.
    if (band             /= "J" .and. band             /= "H") call Error_Report('band ∈ {J,H} is required'            //{introspection:location})
    ! Set zero point for random error polynomial to be M* of the Schechter function fit.
    select case (band)
    case ('J')
       select case (redshiftInterval)
       case (0)
          errorPolynomialZeroPoint=-23.72d0
       case (1)
          errorPolynomialZeroPoint=-23.60d0
       case (2)
          errorPolynomialZeroPoint=-23.42d0
       case (3)
          errorPolynomialZeroPoint=-23.28d0
       end select
    case ('H')
       select case (redshiftInterval)
       case (0)
          errorPolynomialZeroPoint=-24.03d0
       case (1)
          errorPolynomialZeroPoint=-23.94d0
       case (2)
          errorPolynomialZeroPoint=-23.74d0
       case (3)
          errorPolynomialZeroPoint=-23.39d0
       end select
    end select
    ! Set redshift interval.
    select case (redshiftInterval)
    case (0)
       redshiftMinimum=1.5d0
       redshiftMaximum=2.0d0
    case (1)
       redshiftMinimum=2.0d0
       redshiftMaximum=2.5d0
    case (2)
       redshiftMinimum=2.5d0
       redshiftMaximum=3.0d0
    case (3)
       redshiftMinimum=3.0d0
       redshiftMaximum=3.5d0
    end select
    ! Construct the file name.
    write (redshiftLabel,'(f3.1,"-",f3.1)') redshiftMinimum,redshiftMaximum
    fileName=String_Lower_Case(band)//"Band:z"//redshiftLabel//":LuminosityFunctionStefanonMarchesini2013.hdf5"
    ! Build a filter which select galaxies with stellar mass 10³M☉ or greater.
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterStellarMass(massThreshold=1.0d3)"/>
    !!]
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
        &amp;                        OmegaBaryon    = 0.04550d0  &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
      cosmologyFunctionsMatterLambda(                             &amp;
        &amp;                        cosmologyParametersData      &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a full-sky geometry.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorRandomErrorPlynml_">
     <constructor>
      outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
        &amp;                                              randomErrorMinimum              , &amp;
        &amp;                                              randomErrorMaximum              , &amp;
        &amp;                                              errorPolynomialZeroPoint        , &amp;
        &amp;                                              randomErrorPolynomialCoefficient  &amp;
        &amp;                                             )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorGrvtnlLnsng_">
     <constructor>
      outputAnalysisDistributionOperatorGrvtnlLnsng       (                                  &amp;
        &amp;                                              gravitationalLensing_           , &amp;
        &amp;                                              outputTimes_                    , &amp;
        &amp;                                              sizeSourceLensing                 &amp;
        &amp;                                             )
     </constructor>
    </referenceConstruct>
    !!]
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence            )
    allocate(distributionOperatorSequence       %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorSequence          (                                  &amp;
        &amp;                                              distributionOperatorSequence      &amp;
        &amp;                                             )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the object.
    self%outputAnalysisLuminosityFunction=                                                                                                         &
         & outputAnalysisLuminosityFunction(                                                                                                       &
         &                                               var_str('StefanonMarchesini2013')//band//'z'//redshiftInterval                          , &
         &                                               var_str(band//'-band luminosity function for the Stefanon & Marchesini (2013) analysis'), &
         &                                               char   (inputPath(pathTypeDataStatic)//'/observations/luminosityFunctions/'//fileName  ), &
         &                                               galacticFilter_                                                                         , &
         &                                               surveyGeometry_                                                                         , &
         &                                               cosmologyFunctions_                                                                     , &
         &                                               cosmologyFunctionsData                                                                  , &
         &                                               outputAnalysisPropertyOperator_                                                         , &
         &                                               outputAnalysisDistributionOperator_                                                     , &
         &                                               outputTimes_                                                                            , &
         &                                               covarianceBinomialBinsPerDecade                                                         , &
         &                                               covarianceBinomialMassHaloMinimum                                                       , &
         &                                               covarianceBinomialMassHaloMaximum                                                       , &
         &                                  filterName  ='2MASS_'//band                                                                          , &
         &                                  filterType  ='rest'                                                                                    &
         &                                 )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                     />
    <objectDestructor name="galacticFilter_"                                     />
    <objectDestructor name="cosmologyParametersData"                             />
    <objectDestructor name="cosmologyFunctionsData"                              />
    <objectDestructor name="outputAnalysisPropertyOperator_"                     />
    <objectDestructor name="outputAnalysisDistributionOperator_"                 />
    <objectDestructor name="outputAnalysisDistributionOperatorGrvtnlLnsng_"      />
    <objectDestructor name="outputAnalysisDistributionOperatorRandomErrorPlynml_"/>
    !!]
    nullify(distributionOperatorSequence)
    return
  end function luminosityFunctionStefanonMarchesini2013ConstructorInternal

  subroutine luminosityFunctionStefanonMarchesini2013Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLuminosityFunctionStefanonMarchesini2013} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLuminosityFunctionStefanonMarchesini2013), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine luminosityFunctionStefanonMarchesini2013Destructor
  
