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
Implements a stellar mass function output analysis class.
!!}


  !![
  <outputAnalysis name="outputAnalysisLuminosityFunctionSobral2013HiZELS">
   <description>An SDSS H$\alpha$ luminosity function output analysis class for the \cite{sobral_large_2013} analysis.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisLuminosityFunctionHalpha) :: outputAnalysisLuminosityFunctionSobral2013HiZELS
     !!{
     An SDSS H$\alpha$ luminosity function output analysis class for the \cite{sobral_large_2013} analysis.
     !!}
     private
     class           (gravitationalLensingClass), pointer                     :: gravitationalLensing_            => null()
     double precision                           , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     integer                                                                  :: redshiftInterval
     double precision                                                         :: randomErrorMinimum                        , randomErrorMaximum                  , &
          &                                                                      sizeSourceLensing
   contains
     final :: luminosityFunctionSobral2013HiZELSDestructor
  end type outputAnalysisLuminosityFunctionSobral2013HiZELS

  interface outputAnalysisLuminosityFunctionSobral2013HiZELS
     !!{
     Constructors for the {\normalfont \ttfamily luminosityFunctionSobral2013HiZELS} output analysis class.
     !!}
     module procedure luminosityFunctionSobral2013HiZELSConstructorParameters
     module procedure luminosityFunctionSobral2013HiZELSConstructorInternal
  end interface outputAnalysisLuminosityFunctionSobral2013HiZELS

contains

  function luminosityFunctionSobral2013HiZELSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily luminosityFunctionSobral2013HiZELS} output analysis class which takes a parameter set as input.
    !!}
    use :: Gravitational_Lensing         , only : gravitationalLensing           , gravitationalLensingClass
    use :: Input_Parameters              , only : inputParameter                 , inputParameters
    use :: Star_Formation_Rates_Disks    , only : starFormationRateDisksClass
    use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass
    implicit none
    type            (outputAnalysisLuminosityFunctionSobral2013HiZELS)                              :: self
    type            (inputParameters                                 ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                         ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                ), pointer                     :: outputTimes_
    class           (gravitationalLensingClass                       ), pointer                     :: gravitationalLensing_
    class           (starFormationRateDisksClass                     ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                 ), pointer                     :: starFormationRateSpheroids_
    class           (stellarSpectraDustAttenuationClass              ), pointer                     :: stellarSpectraDustAttenuation_
    double precision                                                  , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                         :: covarianceBinomialBinsPerDecade  , redshiftInterval
    double precision                                                                                :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                             randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                             sizeSourceLensing                , depthOpticalISMCoefficient

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
      <name>redshiftInterval</name>
      <source>parameters</source>
      <variable>redshiftInterval</variable>
      <description>The redshift interval (1, 2, 3, or 4) to use.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for SDSS H$\alpha$ luminosities.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for SDSS H$\alpha$ luminosities.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.1d0]</defaultValue>
      <description>The coefficients of the random error polynomial for SDSS H$\alpha$ luminosities.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for SDSS H$\alpha$ luminosities.</description>
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
      <description>The number of bins per decade of halo mass to use when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>depthOpticalISMCoefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>Multiplicative coefficient for optical depth in the ISM.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"            name="cosmologyFunctions_"            source="parameters"/>
    <objectBuilder class="outputTimes"                   name="outputTimes_"                   source="parameters"/>
    <objectBuilder class="gravitationalLensing"          name="gravitationalLensing_"          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"        name="starFormationRateDisks_"        source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"    name="starFormationRateSpheroids_"    source="parameters"/>
    <objectBuilder class="stellarSpectraDustAttenuation" name="stellarSpectraDustAttenuation_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisLuminosityFunctionSobral2013HiZELS(cosmologyFunctions_,gravitationalLensing_,stellarSpectraDustAttenuation_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,redshiftInterval,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,depthOpticalISMCoefficient)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"           />
    <objectDestructor name="outputTimes_"                  />
    <objectDestructor name="gravitationalLensing_"         />
    <objectDestructor name="starFormationRateDisks_"       />
    <objectDestructor name="starFormationRateSpheroids_"   />
    <objectDestructor name="stellarSpectraDustAttenuation_"/>
    !!]
    return
  end function luminosityFunctionSobral2013HiZELSConstructorParameters

  function luminosityFunctionSobral2013HiZELSConstructorInternal(cosmologyFunctions_,gravitationalLensing_,stellarSpectraDustAttenuation_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,redshiftInterval,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,depthOpticalISMCoefficient) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily luminosityFunctionSobral2013HiZELS} output analysis class for internal use.
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
    use :: Star_Formation_Rates_Disks            , only : starFormationRateDisksClass
    use :: Star_Formation_Rates_Spheroids        , only : starFormationRateSpheroidsClass
    use :: String_Handling                       , only : operator(//)
    implicit none
    type            (outputAnalysisLuminosityFunctionSobral2013HiZELS   )                              :: self
    class           (cosmologyFunctionsClass                            ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target       :: outputTimes_
    class           (gravitationalLensingClass                          ), intent(in   ), target       :: gravitationalLensing_
    class           (stellarSpectraDustAttenuationClass                 ), intent(in   ), target       :: stellarSpectraDustAttenuation_
    class           (starFormationRateDisksClass                        ), intent(in   ), target       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target       :: starFormationRateSpheroids_
    integer                                                              , intent(in   )               :: redshiftInterval
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                  , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing                                   , depthOpticalISMCoefficient
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient
    integer                                                              , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )               :: covarianceBinomialMassHaloMinimum                   , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryFullSky                              )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng      )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence         )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                           )               , pointer      :: distributionOperatorSequence
    double precision                                                                    , parameter    :: errorPolynomialZeroPoint                            =40.0d0
    type            (varying_string                                     )                              :: fileName
    !![
    <constructorAssign variables="randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, redshiftInterval, randomErrorMinimum, randomErrorMaximum, sizeSourceLensing, *gravitationalLensing_"/>
    !!]
    
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
      cosmologyFunctionsMatterLambda(                            &amp;
        &amp;                        cosmologyParametersData     &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the survey geometry. Since these are narrow band surveys with narrow redshift ranges we simply use a full sky geometry with matched redshift intervals.
    allocate(surveyGeometry_)
    select case (redshiftInterval)
    case (1)
       !![
       <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=0.401d0-0.010d0,redshiftMaximum=0.401d0+0.010d0,cosmologyFunctions_=cosmologyFunctions_)"/>
       !!]
       fileName       ='hAlphaLuminosityFunctionSobral2013HiZELSZ0.4.hdf5'
    case (2)
       !![
       <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=0.845d0-0.015d0,redshiftMaximum=0.845d0+0.015d0,cosmologyFunctions_=cosmologyFunctions_)"/>
       !!]
       fileName       ='hAlphaLuminosityFunctionSobral2013HiZELSZ0.84.hdf5'
    case (3)
       !![
       <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=1.466d0-0.016d0,redshiftMaximum=1.466d0+0.016d0,cosmologyFunctions_=cosmologyFunctions_)"/>
       !!]
       fileName       ='hAlphaLuminosityFunctionSobral2013HiZELSZ1.47.hdf5'
    case (4)
       !![
       <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=2.237d0-0.023d0,redshiftMaximum=2.237d0+0.023d0,cosmologyFunctions_=cosmologyFunctions_)"/>
       !!]
       fileName       ='hAlphaLuminosityFunctionSobral2013HiZELSZ2.23.hdf5'
    case default
       call Error_Report('redshift interval must be 1, 2, 3, or 4'//{introspection:location})
    end select
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"    constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
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
    self%outputAnalysisLuminosityFunctionHalpha=                                                                                                                                &
         & outputAnalysisLuminosityFunctionHalpha(                                                                                                                              &
         &                                        var_str('Sobral2013HiZELSZ'                                                                             )//redshiftInterval , &
         &                                        var_str('H$\alpha$ luminosity function for the Sobral et al. (2013) HiZELS analysis; redshift interval ')//redshiftInterval , &
         &                                        char   (inputPath(pathTypeDataStatic)//'/observations/luminosityFunctions/'                              //fileName        ), &
         &                                        .false.                                                                                                                     , &
         &                                        depthOpticalISMCoefficient                                                                                                  , &
         &                                        galacticFilter_                                                                                                             , &
         &                                        surveyGeometry_                                                                                                             , &
         &                                        stellarSpectraDustAttenuation_                                                                                              , &
         &                                        cosmologyFunctions_                                                                                                         , &
         &                                        cosmologyFunctionsData                                                                                                      , &
         &                                        outputAnalysisPropertyOperator_                                                                                             , &
         &                                        outputAnalysisDistributionOperator_                                                                                         , &
         &                                        outputTimes_                                                                                                                , &
         &                                        starFormationRateDisks_                                                                                                     , &
         &                                        starFormationRateSpheroids_                                                                                                 , &
         &                                        covarianceBinomialBinsPerDecade                                                                                             , &
         &                                        covarianceBinomialMassHaloMinimum                                                                                           , &
         &                                        covarianceBinomialMassHaloMaximum                                                                                             &
         &                                       )
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
  end function luminosityFunctionSobral2013HiZELSConstructorInternal

  subroutine luminosityFunctionSobral2013HiZELSDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily luminosityFunctionSobral2013HiZELS} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLuminosityFunctionSobral2013HiZELS), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine luminosityFunctionSobral2013HiZELSDestructor
  
