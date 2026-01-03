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
Implements a stellar mass function output analysis class.
!!}

  !![
  <outputAnalysis name="outputAnalysisLuminosityFunctionMonteroDorta2009SDSS">
   <description>An SDSS luminosity function output analysis class for the \cite{montero-dorta_sdss_2009} analysis.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisLuminosityFunction) :: outputAnalysisLuminosityFunctionMonteroDorta2009SDSS
     !!{
     An SDSS luminosity function output analysis class for the \cite{montero-dorta_sdss_2009} analysis.
     !!}
     private
     class           (gravitationalLensingClass), pointer                     :: gravitationalLensing_            => null()
     double precision                           , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                         :: randomErrorMinimum                        , randomErrorMaximum                  , &
          &                                                                      sizeSourceLensing
     character       (len=1                    )                              :: band
   contains
     final :: luminosityFunctionMonteroDorta2009SDSSDestructor
  end type outputAnalysisLuminosityFunctionMonteroDorta2009SDSS

  interface outputAnalysisLuminosityFunctionMonteroDorta2009SDSS
     !!{
     Constructors for the \refClass{outputAnalysisLuminosityFunctionMonteroDorta2009SDSS} output analysis class.
     !!}
     module procedure luminosityFunctionMonteroDorta2009SDSSConstructorParameters
     module procedure luminosityFunctionMonteroDorta2009SDSSConstructorInternal
  end interface outputAnalysisLuminosityFunctionMonteroDorta2009SDSS

contains

  function luminosityFunctionMonteroDorta2009SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunctionMonteroDorta2009SDSS} output analysis class which takes a parameter set as input.
    !!}
    use :: Gravitational_Lensing, only : gravitationalLensing, gravitationalLensingClass
    use :: Input_Parameters     , only : inputParameter      , inputParameters
    implicit none
    type            (outputAnalysisLuminosityFunctionMonteroDorta2009SDSS)                              :: self
    type            (inputParameters                                     ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                             ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                    ), pointer                     :: outputTimes_
    class           (gravitationalLensingClass                           ), pointer                     :: gravitationalLensing_
    double precision                                                      , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                             :: covarianceBinomialBinsPerDecade
    double precision                                                                                    :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                                 randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                                 sizeSourceLensing
    character       (len=1                                               )                              :: band

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                   ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                   ))
    end if
    !![
    <inputParameter>
      <name>band</name>
      <source>parameters</source>
      <description>The band (u, g, r, i, or z) for which the luminosity function should be computed.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>Computed from the quoted 95\% (S/N$\approx$1.64) depth of $r=22.2$ (\href{http://classic.sdss.org/dr7/}{http://classic.sdss.org/dr7/}), and assuming that most galaxies are at the limiting magnitude of $17.77$ for this sample using $\sigma_M=2.5 \log_{10}[1+1/\left\{\hbox{S/N}_\mathrm{lim} 10^{-0.4(m-m_\mathrm{lim})}\right\}]$.</defaultSource>
      <description>The minimum random error for SDSS absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>Computed from the quoted 95\% (S/N$\approx$1.64) depth of $r=22.2$ (\href{http://classic.sdss.org/dr7/}{http://classic.sdss.org/dr7/}), and assuming that most galaxies are at the limiting magnitude of $17.77$ for this sample using $\sigma_M=2.5 \log_{10}[1+1/\left\{\hbox{S/N}_\mathrm{lim} 10^{-0.4(m-m_\mathrm{lim})}\right\}]$.</defaultSource>
      <description>The minimum random error for SDSS absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.01d0]</defaultValue>
      <defaultSource>Computed from the quoted 95\% (S/N$\approx$1.64) depth of $r=22.2$ (\href{http://classic.sdss.org/dr7/}{http://classic.sdss.org/dr7/}), and assuming that most galaxies are at the limiting magnitude of $17.77$ for this sample using $\sigma_M=2.5 \log_{10}[1+1/\left\{\hbox{S/N}_\mathrm{lim} 10^{-0.4(m-m_\mathrm{lim})}\right\}]$.</defaultSource>
      <description>The coefficients of the random error polynomial for SDSS absolute magnitudes.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for SDSS absolute magnitudes.</description>
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
      <description>The number of bins per decade of halo mass to use when constructing SDSS luminosity function covariance matrices for main branch galaxies.</description>
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
    self=outputAnalysisLuminosityFunctionMonteroDorta2009SDSS(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,band)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="outputTimes_"         />
    <objectDestructor name="gravitationalLensing_"/>
    !!]
    return
  end function luminosityFunctionMonteroDorta2009SDSSConstructorParameters

  function luminosityFunctionMonteroDorta2009SDSSConstructorInternal(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,band) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunctionMonteroDorta2009SDSS} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                        , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                      , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryMonteroDorta2009SDSS
    use :: Gravitational_Lensing                 , only : gravitationalLensingClass
    use :: Output_Analysis_Distribution_Operators, only : distributionOperatorList                       , outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorRandomErrorPlynml, outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    implicit none
    type            (outputAnalysisLuminosityFunctionMonteroDorta2009SDSS)                              :: self
    class           (cosmologyFunctionsClass                             ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                                    ), intent(inout), target       :: outputTimes_
    class           (gravitationalLensingClass                           ), intent(in   ), target       :: gravitationalLensing_
    double precision                                                      , intent(in   )               :: randomErrorMinimum                                  , randomErrorMaximum                  , &
         &                                                                                                 sizeSourceLensing
    double precision                                                      , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient
    integer                                                               , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                      , intent(in   )               :: covarianceBinomialMassHaloMinimum                   , covarianceBinomialMassHaloMaximum
    character       (len=1                                               ), intent(in   )               :: band
    type            (galacticFilterStellarMass                           )               , pointer      :: galacticFilter_
    type            (surveyGeometryMonteroDorta2009SDSS                  )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial     )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml )               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng       )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence          )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                           )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                      )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                            )               , pointer      :: distributionOperatorSequence
    double precision                                                                                    :: errorPolynomialZeroPoint
    !![
    <constructorAssign variables="randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient,randomErrorMinimum, randomErrorMaximum, sizeSourceLensing, *gravitationalLensing_"/>
    !!]
    
    ! Validate band, and set zero point for random error polynomial to be M* of the Schechter function fit.
    select case (band)
    case ('u')
       errorPolynomialZeroPoint=-17.72d0
    case ('g')
       errorPolynomialZeroPoint=-19.53d0
    case ('r')
       errorPolynomialZeroPoint=-20.71d0
    case ('i')
       errorPolynomialZeroPoint=-20.93d0
    case ('z')
       errorPolynomialZeroPoint=-21.40d0
    case default
       call Error_Report('band ∈ {u,g,r,i,z} is required'//{introspection:location})
    end select
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
      cosmologyParametersSimple     (                             &amp;
        &amp;                        OmegaMatter    =  0.30000d0, &amp;
        &amp;                        OmegaDarkEnergy=  0.70000d0, &amp;
        &amp;                        HubbleConstant =100.00000d0, &amp;
        &amp;                        temperatureCMB =  2.72548d0, &amp;
        &amp;                        OmegaBaryon    =  0.04550d0  &amp;
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
    ! Build the SDSS survey geometry of Montero-Dorta & Prada (2009).
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryMonteroDorta2009SDSS(band,cosmologyFunctionsData)"/>
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
    self%outputAnalysisLuminosityFunction=                                                                                                                                               &
         & outputAnalysisLuminosityFunction(                                                                                                                                             &
         &                                               var_str('MonteroDorta2009SDSS'//band                                                        )                                 , &
         &                                               var_str(band//'-band luminosity function for the Montero-Dorta & Prada (2009) SDSS analysis')                                 , &
         &                                               char(inputPath(pathTypeDataStatic)//'/observations/luminosityFunctions/'//band//'LuminosityFunctionMonteroDorta2009SDSS.hdf5'), &
         &                                               galacticFilter_                                                                                                               , &
         &                                               surveyGeometry_                                                                                                               , &
         &                                               cosmologyFunctions_                                                                                                           , &
         &                                               cosmologyFunctionsData                                                                                                        , &
         &                                               outputAnalysisPropertyOperator_                                                                                               , &
         &                                               outputAnalysisDistributionOperator_                                                                                           , &
         &                                               outputTimes_                                                                                                                  , &
         &                                               covarianceBinomialBinsPerDecade                                                                                               , &
         &                                               covarianceBinomialMassHaloMinimum                                                                                             , &
         &                                               covarianceBinomialMassHaloMaximum                                                                                             , &
         &                                  filterName  ='SDSS_'//band                                                                                                                 , &
         &                                  filterType  ='observed'                                                                                                                    , &
         &                                  redshiftBand=0.1d0                                                                                                                           &
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
  end function luminosityFunctionMonteroDorta2009SDSSConstructorInternal

  subroutine luminosityFunctionMonteroDorta2009SDSSDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLuminosityFunctionMonteroDorta2009SDSS} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLuminosityFunctionMonteroDorta2009SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine luminosityFunctionMonteroDorta2009SDSSDestructor
  
