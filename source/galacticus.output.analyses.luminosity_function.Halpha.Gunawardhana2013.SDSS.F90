!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a stellar mass function output analysis class.

  use Gravitational_Lensing
  
  !# <outputAnalysis name="outputAnalysisLuminosityFunctionGunawardhana2013SDSS" defaultThreadPrivate="yes">
  !#  <description>An SDSS H$\alpha$ luminosity function output analysis class for the \cite{gunawardhana_galaxy_2013} analysis.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisLuminosityFunctionHalpha) :: outputAnalysisLuminosityFunctionGunawardhana2013SDSS
     !% An SDSS H$\alpha luminosity function output analysis class for the \cite{gunawardhana_galaxy_2013} analysis.
     private
  end type outputAnalysisLuminosityFunctionGunawardhana2013SDSS

  interface outputAnalysisLuminosityFunctionGunawardhana2013SDSS
     !% Constructors for the ``luminosityFunctionGunawardhana2013SDSS'' output analysis class.
     module procedure luminosityFunctionGunawardhana2013SDSSConstructorParameters
     module procedure luminosityFunctionGunawardhana2013SDSSConstructorInternal
  end interface outputAnalysisLuminosityFunctionGunawardhana2013SDSS

contains

  function luminosityFunctionGunawardhana2013SDSSConstructorParameters(parameters) result (self)
    !% Constructor for the ``luminosityFunctionGunawardhana2013SDSS'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisLuminosityFunctionGunawardhana2013SDSS)                              :: self
    type            (inputParameters                                     ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                             ), pointer                     :: cosmologyFunctions_
    class           (gravitationalLensingClass                           ), pointer                     :: gravitationalLensing_
    double precision                                                      , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                             :: covarianceBinomialBinsPerDecade
    double precision                                                                                    :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                                 randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                                 sizeSourceLensing                , depthOpticalISMCoefficient

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
    !# <inputParameter>
    !#   <name>randomErrorMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMinimum</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The minimum random error for SDSS H$\alpha$ luminosities.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMaximum</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The minimum random error for SDSS H$\alpha$ luminosities.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.1d0]</defaultValue>
    !#   <description>The coefficients of the random error polynomial for SDSS H$\alpha$ luminosities.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial for SDSS H$\alpha$ luminosities.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sizeSourceLensing</name>
    !#   <source>parameters</source>
    !#   <variable>sizeSourceLensing</variable>
    !#   <defaultValue>2.0d-3</defaultValue>
    !#   <description>The characteristic source size for gravitational lensing calculations.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialBinsPerDecade</variable>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMinimum</variable>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMaximum</variable>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing SDSS H$\alpha$ luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>depthOpticalISMCoefficient</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Multiplicative coefficient for optical depth in the ISM.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !# <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisLuminosityFunctionGunawardhana2013SDSS(cosmologyFunctions_,gravitationalLensing_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,depthOpticalISMCoefficient)
    !# <inputParametersValidate source="parameters"/>
    return
  end function luminosityFunctionGunawardhana2013SDSSConstructorParameters

  function luminosityFunctionGunawardhana2013SDSSConstructorInternal(cosmologyFunctions_,gravitationalLensing_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing,depthOpticalISMCoefficient) result (self)
    !% Constructor for the ``luminosityFunctionGunawardhana2013SDSS'' output analysis class for internal use.
    use Input_Parameters
    use Galacticus_Paths
    use Output_Analysis_Distribution_Operators
    use Cosmology_Parameters
    use Galacticus_Error
    implicit none
    type            (outputAnalysisLuminosityFunctionGunawardhana2013SDSS)                              :: self
    class           (cosmologyFunctionsClass                             ), intent(in   ), target       :: cosmologyFunctions_
    class           (gravitationalLensingClass                           ), intent(in   ), target       :: gravitationalLensing_
    double precision                                                      , intent(in   )               :: randomErrorMinimum                                  , randomErrorMaximum                  , &
         &                                                                                                 sizeSourceLensing                                   , depthOpticalISMCoefficient
    double precision                                                      , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient
    integer                                                               , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                      , intent(in   )               :: covarianceBinomialMassHaloMinimum                   , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                           )               , pointer      :: galacticFilter_
    type            (surveyGeometryGunawardhana2013SDSS                  )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial     )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml )               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng       )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence          )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                           )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                      )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                            )               , pointer      :: distributionOperatorSequence
    double precision                                                                                    :: errorPolynomialZeroPoint
    
    ! Build a filter which select galaxies with stellar mass 10³M☉ or greater.
    allocate(galacticFilter_)
    galacticFilter_=galacticFilterStellarMass(massThreshold=1.0d3)
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.30000d0, &
         &                                                 OmegaDarkEnergy= 0.70000d0, &
         &                                                 HubbleConstant =70.00000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.04550d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Build the SDSS survey geometry of Gunawardhana et al. (2013).
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryGunawardhana2013SDSS(cosmologyFunctionsData)
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    outputAnalysisPropertyOperator_    =outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    outputAnalysisDistributionOperatorRandomErrorPlynml_ =  outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &
         &                                                                                                       randomErrorMinimum              , &
         &                                                                                                       randomErrorMaximum              , &
         &                                                                                                       errorPolynomialZeroPoint        , &
         &                                                                                                       randomErrorPolynomialCoefficient  &
         &                                                                                                      )
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    outputAnalysisDistributionOperatorGrvtnlLnsng_       =  outputAnalysisDistributionOperatorGrvtnlLnsng       (                                  &
         &                                                                                                       gravitationalLensing_           , &
         &                                                                                                       sizeSourceLensing                 &
         &                                                                                                      )
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence            )
    allocate(distributionOperatorSequence       %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    outputAnalysisDistributionOperator_                  =  outputAnalysisDistributionOperatorSequence          (                                  &
         &                                                                                                       distributionOperatorSequence      &
         &                                                                                                      )
    ! Build the object.
    self%outputAnalysisLuminosityFunctionHalpha=                                                                                                                     &
         & outputAnalysisLuminosityFunctionHalpha(                                                                                                                   &
         &                                  var_str('Gunawardhana2013SDSS'                                                   )                                     , &
         &                                  var_str('Hα luminosity function for the Gunawardhana et al. (2013) SDSS analysis')                                     , &
         &                                  char(galacticusPath(pathTypeDataStatic)//'/observations/luminosityFunctions/hAlphaLuminosityFunctionGunawardhana13SDSS.hdf5'), &
         &                                  .false.                                                                                                                , &
         &                                  depthOpticalISMCoefficient                                                                                             , &
         &                                  galacticFilter_                                                                                                        , &
         &                                  surveyGeometry_                                                                                                        , &
         &                                  cosmologyFunctions_                                                                                                    , &
         &                                  cosmologyFunctionsData                                                                                                 , &
         &                                  outputAnalysisPropertyOperator_                                                                                        , &
         &                                  outputAnalysisDistributionOperator_                                                                                    , &
         &                                  covarianceBinomialBinsPerDecade                                                                                        , &
         &                                  covarianceBinomialMassHaloMinimum                                                                                      , &
         &                                  covarianceBinomialMassHaloMaximum                                                                                        &
         &                                 )
    ! Clean up.
    nullify(surveyGeometry_                                     )
    nullify(galacticFilter_                                     )
    nullify(cosmologyParametersData                             )
    nullify(cosmologyFunctionsData                              )
    nullify(outputAnalysisDistributionOperator_                 )
    nullify(outputAnalysisDistributionOperatorGrvtnlLnsng_      )
    nullify(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    nullify(distributionOperatorSequence                        )
    return
  end function luminosityFunctionGunawardhana2013SDSSConstructorInternal
