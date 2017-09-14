!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  
  !# <outputAnalysis name="outputAnalysisMassFunctionStellarSDSS" defaultThreadPrivate="yes">
  !#  <description>An SDSS stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMassFunctionStellar) :: outputAnalysisMassFunctionStellarSDSS
     !% An SDSS stellar mass function output analysis class.
     private
  end type outputAnalysisMassFunctionStellarSDSS

  interface outputAnalysisMassFunctionStellarSDSS
     !% Constructors for the ``massFunctionStellarSDSS'' output analysis class.
     module procedure massFunctionStellarSDSSConstructorParameters
     module procedure massFunctionStellarSDSSConstructorInternal
  end interface outputAnalysisMassFunctionStellarSDSS

contains

  function massFunctionStellarSDSSConstructorParameters(parameters) result (self)
    !% Constructor for the ``massFunctionStellarSDSS'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS)                              :: self
    type            (inputParameters                      ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass              ), pointer                     :: cosmologyFunctions_
    class           (gravitationalLensingClass            ), pointer                     :: gravitationalLensing_
    double precision                                       , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                              :: covarianceBinomialBinsPerDecade
    double precision                                                                     :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                  randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                  sizeSourceLensing

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
    !#   <defaultValue>0.07d0</defaultValue>
    !#   <description>The minimum random error for SDSS stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMaximum</variable>
    !#   <defaultValue>0.07d0</defaultValue>
    !#   <description>The minimum random error for SDSS stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.07d0]</defaultValue>
    !#   <description>The coefficients of the random error polynomial for SDSS stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial for SDSS stellar masses.</description>
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
    !#   <description>The number of bins per decade of halo mass to use when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMinimum</variable>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMaximum</variable>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !# <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisMassFunctionStellarSDSS(cosmologyFunctions_,gravitationalLensing_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massFunctionStellarSDSSConstructorParameters

  function massFunctionStellarSDSSConstructorInternal(cosmologyFunctions_,gravitationalLensing_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !% Constructor for the ``massFunctionStellarSDSS'' output analysis class for internal use.
    use Input_Parameters
    use Galacticus_Input_Paths
    use Output_Analysis_Distribution_Operators
    use Cosmology_Parameters
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS              )                              :: self
    class           (cosmologyFunctionsClass                            ), intent(in   ), target       :: cosmologyFunctions_
    class           (gravitationalLensingClass                          ), intent(in   ), target       :: gravitationalLensing_
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                         , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient
    integer                                                              , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )               :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryLiWhite2009SDSS                      )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng      )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence         )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                           )               , pointer      :: distributionOperatorSequence
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                            =11.3d+0

    ! Build a filter which select galaxies with stellar mass 10⁶M☉ or greater.
    allocate(galacticFilter_)
    galacticFilter_=galacticFilterStellarMass(massThreshold=1.0d6)
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
    ! Build the SDSS survey geometry of Li & White (2009) with their imposed redshift limits.
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryLiWhite2009SDSS(redshiftMinimum=1.0d-3,redshiftMaximum=5000.0d0,cosmologyFunctions_=cosmologyFunctions_)
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
    allocate(distributionOperatorSequence                )
    allocate(distributionOperatorSequence           %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    outputAnalysisDistributionOperator_                  =  outputAnalysisDistributionOperatorSequence          (                                  &
         &                                                                                                       distributionOperatorSequence      &
         &                                                                                                      )
    ! Build the object.
    self%outputAnalysisMassFunctionStellar=                                                                                                                     &
         & outputAnalysisMassFunctionStellar(                                                                                                                   &
         &                                   var_str('LiWhite2009SDSS'                                              )                                         , &
         &                                   var_str('Stellar mass function for the Li & White (2009) SDSS analysis')                                         , &
         &                                   char(Galacticus_Input_Path()//'/data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.hdf5'), &
         &                                   galacticFilter_                                                                                                  , &
         &                                   surveyGeometry_                                                                                                  , &
         &                                   cosmologyFunctions_                                                                                              , &
         &                                   cosmologyFunctionsData                                                                                           , &
         &                                   outputAnalysisPropertyOperator_                                                                                  , &
         &                                   outputAnalysisDistributionOperator_                                                                              , &
         &                                   covarianceBinomialBinsPerDecade                                                                                  , &
         &                                   covarianceBinomialMassHaloMinimum                                                                                , &
         &                                   covarianceBinomialMassHaloMaximum                                                                                  &
         &                                  )
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
  end function massFunctionStellarSDSSConstructorInternal
