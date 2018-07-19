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

!% Contains a module which implements an ALFALFA HI mass function output analysis class.

  use Gravitational_Lensing
  
  !# <outputAnalysis name="outputAnalysisMassFunctionHIALFALFAMartin2010" defaultThreadPrivate="yes">
  !#  <description>An ALFALFA HI mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMassFunctionHI) :: outputAnalysisMassFunctionHIALFALFAMartin2010
     !% An ALFALFA HI mass function output analysis class.
     private
  end type outputAnalysisMassFunctionHIALFALFAMartin2010

  interface outputAnalysisMassFunctionHIALFALFAMartin2010
     !% Constructors for the ``massFunctionHIALFALFAMartin2010'' output analysis class.
     module procedure massFunctionHIALFALFAMartin2010ConstructorParameters
     module procedure massFunctionHIALFALFAMartin2010ConstructorInternal
  end interface outputAnalysisMassFunctionHIALFALFAMartin2010

contains

  function massFunctionHIALFALFAMartin2010ConstructorParameters(parameters) result (self)
    !% Constructor for the ``massFunctionHIALFALFAMartin2010'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Output_Analysis_Molecular_Ratios
    implicit none
    type            (outputAnalysisMassFunctionHIALFALFAMartin2010)                              :: self
    type            (inputParameters                              ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                      ), pointer                     :: cosmologyFunctions_
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
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial for ALFALFA HI masses.</description>
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
    !#   <description>The number of bins per decade of halo mass to use when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMinimum</variable>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMaximum</variable>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing ALFALFA HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                            source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                           source="parameters"/>
    !# <objectBuilder class="gravitationalLensing"               name="gravitationalLensing_"                          source="parameters"/>
    !# <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperatorRandomError_" source="parameters"/>
    !# <objectBuilder class="outputAnalysisMolecularRatio"       name="outputAnalysisMolecularRatio_"                  source="parameters"/>
    ! Build the object.
    self=outputAnalysisMassFunctionHIALFALFAMartin2010(cosmologyFunctions_,cosmologyParameters_,outputAnalysisDistributionOperatorRandomError_,outputAnalysisMolecularRatio_,gravitationalLensing_,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massFunctionHIALFALFAMartin2010ConstructorParameters

  function massFunctionHIALFALFAMartin2010ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,outputAnalysisDistributionOperatorRandomError_,outputAnalysisMolecularRatio_,gravitationalLensing_,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !% Constructor for the ``massFunctionHIALFALFAMartin2010'' output analysis class for internal use.
    use Input_Parameters
    use Galacticus_Paths
    use Output_Analysis_Distribution_Operators
    use Cosmology_Parameters
    use Output_Analysis_Molecular_Ratios
    implicit none
    type            (outputAnalysisMassFunctionHIALFALFAMartin2010  )                              :: self
    class           (cosmologyFunctionsClass                        ), intent(in   ), target       :: cosmologyFunctions_
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

    ! Build a filter which select galaxies with ISM mass 10⁴M☉ or greater.
    allocate(galacticFilter_)
    galacticFilter_=galacticFilterISMMass(massThreshold=1.0d4)
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
    ! Build the ALFALFA survey geometry of Martin et al. (2010).
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryMartin2010ALFALFA(cosmologyParameters_)
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    outputAnalysisPropertyOperator_    =outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    outputAnalysisDistributionOperatorGrvtnlLnsng_        =  outputAnalysisDistributionOperatorGrvtnlLnsng       (                              &
         &                                                                                                        gravitationalLensing_       , &
         &                                                                                                        sizeSourceLensing             &
         &                                                                                                       )
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence                )
    allocate(distributionOperatorSequence           %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomError_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    outputAnalysisDistributionOperator_                  =  outputAnalysisDistributionOperatorSequence          (                               &
         &                                                                                                       distributionOperatorSequence   &
         &                                                                                                      )
    ! Build the object.
    self%outputAnalysisMassFunctionHI=                                                                                                          &
         & outputAnalysisMassFunctionHI(                                                                                                        &
         &                              var_str('Martin2010ALFALFA'                                             )                             , &
         &                              var_str('HI mass function for the Martin et al. (2010) ALFALFA analysis')                             , &
         &                              char(galacticusPath(pathTypeDataStatic)//'/observations/massFunctionsHI/HI_Mass_Function_ALFALFA_2010.hdf5'), &
         &                              galacticFilter_                                                                                       , &
         &                              surveyGeometry_                                                                                       , &
         &                              cosmologyFunctions_                                                                                   , &
         &                              cosmologyFunctionsData                                                                                , &
         &                              outputAnalysisPropertyOperator_                                                                       , &
         &                              outputAnalysisDistributionOperator_                                                                   , &
         &                              outputAnalysisMolecularRatio_                                                                         , &
         &                              covarianceBinomialBinsPerDecade                                                                       , &
         &                              covarianceBinomialMassHaloMinimum                                                                     , &
         &                              covarianceBinomialMassHaloMaximum                                                                       &
         &                             )
    ! Clean up.
    nullify(surveyGeometry_                               )
    nullify(galacticFilter_                               )
    nullify(cosmologyParametersData                       )
    nullify(cosmologyFunctionsData                        )
    nullify(outputAnalysisDistributionOperator_           )
    nullify(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    nullify(distributionOperatorSequence                  )
    return
  end function massFunctionHIALFALFAMartin2010ConstructorInternal
