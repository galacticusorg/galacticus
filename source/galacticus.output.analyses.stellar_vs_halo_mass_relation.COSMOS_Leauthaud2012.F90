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
  
  !% Contains a module which implements a stellar vs halo mass relation analysis class.
  
  !# <outputAnalysis name="outputAnalysisStellarVsHaloMassRelationLeauthaud2012" defaultThreadPrivate="yes">
  !#  <description>A stellar vs halo mass relation output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisStellarVsHaloMassRelationLeauthaud2012
     !% A stellar vs halo mass relation output analysis class.
     private
  end type outputAnalysisStellarVsHaloMassRelationLeauthaud2012

  interface outputAnalysisStellarVsHaloMassRelationLeauthaud2012
     !% Constructors for the ``stellarVsHaloMassRelationLeauthaud2012'' output analysis class.
     module procedure stellarVsHaloMassRelationLeauthaud2012ConstructorParameters
     module procedure stellarVsHaloMassRelationLeauthaud2012ConstructorInternal
  end interface outputAnalysisStellarVsHaloMassRelationLeauthaud2012

contains

  function stellarVsHaloMassRelationLeauthaud2012ConstructorParameters(parameters) result (self)
    !% Constructor for the ``stellarVsHaloMassRelationLeauthaud2012'' output analysis class which takes a parameter set as input.
    use Cosmology_Functions
    use Input_Parameters
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelationLeauthaud2012)                              :: self
    type            (inputParameters                                     ), intent(inout)               :: parameters
    double precision                                                      , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                             ), pointer                     :: cosmologyFunctions_
    integer                                                                                             :: redshiftInterval
    
    ! Check and read parameters.
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    !# <inputParameter>
    !#   <name>redshiftInterval</name>
    !#   <source>parameters</source>
    !#   <variable>redshiftInterval</variable>
    !#   <description>The redshift interval (1, 2, or 3) to use.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial for stellar vs halo mass relation.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisStellarVsHaloMassRelationLeauthaud2012(redshiftInterval,systematicErrorPolynomialCoefficient,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function stellarVsHaloMassRelationLeauthaud2012ConstructorParameters

  function stellarVsHaloMassRelationLeauthaud2012ConstructorInternal(redshiftInterval,systematicErrorPolynomialCoefficient,cosmologyFunctions_) result (self)
    !% Constructor for the ``stellarVsHaloMassRelationLeauthaud2012'' output analysis class for internal use.
    use ISO_Varying_String
    use Numerical_Constants_Astronomical
    use Galacticus_Output_Times
    use Numerical_Ranges
    use Numerical_Comparison
    use Output_Analysis_Property_Operators
    use Output_Analysis_Property_Extractions
    use Output_Analysis_Distribution_Operators
    use Output_Analysis_Weight_Operators
    use Output_Analysis_Utilities
    use Memory_Management
    use Cosmology_Parameters
    use Cosmology_Functions
    use String_Handling
    use Galacticus_Error
    use Virial_Density_Contrast
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelationLeauthaud2012        )                                :: self
    integer                                                                       , intent(in   )                 :: redshiftInterval
    double precision                                                              , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                                     ), intent(inout), target         :: cosmologyFunctions_
    integer         (c_size_t                                                    ), parameter                     :: massHaloCount                                         =26
    double precision                                                              , parameter                     :: massHaloMinimum                                       = 1.0d10, massHaloMaximum                  =1.0d15
    integer                                                                       , parameter                     :: covarianceBinomialBinsPerDecade                       =10
    double precision                                                              , parameter                     :: covarianceBinomialMassHaloMinimum                     = 1.0d08, covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                              , allocatable  , dimension(:  ) :: massHalo
    double precision                                                              , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterStellarMass                                   ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterHaloIsolated                                  ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterAll                                           ), pointer                       :: galacticFilterAll_
    type            (filterList                                                  ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorIdentity                  ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity                        ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorLog10                         ), pointer                       :: outputAnalysisPropertyOperator_                              , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10                     ), pointer                       :: outputAnalysisPropertyUnoperator_                            , outputAnalysisWeightPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorSequence                      ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc             ), pointer                       :: outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial             ), pointer                       :: outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorFilterHighPass                ), pointer                       :: outputAnalysisWeightPropertyOperatorFilterHighPass_
    type            (outputAnalysisPropertyExtractorMassHalo                     ), pointer                       :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyExtractorMassStellar                  ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (propertyOperatorList                                        ), pointer                       :: propertyOperators_
    type            (cosmologyParametersSimple                                   ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                              ), pointer                       :: cosmologyFunctionsData
    type            (virialDensityContrastFixed                                  ), pointer                       :: virialDensityContrast_
    double precision                                                              , parameter                     :: errorPolynomialZeroPoint                              =11.3d00
    integer         (c_size_t                                                    )                                :: iBin
    double precision                                                                                              :: massStellarLimit                                              , redshiftMinimum                         , &
         &                                                                                                           redshiftMaximum
    type            (surveyGeometryFullSky                                       )                                :: surveyGeometry_
    
    ! Construct survey geometry.
    select case (redshiftInterval)
    case (1)
       redshiftMinimum = 0.22d0
       redshiftMaximum = 0.48d0
       massStellarLimit=10.00d0**8.7d0
    case (2)
       redshiftMinimum = 0.48d0
       redshiftMaximum = 0.74d0
       massStellarLimit=10.00d0**9.3d0
    case (3)
       redshiftMinimum = 0.74d0
       redshiftMaximum = 1.00d0
       massStellarLimit=10.00d0**9.8d0
    case default
       call Galacticus_Error_Report('stellarVsHaloMassRelationLeauthaud2012ConstructorInternal','redshiftInterval ∈ {1,2,3}')
    end select
    surveyGeometry_=surveyGeometryFullSky(redshiftMinimum,redshiftMaximum,cosmologyFunctions_)
    ! Create output time weights.
    call allocateArray(outputWeight,[massHaloCount,Galacticus_Output_Time_Count()])
    outputWeight(1,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,massStellarLimit,allowSingleEpoch=.true.)
    forall(iBin=2:massHaloCount)
       outputWeight(iBin,:)=outputWeight(1,:)
    end forall
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.25800d0, &
         &                                                 OmegaDarkEnergy= 0.74200d0, &
         &                                                 HubbleConstant =72.00000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.04385d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Create bins in halo mass.
    massHalo=Make_Range(log10(massHaloMinimum),log10(massHaloMaximum),int(massHaloCount),rangeType=rangeTypeLinear)
    ! Build a filter which select central galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilterStellarMass_      )
    allocate(galacticFilterHaloIsolated_     )
    allocate(galacticFilterAll_              )
    allocate(filters_                        )
    allocate(filters_                   %next)
    filters_                        %filter_ => galacticFilterHaloIsolated_
    filters_                   %next%filter_ => galacticFilterStellarMass_
    galacticFilterStellarMass_               =  galacticFilterStellarMass  (massThreshold=1.0d7   )
    galacticFilterHaloIsolated_              =  galacticFilterHaloIsolated (                      )
    galacticFilterAll_                       =  galacticFilterAll          (              filters_)
    ! Build identity distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    outputAnalysisDistributionOperator_                    =  outputAnalysisDistributionOperatorIdentity            (                                                                  )
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    outputAnalysisWeightOperator_                          =  outputAnalysisWeightOperatorIdentity                  (                                                                  )
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperator_                       )
    outputAnalysisPropertyOperator_                        =  outputAnalysisPropertyOperatorLog10                   (                                                                  )
    ! Build a sequence (log10, polynomial systematic, anti-log10, cosmological luminosity distance, high-pass filter) of weight property operators.
    allocate(outputAnalysisWeightPropertyOperatorFilterHighPass_)
    outputAnalysisWeightPropertyOperatorFilterHighPass_    =  outputAnalysisPropertyOperatorFilterHighPass          (massStellarLimit                                                  )
    allocate(outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_)
    outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_ =  outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc       (cosmologyFunctions_     ,cosmologyFunctionsData                   )
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    outputAnalysisWeightPropertyOperatorSystmtcPolynomial_ =  outputAnalysisPropertyOperatorSystmtcPolynomial       (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient     )
    allocate(outputAnalysisWeightPropertyOperatorLog10_            )
    outputAnalysisWeightPropertyOperatorLog10_             =  outputAnalysisPropertyOperatorLog10                   (                                                                  )
    allocate(outputAnalysisWeightPropertyOperatorAntiLog10_        )
    outputAnalysisWeightPropertyOperatorAntiLog10_         =  outputAnalysisPropertyOperatorAntiLog10               (                                                                  )
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    allocate(propertyOperators_%next%next%next                     )
    allocate(propertyOperators_%next%next%next%next                )
    propertyOperators_                    %operator_       => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperators_%next               %operator_       => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperators_%next%next          %operator_       => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperators_%next%next%next     %operator_       => outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next%next%next%next%operator_       => outputAnalysisWeightPropertyOperatorFilterHighPass_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    outputAnalysisWeightPropertyOperator_                  =  outputAnalysisPropertyOperatorSequence                (propertyOperators_                                               )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    outputAnalysisPropertyUnoperator_                      =  outputAnalysisPropertyOperatorAntiLog10               (                                                                 )
    ! Create a stellar mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    outputAnalysisWeightPropertyExtractor_                 =  outputAnalysisPropertyExtractorMassStellar            (                                                                 )
    ! Create a halo mass weight property extractor.
    allocate(virialDensityContrast_                                )
    virialDensityContrast_                                 =  virialDensityContrastFixed                            (200.0d0                 ,fixedDensityTypeMean,cosmologyFunctions_)
    allocate(outputAnalysisPropertyExtractor_                      )
    outputAnalysisPropertyExtractor_                       =  outputAnalysisPropertyExtractorMassHalo               (virialDensityContrast_                                           )
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                    &
         &                                                         var_str('stellarHaloMassRelationLeauthaud2012z')//redshiftInterval, &
         &                                                         var_str('Stellar vs. halo mass relation'       )                  , &
         &                                                         var_str('massHalo'                             )                  , &
         &                                                         var_str('Halo mass'                            )                  , &
         &                                                         var_str('M☉'                                  )                  , &
         &                                                         massSolar                                                         , &
         &                                                         var_str('massStellar'                          )                  , &
         &                                                         var_str('Stellar mass'                         )                  , &
         &                                                         var_str('M☉'                                  )                  , &
         &                                                         massSolar                                                         , &
         &                                                         massHalo                                                          , &
         &                                                         0_c_size_t                                                        , &
         &                                                         outputWeight                                                      , &
         &                                                         outputAnalysisPropertyExtractor_                                  , &
         &                                                         outputAnalysisWeightPropertyExtractor_                            , &
         &                                                         outputAnalysisPropertyOperator_                                   , &
         &                                                         outputAnalysisWeightPropertyOperator_                             , &
         &                                                         outputAnalysisPropertyUnoperator_                                 , &
         &                                                         outputAnalysisWeightOperator_                                     , &
         &                                                         outputAnalysisDistributionOperator_                               , &
         &                                                         galacticFilterAll_                                                , &
         &                                                         outputAnalysisCovarianceModelBinomial                             , &
         &                                                         covarianceBinomialBinsPerDecade                                   , &
         &                                                         covarianceBinomialMassHaloMinimum                                 , &
         &                                                         covarianceBinomialMassHaloMaximum                                   &
         &                                                        )
    ! Clean up.
    nullify(galacticFilterAll_                    )
    nullify(galacticFilterStellarMass_            )
    nullify(galacticFilterHaloIsolated_           )
    nullify(filters_                              )
    nullify(outputAnalysisDistributionOperator_   )
    nullify(outputAnalysisWeightOperator_         )
    nullify(outputAnalysisPropertyOperator_       )
    nullify(outputAnalysisPropertyUnoperator_     )
    nullify(outputAnalysisWeightPropertyOperator_ )
    nullify(outputAnalysisWeightPropertyExtractor_)
    nullify(outputAnalysisPropertyExtractor_      )
    nullify(cosmologyParametersData               )
    nullify(cosmologyFunctionsData                )
    nullify(virialDensityContrast_                )
    return
  end function stellarVsHaloMassRelationLeauthaud2012ConstructorInternal

