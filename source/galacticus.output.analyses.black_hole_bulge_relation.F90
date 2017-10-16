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
  
  !% Contains a module which implements a mass-metallicity relation analysis class.
  
  !# <outputAnalysis name="outputAnalysisBlackHoleBulgeRelation" defaultThreadPrivate="yes">
  !#  <description>A mass-metallicity relation output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisBlackHoleBulgeRelation
     !% A mass-metallicity relation output analysis class.
     private
  end type outputAnalysisBlackHoleBulgeRelation

  interface outputAnalysisBlackHoleBulgeRelation
     !% Constructors for the ``blackHoleBulgeRelation'' output analysis class.
     module procedure blackHoleBulgeRelationConstructorParameters
     module procedure blackHoleBulgeRelationConstructorInternal
  end interface outputAnalysisBlackHoleBulgeRelation

contains

  function blackHoleBulgeRelationConstructorParameters(parameters) result (self)
    !% Constructor for the ``blackHoleBulgeRelation'' output analysis class which takes a parameter set as input.
    use Cosmology_Functions
    use Input_Parameters
    implicit none
    type            (outputAnalysisBlackHoleBulgeRelation)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    double precision                                      , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass             ), pointer                     :: cosmologyFunctions_
    double precision                                                                    :: randomErrorMinimum                  , randomErrorMaximum 

    
    ! Check and read parameters.
    allocate(systematicErrorPolynomialCoefficient(max(1,parameters%count('systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorPolynomialCoefficient(max(1,parameters%count(    'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.09d0]</defaultValue>
    !#   <description>The coefficients of the random error polynomial.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMinimum</variable>
    !#   <defaultValue>0.09d0</defaultValue>
    !#   <description>The minimum random error for stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMaximum</variable>
    !#   <defaultValue>0.09d0</defaultValue>
    !#   <description>The minimum random error for stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisBlackHoleBulgeRelation(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function blackHoleBulgeRelationConstructorParameters

  function blackHoleBulgeRelationConstructorInternal(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_) result (self)
    !% Constructor for the ``blackHoleBulgeRelation'' output analysis class for internal use.
    use Memory_Management
    use IO_HDF5
    use Galacticus_Input_Paths  
    use Galacticus_Output_Times
    use Output_Analysis_Property_Operators
    use Output_Analysis_Property_Extractions
    use Output_Analysis_Distribution_Operators
    use Output_Analysis_Weight_Operators
    use Output_Analysis_Utilities
    use Cosmology_Parameters
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Numerical_Ranges
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    type            (outputAnalysisBlackHoleBulgeRelation               )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                    , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: masses
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterSpheroidStellarMass                  ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyExtractorMassStellarSpheroid ), pointer                       :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyExtractorMassBlackHole       ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.3d00
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    double precision                                                     , parameter                     :: massStellarMinimum                                       = 3.0d09, massStellarMaximum                  =3.0d11
    integer         (c_size_t                                           ), parameter                     :: massStellarCount                                         =5
    integer         (c_size_t                                           )                                :: iOutput

    ! Specify mass bins.
    masses=Make_Range(log10(massStellarMinimum),log10(massStellarMaximum),int(massStellarCount),rangeType=rangeTypeLinear)
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[massStellarCount,Galacticus_Output_Time_Count()])
    outputWeight=0.0d0
    do iOutput=1,Galacticus_Output_Time_Count()
       if (Values_Agree(Galacticus_Output_Redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Galacticus_Error_Report('zero redshift output is required'//{introspection:location})
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.30000d0, &
         &                                                 OmegaDarkEnergy= 0.70000d0, &
         &                                                 HubbleConstant =70.50000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.00000d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilter_                                       )
    galacticFilter_                                  =  galacticFilterSpheroidStellarMass                  (massThreshold=1.0d8                                          )
     ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    outputAnalysisWeightOperator_                    =  outputAnalysisWeightOperatorIdentity               (                                                             )
    ! Build luminosity distance, systematic, and log10() property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_ =  outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    (cosmologyFunctions_     ,cosmologyFunctionsData              )
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    outputAnalysisPropertyOperatorSystmtcPolynomial_ =  outputAnalysisPropertyOperatorSystmtcPolynomial    (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    outputAnalysisPropertyOperatorLog10_             =  outputAnalysisPropertyOperatorLog10                (                                                             )
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_               %operator_      => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next          %operator_      => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next     %operator_      => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    outputAnalysisPropertyOperator_                  =  outputAnalysisPropertyOperatorSequence             (propertyOperators_                                           )
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    outputAnalysisDistributionOperator_              =  outputAnalysisDistributionOperatorRandomErrorPlynml(                                  &
         &                                                                                                  randomErrorMinimum              , &
         &                                                                                                  randomErrorMaximum              , &
         &                                                                                                  errorPolynomialZeroPoint        , &
         &                                                                                                  randomErrorPolynomialCoefficient  &
         &                                                                                                 )
    ! Build a metallicity weight property operator.
    allocate(outputAnalysisWeightPropertyOperator_                 )
    outputAnalysisWeightPropertyOperator_            =  outputAnalysisPropertyOperatorLog10                (                                                             )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    outputAnalysisPropertyUnoperator_                =  outputAnalysisPropertyOperatorAntiLog10            (                                                             )
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisPropertyExtractor_                      )
    outputAnalysisPropertyExtractor_                 =  outputAnalysisPropertyExtractorMassStellarSpheroid (                                                             )
    ! Create an ISM metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    outputAnalysisWeightPropertyExtractor_           =  outputAnalysisPropertyExtractorMassBlackHole       (                                                             )
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                &
         &                                                         var_str('blackHoleBulgeRelation'             ), &
         &                                                         var_str('Black hole mass-bulge mass relation'), &
         &                                                         var_str('massStellarSpheroid'                ), &
         &                                                         var_str('Stellar mass of spheroid'           ), &
         &                                                         var_str('M☉'                                 ), &
         &                                                         massSolar                                     , &
         &                                                         var_str('massBlackHole'                      ), &
         &                                                         var_str('Mass of central black hole'         ), &
         &                                                         var_str('M☉'                                 ), &
         &                                                         massSolar                                     , &
         &                                                         masses                                        , &
         &                                                         bufferCount                                   , &
         &                                                         outputWeight                                  , &
         &                                                         outputAnalysisPropertyExtractor_              , &
         &                                                         outputAnalysisWeightPropertyExtractor_        , &
         &                                                         outputAnalysisPropertyOperator_               , &
         &                                                         outputAnalysisWeightPropertyOperator_         , &
         &                                                         outputAnalysisPropertyUnoperator_             , &
         &                                                         outputAnalysisWeightOperator_                 , &
         &                                                         outputAnalysisDistributionOperator_           , &
         &                                                         galacticFilter_                               , &
         &                                                         outputAnalysisCovarianceModelBinomial         , &
         &                                                         covarianceBinomialBinsPerDecade               , &
         &                                                         covarianceBinomialMassHaloMinimum             , &
         &                                                         covarianceBinomialMassHaloMaximum               &
         &                                                        )
    ! Clean up.
    nullify(galacticFilter_                                 )
    nullify(outputAnalysisDistributionOperator_             )
    nullify(outputAnalysisWeightOperator_                   )
    nullify(outputAnalysisPropertyOperator_                 )
    nullify(outputAnalysisPropertyOperatorLog10_            )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    nullify(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    nullify(outputAnalysisPropertyUnoperator_               )
    nullify(outputAnalysisWeightPropertyOperator_           )
    nullify(outputAnalysisWeightPropertyExtractor_          )
    nullify(outputAnalysisPropertyExtractor_                )
    nullify(cosmologyParametersData                         )
    nullify(cosmologyFunctionsData                          )
    nullify(propertyOperators_                              )
    return
  end function blackHoleBulgeRelationConstructorInternal

