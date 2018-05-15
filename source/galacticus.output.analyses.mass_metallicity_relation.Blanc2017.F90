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

  !% Contains a module which implements a mass-metallicity relation analysis class.
  
  !# <outputAnalysis name="outputAnalysisMassMetallicityBlanc2017" defaultThreadPrivate="yes">
  !#  <description>A mass-metallicity relation output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisMassMetallicityBlanc2017
     !% A mass-metallicity relation output analysis class.
     private
  end type outputAnalysisMassMetallicityBlanc2017

  interface outputAnalysisMassMetallicityBlanc2017
     !% Constructors for the ``massMetallicityBlanc2017'' output analysis class.
     module procedure massMetallicityBlanc2017ConstructorParameters
     module procedure massMetallicityBlanc2017ConstructorInternal
  end interface outputAnalysisMassMetallicityBlanc2017

contains

  function massMetallicityBlanc2017ConstructorParameters(parameters) result (self)
    !% Constructor for the ``massMetallicityBlanc2017'' output analysis class which takes a parameter set as input.
    use Cosmology_Functions
    use Input_Parameters
    implicit none
    type            (outputAnalysisMassMetallicityBlanc2017)                              :: self
    type            (inputParameters                       ), intent(inout)               :: parameters
    double precision                                        , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient           , randomErrorPolynomialCoefficient, &
         &                                                                                   metallicitySystematicErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass               ), pointer                     :: cosmologyFunctions_
    double precision                                                                      :: randomErrorMinimum                             , randomErrorMaximum 

    
    ! Check and read parameters.
    allocate(metallicitySystematicErrorPolynomialCoefficient(max(1,parameters%count('metallicitySystematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(           systematicErrorPolynomialCoefficient(max(1,parameters%count(           'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(               randomErrorPolynomialCoefficient(max(1,parameters%count(               'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !# <inputParameter>
    !#   <name>metallicitySystematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>metallicitySystematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the metallicity systematic error polynomial.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
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
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the random error polynomial.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMinimum</variable>
    !#   <defaultValue>0.07d0</defaultValue>
    !#   <description>The minimum random error for stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMaximum</variable>
    !#   <defaultValue>0.07d0</defaultValue>
    !#   <description>The minimum random error for stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisMassMetallicityBlanc2017(metallicitySystematicErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massMetallicityBlanc2017ConstructorParameters

  function massMetallicityBlanc2017ConstructorInternal(metallicitySystematicErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_) result (self)
    !% Constructor for the ``massMetallicityBlanc2017'' output analysis class for internal use.
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
    use Abundances_Structure
    use Atomic_Data
    implicit none
    type            (outputAnalysisMassMetallicityBlanc2017             )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: metallicitySystematicErrorPolynomialCoefficient         , systematicErrorPolynomialCoefficient                          , &
         &                                                                                                  randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum                      =1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: masses
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterStellarMass                          ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterStarFormationRate                    ), pointer                       :: galacticFilterStarFormationRate_
    type            (galacticFilterAll                                  ), pointer                       :: galacticFilter_
    type            (filterList                                         ), pointer                       :: filters_                                       , filter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_                         , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorFilterHighPass       ), pointer                       :: outputAnalysisPropertyOperatorFilterHighPass_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorMetallicity12LogNH   ), pointer                       :: outputAnalysisPropertyOperatorMetallicity12LogNH_
    type            (outputAnalysisPropertyExtractorMassStellar         ), pointer                       :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyExtractorMetallicityISM      ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_        , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_                                      , weightPropertyOperators_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.3d00
    double precision                                                     , parameter                     :: metallicityErrorPolynomialZeroPoint             = 8.8d00
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    integer         (c_size_t                                           )                                :: iBin                                                    , binCount
    type            (surveyGeometryLiWhite2009SDSS                      )                                :: surveyGeometry_
    type            (hdf5Object                                         )                                :: dataFile

    ! Read masses at which fraction was measured.
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(Galacticus_Input_Path())//"data/observations/abundances/massMetallicityRelationBlanc2017.hdf5",readOnly=.true.)
    call dataFile%readDataset("massStellar"                                                                                      ,         masses)
    call dataFile%close      (                                                                                                                   )
    !$ call hdf5Access%unset()
    ! Construct survey geometry. Use a lower redshift limit than actually used by Blanc et al. to ensure that low mass bins have non-zero weight.
    surveyGeometry_=surveyGeometryLiWhite2009SDSS(redshiftMinimum=0.02d0,redshiftMaximum=0.25d0,cosmologyFunctions_=cosmologyFunctions_)
    ! Compute weights that apply to each output redshift.
    binCount=size(masses,kind=c_size_t)
    call allocateArray(outputWeight,[binCount,Galacticus_Output_Time_Count()])
    do iBin=1,binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,masses(iBin))
    end do
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.30000d0, &
         &                                                 OmegaDarkEnergy= 0.70000d0, &
         &                                                 HubbleConstant =70.00000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.00000d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample and with a
    ! selection for the star forming main sequence (based on that proposed by Renzini & Peng; 2015;
    ! http://adsabs.harvard.edu/abs/2015ApJ...801L..29R with a downward shift of 0.6dex to allow for the width of the
    ! distribution).
    allocate(galacticFilterStellarMass_                            )
    galacticFilterStellarMass_                       =  galacticFilterStellarMass                              (massThreshold=1.00d8                                         )
    allocate(galacticFilterStarFormationRate_                                       )
    galacticFilterStarFormationRate_                 =  galacticFilterStarFormationRate                        (                                                              &
         &                                                                                                      logM0        =0.00d0,                                         &
         &                                                                                                      logSFR0      =0.76d0,                                         &
         &                                                                                                      logSFR1      =0.76d0                                          &
         &                                                                                                     )
    allocate(filters_                                              )
    filter_ => filters_
    filter_%filter_ =>galacticFilterStellarMass_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterStarFormationRate_
    allocate(galacticFilter_                                       )
    galacticFilter_                                  =  galacticFilterAll                                      (filters_                                                     )
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    outputAnalysisWeightOperator_                    =  outputAnalysisWeightOperatorIdentity                   (                                                             )
    ! Build luminosity distance, systematic, and log10() property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_ =  outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc        (cosmologyFunctions_     ,cosmologyFunctionsData              )
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    outputAnalysisPropertyOperatorSystmtcPolynomial_ =  outputAnalysisPropertyOperatorSystmtcPolynomial        (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    outputAnalysisPropertyOperatorLog10_             =  outputAnalysisPropertyOperatorLog10                    (                                                             )
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_          %operator_           => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next     %operator_           => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next%operator_           => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    outputAnalysisPropertyOperator_                  =  outputAnalysisPropertyOperatorSequence                 (propertyOperators_                                           )
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    outputAnalysisDistributionOperator_              =  outputAnalysisDistributionOperatorRandomErrorPlynml    (                                  &
         &                                                                                                      randomErrorMinimum              , &
         &                                                                                                      randomErrorMaximum              , &
         &                                                                                                      errorPolynomialZeroPoint        , &
         &                                                                                                      randomErrorPolynomialCoefficient  &
         &                                                                                                     )
    ! Build a metallicity weight property operator.
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    outputAnalysisWeightPropertyOperatorSystmtcPolynomial_=  outputAnalysisPropertyOperatorSystmtcPolynomial        (metallicityErrorPolynomialZeroPoint,metallicitySystematicErrorPolynomialCoefficient)
    allocate(outputAnalysisPropertyOperatorMetallicity12LogNH_     )
    outputAnalysisPropertyOperatorMetallicity12LogNH_     =  outputAnalysisPropertyOperatorMetallicity12LogNH       (                                                          &
         &                                                                                                           Atomic_Mass(shortLabel="O")                               &
         &                                                                                                          )
    allocate(outputAnalysisPropertyOperatorFilterHighPass_         )
    outputAnalysisPropertyOperatorFilterHighPass_         =  outputAnalysisPropertyOperatorFilterHighPass           (0.0d0                                                        )
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    allocate(weightPropertyOperators_%next%next                    )
    weightPropertyOperators_          %operator_ => outputAnalysisPropertyOperatorMetallicity12LogNH_
    weightPropertyOperators_%next     %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    weightPropertyOperators_%next%next%operator_ => outputAnalysisPropertyOperatorFilterHighPass_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    outputAnalysisWeightPropertyOperator_                 =  outputAnalysisPropertyOperatorSequence                 (weightPropertyOperators_                                     )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    outputAnalysisPropertyUnoperator_                     =  outputAnalysisPropertyOperatorAntiLog10                (                                                             )
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisPropertyExtractor_                      )
    outputAnalysisPropertyExtractor_                      =  outputAnalysisPropertyExtractorMassStellar             (                                                             )
    ! Create an ISM metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    outputAnalysisWeightPropertyExtractor_                =  outputAnalysisPropertyExtractorMetallicityISM          (Abundances_Index_From_Name("O")                              )
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                &
         &                                                         var_str('massMetallicityBlanc2017'           ), &
         &                                                         var_str('Mass-metallicity relation'          ), &
         &                                                         var_str('massStellar'                        ), &
         &                                                         var_str('Stellar mass'                       ), &
         &                                                         var_str('M☉'                                 ), &
         &                                                         massSolar                                     , &
         &                                                         var_str('metallicityMean'                    ), &
         &                                                         var_str('Mean metallicity'                   ), &
         &                                                         var_str('dimensionless'                      ), &
         &                                                         0.0d0                                         , &
         &                                                         log10(masses)                                 , &
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
    nullify(galacticFilter_                                       )
    nullify(galacticFilterStellarMass_                            )
    nullify(galacticFilterStarFormationRate_                      )
    nullify(outputAnalysisDistributionOperator_                   )
    nullify(outputAnalysisWeightOperator_                         )
    nullify(outputAnalysisPropertyOperator_                       )
    nullify(outputAnalysisPropertyOperatorLog10_                  )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    nullify(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    nullify(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    nullify(outputAnalysisPropertyOperatorMetallicity12LogNH_     )
    nullify(outputAnalysisPropertyOperatorFilterHighPass_         )
    nullify(outputAnalysisPropertyUnoperator_                     )
    nullify(outputAnalysisWeightPropertyOperator_                 )
    nullify(outputAnalysisWeightPropertyExtractor_                )
    nullify(outputAnalysisPropertyExtractor_                      )
    nullify(cosmologyParametersData                               )
    nullify(cosmologyFunctionsData                                )
    nullify(propertyOperators_                                    )
    nullify(weightPropertyOperators_                              )
    nullify(filter_                                               ) 
    nullify(filters_                                              )    
    return
  end function massMetallicityBlanc2017ConstructorInternal

