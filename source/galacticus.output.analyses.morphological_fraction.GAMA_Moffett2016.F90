!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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
  
  !# <outputAnalysis name="outputAnalysisMorphologicalFractionGAMAMoffett2016">
  !#  <description>A morphological fraction output analysis class for the analysis of \cite{moffett_galaxy_2016}.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisMorphologicalFractionGAMAMoffett2016
     !% A morphological fraction output analysis class.
     private
  end type outputAnalysisMorphologicalFractionGAMAMoffett2016

  interface outputAnalysisMorphologicalFractionGAMAMoffett2016
     !% Constructors for the ``morphologicalFractionGAMAMoffett2016'' output analysis class.
     module procedure morphologicalFractionGAMAMoffett2016ConstructorParameters
     module procedure morphologicalFractionGAMAMoffett2016ConstructorInternal
  end interface outputAnalysisMorphologicalFractionGAMAMoffett2016

contains

  function morphologicalFractionGAMAMoffett2016ConstructorParameters(parameters) result (self)
    !% Constructor for the ``morphologicalFractionGAMAMoffett2016'' output analysis class which takes a parameter set as input.
    use Cosmology_Functions
    use Input_Parameters
    implicit none
    type            (outputAnalysisMorphologicalFractionGAMAMoffett2016)                              :: self
    type            (inputParameters                                   ), intent(inout)               :: parameters
    double precision                                                    , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                           ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                  ), pointer                     :: outputTimes_
    double precision                                                                                  :: ratioEarlyType                      , ratioEarlyTypeError             , &
         &                                                                                               randomErrorMinimum                  , randomErrorMaximum 

    
    ! Check and read parameters.
    allocate(systematicErrorPolynomialCoefficient(max(1,parameters%count('systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorPolynomialCoefficient(max(1,parameters%count(    'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !# <inputParameter>
    !#   <name>ratioEarlyType</name>
    !#   <cardinality>0..1</cardinality>
    !#   <defaultValue>0.5d0</defaultValue>
    !#   <description>The minimum spheroid-to-total ratio for a galaxy to be classified as ``early-type'' when constructing the \gls{gama} early-type fraction function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>ratioEarlyTypeError</name>
    !#   <cardinality>0..1</cardinality>
    !#   <defaultValue>0.3d0</defaultValue>
    !#   <description>The error in spheroid fraction to be used when constructing the \gls{gama} early-type fraction function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
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
    !# <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisMorphologicalFractionGAMAMoffett2016(ratioEarlyType,ratioEarlyTypeError,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    !# <objectDestructor name="outputTimes_"       />
    return
  end function morphologicalFractionGAMAMoffett2016ConstructorParameters

  function morphologicalFractionGAMAMoffett2016ConstructorInternal(ratioEarlyType,ratioEarlyTypeError,systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyFunctions_,outputTimes_) result (self)
    !% Constructor for the ``morphologicalFractionGAMAMoffett2016'' output analysis class for internal use.
    use Memory_Management
    use IO_HDF5
    use Galacticus_Paths  
    use Output_Times
    use Output_Analysis_Property_Operators
    use Output_Analysis_Property_Extractions
    use Output_Analysis_Distribution_Operators
    use Output_Analysis_Weight_Operators
    use Output_Analysis_Utilities
    use Cosmology_Parameters
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    implicit none
    type            (outputAnalysisMorphologicalFractionGAMAMoffett2016   )                                :: self
    double precision                                                       , intent(in   )                 :: ratioEarlyType                                          , ratioEarlyTypeError                     , &
         &                                                                                                    randomErrorMinimum                                      , randomErrorMaximum
    double precision                                                       , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                    , randomErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                              ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                     ), intent(inout), target         :: outputTimes_
    integer                                                                , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                       , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                       , allocatable  , dimension(:  ) :: masses
    double precision                                                       , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterStellarMass                            ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml  ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity                 ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence               ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                  ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10              ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorNormal                 ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyExtractorMassStellar           ), pointer                       :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyExtractorMassStellarMorphology ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc      ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial      ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (cosmologyParametersSimple                            ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                       ), pointer                       :: cosmologyFunctionsData
    type            (propertyOperatorList                                 ), pointer                       :: propertyOperators_
    double precision                                                       , parameter                     :: errorPolynomialZeroPoint                        =11.3d00
    integer         (c_size_t                                             ), parameter                     :: bufferCount                                     =10
    integer         (c_size_t                                             )                                :: iBin                                                    , binCount
    type            (surveyGeometryBaldry2012GAMA                         )                                :: surveyGeometry_
    type            (hdf5Object                                           )                                :: dataFile

    ! Read masses at which fraction was measured.
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(galacticusPath(pathTypeDataStatic))//"observations/morphology/earlyTypeFractionGAMA.hdf5",readOnly=.true.)
    call dataFile%readDataset("mass"                                                                                  ,         masses)
    call dataFile%close      (                                                                                                        )
    !$ call hdf5Access%unset()
    ! Construct survey geometry.
    surveyGeometry_=surveyGeometryBaldry2012GAMA(cosmologyFunctions_)
    ! Compute weights that apply to each output redshift.
    binCount=size(masses,kind=c_size_t)
    call allocateArray(outputWeight,[binCount,outputTimes_%count()])
    do iBin=1,binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,masses(iBin))
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
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilter_                                       )
    galacticFilter_                                  =  galacticFilterStellarMass                           (massThreshold=1.0d8                                                       )
     ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    outputAnalysisWeightOperator_                    =  outputAnalysisWeightOperatorIdentity                (                                                                          )
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_      )
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_ =  outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc     (cosmologyFunctions_     ,cosmologyFunctionsData              ,outputTimes_)
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    outputAnalysisPropertyOperatorSystmtcPolynomial_ =  outputAnalysisPropertyOperatorSystmtcPolynomial     (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient             )
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    outputAnalysisPropertyOperatorLog10_                        =  outputAnalysisPropertyOperatorLog10      (                                                                          )
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    allocate(propertyOperators_%next%next                          )
    propertyOperators_          %operator_           => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next     %operator_           => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%next%operator_           => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    outputAnalysisPropertyOperator_                  =  outputAnalysisPropertyOperatorSequence              (propertyOperators_                                                        )
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    outputAnalysisDistributionOperator_              =  outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &
         &                                                                                                   randomErrorMinimum              , &
         &                                                                                                   randomErrorMaximum              , &
         &                                                                                                   errorPolynomialZeroPoint        , &
         &                                                                                                   randomErrorPolynomialCoefficient  &
         &                                                                                                  )
    ! Build a weight property operator.
    allocate(outputAnalysisWeightPropertyOperator_                 )
    outputAnalysisWeightPropertyOperator_            =  outputAnalysisPropertyOperatorNormal                (                                  &
         &                                                                                                   rangeLower  =ratioEarlyType     , &
         &                                                                                                   rangeUpper  =1.0d0              , &
         &                                                                                                   extentLower =0.0d0              , &
         &                                                                                                   extentUpper =1.0d0              , &
         &                                                                                                   rootVariance=ratioEarlyTypeError  &
         &                                                                                                  )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    outputAnalysisPropertyUnoperator_                =  outputAnalysisPropertyOperatorAntiLog10             (                                                                          )
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisPropertyExtractor_                      )
    outputAnalysisPropertyExtractor_                 =  outputAnalysisPropertyExtractorMassStellar          (                                                                          )
    ! Create a morpology weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    outputAnalysisWeightPropertyExtractor_           =  outputAnalysisPropertyExtractorMassStellarMorphology(                                                                          )
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                 &
         &                                                         var_str('morphologicalFractionGAMAMoffett2016'), &
         &                                                         var_str('Early-type fraction'                 ), &
         &                                                         var_str('massStellar'                         ), &
         &                                                         var_str('Stellar mass'                        ), &
         &                                                         var_str('M☉'                                  ), &
         &                                                         massSolar                                      , &
         &                                                         var_str('earlyTypeFraction'                   ), &
         &                                                         var_str('Early-type fraction'                 ), &
         &                                                         var_str(' '                                   ), &
         &                                                         0.0d0                                          , &
         &                                                         log10(masses)                                  , &
         &                                                         bufferCount                                    , &
         &                                                         outputWeight                                   , &
         &                                                         outputAnalysisPropertyExtractor_               , &
         &                                                         outputAnalysisWeightPropertyExtractor_         , &
         &                                                         outputAnalysisPropertyOperator_                , &
         &                                                         outputAnalysisWeightPropertyOperator_          , &
         &                                                         outputAnalysisPropertyUnoperator_              , &
         &                                                         outputAnalysisWeightOperator_                  , &
         &                                                         outputAnalysisDistributionOperator_            , &
         &                                                         galacticFilter_                                , &
         &                                                         outputTimes_                                   , &
         &                                                         outputAnalysisCovarianceModelBinomial          , &
         &                                                         covarianceBinomialBinsPerDecade                , &
         &                                                         covarianceBinomialMassHaloMinimum              , &
         &                                                         covarianceBinomialMassHaloMaximum                &
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
  end function morphologicalFractionGAMAMoffett2016ConstructorInternal

