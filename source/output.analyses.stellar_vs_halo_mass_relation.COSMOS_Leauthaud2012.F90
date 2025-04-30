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
  Implements a stellar vs halo mass relation analysis class.
  !!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <outputAnalysis name="outputAnalysisStellarVsHaloMassRelationLeauthaud2012">
   <description>A stellar vs halo mass relation output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisStellarVsHaloMassRelationLeauthaud2012
     !!{
     A stellar vs halo mass relation output analysis class.
     !!}
     private
     class           (outputAnalysisClass       ), pointer                     :: outputAnalysis_                        => null()
     class           (cosmologyParametersClass  ), pointer                     :: cosmologyParameters_                   => null()
     class           (cosmologyFunctionsClass   ), pointer                     :: cosmologyFunctions_                    => null()
     class           (darkMatterProfileDMOClass ), pointer                     :: darkMatterProfileDMO_                  => null()
     class           (virialDensityContrastClass), pointer                     :: virialDensityContrast_                 => null()
     class           (outputTimesClass          ), pointer                     :: outputTimes_                           => null()
     logical                                                                   :: computeScatter                                  , likelihoodNormalize
     integer         (c_size_t                  ), allocatable, dimension(:  ) :: likelihoodBins
     integer                                                                   :: redshiftInterval
     double precision                            , allocatable, dimension(:  ) :: systematicErrorPolynomialCoefficient            , systematicErrorMassHaloPolynomialCoefficient, &
          &                                                                       massStellarLogarithmicTarget
     double precision                            , allocatable, dimension(:,:) :: massStellarLogarithmicCovarianceTarget
     type            (varying_string            )                              :: analysisLabel
   contains
     final     ::                  stellarVsHaloMassRelationLeauthaud2012Destructor
     procedure :: analyze       => stellarVsHaloMassRelationLeauthaud2012Analyze
     procedure :: finalize      => stellarVsHaloMassRelationLeauthaud2012Finalize
     procedure :: reduce        => stellarVsHaloMassRelationLeauthaud2012Reduce
     procedure :: logLikelihood => stellarVsHaloMassRelationLeauthaud2012LogLikelihood
  end type outputAnalysisStellarVsHaloMassRelationLeauthaud2012

  interface outputAnalysisStellarVsHaloMassRelationLeauthaud2012
     !!{
     Constructors for the {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis class.
     !!}
     module procedure stellarVsHaloMassRelationLeauthaud2012ConstructorParameters
     module procedure stellarVsHaloMassRelationLeauthaud2012ConstructorInternal
  end interface outputAnalysisStellarVsHaloMassRelationLeauthaud2012

contains

  function stellarVsHaloMassRelationLeauthaud2012ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions     , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    use :: Input_Parameters        , only : inputParameters
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelationLeauthaud2012)                              :: self
    type            (inputParameters                                     ), intent(inout)               :: parameters
    double precision                                                      , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, systematicErrorMassHaloPolynomialCoefficient
    class           (cosmologyParametersClass                            ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                             ), pointer                     :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                           ), pointer                     :: darkMatterProfileDMO_
    class           (virialDensityContrastClass                          ), pointer                     :: virialDensityContrast_
    class           (outputTimesClass                                    ), pointer                     :: outputTimes_
    integer                                                                                             :: redshiftInterval
    logical                                                                                             :: computeScatter                      , likelihoodNormalize
    integer         (c_size_t                                            ), allocatable  , dimension(:) :: likelihoodBins

    ! Check and read parameters.
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    if (parameters%isPresent('systematicErrorMassHaloPolynomialCoefficient')) then
       allocate(systematicErrorMassHaloPolynomialCoefficient(parameters%count('systematicErrorMassHaloPolynomialCoefficient')))
    else
       allocate(systematicErrorMassHaloPolynomialCoefficient(1                                                               ))
    end if
    if (parameters%isPresent('likelihoodBins')) then
       allocate(likelihoodBins(parameters%count('likelihoodBins')))
       !![
       <inputParameter>
	 <name>likelihoodBins</name>
	 <source>parameters</source>
	 <description>If $>0$ then use only the mass bin given by this value in the likelihood calculation.</description>
       </inputParameter>
       !!]
    else
       allocate(likelihoodBins(               0                  ))
    end if
    !![
    <inputParameter>
      <name>redshiftInterval</name>
      <source>parameters</source>
      <variable>redshiftInterval</variable>
      <description>The redshift interval (1, 2, or 3) to use.</description>
    </inputParameter>
    <inputParameter>
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <variable>likelihoodNormalize</variable>
      <defaultValue>.false.</defaultValue>
      <description>If true, then normalize the likelihood to make it a probability density.</description>
    </inputParameter>
    <inputParameter>
      <name>computeScatter</name>
      <source>parameters</source>
      <variable>computeScatter</variable>
      <defaultValue>.false.</defaultValue>
      <description>If true, the scatter in log10(stellar mass) is computed. Otherwise, the mean is computed.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for stellar mass in the stellar vs halo mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorMassHaloPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorMassHaloPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for halo mass in the stellar vs halo mass relation.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisStellarVsHaloMassRelationLeauthaud2012(redshiftInterval,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassHaloPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,outputTimes_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    <objectDestructor name="outputTimes_"          />
    !!]
    return
  end function stellarVsHaloMassRelationLeauthaud2012ConstructorParameters

  function stellarVsHaloMassRelationLeauthaud2012ConstructorInternal(redshiftInterval,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassHaloPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,outputTimes_) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                    , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass                   , cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterHaloIsolated
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                  , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: ISO_Varying_String                    , only : var_str                                    , varying_string
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassHalo              , nodePropertyExtractorMassStellar
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Numerical_Interpolation               , only : gsl_interp_cspline
    use :: Numerical_Ranges                      , only : Make_Range                                 , rangeTypeLinear
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10    , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorFilterHighPass, outputAnalysisPropertyOperatorLog10, &
          &                                               outputAnalysisPropertyOperatorSequence     , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorProperty
    use :: String_Handling                       , only : operator(//)
    use :: Tables                                , only : table                                      , table1DGeneric
    use :: Virial_Density_Contrast               , only : fixedDensityTypeMean                       , virialDensityContrastClass                     , virialDensityContrastFixed
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelationLeauthaud2012)                                :: self
    integer                                                               , intent(in   )                 :: redshiftInterval
    logical                                                               , intent(in   )                 :: computeScatter                                                , likelihoodNormalize
    integer         (c_size_t                                            ), intent(in   ), dimension(:  ) :: likelihoodBins
    double precision                                                      , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                          , systematicErrorMassHaloPolynomialCoefficient
    class           (cosmologyParametersClass                            ), intent(in   ), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                             ), intent(inout), target         :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                           ), intent(inout), target         :: darkMatterProfileDMO_
    class           (virialDensityContrastClass                          ), intent(in   ), target         :: virialDensityContrast_
    class           (outputTimesClass                                    ), intent(inout), target         :: outputTimes_
    integer         (c_size_t                                            ), parameter                     :: massHaloCount                                         =26
    double precision                                                      , allocatable  , dimension(:  ) :: massHalo                                                      , massStellarDataLogarithmic                                         , &
         &                                                                                                   massHaloMeanDataLogarithmic                                   , massHaloLowDataLogarithmic                                         , &
         &                                                                                                   massHaloHighDataLogarithmic                                   , massHaloErrorDataLogarithmic                                       , &
         &                                                                                                   massStellarLogarithmicTarget                                  , massHaloHighData                                                   , &
         &                                                                                                   massHaloMeanData                                              , massHaloLowData                                                    , &
         &                                                                                                   massStellarData
    double precision                                                      , allocatable  , dimension(:,:) :: outputWeight                                                  , massStellarLogarithmicCovarianceTarget
    type            (galacticFilterHaloIsolated                          ), pointer                       :: galacticFilter_
    type            (outputAnalysisDistributionOperatorIdentity          ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorProperty                ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorLog10                 ), pointer                       :: outputAnalysisPropertyOperatorLog10_                         , outputAnalysisWeightPropertyOperatorLog10_                          , &
         &                                                                                                   outputAnalysisWeightPropertyOperatorLog10Second_
    type            (outputAnalysisPropertyOperatorAntiLog10             ), pointer                       :: outputAnalysisPropertyUnoperator_                            , outputAnalysisWeightPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorSequence              ), pointer                       :: outputAnalysisWeightPropertyOperator_                        , outputAnalysisPropertyOperator_                                     , &
         &                                                                                                   outputAnalysisWeightPropertyOperatorNormalized_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc     ), pointer                       :: outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial     ), pointer                       :: outputAnalysisWeightPropertyOperatorSystmtcPolynomial_       , outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorFilterHighPass        ), pointer                       :: outputAnalysisWeightPropertyOperatorFilterHighPass_          , outputAnalysisWeightPropertyOperatorFilterHighPassNormalized_
    type            (nodePropertyExtractorMassHalo                       ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellar                    ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (propertyOperatorList                                ), pointer                       :: propertyOperators_                                           , propertyOperatorsMassHalo_                                          , &
         &                                                                                                   propertyOperatorsNormalized_
    type            (cosmologyParametersSimple                           ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                      ), pointer                       :: cosmologyFunctionsData
    type            (virialDensityContrastFixed                          ), pointer                       :: virialDensityContrastDefinition_
    type            (surveyGeometryFullSky                               ), pointer                       :: surveyGeometry_
    double precision                                                      , parameter                     :: errorPolynomialZeroPoint                              =11.3d0, errorPolynomialMassHaloZeroPoint                             =12.0d0
    double precision                                                      , parameter                     :: covarianceLarge                                       = 1.0d4
    integer         (c_size_t                                            )                                :: iBin
    double precision                                                                                      :: massStellarLimit                                              , redshiftMinimum                                                   , &
         &                                                                                                   redshiftMaximum                                               , massHaloMinimum                                                   , &
         &                                                                                                   massHaloMaximum                                               , widthFilter
    type            (varying_string                                      )                                :: analysisLabel                                                 , weightPropertyLabel                                               , &
         &                                                                                                   weightPropertyDescription                                     , groupRedshiftName
    type            (hdf5Object                                          )                                :: fileData                                                      , groupRedshift
    type            (table1DGeneric                                      )                                :: interpolator
    character       (len=4                                               )                                :: redshiftMinimumLabel                                          , redshiftMaximumLabel
    !![
    <constructorAssign variables="redshiftInterval, likelihoodBins, likelihoodNormalize, computeScatter, systematicErrorPolynomialCoefficient, systematicErrorMassHaloPolynomialCoefficient, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *outputTimes_"/>
    !!]

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
       call Error_Report('redshiftInterval ∈ {1,2,3}'//{introspection:location})
    end select
    widthFilter=0.05d0
    write (redshiftMinimumLabel,'(f4.2)') redshiftMinimum
    write (redshiftMaximumLabel,'(f4.2)') redshiftMaximum
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create output time weights.
    allocate(outputWeight(massHaloCount,outputTimes_%count()))
    outputWeight(1,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,massStellarLimit,allowSingleEpoch=.true.)
    forall(iBin=2:massHaloCount)
       outputWeight(iBin,:)=outputWeight(1,:)
    end forall
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple     (                            &amp;
         &amp;                       OmegaMatter    = 0.25800d0, &amp;
         &amp;                       OmegaDarkEnergy= 0.74200d0, &amp;
         &amp;                       HubbleConstant =72.00000d0, &amp;
         &amp;                       temperatureCMB = 2.72548d0, &amp;
         &amp;                       OmegaBaryon    = 0.04385d0  &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
      cosmologyFunctionsMatterLambda(                            &amp;
         &amp;                       cosmologyParametersData     &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! Read observational data and convert masses to logarithmic.
    !$ call hdf5Access%set()
    call fileData%openFile(char(inputPath(pathTypeDataStatic))//"observations/stellarHaloMassRelation/stellarHaloMassRelation_COSMOS_Leauthaud2012.hdf5",readOnly=.true.)
    groupRedshiftName=var_str('redshiftInterval')//redshiftInterval
    groupRedshift=fileData%openGroup(char(groupRedshiftName))
    call groupRedshift%readDataset('massStellar' ,massStellarData )
    call groupRedshift%readDataset('massHaloMean',massHaloMeanData)
    call groupRedshift%readDataset('massHaloLow' ,massHaloLowData )
    call groupRedshift%readDataset('massHaloHigh',massHaloHighData)
    call groupRedshift%close      (                               )
    call fileData     %close      (                               )
    !$ call hdf5Access%unset()
    ! Create bins in halo mass.
    massHaloMinimum=massHaloMeanData(1                     )
    massHaloMaximum=massHaloMeanData(size(massHaloMeanData))
    massHalo       =Make_Range(log10(massHaloMinimum),log10(massHaloMaximum),int(massHaloCount),rangeType=rangeTypeLinear)
    ! Find a spline fit to the observed data, and compute the uncertainty in logarithm of halo mass.
    allocate(massHaloErrorDataLogarithmic(size(massStellarData)))
    massStellarDataLogarithmic  =log(massStellarData )
    massHaloMeanDataLogarithmic =log(massHaloMeanData)
    massHaloLowDataLogarithmic  =log(massHaloLowData )
    massHaloHighDataLogarithmic =log(massHaloHighData)
    massHaloErrorDataLogarithmic=+0.5d0                         &
         &                       *(                             &
         &                         +massHaloHighDataLogarithmic &
         &                         -massHaloLowDataLogarithmic  &
         &                        )
    call interpolator%create  (massHaloMeanDataLogarithmic ,tableCount=2,interpolationType=gsl_interp_cspline)
    call interpolator%populate(massStellarDataLogarithmic  ,table     =1                                     )
    call interpolator%populate(massHaloErrorDataLogarithmic,table     =2                                     )
    ! Interpolate observational data to model points.
    allocate(massStellarLogarithmicTarget          (massHaloCount              ))
    allocate(massStellarLogarithmicCovarianceTarget(massHaloCount,massHaloCount))
    massStellarLogarithmicCovarianceTarget=0.0d0
    do iBin=1,massHaloCount
       massStellarLogarithmicTarget          (iBin     )=+  interpolator%interpolate        (massHalo(iBin)*log(10.0d0),table=1)
       massStellarLogarithmicCovarianceTarget(iBin,iBin)=+(                                                                      &
            &                                              +interpolator%interpolateGradient(massHalo(iBin)*log(10.0d0),table=1) &
            &                                              *interpolator%interpolate        (massHalo(iBin)*log(10.0d0),table=2) &
            &                                             )**2
    end do
    call interpolator%destroy()
    massStellarLogarithmicTarget          =massStellarLogarithmicTarget          /log(10.0d0)
    massStellarLogarithmicCovarianceTarget=massStellarLogarithmicCovarianceTarget/log(10.0d0)**2
    if (size(self%likelihoodBins) > 0) then
       ! Assume that only a subset of bins of the relation are to be populated. Set the target dataset in all other bins to zero
       ! so that they do not contribute to the likelihood.
       if (any(self%likelihoodBins > massHaloCount .or. self%likelihoodBins < 1_c_size_t)) call Error_Report('likelihoodBins is out of range'//{introspection:location})
       do iBin=1,massHaloCount
          if (.not.any(self%likelihoodBins == iBin)) then
             massStellarLogarithmicTarget          (iBin     )=0.0d0
             massStellarLogarithmicCovarianceTarget(iBin,iBin)=covarianceLarge
          end if
       end do
    end if
    allocate(self%massStellarLogarithmicTarget          ,source=massStellarLogarithmicTarget          )
    allocate(self%massStellarLogarithmicCovarianceTarget,source=massStellarLogarithmicCovarianceTarget)
    ! Build a filter which select central galaxies. Note that no filter on stellar mass is included here. While Leauthaud et
    ! al. (2012) did include a threshold on stellar mass (see their Table 1) in their analysis, what we are fitting to here is
    ! their derived model - they accounted for the effects of the stellar mass threshold when deriving this model.
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterHaloIsolated()"/>
    !!]
    ! Build identity distribution operator.
    allocate   (outputAnalysisDistributionOperator_                          )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"                           constructor="outputAnalysisDistributionOperatorIdentity            (                                                                                                                              )"/>
    !!]
    ! Build a sequence (log10(), polynomial systematic) property operator.
    allocate   (outputAnalysisPropertyOperatorSystmtcPolynomial_             )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_"              constructor="outputAnalysisPropertyOperatorSystmtcPolynomial       (errorPolynomialMassHaloZeroPoint,systematicErrorMassHaloPolynomialCoefficient                                                 )"/>
    !!]
    allocate   (outputAnalysisPropertyOperatorLog10_                         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"                          constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                                                              )"/>
    !!]
    allocate(propertyOperatorsMassHalo_     )
    allocate(propertyOperatorsMassHalo_%next)
    propertyOperatorsMassHalo_     %operator_ => outputAnalysisPropertyOperatorLog10_
    propertyOperatorsMassHalo_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                          )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                               constructor="outputAnalysisPropertyOperatorSequence                (propertyOperatorsMassHalo_                                                                                                    )"/>
    !!]    
    ! Build a sequence (log10, polynomial systematic, anti-log10, cosmological luminosity distance, high-pass filter) of weight property operators.
    allocate   (outputAnalysisWeightPropertyOperatorFilterHighPassNormalized_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorFilterHighPassNormalized_" constructor="outputAnalysisPropertyOperatorFilterHighPass          (log10(massStellarLimit),widthFilter,normalized=.true.                                                                         )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorFilterHighPass_          )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorFilterHighPass_"           constructor="outputAnalysisPropertyOperatorFilterHighPass          (log10(massStellarLimit),widthFilter                                                                                           )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_       )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_"        constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc       (cosmologyFunctions_     ,cosmologyFunctionsData              ,outputTimes_                                                    )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorSystmtcPolynomial_       )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"        constructor="outputAnalysisPropertyOperatorSystmtcPolynomial       (errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient                                                                 )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorLog10_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"                    constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                                                              )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorLog10Second_             )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10Second_"              constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                                                              )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorAntiLog10_               )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorAntiLog10_"                constructor="outputAnalysisPropertyOperatorAntiLog10               (                                                                                                                              )"/>
    !!]
    allocate       (propertyOperators_                                        )
    allocate       (propertyOperators_          %next                         )
    allocate       (propertyOperators_          %next%next                    )
    allocate       (propertyOperators_          %next%next%next               )
    allocate       (propertyOperators_          %next%next%next%next          )
    allocate       (propertyOperatorsNormalized_                              )
    allocate       (propertyOperatorsNormalized_%next                         )
    allocate       (propertyOperatorsNormalized_%next%next                    )
    allocate       (propertyOperatorsNormalized_%next%next%next               )
    allocate       (propertyOperatorsNormalized_%next%next%next%next          )
    allocate       (propertyOperatorsNormalized_%next%next%next%next%next     )
    propertyOperators_                                   %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperators_          %next                    %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperators_          %next%next               %operator_ => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperators_          %next%next%next          %operator_ => outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_          %next%next%next%next     %operator_ => outputAnalysisWeightPropertyOperatorLog10Second_
    propertyOperatorsNormalized_                         %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperatorsNormalized_%next                    %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperatorsNormalized_%next%next               %operator_ => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperatorsNormalized_%next%next%next          %operator_ => outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperatorsNormalized_%next%next%next%next     %operator_ => outputAnalysisWeightPropertyOperatorLog10Second_
    propertyOperatorsNormalized_%next%next%next%next%next%operator_ => outputAnalysisWeightPropertyOperatorFilterHighPassNormalized_
    ! Create a stellar mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                          )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                        constructor="nodePropertyExtractorMassStellar                      (                                                                                                                              )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                         constructor="outputAnalysisPropertyOperatorSequence                (propertyOperators_                                                                                                            )"/>
    !!]
   allocate(outputAnalysisWeightPropertyOperatorNormalized_                  )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorNormalized_"               constructor="outputAnalysisPropertyOperatorSequence                (propertyOperatorsNormalized_                                                                                                  )"/>
    !!]
    ! Build weight operator.
    allocate   (outputAnalysisWeightOperator_                                )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                                 constructor="outputAnalysisWeightOperatorProperty                  (outputAnalysisWeightPropertyExtractor_,outputAnalysisWeightPropertyOperatorNormalized_                                        )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                             constructor="outputAnalysisPropertyOperatorAntiLog10               (                                                                                                                              )"/>
    !!]
    ! Create a halo mass weight property extractor.
    allocate(virialDensityContrastDefinition_                                )
    !![
    <referenceConstruct object="virialDensityContrastDefinition_"                              constructor="virialDensityContrastFixed                            (200.0d0,fixedDensityTypeMean,2.0d0,cosmologyParameters_,cosmologyFunctions_                                                   )"/>
    !!]
    allocate(nodePropertyExtractor_                                          )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                        constructor="nodePropertyExtractorMassHalo                         (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    !!]
    ! Build the object.
    if (computeScatter) then
       analysisLabel            =var_str('stellarHaloMassRelationScatterLeauthaud2012z')//redshiftInterval
       weightPropertyLabel      =var_str('massStellarLog10Scatter'                     )
       weightPropertyDescription=var_str('σ_{log₁₀(Stellar mass/M☉)}'                  )
       allocate(outputAnalysisScatterFunction1D :: self%outputAnalysis_)
       ! For the scatter we need to set an appropriate target and covariance for likelihood calculation. This measurement is based
       ! on the constraint on σ_{log₁₀L}=0.16±0.04 from More et al. (2009; MNRAS; 392; 801) for SDSS galaxies.
       do iBin=1,massHaloCount
          if (size(self%likelihoodBins) == 0 .or. any(self%likelihoodBins == iBin)) then
             massStellarLogarithmicTarget          (iBin     )=0.16d0
             massStellarLogarithmicCovarianceTarget(iBin,iBin)=0.04d0**2
          end if
       end do
    else
       analysisLabel            =var_str('stellarHaloMassRelationLeauthaud2012z'       )//redshiftInterval
       weightPropertyLabel      =var_str('massStellarLog10'                            )
       weightPropertyDescription=var_str('⟨log₁₀(Stellar mass/M☉)⟩'                    )
       allocate(outputAnalysisMeanFunction1D    :: self%outputAnalysis_)
    end if
    self%analysisLabel=analysisLabel
    ! Construct the analysis objects. Note that we use a "Poisson" model for the covariance here, not the "binomial" model which
    ! is appropriate to counting analyses (e.g. mass functions), but not to this type of mean or scatter analysis.
    select type (outputAnalysis_ => self%outputAnalysis_)
    type is (outputAnalysisScatterFunction1D)
       !![
       <referenceConstruct isResult="yes" object="outputAnalysis_">
        <constructor>
         outputAnalysisScatterFunction1D(                                                                                                                                                 &amp;
          &amp;                                                  analysisLabel                                                                                                          , &amp;
          &amp;                                                  var_str('Scatter in stellar-halo mass relation for $')//redshiftMinimumLabel//' \le z &lt; '//redshiftMaximumLabel//'$', &amp;
          &amp;                                                  var_str('massHalo'                                    )                                                                , &amp;
          &amp;                                                  var_str('Halo mass'                                   )                                                                , &amp;
          &amp;                                                  var_str('M☉'                                          )                                                                , &amp;
          &amp;                                                  massSolar                                                                                                              , &amp;
          &amp;                                                  weightPropertyLabel                                                                                                    , &amp;
          &amp;                                                  weightPropertyDescription                                                                                              , &amp;
          &amp;                                                  var_str(' '                                           )                                                                , &amp;
          &amp;                                                  0.0d0                                                                                                                  , &amp;
          &amp;                                                  massHalo                                                                                                               , &amp;
          &amp;                                                  0_c_size_t                                                                                                             , &amp;
          &amp;                                                  outputWeight                                                                                                           , &amp;
          &amp;                                                  nodePropertyExtractor_                                                                                                 , &amp;
          &amp;                                                  outputAnalysisWeightPropertyExtractor_                                                                                 , &amp;
          &amp;                                                  outputAnalysisPropertyOperator_                                                                                        , &amp;
          &amp;                                                  outputAnalysisWeightPropertyOperator_                                                                                  , &amp;
          &amp;                                                  outputAnalysisPropertyUnoperator_                                                                                      , &amp;
          &amp;                                                  outputAnalysisWeightOperator_                                                                                          , &amp;
          &amp;                                                  outputAnalysisDistributionOperator_                                                                                    , &amp;
          &amp;                                                  galacticFilter_                                                                                                        , &amp;
          &amp;                                                  outputTimes_                                                                                                           , &amp;
          &amp;                                                  outputAnalysisCovarianceModelPoisson                                                                                   , &amp;
          &amp;                          likelihoodNormalize    =likelihoodNormalize                                                                                                    , &amp;
          &amp;                          xAxisLabel             =var_str('M_\mathrm{halo}/\mathrm{M}_\odot'            )                                                                , &amp;
          &amp;                          yAxisLabel             =var_str('\sigma_{\log_{10}(M_\star/\mathrm{M}_\odot)}')                                                                , &amp;
          &amp;                          xAxisIsLog             =.true.                                                                                                                 , &amp;
          &amp;                          yAxisIsLog             =.false.                                                                                                                , &amp;
          &amp;                          targetLabel            =var_str('More et al. (2009)'                          )                                                                , &amp;
          &amp;                          scatterValueTarget     =massStellarLogarithmicTarget                                                                                           , &amp;
          &amp;                          scatterCovarianceTarget=massStellarLogarithmicCovarianceTarget                                                                                   &amp;
          &amp;                         )
        </constructor>
       </referenceConstruct>
       !!]
    type is (outputAnalysisMeanFunction1D   )
       !![
       <referenceConstruct isResult="yes" object="outputAnalysis_">
        <constructor>
         outputAnalysisMeanFunction1D   (                                                                                                                                   &amp;
          &amp;                                               analysisLabel                                                                                               , &amp;
          &amp;                                               var_str('Stellar-halo mass relation for $')//redshiftMinimumLabel//' \le z &lt; '//redshiftMaximumLabel//'$', &amp;
          &amp;                                               var_str('massHalo'                             )                                                            , &amp;
          &amp;                                               var_str('Halo mass'                            )                                                            , &amp;
          &amp;                                               var_str('M☉'                                   )                                                            , &amp;
          &amp;                                               massSolar                                                                                                   , &amp;
          &amp;                                               weightPropertyLabel                                                                                         , &amp;
          &amp;                                               weightPropertyDescription                                                                                   , &amp;
          &amp;                                               var_str(' '                                    )                                                            , &amp;
          &amp;                                               0.0d0                                                                                                       , &amp;
          &amp;                                               massHalo                                                                                                    , &amp;
          &amp;                                               0_c_size_t                                                                                                  , &amp;
          &amp;                                               outputWeight                                                                                                , &amp;
          &amp;                                               nodePropertyExtractor_                                                                                      , &amp;
          &amp;                                               outputAnalysisWeightPropertyExtractor_                                                                      , &amp;
          &amp;                                               outputAnalysisPropertyOperator_                                                                             , &amp;
          &amp;                                               outputAnalysisWeightPropertyOperator_                                                                       , &amp;
          &amp;                                               outputAnalysisPropertyUnoperator_                                                                           , &amp;
          &amp;                                               outputAnalysisWeightOperator_                                                                               , &amp;
          &amp;                                               outputAnalysisDistributionOperator_                                                                         , &amp;
          &amp;                                               galacticFilter_                                                                                             , &amp;
          &amp;                                               outputTimes_                                                                                                , &amp;
          &amp;                                               outputAnalysisCovarianceModelPoisson                                                                        , &amp;
          &amp;                          likelihoodNormalize =likelihoodNormalize                                                                                         , &amp;
          &amp;                          xAxisLabel          =var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'   )                                                            , &amp;
          &amp;                          yAxisLabel          =var_str('$\log_{10}(M_\star/\mathrm{M}_\odot)$')                                                            , &amp;
          &amp;                          xAxisIsLog          =.true.                                                                                                      , &amp;
          &amp;                          yAxisIsLog          =.false.                                                                                                     , &amp;
          &amp;                          targetLabel         =var_str('Leauthaud et al. (2012)'              )                                                            , &amp;
          &amp;                          meanValueTarget     =massStellarLogarithmicTarget                                                                                , &amp;
          &amp;                          meanCovarianceTarget=massStellarLogarithmicCovarianceTarget                                                                        &amp;
          &amp;                         )
        </constructor>
       </referenceConstruct>
       !!]
    end select
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                              />
    <objectDestructor name="cosmologyParametersData"                                      />
    <objectDestructor name="cosmologyFunctionsData"                                       />
    <objectDestructor name="galacticFilter_"                                              />
    <objectDestructor name="outputAnalysisDistributionOperator_"                          />
    <objectDestructor name="outputAnalysisWeightOperator_"                                />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                         />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"             />
    <objectDestructor name="outputAnalysisPropertyOperator_"                              />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorNormalized_"              />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorFilterHighPass_"          />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorFilterHighPassNormalized_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_"       />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"       />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"                   />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10Second_"             />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorAntiLog10_"               />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                        />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                            />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                       />
    <objectDestructor name="virialDensityContrastDefinition_"                             />
    <objectDestructor name="nodePropertyExtractor_"                                       />
    !!]
    nullify(propertyOperatorsMassHalo_  )
    nullify(propertyOperators_          )
    nullify(propertyOperatorsNormalized_)
    return
  end function stellarVsHaloMassRelationLeauthaud2012ConstructorInternal

  subroutine stellarVsHaloMassRelationLeauthaud2012Analyze(self,node,iOutput)
    !!{
    Implement a stellarVsHaloMassRelationLeauthaud2012 output analysis.
    !!}
    implicit none
    class  (outputAnalysisStellarVsHaloMassRelationLeauthaud2012), intent(inout) :: self
    type   (treeNode                                            ), intent(inout) :: node
    integer(c_size_t                                            ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine stellarVsHaloMassRelationLeauthaud2012Analyze

  subroutine stellarVsHaloMassRelationLeauthaud2012Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis class.
    !!}
    implicit none
    type(outputAnalysisStellarVsHaloMassRelationLeauthaud2012), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"       />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%virialDensityContrast_"/>
    <objectDestructor name="self%outputTimes_"          />
    !!]
    return
  end subroutine stellarVsHaloMassRelationLeauthaud2012Destructor

  subroutine stellarVsHaloMassRelationLeauthaud2012Reduce(self,reduced)
    !!{
    Implement reduction for the {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisStellarVsHaloMassRelationLeauthaud2012), intent(inout) :: self
    class(outputAnalysisClass                                 ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisStellarVsHaloMassRelationLeauthaud2012)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine stellarVsHaloMassRelationLeauthaud2012Reduce

  subroutine stellarVsHaloMassRelationLeauthaud2012Finalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisStellarVsHaloMassRelationLeauthaud2012), intent(inout)           :: self
    type (varying_string                                      ), intent(in   ), optional :: groupName
    type (hdf5Object                                          )               , target   :: analysesGroup, subGroup
    type (hdf5Object                                          )               , pointer  :: inGroup
    type (hdf5Object                                          )                          :: analysisGroup

    call self%outputAnalysis_%finalize(groupName)
    ! Overwrite the log-likelihood - this allows us to handle cases where the model is zero everywhere.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup(char(self%analysisLabel))
    call    analysisGroup%writeAttribute(self%logLikelihood(),'logLikelihood')
    call    analysisGroup%close         (                                    )
    if (present(groupName)) &
         & call subGroup %close         (                                    )
    call    analysesGroup%close         (                                    )
    !$ call hdf5Access%unset()
    return
  end subroutine stellarVsHaloMassRelationLeauthaud2012Finalize

  double precision function stellarVsHaloMassRelationLeauthaud2012LogLikelihood(self) result(logLikelihood)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily stellarVsHaloMassRelationLeauthaud2012} output analysis.
    !!}
    use :: Error                       , only : Error_Report
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisStellarVsHaloMassRelationLeauthaud2012), intent(inout)                 :: self
    double precision                                                      , parameter                     :: massStellarLogarithmicTiny              =1.0d-3
    double precision                                                      , allocatable  , dimension(:,:) :: massStellarLogarithmicCovarianceCombined       , massStellarLogarithmicCovarianceCombinedSelected, &
         &                                                                                                   massStellarLogarithmicCovariance               , massStellarLogarithmicCovarianceTarget
    double precision                                                      , allocatable  , dimension(:  ) :: massStellarLogarithmicDifference               , massStellarLogarithmicDifferenceSelected        , &
         &                                                                                                   massStellarLogarithmic                         , massStellarLogarithmicTarget
    type            (vector                                              )                                :: residual
    type            (matrix                                              )                                :: covariance
    integer                                                                                               :: i                                              , j                                               , &
         &                                                                                                   status

    select type (outputAnalysis_ => self%outputAnalysis_)
    class is (outputAnalysisMeanFunction1D   )
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(   meanValue=massStellarLogarithmic,   meanCovariance=massStellarLogarithmicCovariance)
       allocate(massStellarLogarithmicTarget          ,source=self%massStellarLogarithmicTarget          )
       allocate(massStellarLogarithmicCovarianceTarget,source=self%massStellarLogarithmicCovarianceTarget)
    class is (outputAnalysisScatterFunction1D)
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(scatterValue=massStellarLogarithmic,scatterCovariance=massStellarLogarithmicCovariance)
       allocate(massStellarLogarithmicTarget          ,mold  =self%massStellarLogarithmicTarget          )
       allocate(massStellarLogarithmicCovarianceTarget,mold  =self%massStellarLogarithmicCovarianceTarget)
       massStellarLogarithmicCovarianceTarget=0.0d0
       do i=1,size(massStellarLogarithmicTarget)
          massStellarLogarithmicTarget          (i  )=0.16d0
          massStellarLogarithmicCovarianceTarget(i,i)=0.04d0**2
       end do
    class default
       logLikelihood=+outputAnalysis_%logLikelihood()
       return
    end select
    if     (                                                                                                                       &
         &   (size(self%likelihoodBins) == 0 .and. any(massStellarLogarithmic                      <= massStellarLogarithmicTiny)) &
         &  .or.                                                                                                                   &
         &                                         any(massStellarLogarithmic(self%likelihoodBins) <= massStellarLogarithmicTiny)  &
         & ) then
       ! If any active bins contain zero galaxies, judge this model to be improbable.
       logLikelihood=                     logImprobable
    else
       ! Compute difference with the target dataset.
       allocate(massStellarLogarithmicDifference        ,mold=massStellarLogarithmic          )
       allocate(massStellarLogarithmicCovarianceCombined,mold=massStellarLogarithmicCovariance)
       massStellarLogarithmicDifference        =+massStellarLogarithmic          -massStellarLogarithmicTarget
       massStellarLogarithmicCovarianceCombined=+massStellarLogarithmicCovariance+massStellarLogarithmicCovarianceTarget
       ! Construct a reduced set of bins.
       if (size(self%likelihoodBins) > 0) then
          allocate(massStellarLogarithmicDifferenceSelected        (size(self%likelihoodBins)                          ))
          allocate(massStellarLogarithmicCovarianceCombinedSelected(size(self%likelihoodBins),size(self%likelihoodBins)))
          do i=1,size(self%likelihoodBins)
             massStellarLogarithmicDifferenceSelected           (i  )=massStellarLogarithmicDifference        (self%likelihoodBins(i)                       )
             do j=1,size(self%likelihoodBins)
                massStellarLogarithmicCovarianceCombinedSelected(i,j)=massStellarLogarithmicCovarianceCombined(self%likelihoodBins(i),self%likelihoodBins(j))
             end do
          end do
       else
          allocate(massStellarLogarithmicDifferenceSelected        ,source=massStellarLogarithmicDifference        )
          allocate(massStellarLogarithmicCovarianceCombinedSelected,source=massStellarLogarithmicCovarianceCombined)
       end if
       ! Construct residual vector and covariance matrix.
       residual  =vector(massStellarLogarithmicDifferenceSelected        )
       covariance=matrix(massStellarLogarithmicCovarianceCombinedSelected)
       ! Compute the log-likelihood.
       logLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
       if (status == GSL_Success) then
          if (self%likelihoodNormalize)                                                                  &
               & logLikelihood=+logLikelihood                                                            &
               &               -0.5d0*covariance%logarithmicDeterminant()                                &
               &               -0.5d0*dble(size(massStellarLogarithmicDifferenceSelected))*log(2.0d0*Pi)
       else
          logLikelihood       =+logImprobable
       end if
    end if
    return
  end function stellarVsHaloMassRelationLeauthaud2012LogLikelihood
