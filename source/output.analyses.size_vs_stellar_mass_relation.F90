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
  Implements a galaxy size vs stellar mass relation analysis class.
  !!}

  use, intrinsic :: ISO_C_Binding                             , only : c_size_t
  use            :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use            :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use            :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <outputAnalysis name="outputAnalysisSizeVsStellarMassRelation">
    <description>
      A size vs. stellar mass relation output analysis class. Target data is read from an \gls{hdf5} file specified by the
      {\normalfont \ttfamily [fileNameTarget]} parameter. This file must contain one or more groups named {\normalfont \ttfamily
      sampleN} where {\normalfont \ttfamily N} is an integer. Each such group specifies the galaxy size---stellar mass relation
      for one sample (a combination of redshift interval and any selection criteria), and must contain the following datasets and
      attributes:
      \begin{itemize}
       \item dataset {\normalfont \ttfamily massStellar}: stellar mass in units of $\mathrm{M}_\odot$;
       \item dataset {\normalfont \ttfamily radiusEffective}: effective radius in units of Mpc;
       \item dataset {\normalfont \ttfamily radiusEffectiveError}: uncertainty in effective radius in units of Mpc;
       \item dataset {\normalfont \ttfamily radiusEffectiveScatter}: scatter in effective radius in units of dex;
       \item dataset {\normalfont \ttfamily radiusEffectiveScatterError}: uncertainty in scatter in effective radius in units of dex;
       \item attribute {\normalfont \ttfamily redshiftMinimum}: the minimum redshift associated with this sample;
       \item attribute {\normalfont \ttfamily redshiftMaximum}: the maximum redshift associated with this sample.
      \end{itemize}
      While not required, it is recommended that each of these datasets has attributes {\normalfont \ttfamily description} and
      {\normalfont \ttfamily unitsInSI} that provide a description of the dataset, and the multiplicative factor needed to convert
      them to SI standard units, respectively.

      Additionally, the file must contain a {\normalfont \ttfamily cosmology} group that specifies the cosmological model assumed
      in constructing the dataset, and which has attributes:
      \begin{itemize}
       \item {\normalfont \ttfamily OmegaMatter}: the matter density in units of the critical density, $\Omega_\mathrm{M}$;
       \item {\normalfont \ttfamily OmegaDarkEnergy}: the dark energy density in units of the critical density, $\Omega_\Lambda$;
       \item {\normalfont \ttfamily OmegaBaryon}: the baryon density in units of the critical density, $\Omega_\mathrm{b}$;
       \item {\normalfont \ttfamily HubbleConstant}: the Hubble constant in units of km/s/Mpc.
      \end{itemize}

      Each {\normalfont \ttfamily sampleN} group must have an attribute {\normalfont \ttfamily selection} which specifies the
      selection criterion used in constructing the dataset. Allowed values are:      
      \begin{itemize}
       \item {\normalfont \ttfamily `none'}: no selection criterion will be applied;
       \item {\normalfont \ttfamily `star forming'}: only galaxies on or above the star forming main sequence are included;
       \item {\normalfont \ttfamily `quiescent}: only galaxies below the star forming main sequence are included.
      \end{itemize}
      For the {\normalfont ``\ttfamily star forming}'' and ``{\normalfont \ttfamily quiescent}'' options, a dataset {\normalfont
      \ttfamily mainSequenceSFR} must be specified in the {\normalfont \ttfamily sampleN} group which specifies the mean
      (of the logarithm of star formation rate in units of $\mathrm{M}_\odot/\hbox{yr}^{-1}$) of the star forming main sequence at
      the center of each bin, and an attribute {\normalfont \ttfamily offsetMainSequenceSFR} which specifies an offset below the
      mean of the star forming main sequence below which galaxies are considered to be quiescent. That is, a galaxy will be
      classified as quiescent if
      \begin{equation}
      \log_{10} ( \dot{\phi} / \mathrm{M}_\odot \hbox{yr}^{-1}) &lt; \log_{10} ( \dot{M}_{\star,ms} / \mathrm{M}_\odot \hbox{yr}^{-1}) - \Delta_{\star},
      \end{equation}
      where $\log_{10} ( \dot{\phi} / \mathrm{M}_\odot \hbox{yr}^{-1})$ is the star formation rate in the galaxy, $\log_{10} (
      \dot{M}_{\star,ms} / \mathrm{M}_\odot \hbox{yr}^{-1})$ is the mean of the star forming main sequence (as specified by
      {\normalfont \ttfamily mainSequenceSFR} and interpolated to the stellar mass of the galaxy, and $\Delta_{\star}$ is
      {\normalfont \ttfamily offsetMainSequenceSFR},
								  
      Lastly, the file must have two attributes used to identify and level the dataset:
      \begin{itemize}
       \item {\normalfont \ttfamily label}: a space-free label that will be appended to the analysis group in the output, e.g. {\normalfont \ttfamily vanDerWel2014};
       \item {\normalfont \ttfamily reference}: a reference for the dataset suitable for inclusion in figures, e.g. {\normalfont \ttfamily van der Wel et al. (2014)}.
      \end{itemize}
    </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSizeVsStellarMassRelation
     !!{
     A stellar vs halo mass relation output analysis class.
     !!}
     private
     class           (starFormationRateDisksClass              ), pointer                   :: starFormationRateDisks_                      => null()
     class           (starFormationRateSpheroidsClass          ), pointer                   :: starFormationRateSpheroids_                  => null()
     class           (starFormationRateNuclearStarClustersClass), pointer                   :: starFormationRateNuclearStarClusters_        => null()
     class           (outputAnalysisClass                     ), pointer                     :: outputAnalysis_                             => null()
     class           (cosmologyParametersClass                ), pointer                     :: cosmologyParameters_                        => null()
     class           (cosmologyFunctionsClass                 ), pointer                     :: cosmologyFunctions_                         => null()
     class           (outputTimesClass                        ), pointer                     :: outputTimes_                                => null()
     logical                                                                                 :: computeScatter                                       , likelihoodNormalize                            , &
          &                                                                                     likelihoodBinsAutomatic
     integer         (c_size_t                                ), allocatable, dimension(:  ) :: likelihoodBins
     integer                                                                                 :: sample
     double precision                                          , allocatable, dimension(:  ) :: radiusEffectiveLogarithmicTarget                     , radiusEffectiveScatterTarget                   , &
          &                                                                                     systematicErrorPolynomialCoefficient                 , systematicErrorMassStellarPolynomialCoefficient, &
          &                                                                                     randomErrorMassStellarPolynomialCoefficient                                                           
     double precision                                          , allocatable, dimension(:,:) :: radiusEffectiveLogarithmicCovarianceTarget           , radiusEffectiveScatterCovarianceTarget
     double precision                                                                        :: randomErrorMassStellarMinimum                        , randomErrorMassStellarMaximum
     type            (varying_string                          )                              :: analysisLabel                                        , fileNameTarget
   contains
     final     ::                  sizeVsStellarMassRelationDestructor
     procedure :: analyze       => sizeVsStellarMassRelationAnalyze
     procedure :: finalize      => sizeVsStellarMassRelationFinalize
     procedure :: reduce        => sizeVsStellarMassRelationReduce
     procedure :: logLikelihood => sizeVsStellarMassRelationLogLikelihood
  end type outputAnalysisSizeVsStellarMassRelation

  interface outputAnalysisSizeVsStellarMassRelation
     !!{
     Constructors for the \refClass{outputAnalysisSizeVsStellarMassRelation} output analysis class.
     !!}
     module procedure sizeVsStellarMassRelationConstructorParameters
     module procedure sizeVsStellarMassRelationConstructorInternal
  end interface outputAnalysisSizeVsStellarMassRelation

contains

  function sizeVsStellarMassRelationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSizeVsStellarMassRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions    , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters   , only : cosmologyParametersClass
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    use :: Input_Parameters       , only : inputParameters
    implicit none
    type            (outputAnalysisSizeVsStellarMassRelation  )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    double precision                                           , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient       , systematicErrorMassStellarPolynomialCoefficient, &
         &                                                                                      randomErrorMassStellarPolynomialCoefficient
    class           (cosmologyParametersClass                 ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                  ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass              ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), pointer                     :: starFormationRateNuclearStarClusters_
    integer                                                                                  :: sample
    logical                                                                                  :: computeScatter                             , likelihoodNormalize                            , &
         &                                                                                      likelihoodBinsAutomatic
    integer         (c_size_t                                 ), allocatable  , dimension(:) :: likelihoodBins
    double precision                                                                         :: randomErrorMassStellarMinimum              , randomErrorMassStellarMaximum
    type            (varying_string                           )                              :: fileNameTarget                             , likelihoodBinsText

    ! Check and read parameters.
    allocate(           systematicErrorPolynomialCoefficient(max(1,parameters%count(           'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(systematicErrorMassStellarPolynomialCoefficient(max(1,parameters%count('systematicErrorMassStellarPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorMassStellarPolynomialCoefficient(max(1,parameters%count(    'randomErrorMassStellarPolynomialCoefficient',zeroIfNotPresent=.true.))))
    if (parameters%isPresent('likelihoodBins')) then
       call parameters%value('likelihoodBins',likelihoodBinsText)
       if (likelihoodBinsText == "auto") then
          likelihoodBinsAutomatic=.true.
          allocate(likelihoodBins(               0                  ))
       else
          allocate(likelihoodBins(parameters%count('likelihoodBins')))
          likelihoodBinsAutomatic=.false.
          !![
	  <inputParameter>
	    <name>likelihoodBins</name>
	    <source>parameters</source>
	    <description>
	      Controls which bins in the effective radius--stellar mass relation will be used in computing the likelihood:
	      \begin{itemize}
	      \item \emph{not present}: all bins are included in the likelihood calculation;
	      \item \emph{list of integers}: use only the mass bin(s) given in this list in the likelihood calculation;
	      \item {\normalfont \ttfamily auto}: use only bins which have a non-zero number of halos contributing to them in the likelihood calculation.
	      \end{itemize}
	    </description>
	  </inputParameter>
          !!]
        end if
    else
       allocate   (likelihoodBins(               0                  ))
       likelihoodBinsAutomatic   =.false.
    end if
    !![
    <inputParameter>
      <name>fileNameTarget</name>
      <source>parameters</source>
      <description>The name of the file containing the target data.</description>
    </inputParameter>
    <inputParameter>
      <name>sample</name>
      <source>parameters</source>
      <defaultValue>1</defaultValue>
      <description>The sample to use.</description>
    </inputParameter>
    <inputParameter>
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, then normalize the likelihood to make it a probability density.</description>
    </inputParameter>
    <inputParameter>
      <name>computeScatter</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the scatter in log10(radius effective) is computed. Otherwise, the mean is computed.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for effective radius in the effective radius vs stellar mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorMassStellarPolynomialCoefficient</name>
      <source>parameters</source>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for stellar mass in the effective radius vs stellar mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMassStellarPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorMassStellarPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the random error polynomial for stellar mass in the effective radius vs stellar mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMassStellarMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMassStellarMinimum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMassStellarMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMassStellarMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                  name="cosmologyParameters_"                  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=outputAnalysisSizeVsStellarMassRelation(fileNameTarget,sample,likelihoodBinsAutomatic,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassStellarPolynomialCoefficient,randomErrorMassStellarPolynomialCoefficient,randomErrorMassStellarMinimum,randomErrorMassStellarMaximum,cosmologyParameters_,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"                 />
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function sizeVsStellarMassRelationConstructorParameters

  function sizeVsStellarMassRelationConstructorInternal(fileNameTarget,sample,likelihoodBinsAutomatic,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassStellarPolynomialCoefficient,randomErrorMassStellarPolynomialCoefficient,randomErrorMassStellarMinimum,randomErrorMassStellarMaximum,cosmologyParameters_,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSizeVsStellarMassRelation} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                                       , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass                                      , cosmologyParametersSimple
    use :: Galactic_Filters                      , only : filterList                                                    , galacticFilterAll                              , galacticFilterStellarMass                     , galacticFilterRadiusEffective      , &
         &                                                galacticFilterAlways                                          , galacticFilterStarFormationRateNonParametric   , filterTypeLowPass                             , filterTypeHighPass                 , &
         &                                                enumerationFilterTypeType
    use :: Error                                 , only : Error_Report
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: ISO_Varying_String                    , only : var_str                                                       , varying_string
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassStellar                              , nodePropertyExtractorRadiusEffectiveStellar
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Numerical_Constants_Prefixes          , only : giga
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10                       , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, propertyOperatorList                          , outputAnalysisPropertyOperatorLog10, &
          &                                               outputAnalysisPropertyOperatorSequence                        , outputAnalysisPropertyOperatorSystmtcPolynomial, outputAnalysisPropertyOperatorCsmlgyAnglrDstnc
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: String_Handling                       , only : operator(//)
    implicit none
    type            (outputAnalysisSizeVsStellarMassRelation            )                                :: self
    type            (varying_string                                     ), intent(in   )                 :: fileNameTarget
    integer                                                              , intent(in   )                 :: sample
    logical                                                              , intent(in   )                 :: computeScatter                                                , likelihoodNormalize                                   , &
         &                                                                                                  likelihoodBinsAutomatic
    integer         (c_size_t                                           ), intent(in   ), dimension(:  ) :: likelihoodBins
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                          , systematicErrorMassStellarPolynomialCoefficient       , &
         &                                                                                                  randomErrorMassStellarPolynomialCoefficient
    double precision                                                     , intent(in   )                 :: randomErrorMassStellarMinimum                                 , randomErrorMassStellarMaximum
    class           (cosmologyParametersClass                           ), intent(in   ), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    class           (starFormationRateDisksClass                        ), intent(in   ), target         :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target         :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target         :: starFormationRateNuclearStarClusters_
    double precision                                                     , allocatable  , dimension(:  ) :: massStellar                                                   , radiusEffectiveLogarithmicTarget                      , &
         &                                                                                                  radiusEffectiveTarget                                         , radiusEffectiveScatterTarget                          , &
         &                                                                                                  radiusEffectiveErrorTarget                                    , radiusEffectiveScatterErrorTarget                     , &
         &                                                                                                  massStellarLogarithmic                                        , rateStarFormationMainSequence
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                                  , radiusEffectiveLogarithmicCovarianceTarget            , &
         &                                                                                                  radiusEffectiveScatterCovarianceTarget
    type            (galacticFilterStellarMass                          ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterRadiusEffective                      ), pointer                       :: galacticFilterRadiusEffective_
    class           (galacticFilterClass                                ), pointer                       :: galacticFilterSelection_
    type            (galacticFilterAll                                  ), pointer                       :: galacticFilterAll_
    type            (filterList                                         ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_                         , outputAnalysisWeightPropertyOperatorLog10_             , &
         &                                                                                                  outputAnalysisWeightPropertyOperatorLog10Second_             , outputAnalysisPropertyOperatorLog10Second_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_                            , outputAnalysisWeightPropertyOperatorAntiLog10_         , &
         &                                                                                                  outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisWeightPropertyOperator_                        , outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc     ), pointer                       :: outputAnalysisWeightPropertyOperatorCsmlgyAnglrDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisWeightPropertyOperatorSystmtcPolynomial_       , outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (nodePropertyExtractorMassStellar                   ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorRadiusEffectiveStellar        ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_                                           , propertyOperatorsMassStellar_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersTarget
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsTarget
    type            (surveyGeometryFullSky                              ), pointer                       :: surveyGeometry_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                              =10.0d0
    double precision                                                     , parameter                     :: covarianceLarge                                       = 1.0d4
    double precision                                                     , parameter                     :: fractionMassStellarLimit                              = 1.0d-2 ! A minimal fraction of the stellar mass of the lowest mass bin to consider (to avoid attempting to analyze galaxies with non-positive masses).
    double precision                                                     , parameter                     :: radiusEffectiveLimit                                  = 1.0d-6 ! A minimal effective radius to consider (to avoid attempting to analyze galaxies with non-positive radii).
    integer         (c_size_t                                           )                                :: iBin                                                          , massStellarCount
    double precision                                                                                     :: redshiftMinimum                                               , redshiftMaximum                                       , &
         &                                                                                                  offsetMainSequenceSFR                                         , HubbleConstantTarget                                  , &
         &                                                                                                  OmegaMatterTarget                                             , OmegaDarkEnergyTarget                                 , &
         &                                                                                                  OmegaBaryonTarget
    type            (varying_string                                     )                                :: analysisLabel                                                 , weightPropertyLabel                                   , &
         &                                                                                                  weightPropertyDescription                                     , groupSampleName                                       , &
         &                                                                                                  selection                                                     , referenceTarget                                       , &
         &                                                                                                  labelTarget
    type            (hdf5Object                                         )                                :: fileTarget                                                    , groupSample                                           , &
         &                                                                                                  groupCosmology
    character       (len=4                                              )                                :: redshiftMinimumLabel                                          , redshiftMaximumLabel
    type            (enumerationFilterTypeType                          )                                :: filterType
    !![
    <constructorAssign variables="fileNameTarget, sample, likelihoodBins, likelihoodNormalize, computeScatter, systematicErrorPolynomialCoefficient, systematicErrorMassStellarPolynomialCoefficient, randomErrorMassStellarPolynomialCoefficient, randomErrorMassStellarMinimum, randomErrorMassStellarMaximum, *cosmologyParameters_, *cosmologyFunctions_, *outputTimes_, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]

    ! Open the target data file and read basic information.
    !$ call hdf5Access%set()
    call fileTarget%openFile(self%fileNameTarget,readOnly=.true.)
    ! Find the requested sample.
    groupSampleName=var_str('sample')//sample
    if (.not.fileTarget%hasGroup(char(groupSampleName))) call Error_Report(var_str('redshift interval ')//sample//' is not present in `'//self%fileNameTarget//'`'//{introspection:location})
    groupSample=fileTarget%openGroup(char(groupSampleName))
    ! Read the redshift range and target data.
    call groupSample%readAttribute   ('redshiftMinimum'            ,redshiftMinimum                  )
    call groupSample%readAttribute   ('redshiftMaximum'            ,redshiftMaximum                  )
    call groupSample%readDataset     ('massStellar'                ,massStellar                      )
    call groupSample%readDataset     ('radiusEffective'            ,radiusEffectiveTarget            )
    call groupSample%readDataset     ('radiusEffectiveError'       ,radiusEffectiveErrorTarget       )
    call groupSample%readDataset     ('radiusEffectiveScatter'     ,radiusEffectiveScatterTarget     )
    call groupSample%readDataset     ('radiusEffectiveScatterError',radiusEffectiveScatterErrorTarget)
    call groupSample%readAttribute   ('selection'                  ,selection                        )
    if (selection == "star forming" .or. selection == "quiescent") then
       call groupSample%readDataset  ('mainSequenceSFR'            ,rateStarFormationMainSequence    )
       call groupSample%readAttribute('offsetMainSequenceSFR'      ,offsetMainSequenceSFR            )
    end if
    call groupSample%close           (                                                               )
    ! Get the cosmological parameters used in analyzing the target data.
    groupCosmology=fileTarget%openGroup('cosmology')
    call groupCosmology%readAttribute('OmegaMatter'    ,OmegaMatterTarget    )
    call groupCosmology%readAttribute('OmegaDarkEnergy',OmegaDarkEnergyTarget)
    call groupCosmology%readAttribute('OmegaBaryon'    ,OmegaBaryonTarget    )
    call groupCosmology%readAttribute('HubbleConstant' ,HubbleConstantTarget )
    call groupCosmology%close()
    ! Get the analysis label and target dataset reference.
    call fileTarget%readAttribute('label'    ,labelTarget    )
    call fileTarget%readAttribute('reference',referenceTarget)
    ! Close the target data file.
    call fileTarget%close()
    !$ call hdf5Access%unset()
    ! Construct survey geometry. A fully-sky geometry is used here as only the redshift range is important.
    write (redshiftMinimumLabel,'(f4.2)') redshiftMinimum
    write (redshiftMaximumLabel,'(f4.2)') redshiftMaximum
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create output time weights.
    massStellarCount=size(massStellar)
    allocate(outputWeight(massStellarCount,outputTimes_%count()))
    outputWeight(1,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,radiusEffectiveLimit,allowSingleEpoch=.true.)
    forall(iBin=2:massStellarCount)
       outputWeight(iBin,:)=outputWeight(1,:)
    end forall
    ! Convert masses to logarithmic form, and construct covariances.
    allocate(radiusEffectiveLogarithmicCovarianceTarget(massStellarCount,massStellarCount))
    allocate(radiusEffectiveScatterCovarianceTarget    (massStellarCount,massStellarCount))
    massStellarLogarithmic                    =log10(massStellar         )
    radiusEffectiveLogarithmicTarget          =log10(radiusEffectiveTarget)
    radiusEffectiveLogarithmicCovarianceTarget=0.0d0
    radiusEffectiveScatterCovarianceTarget    =0.0d0
    do iBin=1,massStellarCount
       radiusEffectiveLogarithmicCovarianceTarget(iBin,iBin)=+(                                          &
            &                                                  +radiusEffectiveErrorTarget        (iBin) &
            &                                                  /radiusEffectiveTarget             (iBin) &
            &                                                  /log(10.0d0)                              &
            &                                                 )**2
       radiusEffectiveScatterCovarianceTarget    (iBin,iBin)=+(                                          &
            &                                                   +radiusEffectiveScatterErrorTarget(iBin) &
            &                                                 )**2
    end do
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersTarget)
    allocate(cosmologyFunctionsTarget )
    !![
    <referenceConstruct object="cosmologyParametersTarget">
     <constructor>
      cosmologyParametersSimple     (                                                &amp;
         &amp;                       OmegaMatter         =OmegaMatterTarget        , &amp;
         &amp;                       OmegaDarkEnergy     =OmegaDarkEnergyTarget    , &amp;
         &amp;                       OmegaBaryon         =OmegaBaryonTarget        , &amp;
         &amp;                       HubbleConstant      =HubbleConstantTarget     , &amp;
         &amp;                       temperatureCMB      =2.72548d0                  &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsTarget">
     <constructor>
      cosmologyFunctionsMatterLambda(                                                &amp;
         &amp;                       cosmologyParameters_=cosmologyParametersTarget  &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! If only a subset of bins of the relation are to be populated adjust target data now.
    if (size(self%likelihoodBins) > 0) then
       ! Set the target dataset in all other bins to zero so that they do not contribute to the likelihood.
       if (any(self%likelihoodBins > massStellarCount .or. self%likelihoodBins < 1_c_size_t)) call Error_Report('likelihoodBins is out of range'//{introspection:location})
       do iBin=1,massStellarCount
          if (.not.any(self%likelihoodBins == iBin)) then
             radiusEffectiveLogarithmicTarget          (iBin     )=0.0d0
             radiusEffectiveScatterTarget              (iBin     )=0.0d0          
             radiusEffectiveLogarithmicCovarianceTarget(iBin,iBin)=covarianceLarge
             radiusEffectiveScatterCovarianceTarget    (iBin,iBin)=covarianceLarge
          end if
       end do
    end if
    allocate(self%radiusEffectiveLogarithmicTarget          ,source=radiusEffectiveLogarithmicTarget          )
    allocate(self%radiusEffectiveScatterTarget              ,source=radiusEffectiveScatterTarget              )
    allocate(self%radiusEffectiveLogarithmicCovarianceTarget,source=radiusEffectiveLogarithmicCovarianceTarget)
    allocate(self%radiusEffectiveScatterCovarianceTarget    ,source=radiusEffectiveScatterCovarianceTarget    )
    ! Build a filter which galaxies with stellar mass and effective radius above some coarse lower limit suitable for this sample,
    ! and with the specified selection criteria imposed.
    allocate(galacticFilterStellarMass_           )
    allocate(galacticFilterRadiusEffective_       )
    allocate(galacticFilterAll_                   )
    allocate(filters_                             )
    allocate(filters_                   %next     )
    allocate(filters_                   %next%next)
    !![
    <referenceConstruct object="galacticFilterStellarMass_"     constructor="galacticFilterStellarMass    (  massThreshold=fractionMassStellarLimit*massStellar(1))"/>
    <referenceConstruct object="galacticFilterRadiusEffective_" constructor="galacticFilterRadiusEffective(radiusThreshold=radiusEffectiveLimit                   )"/>
    !!]
    if      (                             &
         &    selection == 'none'         &
         &  ) then
       allocate(galacticFilterAlways :: galacticFilterSelection_)
       select type (galacticFilterSelection_)
       type is (galacticFilterAlways)
          !![
	  <referenceConstruct object="galacticFilterSelection_" constructor="galacticFilterAlways()"/>
          !!]
       end select
    else if (                             &
         &    selection == 'star forming' &
         &   .or.                         &
         &    selection == 'quiescent'    &
         &  ) then
       if (selection == 'star forming') then
          filterType=filterTypeHighPass
       else
          filterType=filterTypeLowPass
       end if
       allocate(galacticFilterStarFormationRateNonParametric :: galacticFilterSelection_)
       select type (galacticFilterSelection_)
       type is (galacticFilterStarFormationRateNonParametric)
          !![
	  <referenceConstruct object="galacticFilterSelection_" constructor="galacticFilterStarFormationRateNonParametric(filterType,massStellar,giga*10.0d0**(rateStarFormationMainSequence-offsetMainSequenceSFR),starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)"/>
          !!]
       end select
    else
       call Error_Report("unknown selection '"//selection//"'"//{introspection:location})
    end if
    filters_                             %filter_ => galacticFilterStellarMass_
    filters_                   %next     %filter_ => galacticFilterRadiusEffective_
    filters_                   %next%next%filter_ => galacticFilterSelection_
    !![
    <referenceConstruct object="galacticFilterAll_" constructor="galacticFilterAll(filters_)"/>
    !!]
    ! Build an identity distribution operator.
    allocate   (outputAnalysisDistributionOperator_                          )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
     outputAnalysisDistributionOperatorRandomErrorPlynml(                                             &amp;
        &amp;                                            randomErrorMassStellarMinimum              , &amp;
        &amp;                                            randomErrorMassStellarMaximum              , &amp;
        &amp;                                            errorPolynomialZeroPoint                   , &amp;
        &amp;                                            randomErrorMassStellarPolynomialCoefficient  &amp;
        &amp;                                           )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a sequence (log10(), polynomial systematic anti-log10, cosmological luminosity distance, log10) property operator.
    allocate   (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_              )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_     ,cosmologyFunctionsTarget                      ,outputTimes_)"/>
    !!]
    allocate   (outputAnalysisPropertyOperatorSystmtcPolynomial_             )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorMassStellarPolynomialCoefficient            )"/>
    !!]
    allocate   (outputAnalysisPropertyOperatorLog10_                         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                                    )"/>
    !!]
    allocate   (outputAnalysisPropertyOperatorLog10Second_                   )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10Second_"       constructor="outputAnalysisPropertyOperatorLog10            (                                                                                    )"/>
    !!]
    allocate   (outputAnalysisPropertyOperatorAntiLog10_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorAntiLog10_"         constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                    )"/>
    !!]
    allocate(propertyOperatorsMassStellar_                    )
    allocate(propertyOperatorsMassStellar_%next               )
    allocate(propertyOperatorsMassStellar_%next%next          )
    allocate(propertyOperatorsMassStellar_%next%next%next     )
    allocate(propertyOperatorsMassStellar_%next%next%next%next)
    propertyOperatorsMassStellar_                    %operator_ => outputAnalysisPropertyOperatorLog10_
    propertyOperatorsMassStellar_%next               %operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    propertyOperatorsMassStellar_%next%next          %operator_ => outputAnalysisPropertyOperatorAntiLog10_
    propertyOperatorsMassStellar_%next%next%next     %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperatorsMassStellar_%next%next%next%next%operator_ => outputAnalysisPropertyOperatorLog10Second_
    allocate(outputAnalysisPropertyOperator_                          )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                       constructor="outputAnalysisPropertyOperatorSequence          (propertyOperatorsMassStellar_                                               )"/>
    !!]
    ! Build a sequence (log10, polynomial systematic, anti-log10, cosmological angular distance, log10) of weight property operators.
    allocate   (outputAnalysisWeightPropertyOperatorCsmlgyAnglrDstnc_        )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorCsmlgyAnglrDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc  (cosmologyFunctions_     ,cosmologyFunctionsTarget              ,outputTimes_)"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorSystmtcPolynomial_       )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient               )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorLog10_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"            constructor="outputAnalysisPropertyOperatorLog10             (                                                                            )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorLog10Second_             )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10Second_"      constructor="outputAnalysisPropertyOperatorLog10             (                                                                            )"/>
    !!]
    allocate   (outputAnalysisWeightPropertyOperatorAntiLog10_               )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorAntiLog10_"        constructor="outputAnalysisPropertyOperatorAntiLog10         (                                                                            )"/>
    !!]
    allocate(propertyOperators_                              )
    allocate(propertyOperators_          %next               )
    allocate(propertyOperators_          %next%next          )
    allocate(propertyOperators_          %next%next%next     )
    allocate(propertyOperators_          %next%next%next%next)
    propertyOperators_                              %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperators_          %next               %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperators_          %next%next          %operator_ => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperators_          %next%next%next     %operator_ => outputAnalysisWeightPropertyOperatorCsmlgyAnglrDstnc_
    propertyOperators_          %next%next%next%next%operator_ => outputAnalysisWeightPropertyOperatorLog10Second_
    ! Create a stellar mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                          )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_" constructor="nodePropertyExtractorRadiusEffectiveStellar(                  )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"  constructor="outputAnalysisPropertyOperatorSequence     (propertyOperators_)"/>
    !!]
    ! Build weight operator.
    allocate   (outputAnalysisWeightOperator_                                )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"          constructor="outputAnalysisWeightOperatorIdentity       (                  )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"      constructor="outputAnalysisPropertyOperatorAntiLog10    (                  )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                                          )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                 constructor="nodePropertyExtractorMassStellar           (                  )"/>
    !!]
    ! Build the object.
    if (computeScatter) then
       analysisLabel            =var_str('stellarHaloMassRelationScatter' )//labelTarget//'Sample'//sample
       weightPropertyLabel      =var_str('radiusEffectiveLog10Scatter'    )
       weightPropertyDescription=var_str('σ_{log₁₀(Effective radius/Mpc)}')
       allocate(outputAnalysisScatterFunction1D :: self%outputAnalysis_)
    else
       analysisLabel            =var_str('stellarHaloMassRelation'        )//labelTarget//'Sample'//sample
       weightPropertyLabel      =var_str('radiusEffectiveLog10'           )
       weightPropertyDescription=var_str('⟨log₁₀(Effective radius/Mpc)⟩'  )
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
         outputAnalysisScatterFunction1D(                                                                                                                                                             &amp;
          &amp;                                                  analysisLabel                                                                                                                      , &amp;
          &amp;                                                  var_str('Scatter in effective radius-stellar mass relation for $')//redshiftMinimumLabel//' \le z &lt; '//redshiftMaximumLabel//'$', &amp;
          &amp;                                                  var_str('massStellar'                                            )                                                                 , &amp;
          &amp;                                                  var_str('Stellar mass'                                           )                                                                 , &amp;
          &amp;                                                  var_str('M☉'                                                     )                                                                 , &amp;
          &amp;                                                  massSolar                                                                                                                          , &amp;
          &amp;                                                  weightPropertyLabel                                                                                                                , &amp;
          &amp;                                                  weightPropertyDescription                                                                                                          , &amp;
          &amp;                                                  var_str(' '                                                      )                                                                 , &amp;
          &amp;                                                  0.0d0                                                                                                                              , &amp;
          &amp;                                                  massStellarLogarithmic                                                                                                             , &amp;
          &amp;                                                  0_c_size_t                                                                                                                         , &amp;
          &amp;                                                  outputWeight                                                                                                                       , &amp;
          &amp;                                                  nodePropertyExtractor_                                                                                                             , &amp;
          &amp;                                                  outputAnalysisWeightPropertyExtractor_                                                                                             , &amp;
          &amp;                                                  outputAnalysisPropertyOperator_                                                                                                    , &amp;
          &amp;                                                  outputAnalysisWeightPropertyOperator_                                                                                              , &amp;
          &amp;                                                  outputAnalysisPropertyUnoperator_                                                                                                  , &amp;
          &amp;                                                  outputAnalysisWeightOperator_                                                                                                      , &amp;
          &amp;                                                  outputAnalysisDistributionOperator_                                                                                                , &amp;
          &amp;                                                  galacticFilterAll_                                                                                                                 , &amp;
          &amp;                                                  outputTimes_                                                                                                                       , &amp;
          &amp;                                                  outputAnalysisCovarianceModelPoisson                                                                                               , &amp;
          &amp;                          likelihoodNormalize    =likelihoodNormalize                                                                                                                , &amp;
          &amp;                          xAxisLabel             =var_str('$M_\star/\mathrm{M}_\odot$'                       )                                                                       , &amp;
          &amp;                          yAxisLabel             =var_str('$\sigma_{\log_{10}(R_\mathrm{eff}/\mathrm{Mpc})}$')                                                                       , &amp;
          &amp;                          xAxisIsLog             =.true.                                                                                                                             , &amp;
          &amp;                          yAxisIsLog             =.false.                                                                                                                            , &amp;
          &amp;                          targetLabel            =referenceTarget                                                                                                                    , &amp;
          &amp;                          scatterValueTarget     =radiusEffectiveScatterTarget                                                                                                       , &amp;
          &amp;                          scatterCovarianceTarget=radiusEffectiveScatterCovarianceTarget                                                                                               &amp;
          &amp;                         )
        </constructor>
       </referenceConstruct>
       !!]
    type is (outputAnalysisMeanFunction1D   )
       !![
       <referenceConstruct isResult="yes" object="outputAnalysis_">
        <constructor>
         outputAnalysisMeanFunction1D   (                                                                                                                                               &amp;
          &amp;                                               analysisLabel                                                                                                           , &amp;
          &amp;                                               var_str('Effective radius-stellar mass relation for $')//redshiftMinimumLabel//' \le z &lt; '//redshiftMaximumLabel//'$', &amp;
          &amp;                                               var_str('massStellar'                                 )                                                                 , &amp;
          &amp;                                               var_str('Stellar mass'                                )                                                                 , &amp;
          &amp;                                               var_str('M☉'                                          )                                                                 , &amp;
          &amp;                                               massSolar                                                                                                               , &amp;
          &amp;                                               weightPropertyLabel                                                                                                     , &amp;
          &amp;                                               weightPropertyDescription                                                                                               , &amp;
          &amp;                                               var_str(' '                                           )                                                                 , &amp;
          &amp;                                               0.0d0                                                                                                                   , &amp;
          &amp;                                               massStellarLogarithmic                                                                                                  , &amp;
          &amp;                                               0_c_size_t                                                                                                              , &amp;
          &amp;                                               outputWeight                                                                                                            , &amp;
          &amp;                                               nodePropertyExtractor_                                                                                                  , &amp;
          &amp;                                               outputAnalysisWeightPropertyExtractor_                                                                                  , &amp;
          &amp;                                               outputAnalysisPropertyOperator_                                                                                         , &amp;
          &amp;                                               outputAnalysisWeightPropertyOperator_                                                                                   , &amp;
          &amp;                                               outputAnalysisPropertyUnoperator_                                                                                       , &amp;
          &amp;                                               outputAnalysisWeightOperator_                                                                                           , &amp;
          &amp;                                               outputAnalysisDistributionOperator_                                                                                     , &amp;
          &amp;                                               galacticFilterAll_                                                                                                      , &amp;
          &amp;                                               outputTimes_                                                                                                            , &amp;
          &amp;                                               outputAnalysisCovarianceModelPoisson                                                                                    , &amp;
          &amp;                          likelihoodNormalize =likelihoodNormalize                                                                                                     , &amp;
          &amp;                          xAxisLabel          =var_str('$M_\star\mathrm{M}_\odot$'               )                                                                     , &amp;
          &amp;                          yAxisLabel          =var_str('$\log_{10}(R_\mathrm{eff}/\mathrm{Mpc})$')                                                                     , &amp;
          &amp;                          xAxisIsLog          =.true.                                                                                                                  , &amp;
          &amp;                          yAxisIsLog          =.false.                                                                                                                 , &amp;
          &amp;                          targetLabel         =referenceTarget                                                                                                         , &amp;
          &amp;                          meanValueTarget     =radiusEffectiveLogarithmicTarget                                                                                        , &amp;
          &amp;                          meanCovarianceTarget=radiusEffectiveLogarithmicCovarianceTarget                                                                                &amp;
          &amp;                         )
        </constructor>
       </referenceConstruct>
       !!]
    end select
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                       />
    <objectDestructor name="cosmologyParametersTarget"                             />
    <objectDestructor name="cosmologyFunctionsTarget"                              />
    <objectDestructor name="galacticFilterStellarMass_"                            />
    <objectDestructor name="galacticFilterSelection_"                              />
    <objectDestructor name="galacticFilterRadiusEffective_"                        />
    <objectDestructor name="galacticFilterAll_"                                    />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10Second_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorAntiLog10_"              />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"      />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10Second_"      />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorAntiLog10_"        />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorCsmlgyAnglrDstnc_" />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="nodePropertyExtractor_"                                />
    !!]
    nullify(propertyOperatorsMassStellar_)
    nullify(propertyOperators_           )
    nullify(filters_                     )
    return
  end function sizeVsStellarMassRelationConstructorInternal

  subroutine sizeVsStellarMassRelationAnalyze(self,node,iOutput)
    !!{
    Implement a sizeVsStellarMassRelation output analysis.
    !!}
    implicit none
    class  (outputAnalysisSizeVsStellarMassRelation), intent(inout) :: self
    type   (treeNode                               ), intent(inout) :: node
    integer(c_size_t                               ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine sizeVsStellarMassRelationAnalyze

  subroutine sizeVsStellarMassRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSizeVsStellarMassRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisSizeVsStellarMassRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"                      />
    <objectDestructor name="self%cosmologyParameters_"                 />
    <objectDestructor name="self%cosmologyFunctions_"                  />
    <objectDestructor name="self%outputTimes_"                         />
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine sizeVsStellarMassRelationDestructor

  subroutine sizeVsStellarMassRelationReduce(self,reduced)
    !!{
    Implement reduction for the {\normalfont \ttfamily sizeVsStellarMassRelation} output analysis class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisSizeVsStellarMassRelation), intent(inout) :: self
    class(outputAnalysisClass                    ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisSizeVsStellarMassRelation)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine sizeVsStellarMassRelationReduce

  subroutine sizeVsStellarMassRelationFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily sizeVsStellarMassRelation} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisSizeVsStellarMassRelation), intent(inout)           :: self
    type (varying_string                         ), intent(in   ), optional :: groupName
    type (hdf5Object                             )               , target   :: analysesGroup, subGroup
    type (hdf5Object                             )               , pointer  :: inGroup
    type (hdf5Object                             )                          :: analysisGroup

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
  end subroutine sizeVsStellarMassRelationFinalize

  double precision function sizeVsStellarMassRelationLogLikelihood(self) result(logLikelihood)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily sizeVsStellarMassRelation} output analysis.
    !!}
    use :: Error                       , only : Error_Report
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisSizeVsStellarMassRelation), intent(inout)                 :: self
    double precision                                         , parameter                     :: radiusEffectiveLogarithmicTiny              =-6.0d0
    double precision                                         , allocatable  , dimension(:,:) :: radiusEffectiveLogarithmicCovarianceCombined       , radiusEffectiveLogarithmicCovarianceCombinedSelected, &
         &                                                                                      radiusEffectiveLogarithmicCovariance               , radiusEffectiveLogarithmicCovarianceTarget
    double precision                                         , allocatable  , dimension(:  ) :: radiusEffectiveLogarithmicDifference               , radiusEffectiveLogarithmicDifferenceSelected        , &
         &                                                                                      radiusEffectiveLogarithmic                         , radiusEffectiveLogarithmicTarget
    integer         (c_size_t                               ), allocatable  , dimension(:  ) :: likelihoodBins
    type            (vector                                 )                                :: residual
    type            (matrix                                 )                                :: covariance
    integer                                                                                  :: i                                                  , j                                                   , &
         &                                                                                      status

    select type (outputAnalysis_ => self%outputAnalysis_)
    class is (outputAnalysisMeanFunction1D   )
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(   meanValue=radiusEffectiveLogarithmic,   meanCovariance=radiusEffectiveLogarithmicCovariance)
       allocate(radiusEffectiveLogarithmicTarget          ,source=self%radiusEffectiveLogarithmicTarget          )
       allocate(radiusEffectiveLogarithmicCovarianceTarget,source=self%radiusEffectiveLogarithmicCovarianceTarget)
    class is (outputAnalysisScatterFunction1D)
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(scatterValue=radiusEffectiveLogarithmic,scatterCovariance=radiusEffectiveLogarithmicCovariance)
       allocate(radiusEffectiveLogarithmicTarget          ,source=self%radiusEffectiveScatterTarget              )
       allocate(radiusEffectiveLogarithmicCovarianceTarget,source=self%radiusEffectiveScatterCovarianceTarget    )
    class default
       logLikelihood=+outputAnalysis_%logLikelihood()
       return
    end select
    ! Determine which bins to use in the likelihood analysis.
    if (self%likelihoodBinsAutomatic) then
       j=0
       do i=1,size(radiusEffectiveLogarithmic)
          if (radiusEffectiveLogarithmic(i) /= 0.0d0) j=j+1
       end do
       allocate(likelihoodBins(j))
       j=0
       do i=1,size(radiusEffectiveLogarithmic)
          if (radiusEffectiveLogarithmic(i) /= 0.0d0) then
             j=j+1
             likelihoodBins(j)=i
          end if
       end do
    else
       allocate(likelihoodBins,source=self%likelihoodBins)
    end if
    if     (                                                                                                                     &
         &   (size(likelihoodBins) == 0 .and. any(radiusEffectiveLogarithmic                 <= radiusEffectiveLogarithmicTiny)) &
         &  .or.                                                                                                                 &
         &                                    any(radiusEffectiveLogarithmic(likelihoodBins) <= radiusEffectiveLogarithmicTiny)  &
         & ) then
       ! If any active bins contain zero galaxies, judge this model to be improbable.
       logLikelihood=                     logImprobable
    else
       ! Compute difference with the target dataset.
       allocate(radiusEffectiveLogarithmicDifference        ,mold=radiusEffectiveLogarithmic          )
       allocate(radiusEffectiveLogarithmicCovarianceCombined,mold=radiusEffectiveLogarithmicCovariance)
       radiusEffectiveLogarithmicDifference        =+radiusEffectiveLogarithmic          -radiusEffectiveLogarithmicTarget
       radiusEffectiveLogarithmicCovarianceCombined=+radiusEffectiveLogarithmicCovariance+radiusEffectiveLogarithmicCovarianceTarget
       ! Construct a reduced set of bins.
       if (size(likelihoodBins) > 0) then
          allocate(radiusEffectiveLogarithmicDifferenceSelected        (size(likelihoodBins)                     ))
          allocate(radiusEffectiveLogarithmicCovarianceCombinedSelected(size(likelihoodBins),size(likelihoodBins)))
          do i=1,size(likelihoodBins)
             radiusEffectiveLogarithmicDifferenceSelected           (i  )=radiusEffectiveLogarithmicDifference        (likelihoodBins(i)                  )
             do j=1,size(likelihoodBins)
                radiusEffectiveLogarithmicCovarianceCombinedSelected(i,j)=radiusEffectiveLogarithmicCovarianceCombined(likelihoodBins(i),likelihoodBins(j))
             end do
          end do
       else
          allocate(radiusEffectiveLogarithmicDifferenceSelected        ,source=radiusEffectiveLogarithmicDifference        )
          allocate(radiusEffectiveLogarithmicCovarianceCombinedSelected,source=radiusEffectiveLogarithmicCovarianceCombined)
       end if
       ! Construct residual vector and covariance matrix.
       residual  =vector(radiusEffectiveLogarithmicDifferenceSelected        )
       covariance=matrix(radiusEffectiveLogarithmicCovarianceCombinedSelected)
       ! Compute the log-likelihood.
       logLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
       if (status == GSL_Success) then
          if (self%likelihoodNormalize)                                                                      &
               & logLikelihood=+logLikelihood                                                                &
               &               -0.5d0*covariance%logarithmicDeterminant()                                    &
               &               -0.5d0*dble(size(radiusEffectiveLogarithmicDifferenceSelected))*log(2.0d0*Pi)
       else
          logLikelihood       =+logImprobable
       end if
    end if
    return
  end function sizeVsStellarMassRelationLogLikelihood
