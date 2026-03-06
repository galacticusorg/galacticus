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
  Implements a stellar vs halo mass relation analysis class.
  !!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <outputAnalysis name="outputAnalysisStellarVsHaloMassRelation">
    <description>
      A stellar vs. halo mass relation output analysis class. Target data is read from an \gls{hdf5} file specified by the
      {\normalfont \ttfamily [fileNameTarget]} parameter. This file must contain one or more groups named {\normalfont \ttfamily
      redshiftIntervalN} where {\normalfont \ttfamily N} is an integer. Each such group specifies the stellar mass--halo mass
      relation in one redshift interval, and must contain the following datasets and attributes:
      \begin{itemize}
       \item dataset {\normalfont \ttfamily massHalo}: halo mass in units of $\mathrm{M}_\odot$;
       \item dataset {\normalfont \ttfamily massStellar}: stellar mass in units of $\mathrm{M}_\odot$;
       \item dataset {\normalfont \ttfamily massStellarError}: uncertainty in stellar mass in units of $\mathrm{M}_\odot$;
       \item dataset {\normalfont \ttfamily massStellarScatter}: scatter in stellar mass in units of dex;
       \item dataset {\normalfont \ttfamily massStellarScatterError}: uncertainty in scatter in stellar mass in units of dex;
       \item attribute {\normalfont \ttfamily redshiftMinimum}: the minimum redshift associated with this redshift interval;
       \item attribute {\normalfont \ttfamily redshiftMaximum}: the maximum redshift associated with this redshift interval.
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

      The file must have an attribute {\normalfont \ttfamily haloMassDefinition} which specifies the halo mass definition assumed
      in constructing the dataset. Allowed values are:
      \begin{itemize}
       \item {\normalfont \ttfamily `spherical collapse'} or {\normalfont \ttfamily `virial'}: halos are defined as have mean
             density contrasts given by spherical collapse calculations, e.g. \cite{percival_cosmological_2005};
       \item {\normalfont \ttfamily `Bryan \&amp; Norman (1998)'}: halos are defined as have mean density contrasts given by the
             fitting formula of \cite{bryan_statistical_1998};
       \item {\normalfont \ttfamily `X * mean density'}: halos are defined as having mean densities equal to ;{\normalfont
             \ttfamily X} times the mean density of the universe;
       \item {\normalfont \ttfamily `X * critical density'}: halos are defined as having mean densities equal to ;{\normalfont
             \ttfamily X} times the critical density of the universe;
      \end{itemize}
      Lastly, the file must have two attributes used to identify and level the dataset:
      \begin{itemize}
       \item {\normalfont \ttfamily label}: a space-free label that will be appended to the analysis group in the output, e.g. {\normalfont \ttfamily Leauthaud2012};
       \item {\normalfont \ttfamily reference}: a reference for the dataset suitable for inclusion in figures, e.g. {\normalfont \ttfamily Leauthaud et al. (2012)}.
      \end{itemize}
    </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisStellarVsHaloMassRelation
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
     logical                                                                   :: computeScatter                                  , likelihoodNormalize                         , &
          &                                                                       likelihoodBinsAutomatic
     integer         (c_size_t                  ), allocatable, dimension(:  ) :: likelihoodBins
     integer                                                                   :: redshiftInterval
     double precision                            , allocatable, dimension(:  ) :: systematicErrorPolynomialCoefficient            , systematicErrorMassHaloPolynomialCoefficient, &
          &                                                                       massStellarLogarithmicTarget                    , massStellarScatterTarget
     double precision                            , allocatable, dimension(:,:) :: massStellarLogarithmicCovarianceTarget          , massStellarScatterCovarianceTarget
     type            (varying_string            )                              :: analysisLabel                                   , fileNameTarget
   contains
     final     ::                  stellarVsHaloMassRelationDestructor
     procedure :: analyze       => stellarVsHaloMassRelationAnalyze
     procedure :: finalize      => stellarVsHaloMassRelationFinalize
     procedure :: reduce        => stellarVsHaloMassRelationReduce
     procedure :: logLikelihood => stellarVsHaloMassRelationLogLikelihood
  end type outputAnalysisStellarVsHaloMassRelation

  interface outputAnalysisStellarVsHaloMassRelation
     !!{
     Constructors for the \refClass{outputAnalysisStellarVsHaloMassRelation} output analysis class.
     !!}
     module procedure stellarVsHaloMassRelationConstructorParameters
     module procedure stellarVsHaloMassRelationConstructorInternal
  end interface outputAnalysisStellarVsHaloMassRelation

contains

  function stellarVsHaloMassRelationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStellarVsHaloMassRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions     , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    use :: Input_Parameters        , only : inputParameters
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelation)                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    double precision                                         , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, systematicErrorMassHaloPolynomialCoefficient
    class           (cosmologyParametersClass               ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass              ), pointer                     :: darkMatterProfileDMO_
    class           (virialDensityContrastClass             ), pointer                     :: virialDensityContrast_
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    integer                                                                                :: redshiftInterval
    logical                                                                                :: computeScatter                      , likelihoodNormalize                         , &
         &                                                                                    likelihoodBinsAutomatic
    integer         (c_size_t                               ), allocatable  , dimension(:) :: likelihoodBins
    type            (varying_string                         )                              :: fileNameTarget                      , likelihoodBinsText

    ! Check and read parameters.
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient        (parameters%count('systematicErrorPolynomialCoefficient'        )))
    else
       allocate(systematicErrorPolynomialCoefficient        (               1                                                ))
    end if
    if (parameters%isPresent('systematicErrorMassHaloPolynomialCoefficient')) then
       allocate(systematicErrorMassHaloPolynomialCoefficient(parameters%count('systematicErrorMassHaloPolynomialCoefficient')))
    else
       allocate(systematicErrorMassHaloPolynomialCoefficient(               1                                                ))
    end if
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
	      Controls which bins in the stellar mass--halo mass relation will be used in computing the likelihood:
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
      <name>redshiftInterval</name>
      <source>parameters</source>
      <defaultValue>1</defaultValue>
      <description>The redshift interval to use.</description>
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
      <description>If true, the scatter in log10(stellar mass) is computed. Otherwise, the mean is computed.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for stellar mass in the stellar vs halo mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorMassHaloPolynomialCoefficient</name>
      <source>parameters</source>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for halo mass in the stellar vs halo mass relation.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    !!]
    self=outputAnalysisStellarVsHaloMassRelation(fileNameTarget,redshiftInterval,likelihoodBinsAutomatic,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassHaloPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,outputTimes_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    <objectDestructor name="outputTimes_"          />
    !!]
    return
  end function stellarVsHaloMassRelationConstructorParameters

  function stellarVsHaloMassRelationConstructorInternal(fileNameTarget,redshiftInterval,likelihoodBinsAutomatic,likelihoodBins,likelihoodNormalize,computeScatter,systematicErrorPolynomialCoefficient,systematicErrorMassHaloPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStellarVsHaloMassRelation} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                                       , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass                                      , cosmologyParametersSimple
    use :: Galactic_Filters                      , only : filterList                                                    , galacticFilterAll                              , galacticFilterHaloIsolated, galacticFilterStellarMass
    use :: Error                                 , only : Error_Report
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: ISO_Varying_String                    , only : var_str                                                       , varying_string
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassHalo                                 , nodePropertyExtractorMassStellar
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10                       , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, propertyOperatorList      , outputAnalysisPropertyOperatorLog10, &
          &                                               outputAnalysisPropertyOperatorSequence                        , outputAnalysisPropertyOperatorSystmtcPolynomial, 
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: String_Handling                       , only : operator(//)
    use :: Virial_Density_Contrast               , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastClass                     , virialDensityContrastFixed, virialDensityContrastBryanNorman1998, &
         &                                                enumerationFixedDensityTypeType
    implicit none
    type            (outputAnalysisStellarVsHaloMassRelation        )                                :: self
    type            (varying_string                                 ), intent(in   )                 :: fileNameTarget
    integer                                                          , intent(in   )                 :: redshiftInterval
    logical                                                          , intent(in   )                 :: computeScatter                                                , likelihoodNormalize                                                , &
         &                                                                                              likelihoodBinsAutomatic
    integer         (c_size_t                                       ), intent(in   ), dimension(:  ) :: likelihoodBins
    double precision                                                 , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                          , systematicErrorMassHaloPolynomialCoefficient
    class           (cosmologyParametersClass                       ), intent(in   ), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                        ), intent(inout), target         :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                      ), intent(inout), target         :: darkMatterProfileDMO_
    class           (virialDensityContrastClass                     ), intent(in   ), target         :: virialDensityContrast_
    class           (outputTimesClass                               ), intent(inout), target         :: outputTimes_
    double precision                                                 , allocatable  , dimension(:  ) :: massHalo                                                      , massStellarLogarithmicTarget                                       , &
         &                                                                                              massStellarTarget                                             , massStellarScatterTarget                                           , &
         &                                                                                              massStellarErrorTarget                                        , massStellarScatterErrorTarget                                      , &
         &                                                                                              massHaloLogarithmic
    double precision                                                 , allocatable  , dimension(:,:) :: outputWeight                                                  , massStellarLogarithmicCovarianceTarget                             , &
         &                                                                                              massStellarScatterCovarianceTarget
    type            (galacticFilterStellarMass                      ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterHaloIsolated                     ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterAll                              ), pointer                       :: galacticFilterAll_
    type            (filterList                                     ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorIdentity     ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity           ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorLog10            ), pointer                       :: outputAnalysisPropertyOperatorLog10_                         , outputAnalysisWeightPropertyOperatorLog10_                          , &
         &                                                                                              outputAnalysisWeightPropertyOperatorLog10Second_
    type            (outputAnalysisPropertyOperatorAntiLog10        ), pointer                       :: outputAnalysisPropertyUnoperator_                            , outputAnalysisWeightPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorSequence         ), pointer                       :: outputAnalysisWeightPropertyOperator_                        , outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc), pointer                       :: outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial), pointer                       :: outputAnalysisWeightPropertyOperatorSystmtcPolynomial_       , outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (nodePropertyExtractorMassHalo                  ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellar               ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (propertyOperatorList                           ), pointer                       :: propertyOperators_                                           , propertyOperatorsMassHalo_
    type            (cosmologyParametersSimple                      ), pointer                       :: cosmologyParametersTarget
    type            (cosmologyFunctionsMatterLambda                 ), pointer                       :: cosmologyFunctionsTarget
    class           (virialDensityContrastClass                     ), pointer                       :: virialDensityContrastDefinition_
    type            (surveyGeometryFullSky                          ), pointer                       :: surveyGeometry_
    double precision                                                 , parameter                     :: errorPolynomialZeroPoint                              =11.3d0, errorPolynomialMassHaloZeroPoint                             =12.0d0
    double precision                                                 , parameter                     :: covarianceLarge                                       = 1.0d4
    double precision                                                 , parameter                     :: massStellarLimit                                      = 1.0d0 ! A minimal stellar mass to consider (to avoid attempting to analyze galaxies with non-positive masses).
    integer         (c_size_t                                       )                                :: iBin                                                          , massHaloCount
    double precision                                                                                 :: redshiftMinimum                                               , redshiftMaximum                                                   , &
         &                                                                                              densityValue                                                  , HubbleConstantTarget                                              , &
         &                                                                                              OmegaMatterTarget                                             , OmegaDarkEnergyTarget                                             , &
         &                                                                                              OmegaBaryonTarget
    type            (varying_string                                 )                                :: analysisLabel                                                 , weightPropertyLabel                                               , &
         &                                                                                              weightPropertyDescription                                     , groupRedshiftName                                                 , &
         &                                                                                              haloMassDefinition                                            , referenceTarget                                                   , &
         &                                                                                              labelTarget
    type            (hdf5Object                                     )                                :: fileTarget                                                    , groupRedshift                                                     , &
         &                                                                                              groupCosmology
    character       (len=4                                          )                                :: redshiftMinimumLabel                                          , redshiftMaximumLabel
    type(enumerationFixedDensityTypeType) :: densityType
    !![
    <constructorAssign variables="fileNameTarget, redshiftInterval, likelihoodBins, likelihoodBinsAutomatic, likelihoodNormalize, computeScatter, systematicErrorPolynomialCoefficient, systematicErrorMassHaloPolynomialCoefficient, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *outputTimes_"/>
    !!]

    ! Open the target data file and read basic information.
    !$ call hdf5Access%set()
    block
      fileTarget=hdf5Object(self%fileNameTarget,readOnly=.true.)
      ! Find the requested redshift interval.
      groupRedshiftName=var_str('redshiftInterval')//redshiftInterval
      if (.not.fileTarget%hasGroup(char(groupRedshiftName))) call Error_Report(var_str('redshift interval ')//redshiftInterval//' is not present in `'//self%fileNameTarget//'`'//{introspection:location})
      groupRedshift=fileTarget%openGroup(char(groupRedshiftName))
      ! Read the redshift range and target data.
      call groupRedshift%readAttribute('redshiftMinimum'        ,redshiftMinimum              )
      call groupRedshift%readAttribute('redshiftMaximum'        ,redshiftMaximum              )
      call groupRedshift%readDataset  ('massHalo'               ,massHalo                     )
      call groupRedshift%readDataset  ('massStellar'            ,massStellarTarget            )
      call groupRedshift%readDataset  ('massStellarError'       ,massStellarErrorTarget       )
      call groupRedshift%readDataset  ('massStellarScatter'     ,massStellarScatterTarget     )
      call groupRedshift%readDataset  ('massStellarScatterError',massStellarScatterErrorTarget)
      ! Get the cosmological parameters used in analyzing the target data.
      groupCosmology=fileTarget%openGroup('cosmology')
      call groupCosmology%readAttribute('OmegaMatter'    ,OmegaMatterTarget    )
      call groupCosmology%readAttribute('OmegaDarkEnergy',OmegaDarkEnergyTarget)
      call groupCosmology%readAttribute('OmegaBaryon'    ,OmegaBaryonTarget    )
      call groupCosmology%readAttribute('HubbleConstant' ,HubbleConstantTarget )
      ! Get the halo mass definition.
      call fileTarget%readAttribute('haloMassDefinition',haloMassDefinition)
      ! Get the analysis label and target dataset reference.
      call fileTarget%readAttribute('label'    ,labelTarget    )
      call fileTarget%readAttribute('reference',referenceTarget)
    end block
    !$ call hdf5Access%unset()
    ! Construct survey geometry. A fully-sky geometry is used here as only the redshift range is important.
    write (redshiftMinimumLabel,'(f4.2)') redshiftMinimum
    write (redshiftMaximumLabel,'(f4.2)') redshiftMaximum
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create output time weights.
    massHaloCount=size(massHalo)
    allocate(outputWeight(massHaloCount,outputTimes_%count()))
    outputWeight(1,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,massStellarLimit,allowSingleEpoch=.true.)
    forall(iBin=2:massHaloCount)
       outputWeight(iBin,:)=outputWeight(1,:)
    end forall
    ! Convert masses to logarithmic form, and construct covariances.
    allocate(massStellarLogarithmicCovarianceTarget(massHaloCount,massHaloCount))
    allocate(massStellarScatterCovarianceTarget    (massHaloCount,massHaloCount))
    massHaloLogarithmic                   =log10(massHalo         )
    massStellarLogarithmicTarget          =log10(massStellarTarget)
    massStellarLogarithmicCovarianceTarget=0.0d0
    massStellarScatterCovarianceTarget    =0.0d0
    do iBin=1,massHaloCount
       massStellarLogarithmicCovarianceTarget(iBin,iBin)=+(                                      &
            &                                              +massStellarErrorTarget        (iBin) &
            &                                              /massStellarTarget             (iBin) &
            &                                              /log(10.0d0)                          &
            &                                             )**2
       massStellarScatterCovarianceTarget    (iBin,iBin)=+(                                      &
            &                                               +massStellarScatterErrorTarget(iBin) &
            &                                             )**2
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
       if (any(self%likelihoodBins > massHaloCount .or. self%likelihoodBins < 1_c_size_t)) call Error_Report('likelihoodBins is out of range'//{introspection:location})
       do iBin=1,massHaloCount
          if (.not.any(self%likelihoodBins == iBin)) then
             massStellarLogarithmicTarget          (iBin     )=0.0d0
             massStellarScatterTarget              (iBin     )=0.0d0          
             massStellarLogarithmicCovarianceTarget(iBin,iBin)=covarianceLarge
             massStellarScatterCovarianceTarget    (iBin,iBin)=covarianceLarge
          end if
       end do
    end if
    allocate(self%massStellarLogarithmicTarget          ,source=massStellarLogarithmicTarget          )
    allocate(self%massStellarScatterTarget              ,source=massStellarScatterTarget              )
    allocate(self%massStellarLogarithmicCovarianceTarget,source=massStellarLogarithmicCovarianceTarget)
    allocate(self%massStellarScatterCovarianceTarget    ,source=massStellarScatterCovarianceTarget    )
   ! Build a filter which selects central galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilterStellarMass_      )
    allocate(galacticFilterHaloIsolated_     )
    allocate(galacticFilterAll_              )
    allocate(filters_                        )
    allocate(filters_                   %next)
    filters_                        %filter_ => galacticFilterHaloIsolated_
    filters_                   %next%filter_ => galacticFilterStellarMass_
    !![
    <referenceConstruct object="galacticFilterStellarMass_"  constructor="galacticFilterStellarMass  (massThreshold=massStellarLimit)"/>
    <referenceConstruct object="galacticFilterHaloIsolated_" constructor="galacticFilterHaloIsolated (                              )"/>
    <referenceConstruct object="galacticFilterAll_"          constructor="galacticFilterAll          (              filters_        )"/>
    !!]
    ! Build an identity distribution operator.
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
    ! Build a sequence (log10, polynomial systematic, anti-log10, cosmological luminosity distance) of weight property operators.
    allocate   (outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_       )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_"        constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc       (cosmologyFunctions_     ,cosmologyFunctionsTarget              ,outputTimes_                                                    )"/>
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
    allocate(propertyOperators_                              )
    allocate(propertyOperators_          %next               )
    allocate(propertyOperators_          %next%next          )
    allocate(propertyOperators_          %next%next%next     )
    allocate(propertyOperators_          %next%next%next%next)
    propertyOperators_                              %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperators_          %next               %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperators_          %next%next          %operator_ => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperators_          %next%next%next     %operator_ => outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_          %next%next%next%next%operator_ => outputAnalysisWeightPropertyOperatorLog10Second_
    ! Create a stellar mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                          )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                        constructor="nodePropertyExtractorMassStellar                              (                                                                                                                              )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                         constructor="outputAnalysisPropertyOperatorSequence                        (propertyOperators_                                                                                                            )"/>
    !!]
    ! Build weight operator.
    allocate   (outputAnalysisWeightOperator_                                )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                                 constructor="outputAnalysisWeightOperatorIdentity                          (                                                                                                                              )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                             constructor="outputAnalysisPropertyOperatorAntiLog10                       (                                                                                                                              )"/>
    !!]
    ! Create a halo mass property extractor.
    if (haloMassDefinition == "spherical collapse" .or. haloMassDefinition == "virial") then
       allocate(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt :: virialDensityContrastDefinition_                                )
       select type (virialDensityContrastDefinition_)
       type is (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)
          !![
	  <referenceConstruct object="virialDensityContrastDefinition_"                        constructor="virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(tableStore=.true.,cosmologyFunctions_=cosmologyFunctions_                                                                     )"/>
          !!]
       end select
    else if (haloMassDefinition == "Bryan & Norman (1998)") then
       allocate(virialDensityContrastBryanNorman1998 :: virialDensityContrastDefinition_                                )
       select type (virialDensityContrastDefinition_)
       type is (virialDensityContrastBryanNorman1998)
          !![
	  <referenceConstruct object="virialDensityContrastDefinition_"                        constructor="virialDensityContrastBryanNorman1998                          (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_           )"/>
          !!]
       end select
    else if (index(haloMassDefinition,"*") > 0) then
       densityValue=densityValueParse(char(extract(haloMassDefinition,1,index(haloMassDefinition,"*")-1)))
       densityType =densityTypeParse (char(extract(haloMassDefinition  ,index(haloMassDefinition,"*")+1)))
       allocate(virialDensityContrastFixed :: virialDensityContrastDefinition_                                )
       select type (virialDensityContrastDefinition_)
       type is (virialDensityContrastFixed)
          !![
	  <referenceConstruct object="virialDensityContrastDefinition_"                        constructor="virialDensityContrastFixed                                    (densityValue,densityType,2.0d0,cosmologyParameters_,cosmologyFunctions_                                                       )"/>
          !!]
       end select
    else
       call Error_Report("unknown halo mass definition '"//char(haloMassDefinition)//"'"//{introspection:location})
    end if
    allocate(nodePropertyExtractor_                                          )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                        constructor="nodePropertyExtractorMassHalo                                 (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    !!]
    ! Build the object.
    if (computeScatter) then
       analysisLabel            =var_str('stellarHaloMassRelationScatter')//labelTarget//'z'//redshiftInterval
       weightPropertyLabel      =var_str('massStellarLog10Scatter'   )
       weightPropertyDescription=var_str('σ_{log₁₀(Stellar mass/M☉)}')
       allocate(outputAnalysisScatterFunction1D :: self%outputAnalysis_)
    else
       analysisLabel            =var_str('stellarHaloMassRelation'       )//labelTarget//'z'//redshiftInterval
       weightPropertyLabel      =var_str('massStellarLog10'        )
       weightPropertyDescription=var_str('⟨log₁₀(Stellar mass/M☉)⟩')
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
          &amp;                                                  massHaloLogarithmic                                                                                                    , &amp;
          &amp;                                                  0_c_size_t                                                                                                             , &amp;
          &amp;                                                  outputWeight                                                                                                           , &amp;
          &amp;                                                  nodePropertyExtractor_                                                                                                 , &amp;
          &amp;                                                  outputAnalysisWeightPropertyExtractor_                                                                                 , &amp;
          &amp;                                                  outputAnalysisPropertyOperator_                                                                                        , &amp;
          &amp;                                                  outputAnalysisWeightPropertyOperator_                                                                                  , &amp;
          &amp;                                                  outputAnalysisPropertyUnoperator_                                                                                      , &amp;
          &amp;                                                  outputAnalysisWeightOperator_                                                                                          , &amp;
          &amp;                                                  outputAnalysisDistributionOperator_                                                                                    , &amp;
          &amp;                                                  galacticFilterAll_                                                                                                     , &amp;
          &amp;                                                  outputTimes_                                                                                                           , &amp;
          &amp;                                                  outputAnalysisCovarianceModelPoisson                                                                                   , &amp;
          &amp;                          likelihoodNormalize    =likelihoodNormalize                                                                                                    , &amp;
          &amp;                          xAxisLabel             =var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'            )                                                                , &amp;
          &amp;                          yAxisLabel             =var_str('$\sigma_{\log_{10}(M_\star/\mathrm{M}_\odot)}$')                                                                , &amp;
          &amp;                          xAxisIsLog             =.true.                                                                                                                 , &amp;
          &amp;                          yAxisIsLog             =.false.                                                                                                                , &amp;
          &amp;                          targetLabel            =referenceTarget                                                                                                        , &amp;
          &amp;                          scatterValueTarget     =massStellarScatterTarget                                                                                               , &amp;
          &amp;                          scatterCovarianceTarget=massStellarScatterCovarianceTarget                                                                                       &amp;
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
          &amp;                                               massHaloLogarithmic                                                                                         , &amp;
          &amp;                                               0_c_size_t                                                                                                  , &amp;
          &amp;                                               outputWeight                                                                                                , &amp;
          &amp;                                               nodePropertyExtractor_                                                                                      , &amp;
          &amp;                                               outputAnalysisWeightPropertyExtractor_                                                                      , &amp;
          &amp;                                               outputAnalysisPropertyOperator_                                                                             , &amp;
          &amp;                                               outputAnalysisWeightPropertyOperator_                                                                       , &amp;
          &amp;                                               outputAnalysisPropertyUnoperator_                                                                           , &amp;
          &amp;                                               outputAnalysisWeightOperator_                                                                               , &amp;
          &amp;                                               outputAnalysisDistributionOperator_                                                                         , &amp;
          &amp;                                               galacticFilterAll_                                                                                          , &amp;
          &amp;                                               outputTimes_                                                                                                , &amp;
          &amp;                                               outputAnalysisCovarianceModelPoisson                                                                        , &amp;
          &amp;                          likelihoodNormalize =likelihoodNormalize                                                                                         , &amp;
          &amp;                          xAxisLabel          =var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'   )                                                            , &amp;
          &amp;                          yAxisLabel          =var_str('$\log_{10}(M_\star/\mathrm{M}_\odot)$')                                                            , &amp;
          &amp;                          xAxisIsLog          =.true.                                                                                                      , &amp;
          &amp;                          yAxisIsLog          =.false.                                                                                                     , &amp;
          &amp;                          targetLabel         =referenceTarget                                                                                             , &amp;
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
    <objectDestructor name="cosmologyParametersTarget"                                    />
    <objectDestructor name="cosmologyFunctionsTarget"                                     />
    <objectDestructor name="galacticFilterStellarMass_"                                   />
    <objectDestructor name="galacticFilterHaloIsolated_"                                  />
    <objectDestructor name="galacticFilterAll_"                                           />
    <objectDestructor name="outputAnalysisDistributionOperator_"                          />
    <objectDestructor name="outputAnalysisWeightOperator_"                                />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                         />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"             />
    <objectDestructor name="outputAnalysisPropertyOperator_"                              />
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
    nullify(propertyOperatorsMassHalo_)
    nullify(propertyOperators_        )
    nullify(filters_                  )
    return

  contains
    
    double precision function densityValueParse(densityValueText) result(densityValue)
      !!{
      Parse and return the density multiplier value in a halo mass definition.
      !!}
      implicit none
      character(len=*), intent(in   ) :: densityValueText
      integer :: status
      
      read (densityValueText,*,iostat=status) densityValue
      if (status /= 0) call Error_Report("unable to parse density multiplier '"//trim(adjustl(densityValueText))//"' in halo mass definition"//{introspection:location})
      return
    end function densityValueParse
  
    function densityTypeParse(densityTypeText) result(densityType)
      !!{
      Parse and return the density type in a halo mass definition.
      !!}
      use :: Virial_Density_Contrast, only : fixedDensityTypeMean,fixedDensityTypeCritical
      implicit none
      type     (enumerationFixedDensityTypeType)                :: densityType
      character(len=*                          ), intent(in   ) :: densityTypeText

      select case (trim(adjustl(densityTypeText)))
      case ("mean density"    )
         densityType=fixedDensityTypeMean
      case ("critical density")
         densityType=fixedDensityTypeCritical
      case default
         densityType=fixedDensityTypeMean
         call Error_Report("unable to parse density type '"//trim(adjustl(densityTypeText))//"' in halo mass definition"//{introspection:location})
      end select
      return
    end function densityTypeParse
  
  end function stellarVsHaloMassRelationConstructorInternal

  subroutine stellarVsHaloMassRelationAnalyze(self,node,iOutput)
    !!{
    Implement a stellarVsHaloMassRelation output analysis.
    !!}
    implicit none
    class  (outputAnalysisStellarVsHaloMassRelation), intent(inout) :: self
    type   (treeNode                               ), intent(inout) :: node
    integer(c_size_t                               ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine stellarVsHaloMassRelationAnalyze

  subroutine stellarVsHaloMassRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisStellarVsHaloMassRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisStellarVsHaloMassRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"       />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%virialDensityContrast_"/>
    <objectDestructor name="self%outputTimes_"          />
    !!]
    return
  end subroutine stellarVsHaloMassRelationDestructor

  subroutine stellarVsHaloMassRelationReduce(self,reduced)
    !!{
    Implement reduction for the {\normalfont \ttfamily stellarVsHaloMassRelation} output analysis class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisStellarVsHaloMassRelation), intent(inout) :: self
    class(outputAnalysisClass                    ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisStellarVsHaloMassRelation)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine stellarVsHaloMassRelationReduce

  subroutine stellarVsHaloMassRelationFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily stellarVsHaloMassRelation} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisStellarVsHaloMassRelation), intent(inout)           :: self
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
    !$ call hdf5Access%unset()
    return
  end subroutine stellarVsHaloMassRelationFinalize

  double precision function stellarVsHaloMassRelationLogLikelihood(self) result(logLikelihood)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily stellarVsHaloMassRelation} output analysis.
    !!}
    use :: Error                       , only : Error_Report
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisStellarVsHaloMassRelation), intent(inout)                 :: self
    double precision                                         , parameter                     :: massStellarLogarithmicTiny              =1.0d-3
    double precision                                         , allocatable  , dimension(:,:) :: massStellarLogarithmicCovarianceCombined       , massStellarLogarithmicCovarianceCombinedSelected, &
         &                                                                                      massStellarLogarithmicCovariance               , massStellarLogarithmicCovarianceTarget
    double precision                                         , allocatable  , dimension(:  ) :: massStellarLogarithmicDifference               , massStellarLogarithmicDifferenceSelected        , &
         &                                                                                      massStellarLogarithmic                         , massStellarLogarithmicTarget
    integer         (c_size_t                               ), allocatable  , dimension(:  ) :: likelihoodBins
    type            (vector                                 )                                :: residual
    type            (matrix                                 )                                :: covariance
    integer                                                                                  :: i                                              , j                                               , &
         &                                                                                      status

    select type (outputAnalysis_ => self%outputAnalysis_)
    class is (outputAnalysisMeanFunction1D   )
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(   meanValue=massStellarLogarithmic,   meanCovariance=massStellarLogarithmicCovariance)
       allocate(massStellarLogarithmicTarget          ,source=self%massStellarLogarithmicTarget          )
       allocate(massStellarLogarithmicCovarianceTarget,source=self%massStellarLogarithmicCovarianceTarget)
    class is (outputAnalysisScatterFunction1D)
       ! Retrieve the results of the analysis.
       call outputAnalysis_%results(scatterValue=massStellarLogarithmic,scatterCovariance=massStellarLogarithmicCovariance)
       allocate(massStellarLogarithmicTarget          ,source=self%massStellarScatterTarget              )
       allocate(massStellarLogarithmicCovarianceTarget,source=self%massStellarScatterCovarianceTarget    )
    class default
       logLikelihood=+outputAnalysis_%logLikelihood()
       return
    end select
    ! Determine which bins to use in the likelihood analysis.
    if (self%likelihoodBinsAutomatic) then
       j=0
       do i=1,size(massStellarLogarithmic)
          if (massStellarLogarithmic(i) /= 0.0d0) j=j+1
       end do
       allocate(likelihoodBins(j))
       j=0
       do i=1,size(massStellarLogarithmic)
          if (massStellarLogarithmic(i) /= 0.0d0) then
             j=j+1
             likelihoodBins(j)=i
          end if
       end do
    else
       allocate(likelihoodBins,source=self%likelihoodBins)
    end if
    if     (                                                                                                             &
         &   (size(likelihoodBins) == 0 .and. any(massStellarLogarithmic                 <= massStellarLogarithmicTiny)) &
         &  .or.                                                                                                         &
         &                                    any(massStellarLogarithmic(likelihoodBins) <= massStellarLogarithmicTiny)  &
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
       if (size(likelihoodBins) > 0) then
          allocate(massStellarLogarithmicDifferenceSelected        (size(likelihoodBins)                     ))
          allocate(massStellarLogarithmicCovarianceCombinedSelected(size(likelihoodBins),size(likelihoodBins)))
          do i=1,size(likelihoodBins)
             massStellarLogarithmicDifferenceSelected           (i  )=massStellarLogarithmicDifference        (likelihoodBins(i)                  )
             do j=1,size(likelihoodBins)
                massStellarLogarithmicCovarianceCombinedSelected(i,j)=massStellarLogarithmicCovarianceCombined(likelihoodBins(i),likelihoodBins(j))
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
  end function stellarVsHaloMassRelationLogLikelihood
