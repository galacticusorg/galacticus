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
  Implements an output analysis class for the star forming main sequence measurements of \cite{wagner_evolution_2016}.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
   
  ! Enumerations of analyses.
  !![
  <enumeration>
   <name>wagner2016SSFRRedshiftRange</name>
   <description>Specifies the redshift range for the \cite{wagner_evolution_2016} analysis</description>
   <validator>yes</validator>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="low" />
   <entry label="high"/>
  </enumeration>
  <enumeration>
   <name>wagner2016SSFRGalaxyType</name>
   <description>Specifies the galaxy type for the \cite{wagner_evolution_2016} analysis</description>
   <validator>yes</validator>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="quiescent"  />
   <entry label="starForming"/>
  </enumeration>
  !!]
  
  !![
  <outputAnalysis name="outputAnalysisStarFormingMainSequenceWagner2016">
    <description>An output analysis class for the star forming main sequence measurements of \cite{wagner_evolution_2016}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisStarFormingMainSequence) :: outputAnalysisStarFormingMainSequenceWagner2016
     !!{
     An output analysis class for the star forming main sequence measurements of \cite{wagner_evolution_2016}.
     !!}
     private
     class           (cosmologyParametersClass                  ), pointer                     :: cosmologyParameters_             => null()
     class           (virialDensityContrastClass                ), pointer                     :: virialDensityContrast_           => null()
     class           (darkMatterProfileDMOClass                 ), pointer                     :: darkMatterProfileDMO_            => null()
     double precision                                            , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                                          :: randomErrorMinimum                        , randomErrorMaximum
     type            (enumerationWagner2016SSFRRedshiftRangeType)                              :: redshiftRange
     type            (enumerationWagner2016SSFRGalaxyTypeType   )                              :: galaxyType
   contains
     final :: starFormingMainSequenceWagner2016Destructor
  end type outputAnalysisStarFormingMainSequenceWagner2016

  interface outputAnalysisStarFormingMainSequenceWagner2016
     !!{
     Constructors for the \refClass{outputAnalysisStarFormingMainSequenceWagner2016} output analysis class.
     !!}
     module procedure starFormingMainSequenceWagner2016ConstructorParameters
     module procedure starFormingMainSequenceWagner2016ConstructorInternal
  end interface outputAnalysisStarFormingMainSequenceWagner2016

contains

  function starFormingMainSequenceWagner2016ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStarFormingMainSequenceWagner2016} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters   , only : cosmologyParameters       , cosmologyParametersClass
    use :: Cosmology_Functions    , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    use :: Input_Parameters       , only : inputParameter            , inputParameters
    implicit none
    type            (outputAnalysisStarFormingMainSequenceWagner2016)                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    class           (cosmologyParametersClass                       ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                        ), pointer                     :: cosmologyFunctions_
    class           (virialDensityContrastClass                     ), pointer                     :: virialDensityContrast_
    class           (darkMatterProfileDMOClass                      ), pointer                     :: darkMatterProfileDMO_
    class           (outputTimesClass                               ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass                    ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass      ), pointer                     :: starFormationRateNuclearStarClusters_
    double precision                                                 , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient, &
         &                                                                                            weightSystematicErrorPolynomialCoefficient
    double precision                                                                               :: randomErrorMinimum                        , randomErrorMaximum
    type            (varying_string                                 )                              :: redshiftRange                             , galaxyType
    
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient'      )) then
       allocate(          randomErrorPolynomialCoefficient(parameters%count(          'randomErrorPolynomialCoefficient')))
    else
       allocate(          randomErrorPolynomialCoefficient(1                                                   ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient'      )) then
       allocate(      systematicErrorPolynomialCoefficient(parameters%count(      'systematicErrorPolynomialCoefficient')))
    else
       allocate(      systematicErrorPolynomialCoefficient(1                                                   ))
    end if
    if (parameters%isPresent('weightSystematicErrorPolynomialCoefficient')) then
       allocate(weightSystematicErrorPolynomialCoefficient(parameters%count('weightSystematicErrorPolynomialCoefficient')))
    else
       allocate(weightSystematicErrorPolynomialCoefficient(1                                                   ))
    end if
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.07d0]</defaultValue>
      <description>The coefficients of the random error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftRange</name>
      <source>parameters</source>
      <description>The redshift range (``{\normalfont \ttfamily low}'' or ``{\normalfont \ttfamily high}'') for this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>galaxyType</name>
      <source>parameters</source>
      <description>The galaxy type (``{\normalfont \ttfamily quiescent}'' or ``{\normalfont \ttfamily starForming}'') for this analysis.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                  name="cosmologyParameters_"                  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"                 name="darkMatterProfileDMO_"                 source="parameters"/>
    <objectBuilder class="virialDensityContrast"                name="virialDensityContrast_"                source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=outputAnalysisStarFormingMainSequenceWagner2016(enumerationWagner2016SSFRRedshiftRangeEncode(char(redshiftRange),includesPrefix=.false.),enumerationWagner2016SSFRGalaxyTypeEncode(char(galaxyType),includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"                 />
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="virialDensityContrast_"               />
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="darkMatterProfileDMO_"                />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormingMainSequenceWagner2016ConstructorParameters

  function starFormingMainSequenceWagner2016ConstructorInternal(redshiftRange,galaxyType,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisStarFormingMainSequenceWagner2016} output analysis class.
    !!}
    use :: Error                                 , only : Error_Report
    use :: Cosmology_Functions                   , only : cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Input_Paths                           , only : inputPath                                          , pathTypeDataStatic
    use :: Output_Times                          , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors     , only : nbodyHaloMassErrorClass
    use :: String_Handling                       , only : operator(//)
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: Galactic_Filters                      , only : filterList                                         , galacticFilterAll            , galacticFilterClass       , galacticFilterHaloNotIsolated  , &
          &                                               galacticFilterHighPass                             , galacticFilterNot            , galacticFilterNull        , galacticFilterStarFormationRate, &
          &                                               galacticFilterStellarMass
    use :: Node_Property_Extractors              , only : nodePropertyExtractorHostNode                      , nodePropertyExtractorMassHalo
    use :: Virial_Density_Contrast               , only : fixedDensityTypeCritical                           , virialDensityContrastClass   , virialDensityContrastFixed
    implicit none
    type            (outputAnalysisStarFormingMainSequenceWagner2016    )                              :: self
    type            (enumerationWagner2016SSFRRedshiftRangeType         ), intent(in   )               :: redshiftRange
    type            (enumerationWagner2016SSFRGalaxyTypeType            ), intent(in   )               :: galaxyType
    double precision                                                     , intent(in   )               :: randomErrorMinimum                               , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                 , systematicErrorPolynomialCoefficient  , &
         &                                                                                                weightSystematicErrorPolynomialCoefficient
    class           (cosmologyParametersClass                           ), intent(inout), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                            ), intent(inout), target       :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target       :: outputTimes_
    class           (starFormationRateDisksClass                        ), intent(in   ), target       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target       :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target       :: starFormationRateNuclearStarClusters_
    class           (virialDensityContrastClass                         ), intent(in   ), target       :: virialDensityContrast_
    class           (darkMatterProfileDMOClass                          ), intent(inout), target       :: darkMatterProfileDMO_
    type            (galacticFilterHaloNotIsolated                      )               , pointer      :: galacticFilterIsSubhalo_
    type            (galacticFilterHighPass                             )               , pointer      :: galacticFilterHostHaloMass_
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilterStellarMass_
    type            (galacticFilterStarFormationRate                    )               , pointer      :: galacticFilterStarFormationRate_                 , galacticFilterStarFormationRateNonZero_
    type            (galacticFilterAll                                  )               , pointer      :: galacticFilter_
    class           (galacticFilterClass                                )               , pointer      :: galacticFilterGalaxyType_
    type            (virialDensityContrastFixed                         )               , pointer      :: virialDensityContrastDefinition_
    type            (nodePropertyExtractorHostNode                      )               , pointer      :: nodePropertyExtractorHost_
    type            (nodePropertyExtractorMassHalo                      )               , pointer      :: nodePropertyExtractorHostMass_
    type            (filterList                                         )               , pointer      :: filters_
    type            (surveyGeometryFullSky                              )               , pointer      :: surveyGeometry_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_                  , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperator_
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                  =11.0d0
    double precision                                                     , parameter                   :: errorPolynomialZeroPointWeight            = 0.0d0
    double precision                                                                                   :: redshiftMinimum                                  , redshiftMaximum                        , &
         &                                                                                                massHostThreshold
    type            (varying_string                                     )                              :: fileName                                         , label                                  , &
         &                                                                                                description
    !![
    <constructorAssign variables="redshiftRange, galaxyType, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, *cosmologyParameters_, *darkMatterProfileDMO_, *virialDensityContrast_"/>
    !!]

    ! Construct file name and label for the analysis.
    fileName   =inputPath(pathTypeDataStatic)//'observations/starFormationRate/'
    label      ='Wagner2016'
    description='Mean sSFR sequence from Wagner et al. (2016) for '
    !! Determine galaxy type.
    select case (galaxyType%ID)
    case (wagner2016SSFRGalaxyTypeQuiescent  %ID)
       fileName   =fileName   //'quiescent'
       label      =label      //'Quiescent'
       description=description//'quiescent'
    case (wagner2016SSFRGalaxyTypeStarForming%ID)
       fileName   =fileName   //'starForming'
       label      =label      //'StarForming'
       description=description//'star forming'
   case default
       call Error_Report('unrecognized galaxy type'   //{introspection:location})
    end select
    !! Add middle part of file name and description.
    fileName   =fileName   //'MainSequenceWagner2016_z'
    description=description//' galaxies with '
    !! Determine redshift range properties. Host halo mass threshold is judged approximately from Figure 1 of Wagner et al. (2016).
    select case (redshiftRange%ID)
    case (wagner2016SSFRRedshiftRangeLow %ID)
       redshiftMinimum  =0.15d0
       redshiftMaximum  =0.80d0
       massHostThreshold=10.0d0**14.8d0
       fileName         =fileName   // '0.15_0.80'
       label            =label      //'Z0.15_0.80'
       description      =description//'$0.15 < z < 0.80$'
    case (wagner2016SSFRRedshiftRangeHigh%ID)
       redshiftMinimum  =0.80d0
       redshiftMaximum  =1.50d0
       massHostThreshold=10.0d0**14.3d0
       fileName         =fileName   // '0.80_1.50'
       label            =label      //'Z0.80_1.50'
       description      =description//'$0.80 < z < 1.50$'
    case default
       call Error_Report('unrecognized redshift range'//{introspection:location})
    end select
    !! Add final part of file name.
    fileName=fileName//'.hdf5'
    ! Build a filter which selects satellite galaxies in hosts of the relevant mass, above a stellar mass threshold, and either quiescent or star-forming.
    allocate(galacticFilterIsSubhalo_)
    !![
    <referenceConstruct object="galacticFilterIsSubhalo_"                constructor="galacticFilterHaloNotIsolated  (                                                                                                                                                                                                                            )"/>
    !!]
    allocate(galacticFilterStellarMass_)
    !![
    <referenceConstruct object="galacticFilterStellarMass_"              constructor="galacticFilterStellarMass      (massThreshold=1.0d9                                                                                                                                                                                                         )"/>
    !!]
    allocate(galacticFilterStarFormationRateNonZero_)
    !![
    <referenceConstruct object="galacticFilterStarFormationRateNonZero_" constructor="galacticFilterStarFormationRate(logSFR0=-1.0d1,logSFR1=0.0d0,logM0=0.0d0,starFormationRateDisks_=starFormationRateDisks_,starFormationRateSpheroids_=starFormationRateSpheroids_,starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_)"/>
    !!]
    allocate(galacticFilterStarFormationRate_)
    !![
    <referenceConstruct object="galacticFilterStarFormationRate_"        constructor="galacticFilterStarFormationRate(logSFR0=-1.0d0,logSFR1=1.0d0,logM0=0.0d0,starFormationRateDisks_=starFormationRateDisks_,starFormationRateSpheroids_=starFormationRateSpheroids_,starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_)"/>
    !!]
    select case (galaxyType%ID)
    case (wagner2016SSFRGalaxyTypeQuiescent  %ID)
       allocate(galacticFilterNot  :: galacticFilterGalaxyType_)
       select type (galacticFilterGalaxyType_)
       type is (galacticFilterNot )
          !![
          <referenceConstruct object="galacticFilterGalaxyType_" constructor="galacticFilterNot (galacticFilterStarFormationRate_)" />
          !!]
       end select
    case (wagner2016SSFRGalaxyTypeStarForming%ID)
       allocate(galacticFilterNull :: galacticFilterGalaxyType_)
       select type (galacticFilterGalaxyType_)
       type is (galacticFilterNull)
          !![
          <referenceConstruct object="galacticFilterGalaxyType_" constructor="galacticFilterNull(galacticFilterStarFormationRate_)"/>
          !!]
       end select
    case default
       call Error_Report('unrecognized galaxy type'   //{introspection:location})
    end select    
    allocate(virialDensityContrastDefinition_)
    !![
    <referenceConstruct object="virialDensityContrastDefinition_" constructor="virialDensityContrastFixed     (densityContrastValue=200.0d0,densityType=fixedDensityTypeCritical,turnAroundOverVirialRadius=2.0d0,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    allocate(nodePropertyExtractorHostMass_)
    !![
    <referenceConstruct object="nodePropertyExtractorHostMass_"   constructor="nodePropertyExtractorMassHalo  (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                                                      )"/>
    !!]
    allocate(nodePropertyExtractorHost_)
    !![
    <referenceConstruct object="nodePropertyExtractorHost_"       constructor="nodePropertyExtractorHostNode  (nodePropertyExtractor_=nodePropertyExtractorHostMass_                                                                                                                               )"/>
    !!]
    allocate(galacticFilterHostHaloMass_)
    !![
    <referenceConstruct object="galacticFilterHostHaloMass_"      constructor="galacticFilterHighPass         (threshold=massHostThreshold,nodePropertyExtractor_=nodePropertyExtractorHost_                                                                                                       )"/>
    !!]
    allocate(galacticFilter_                    )
    allocate(filters_                           )
    allocate(filters_       %next               )
    allocate(filters_       %next%next          )
    allocate(filters_       %next%next%next     )
    allocate(filters_       %next%next%next%next)
    filters_                    %filter_ => galacticFilterIsSubhalo_
    filters_%next               %filter_ => galacticFilterStellarMass_
    filters_%next%next          %filter_ => galacticFilterStarFormationRateNonZero_
    filters_%next%next%next     %filter_ => galacticFilterGalaxyType_
    filters_%next%next%next%next%filter_ => galacticFilterHostHaloMass_
    !![
    <referenceConstruct object="galacticFilter_"                  constructor="galacticFilterAll              (filters_                                                                                                                                                                       )"/>
    !!]
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
      <constructor>
	cosmologyParametersSimple(                            &amp;
        &amp;                     OmegaMatter    = 0.27200d0, &amp;
        &amp;                     OmegaDarkEnergy= 0.72800d0, &amp;
        &amp;                     HubbleConstant =70.40000d0, &amp;
        &amp;                     temperatureCMB = 2.72548d0, &amp;
        &amp;                     OmegaBaryon    = 0.00000d0  &amp;
        &amp;                    )
      </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
      <constructor>
	cosmologyFunctionsMatterLambda(                        &amp;
        &amp;                          cosmologyParametersData &amp;
        &amp;                         )
      </constructor>
    </referenceConstruct>
    !!]
    ! Build the survey geometry.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_"                       constructor="surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"       constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient                )"/>
    !!]
    !! Systematic error model.
    allocate(outputAnalysisWeightPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,weightSystematicErrorPolynomialCoefficient          )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorRandomErrorPlynml(                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    self%outputAnalysisStarFormingMainSequence=                                         &
         & outputAnalysisStarFormingMainSequence(                                       &
         &                                       char(fileName)                       , &
         &                                       label                                , &
         &                                       description                          , &
         &                                       galacticFilter_                      , &
         &                                       surveyGeometry_                      , &
         &                                       cosmologyFunctions_                  , &
         &                                       cosmologyFunctionsData               , &
         &                                       outputTimes_                         , &
         &                                       outputAnalysisPropertyOperator_      , &
         &                                       outputAnalysisDistributionOperator_  , &
         &                                       outputAnalysisWeightPropertyOperator_, &
         &                                       starFormationRateDisks_              , &
         &                                       starFormationRateSpheroids_          , &
         &                                       starFormationRateNuclearStarClusters_  &
         &                                      )
    !![
    <objectDestructor name="galacticFilterIsSubhalo_"               />
    <objectDestructor name="galacticFilterStellarMass_"             />
    <objectDestructor name="galacticFilterStarFormationRate_"       />
    <objectDestructor name="galacticFilterStarFormationRateNonZero_"/>
    <objectDestructor name="galacticFilterGalaxyType_"              />
    <objectDestructor name="galacticFilter_"                        />
    <objectDestructor name="cosmologyParametersData"                />
    <objectDestructor name="cosmologyFunctionsData"                 />
    <objectDestructor name="surveyGeometry_"                        />
    <objectDestructor name="outputAnalysisPropertyOperator_"        />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"  />
    <objectDestructor name="outputAnalysisDistributionOperator_"    />
    !!]
    nullify(filters_)    
    return
  end function starFormingMainSequenceWagner2016ConstructorInternal

  subroutine starFormingMainSequenceWagner2016Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisStarFormingMainSequenceWagner2016} output analysis class.
    !!}
    implicit none
    type(outputAnalysisStarFormingMainSequenceWagner2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine starFormingMainSequenceWagner2016Destructor
