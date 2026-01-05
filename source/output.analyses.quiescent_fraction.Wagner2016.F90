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
  Implements an output analysis class for the quiescent fraction measurements of \cite{wagner_evolution_2016}.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  ! Enumerations of analyses.
  !![
  <enumeration>
   <name>wagner2016QuiescentRedshiftRange</name>
   <description>Specifies the redshift range for the \cite{wagner_evolution_2016} analysis.</description>
   <validator>yes</validator>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="low" />
   <entry label="mid" />
   <entry label="high"/>
  </enumeration>
  !!]
  
  !![
  <outputAnalysis name="outputAnalysisQuiescentFractionWagner2016">
    <description>An output analysis class for the quiescent fraction measurements of \cite{wagner_evolution_2016}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisQuiescentFraction) :: outputAnalysisQuiescentFractionWagner2016
     !!{
     An output analysis class for the quiescent fraction measurements of \cite{wagner_evolution_2016}.
     !!}
     private
     class           (cosmologyParametersClass                       ), pointer                     :: cosmologyParameters_                       => null()
     class           (virialDensityContrastClass                     ), pointer                     :: virialDensityContrast_                     => null()
     class           (darkMatterProfileDMOClass                      ), pointer                     :: darkMatterProfileDMO_                      => null()
     double precision                                                 , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient, &
         &                                                                                             weightSystematicErrorPolynomialCoefficient
     double precision                                                                               :: randomErrorMinimum                                  , randomErrorMaximum
     type            (enumerationWagner2016QuiescentRedshiftRangeType)                              :: redshiftRange
   contains
     final :: quiescentFractionWagner2016Destructor
  end type outputAnalysisQuiescentFractionWagner2016

  interface outputAnalysisQuiescentFractionWagner2016
     !!{
     Constructors for the \refClass{outputAnalysisQuiescentFractionWagner2016} output analysis class.
     !!}
     module procedure quiescentFractionWagner2016ConstructorParameters
     module procedure quiescentFractionWagner2016ConstructorInternal
  end interface outputAnalysisQuiescentFractionWagner2016

contains

  function quiescentFractionWagner2016ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisQuiescentFractionWagner2016} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters   , only : cosmologyParameters       , cosmologyParametersClass
    use :: Cosmology_Functions    , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    use :: Input_Parameters       , only : inputParameter            , inputParameters
    implicit none
    type            (outputAnalysisQuiescentFractionWagner2016)                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (cosmologyParametersClass                 ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                  ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass              ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), pointer                     :: starFormationRateNuclearStarClusters_
    class           (virialDensityContrastClass               ), pointer                     :: virialDensityContrast_
    class           (darkMatterProfileDMOClass                ), pointer                     :: darkMatterProfileDMO_
    double precision                                           , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient, &
         &                                                                                      weightSystematicErrorPolynomialCoefficient
    double precision                                                                         :: randomErrorMinimum                        , randomErrorMaximum
    type            (varying_string                           )                              :: redshiftRange
    
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
      <name>weightSystematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>weightSystematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for specific star formation rates.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftRange</name>
      <source>parameters</source>
      <description>The redshift range (``{\normalfont \ttfamily low}'' or ``{\normalfont \ttfamily high}'') for this analysis.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                  name="cosmologyParameters_"                  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="virialDensityContrast"                name="virialDensityContrast_"                source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"                 name="darkMatterProfileDMO_"                 source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=outputAnalysisQuiescentFractionWagner2016(enumerationWagner2016QuiescentRedshiftRangeEncode(char(redshiftRange),includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
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
  end function quiescentFractionWagner2016ConstructorParameters

  function quiescentFractionWagner2016ConstructorInternal(redshiftRange,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_, starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisQuiescentFractionWagner2016} output analysis class.
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
    use :: Galactic_Filters                      , only : filterList                                         , galacticFilterAll            , galacticFilterHaloNotIsolated, galacticFilterHighPass, &
          &                                               galacticFilterStarFormationRate                    , galacticFilterStellarMass
    use :: Node_Property_Extractors              , only : nodePropertyExtractorHostNode                      , nodePropertyExtractorMassHalo
    use :: Virial_Density_Contrast               , only : fixedDensityTypeCritical                           , virialDensityContrastClass   , virialDensityContrastFixed
    implicit none
    type            (outputAnalysisQuiescentFractionWagner2016          )                              :: self
    type            (enumerationWagner2016QuiescentRedshiftRangeType    ), intent(in   )               :: redshiftRange
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                    , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                      , systematicErrorPolynomialCoefficient, &
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
    type            (galacticFilterStarFormationRate                    )               , pointer      :: galacticFilterStarFormationRate_
    type            (galacticFilterAll                                  )               , pointer      :: galacticFilter_
    type            (virialDensityContrastFixed                         )               , pointer      :: virialDensityContrastDefinition_
    type            (nodePropertyExtractorHostNode                      )               , pointer      :: nodePropertyExtractorHost_
    type            (nodePropertyExtractorMassHalo                      )               , pointer      :: nodePropertyExtractorHostMass_
    type            (filterList                                         )               , pointer      :: filters_
    type            (surveyGeometryFullSky                              )               , pointer      :: surveyGeometry_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_                       , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperator_
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                     =+11.00d0
    double precision                                                     , parameter                   :: errorPolynomialZeroPointWeight               = +0.00d0
    double precision                                                     , parameter                   :: starFormationRateSpecificQuiescentLogarithmic= -1.00d0
    double precision                                                     , parameter                   :: starFormationRateSpecificLogarithmicError    = +0.15d0
    double precision                                                                                   :: redshiftMinimum                                       , redshiftMaximum                      , &
         &                                                                                                massHostThreshold
    type            (varying_string                                     )                              :: fileName                                              , label                                , &
         &                                                                                                description
    !![
    <constructorAssign variables="redshiftRange, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, *cosmologyParameters_, *darkMatterProfileDMO_, *virialDensityContrast_"/>
    !!]
    
    ! Construct file name and label for the analysis.
    fileName   =inputPath(pathTypeDataStatic)//'observations/starFormationRate/quiescentFractionWagner2016_z'
    label      ='Wagner2016'
    description='Quiescent fraction from Wagner et al. (2016) for galaxies with '
    !! Determine redshift range properties. Host halo mass threshold is judged approximately from Figure 1 of Wagner et al. (2016).
    select case (redshiftRange%ID)
    case (wagner2016QuiescentRedshiftRangeLow %ID)
       redshiftMinimum  =0.15d0
       redshiftMaximum  =0.41d0
       massHostThreshold=10.0d0**14.8d0
       fileName         =fileName   // '0.15_0.41'
       label            =label      //'Z0.15_0.41'
       description      =description//'$0.15 < z < 0.41$'
    case (wagner2016QuiescentRedshiftRangeMid %ID)
       redshiftMinimum  =0.41d0
       redshiftMaximum  =0.80d0
       massHostThreshold=10.0d0**14.8d0
       fileName         =fileName   // '0.41_0.80'
       label            =label      //'Z0.41_0.80'
       description      =description//'$0.41 < z < 0.80$'
    case (wagner2016QuiescentRedshiftRangeHigh%ID)
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
    <referenceConstruct object="galacticFilterIsSubhalo_"         constructor="galacticFilterHaloNotIsolated  (                                                                                                                                                                                                                            )"/>
    !!]
    allocate(galacticFilterStellarMass_)
    !![
    <referenceConstruct object="galacticFilterStellarMass_"       constructor="galacticFilterStellarMass      (massThreshold=1.0d9                                                                                                                                                                                                         )"/>
    !!]
    allocate(galacticFilterStarFormationRate_)
    !![
    <referenceConstruct object="galacticFilterStarFormationRate_" constructor="galacticFilterStarFormationRate(logSFR0=-1.0d0,logSFR1=1.0d0,logM0=0.0d0,starFormationRateDisks_=starFormationRateDisks_,starFormationRateSpheroids_=starFormationRateSpheroids_,starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_)"/>
    !!]
    allocate(virialDensityContrastDefinition_)
    !![
    <referenceConstruct object="virialDensityContrastDefinition_" constructor="virialDensityContrastFixed     (densityContrastValue=200.0d0,densityType=fixedDensityTypeCritical,turnAroundOverVirialRadius=2.0d0,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_                                        )"/>
    !!]
    allocate(nodePropertyExtractorHostMass_)
    !![
    <referenceConstruct object="nodePropertyExtractorHostMass_"   constructor="nodePropertyExtractorMassHalo  (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                                                                                              )"/>
    !!]
    allocate(nodePropertyExtractorHost_)
    !![
    <referenceConstruct object="nodePropertyExtractorHost_"       constructor="nodePropertyExtractorHostNode  (nodePropertyExtractor_=nodePropertyExtractorHostMass_                                                                                                                                                                       )"/>
    !!]
    allocate(galacticFilterHostHaloMass_)
    !![
    <referenceConstruct object="galacticFilterHostHaloMass_"      constructor="galacticFilterHighPass         (threshold=massHostThreshold,nodePropertyExtractor_=nodePropertyExtractorHost_                                                                                                                                               )"/>
    !!]
    allocate(galacticFilter_          )
    allocate(filters_                 )
    allocate(filters_       %next     )
    allocate(filters_       %next%next)
    filters_          %filter_ => galacticFilterIsSubhalo_
    filters_%next     %filter_ => galacticFilterStellarMass_
    filters_%next%next%filter_ => galacticFilterHostHaloMass_
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
    self%outputAnalysisQuiescentFraction=                                                 &
         & outputAnalysisQuiescentFraction(                                               &
         &                                 char(fileName)                               , &
         &                                 label                                        , &
         &                                 description                                  , &
         &                                 starFormationRateSpecificQuiescentLogarithmic, &
         &                                 starFormationRateSpecificLogarithmicError    , &
         &                                 galacticFilter_                              , &
         &                                 surveyGeometry_                              , &
         &                                 cosmologyFunctions_                          , &
         &                                 cosmologyFunctionsData                       , &
         &                                 outputTimes_                                 , &
         &                                 outputAnalysisPropertyOperator_              , &
         &                                 outputAnalysisDistributionOperator_          , &
         &                                 outputAnalysisWeightPropertyOperator_        , &
         &                                 starFormationRateDisks_                      , &
         &                                 starFormationRateSpheroids_                  , &
         &                                 starFormationRateNuclearStarClusters_          &
         &                                )
    !![
    <objectDestructor name="galacticFilterIsSubhalo_"             />
    <objectDestructor name="galacticFilterStellarMass_"           />
    <objectDestructor name="galacticFilterStarFormationRate_"     />
    <objectDestructor name="galacticFilter_"                      />
    <objectDestructor name="cosmologyParametersData"              />
    <objectDestructor name="cosmologyFunctionsData"               />
    <objectDestructor name="surveyGeometry_"                      />
    <objectDestructor name="outputAnalysisPropertyOperator_"      />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"/>
    <objectDestructor name="outputAnalysisDistributionOperator_"  />
    !!]
    nullify(filters_)    
    return
  end function quiescentFractionWagner2016ConstructorInternal

  subroutine quiescentFractionWagner2016Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisQuiescentFractionWagner2016} output analysis class.
    !!}
    implicit none
    type(outputAnalysisQuiescentFractionWagner2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine quiescentFractionWagner2016Destructor
