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
  Implements an output analysis class for the star forming main sequence measurements of \cite{schreiber_herschel_2015}.
  !!}

  !![
  <outputAnalysis name="outputAnalysisStarFormingMainSequenceSchreiber2015">
    <description>An output analysis class for the star forming main sequence measurements of \cite{schreiber_herschel_2015}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisStarFormingMainSequence) :: outputAnalysisStarFormingMainSequenceSchreiber2015
     !!{
     An output analysis class for the star forming main sequence measurements of \cite{schreiber_herschel_2015}.
     !!}
     private
     class           (cosmologyParametersClass), pointer                     :: cosmologyParameters_                       => null()
     double precision                          , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient                    , systematicErrorPolynomialCoefficient, &
         &                                                                      weightSystematicErrorPolynomialCoefficient
     double precision                                                        :: randomErrorMinimum                                  , randomErrorMaximum
     integer                                                                 :: redshiftIndex
   contains
    final :: starFormingMainSequenceSchreiber2015Destructor
  end type outputAnalysisStarFormingMainSequenceSchreiber2015

  interface outputAnalysisStarFormingMainSequenceSchreiber2015
     !!{
     Constructors for the {\normalfont \ttfamily starFormingMainSequenceSchreiber2015} output analysis class.
     !!}
     module procedure starFormingMainSequenceSchreiber2015ConstructorParameters
     module procedure starFormingMainSequenceSchreiber2015ConstructorInternal
  end interface outputAnalysisStarFormingMainSequenceSchreiber2015

contains

  function starFormingMainSequenceSchreiber2015ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormingMainSequenceSchreiber2015} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParameters   , cosmologyParametersClass
    use :: Cosmology_Functions , only : cosmologyFunctions    , cosmologyFunctionsClass
    use :: Input_Parameters    , only : inputParameter        , inputParameters
    implicit none
    type            (outputAnalysisStarFormingMainSequenceSchreiber2015)                              :: self
    type            (inputParameters                                   ), intent(inout)               :: parameters
    class           (cosmologyParametersClass                          ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                  ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass                       ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                   ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass         ), pointer                     :: starFormationRateNuclearStarClusters_
    double precision                                                    , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient, &
         &                                                                                               weightSystematicErrorPolynomialCoefficient
    double precision                                                                                  :: randomErrorMinimum                        , randomErrorMaximum
    integer                                                                                           :: redshiftIndex
    
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
      <name>redshiftIndex</name>
      <source>parameters</source>
      <description>The redshift index (1-6) for this analysis.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                  name="cosmologyParameters_"                  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=outputAnalysisStarFormingMainSequenceSchreiber2015(redshiftIndex,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
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
  end function starFormingMainSequenceSchreiber2015ConstructorParameters

  function starFormingMainSequenceSchreiber2015ConstructorInternal(redshiftIndex,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,weightSystematicErrorPolynomialCoefficient,cosmologyParameters_,cosmologyFunctions_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormingMainSequenceSchreiber2015} output analysis class.
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
    use :: Galactic_Filters                      , only : filterList                                         , galacticFilterAll , galacticFilterStarFormationRate, galacticFilterStellarMass
    implicit none
    type            (outputAnalysisStarFormingMainSequenceSchreiber2015 )                              :: self
    integer                                                              , intent(in   )               :: redshiftIndex
    double precision                                                     , intent(in   )               :: randomErrorMinimum                               , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                 , systematicErrorPolynomialCoefficient  , &
         &                                                                                                weightSystematicErrorPolynomialCoefficient
    class           (cosmologyParametersClass                           ), intent(inout), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                            ), intent(inout), target       :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target       :: outputTimes_
    class           (starFormationRateDisksClass                        ), intent(in   ), target       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target       :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target       :: starFormationRateNuclearStarClusters_
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilterStellarMass_
    type            (galacticFilterStarFormationRate                    )               , pointer      :: galacticFilterStarFormationRate_
    type            (galacticFilterAll                                  )               , pointer      :: galacticFilter_
    type            (filterList                                         )               , pointer      :: filters_
    type            (surveyGeometryFullSky                              )               , pointer      :: surveyGeometry_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_                  , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperator_
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                  =11.0d0
    double precision                                                     , parameter                   :: errorPolynomialZeroPointWeight            = 0.0d0
    double precision                                                                                   :: redshiftMinimum                                  , redshiftMaximum                        , &
         &                                                                                                logSFR0
    type            (varying_string                                     )                              :: fileName                                         , label                                  , &
         &                                                                                                description
    !![
    <constructorAssign variables="redshiftIndex, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, weightSystematicErrorPolynomialCoefficient, *cosmologyParameters_"/>
    !!]
    
    ! Construct file name and label for the analysis.
    fileName   =inputPath(pathTypeDataStatic)//'observations/starFormationRate/starFormingMainSequenceSchreiber2015_z'
    label      ='Schreiber2015'
    description='Mean sSFR sequence from Schreiber et al. (2015) for galaxies with '
    !! Determine redshift range properties. Specific star formation rate thresholds wre judged approximately from Figure A1 of
    !! Schreiber et al. (2015).
    select case (redshiftIndex)
    case (1)
       logSFR0          =-1.15d0
       redshiftMinimum  =+0.30d0
       redshiftMaximum  =+0.70d0
       fileName         =fileName   // '0.3_0.7'
       label            =label      //'Z0.3_0.7'
       description      =description//'$0.3 < z < 0.7$'
    case (2)
       logSFR0          =-0.80d0
       redshiftMinimum  =+0.70d0
       redshiftMaximum  =+1.20d0
       fileName         =fileName   // '0.7_1.2'
       label            =label      //'Z0.7_1.2'
       description      =description//'$0.7 < z < 1.2$'
    case (3)
       logSFR0          =-0.55d0
       redshiftMinimum  =+1.20d0
       redshiftMaximum  =+1.80d0
       fileName         =fileName   // '1.2_1.8'
       label            =label      //'Z1.2_1.8'
       description      =description//'$1.2 < z < 1.8$'
    case (4)
       logSFR0          =-0.35d0
       redshiftMinimum  =+1.80d0
       redshiftMaximum  =+2.50d0
       fileName         =fileName   // '1.8_2.5'
       label            =label      //'Z1.8_2.6'
       description      =description//'$1.8 < z < 2.5$'
    case (5)
       logSFR0          =-0.15d0
       redshiftMinimum  =+2.50d0
       redshiftMaximum  =+3.50d0
       fileName         =fileName   // '2.5_3.5'
       label            =label      //'Z2.5_3.5'
       description      =description//'$2.5 < z < 3.5$'
    case (6)
       logSFR0          =+0.05d0
       redshiftMinimum  =+3.50d0
       redshiftMaximum  =+5.00d0
       fileName         =fileName   // '3.5_5.0'
       label            =label      //'Z3.5_5.0'
       description      =description//'$3.5 < z < 5.0$'
    case default
       call Error_Report('redshift index out of range'//{introspection:location})
    end select
    !! Add final part of file name.
    fileName=fileName//'.hdf5'
    ! Build a filter which selects galaxies above a stellar mass threshold, and star-forming.
    allocate(galacticFilterStellarMass_)
    !![
    <referenceConstruct object="galacticFilterStellarMass_"              constructor="galacticFilterStellarMass      (massThreshold=1.0d8                                                                                                                                                                                                         )"/>
    !!]
    allocate(galacticFilterStarFormationRate_)
    !![
    <referenceConstruct object="galacticFilterStarFormationRate_"        constructor="galacticFilterStarFormationRate(logSFR0=-1.0d0,logSFR1=1.0d0,logM0=0.0d0,starFormationRateDisks_=starFormationRateDisks_,starFormationRateSpheroids_=starFormationRateSpheroids_,starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_)"/>
    !!]
    allocate(galacticFilter_                    )
    allocate(filters_                           )
    allocate(filters_       %next               )
    filters_               %filter_ => galacticFilterStellarMass_
    filters_%next          %filter_ => galacticFilterStarFormationRate_
    !![
    <referenceConstruct object="galacticFilter_"                  constructor="galacticFilterAll              (filters_                                                                                                                                                                                                                           )"/>
    !!]
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
      <constructor>
	cosmologyParametersSimple(                            &amp;
        &amp;                     OmegaMatter    = 0.30000d0, &amp;
        &amp;                     OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                     HubbleConstant =70.00000d0, &amp;
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
    ! Build the survey geometry. A more elaborate model could be used here accounting for the different fields and depths used by Schreiber et al. (2015).
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
  end function starFormingMainSequenceSchreiber2015ConstructorInternal

  subroutine starFormingMainSequenceSchreiber2015Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormingMainSequenceSchreiber2015} output analysis class.
    !!}
    implicit none
    type(outputAnalysisStarFormingMainSequenceSchreiber2015), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    !!]
    return
  end subroutine starFormingMainSequenceSchreiber2015Destructor
