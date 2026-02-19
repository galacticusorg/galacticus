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

!+    Contributions to this file made by Sachi Weerasooriya

  !!{
  Implements an output analysis class for the \cite{robotham_galaxy_2011} star formation rate function.
  !!}


  !![
  <outputAnalysis name="outputAnalysisStarFormationRateFunctionRobotham2011">
   <description>A \cite{robotham_galaxy_2011} star formation rate function output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisStarFormationRateFunction) :: outputAnalysisStarFormationRateFunctionRobotham2011
     !!{
     A \cite{robotham_galaxy_2011} stellar mass function output analysis class.
     !!}
     private
     class           (gravitationalLensingClass), pointer                   :: gravitationalLensing_            => null()
     double precision                           , allocatable, dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                       :: randomErrorMinimum                        , randomErrorMaximum                  , &
          &                                                                    sizeSourceLensing
   contains
     final :: starFormationRateFunctionRobotham2011Destructor
  end type outputAnalysisStarFormationRateFunctionRobotham2011

  interface outputAnalysisStarFormationRateFunctionRobotham2011
     !!{
     Constructors for the \refClass{outputAnalysisStarFormationRateFunctionRobotham2011} output analysis class.
     !!}
     module procedure starFormationRateFunctionRobotham2011ConstructorParameters
     module procedure starFormationRateFunctionRobotham2011ConstructorInternal
  end interface outputAnalysisStarFormationRateFunctionRobotham2011

contains

  function starFormationRateFunctionRobotham2011ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStarFormationRateFunctionRobotham2011} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisStarFormationRateFunctionRobotham2011)                              :: self
    type            (inputParameters                                    ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                            ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                   ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass                        ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), pointer                     :: starFormationRateNuclearStarClusters_
    class           (gravitationalLensingClass                          ), pointer                     :: gravitationalLensing_
    double precision                                                     , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient     , systematicErrorPolynomialCoefficient
    integer                                                                                            :: covarianceBinomialBinsPerDecade
    double precision                                                                                   :: covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum   , &
         &                                                                                                randomErrorMinimum                   , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                       ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.126d0</defaultValue>
      <description>The minimum random error for \cite{robotham_galaxy_2011} star formation rates.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.126d0</defaultValue>
      <description>The minimum random error for \cite{robotham_galaxy_2011} star formation rates.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.126d0]</defaultValue>
      <description>The coefficients of the random error polynomial for \cite{robotham_galaxy_2011} star formation rates.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for \cite{robotham_galaxy_2011} star formation rates.</description>
    </inputParameter>
    <inputParameter>
      <name>sizeSourceLensing</name>
      <source>parameters</source>
      <variable>sizeSourceLensing</variable>
      <defaultValue>2.0d-3</defaultValue>
      <description>The characteristic source size for gravitational lensing calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing \cite{robotham_galaxy_2011} star formation rate function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing \cite{robotham_galaxy_2011} star formation rate function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing \cite{robotham_galaxy_2011} star formation rate function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"/>
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="gravitationalLensing"                 name="gravitationalLensing_"                 source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisStarFormationRateFunctionRobotham2011(cosmologyFunctions_,gravitationalLensing_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="gravitationalLensing_"                />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormationRateFunctionRobotham2011ConstructorParameters

  function starFormationRateFunctionRobotham2011ConstructorInternal(cosmologyFunctions_,gravitationalLensing_,outputTimes_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStarFormationRateFunctionRobotham2011} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                        , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass
    use :: Input_Paths                           , only : inputPath                                      , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryFullSky
    use :: Gravitational_Lensing                 , only : gravitationalLensingClass
    use :: Output_Analysis_Distribution_Operators, only : distributionOperatorList                       , outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorRandomErrorPlynml, outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    implicit none
    type            (outputAnalysisStarFormationRateFunctionRobotham2011)                              :: self
    class           (cosmologyFunctionsClass                            ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target       :: outputTimes_
    class           (gravitationalLensingClass                          ), intent(in   ), target       :: gravitationalLensing_
    class           (starFormationRateDisksClass                        ), intent(in   ), target       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                    ), intent(in   ), target       :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target       :: starFormationRateNuclearStarClusters_
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                          , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                            , systematicErrorPolynomialCoefficient
    integer                                                              , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )               :: covarianceBinomialMassHaloMinimum                           , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryFullSky                              )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng      )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence         )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                           )               , pointer      :: distributionOperatorSequence
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                            =11.3d+0
    !![
    <constructorAssign variables="randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, sizeSourceLensing, *gravitationalLensing_"/>
    !!]

    ! Build a filter which selects galaxies with stellar mass greater than 10^5.25.
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterStellarMass(massThreshold=10.0d0**5.25d0)"/>
    !!]
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple(                             &amp;
        &amp;                   OmegaMatter    =  0.25000d0, &amp;
        &amp;                   OmegaDarkEnergy=  0.75000d0, &amp;
        &amp;                   HubbleConstant =100.00000d0, &amp;
        &amp;                   temperatureCMB =  2.72548d0, &amp;
        &amp;                   OmegaBaryon    =  0.04550d0  &amp;
        &amp;                  )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
    <constructor>
      cosmologyFunctionsMatterLambda(                        &amp;
        &amp;                        cosmologyParametersData &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the survey geometry of GAMA with their imposed redshift limits.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_"                 constructor="surveyGeometryFullSky(redshiftMinimum=0.0d0,redshiftMaximum=0.5d0,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorRandomErrorPlynml_">
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
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorGrvtnlLnsng_">
     <constructor>
      outputAnalysisDistributionOperatorGrvtnlLnsng(                       &amp;
        &amp;                                       gravitationalLensing_, &amp;
        &amp;                                       outputTimes_         , &amp;
        &amp;                                       sizeSourceLensing      &amp;
        &amp;                                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence            )
    allocate(distributionOperatorSequence       %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence     %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence%next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorSequence(                             &amp;
        &amp;                                    distributionOperatorSequence &amp;
        &amp;                                   )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the object.
    self%outputAnalysisStarFormationRateFunction=                                                                                                                                  &
         & outputAnalysisStarFormationRateFunction(                                                                                                                                &
         &                                         var_str('Robotham2011'                                                                                                       ), &
         &                                         var_str('Star formation rate function for the Robotham et al. 2011'                                                          ), &
         &                                         char(inputPath(pathTypeDataStatic)//'/observations/starFormationRateFunctions/Star_Formation_Rate_Function_Robotham2011.hdf5'), &
         &                                         galacticFilter_                                                                                                               , &
         &                                         surveyGeometry_                                                                                                               , &
         &                                         cosmologyFunctions_                                                                                                           , &
         &                                         cosmologyFunctionsData                                                                                                        , &
         &                                         outputAnalysisPropertyOperator_                                                                                               , &
         &                                         outputAnalysisDistributionOperator_                                                                                           , &
         &                                         outputTimes_                                                                                                                  , &
         &                                         starFormationRateDisks_                                                                                                       , &
         &                                         starFormationRateSpheroids_                                                                                                   , &
         &                                         starFormationRateNuclearStarClusters_                                                                                         , &
         &                                         covarianceBinomialBinsPerDecade                                                                                               , &
         &                                         covarianceBinomialMassHaloMinimum                                                                                             , &
         &                                         covarianceBinomialMassHaloMaximum                                                                                               &
         &                                        )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                     />
    <objectDestructor name="galacticFilter_"                                     />
    <objectDestructor name="cosmologyParametersData"                             />
    <objectDestructor name="cosmologyFunctionsData"                              />
    <objectDestructor name="outputAnalysisPropertyOperator_"                     />
    <objectDestructor name="outputAnalysisDistributionOperator_"                 />
    <objectDestructor name="outputAnalysisDistributionOperatorGrvtnlLnsng_"      />
    <objectDestructor name="outputAnalysisDistributionOperatorRandomErrorPlynml_"/>
    !!]
    nullify(distributionOperatorSequence)
    return
  end function starFormationRateFunctionRobotham2011ConstructorInternal

  subroutine starFormationRateFunctionRobotham2011Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisStarFormationRateFunctionRobotham2011} output analysis class.
    !!}
    implicit none
    type(outputAnalysisStarFormationRateFunctionRobotham2011), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine starFormationRateFunctionRobotham2011Destructor
