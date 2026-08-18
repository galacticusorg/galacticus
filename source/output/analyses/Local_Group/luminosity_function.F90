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

!+ Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements an output analysis class that computes the :math:`V`-band luminosity function of Local Group satellite galaxies.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupLuminosityFunction" docformat="rst">
   <description>
   Computes the :math:`V`-band luminosity function of Local Group satellite galaxies, normalized to the number of host
   galaxies, for comparison with the sample of satellites identified by the DELVE Milky Way Census I
   :cite:p:`tan_delve_2025`.

   Absolute magnitude---not stellar mass---is the quantity in which that census, and the observational selection function
   which accompanies it, are defined, and for the faintest satellites the stellar masses reported in the Local Group database
   are themselves inferred from the absolute magnitude under an assumed mass-to-light ratio. This analysis therefore provides
   a more direct comparison than the corresponding
   :galacticus-class:`outputAnalysisLocalGroupStellarMassFunction` analysis.

   Model satellites are weighted by the probability that they would be detected in the census (see the
   :galacticus-class:`outputAnalysisWeightOperatorLocalGroupDetection` weight operator), rather than being subjected to a
   sharp detectability cut.
   </description>
   <deepCopy>
    <functionClass variables="volumeFunctionSatellites, volumeFunctionCentrals"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="volumeFunctionSatellites, volumeFunctionCentrals"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupLuminosityFunction
     !!{RST
     An output analysis class for the luminosity function of Local Group satellite galaxies.
     !!}
     private
     type            (outputAnalysisVolumeFunction1D), pointer                     :: volumeFunctionSatellites          => null(), volumeFunctionCentrals               => null()
     class           (outputTimesClass              ), pointer                     :: outputTimes_                      => null()
     double precision                                , allocatable, dimension(:  ) :: randomErrorPolynomialCoefficient           , systematicErrorPolynomialCoefficient
     double precision                                , allocatable, dimension(:  ) :: magnitudes                                 , luminosityFunction                            , &
          &                                                                           luminosityFunctionTarget
     double precision                                , allocatable, dimension(:,:) :: covariance
     logical                                                                       :: finalized
     integer                                                                       :: covarianceBinomialBinsPerDecade
     double precision                                                              :: covarianceBinomialMassHaloMinimum          , covarianceBinomialMassHaloMaximum             , &
          &                                                                           randomErrorMinimum                         , randomErrorMaximum                            , &
          &                                                                           negativeBinomialScatterFractional          , logLikelihoodZero                             , &
          &                                                                           countFailures                              , radiusHalfLightRatio
     type            (enumerationPositionTypeType   )                              :: positionType
   contains
     !![
     <methods docformat="rst">
       <method description="Finalize analysis." method="finalizeAnalysis" />
     </methods>
     !!]
     final     ::                     localGroupLuminosityFunctionDestructor
     procedure :: analyze          => localGroupLuminosityFunctionAnalyze
     procedure :: finalize         => localGroupLuminosityFunctionFinalize
     procedure :: finalizeAnalysis => localGroupLuminosityFunctionFinalizeAnalysis
     procedure :: reduce           => localGroupLuminosityFunctionReduce
     procedure :: logLikelihood    => localGroupLuminosityFunctionLogLikelihood
  end type outputAnalysisLocalGroupLuminosityFunction

  interface outputAnalysisLocalGroupLuminosityFunction
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisLocalGroupLuminosityFunction` output analysis class.
     !!}
     module procedure localGroupLuminosityFunctionConstructorParameters
     module procedure localGroupLuminosityFunctionConstructorInternal
  end interface outputAnalysisLocalGroupLuminosityFunction

contains

  function localGroupLuminosityFunctionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisLocalGroupLuminosityFunction` output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupLuminosityFunction)                              :: self
    type            (inputParameters                           ), intent(inout)               :: parameters
    class           (outputTimesClass                          ), pointer                     :: outputTimes_
    double precision                                            , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                   :: covarianceBinomialBinsPerDecade
    double precision                                                                          :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                       randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                       negativeBinomialScatterFractional, logLikelihoodZero
    double precision                                                                          :: radiusHalfLightRatio
    type            (varying_string                            )                              :: positionType

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
    <inputParameter docformat="rst">
      <name>negativeBinomialScatterFractional</name>
      <source>parameters</source>
      <defaultValue>0.18d0</defaultValue>
      <defaultSource>
      :cite:p:`boylan-kolchin_theres_2010`
      </defaultSource>
      <description>
      The fractional scatter (relative to the Poisson scatter) in the negative binomial distribution used in likelihood calculations.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.2d0</defaultValue>
      <description>
      The minimum random error for absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.2d0</defaultValue>
      <description>
      The maximum random error for absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.2d0]</defaultValue>
      <description>
      The coefficients of the random error polynomial for absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>
      The coefficients of the systematic error polynomial for absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>
      The number of bins per decade of halo mass to use when constructing Local Group luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The minimum halo mass to consider when constructing Local Group luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>
      The maximum halo mass to consider when constructing Local Group luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusHalfLightRatio</name>
      <source>parameters</source>
      <defaultValue>0.75d0</defaultValue>
      <defaultSource>
      Appropriate for a Plummer profile, for which the projected half-light radius is :math:`0.766` times the three-dimensional
      half-mass radius.
      </defaultSource>
      <description>
      The ratio of the projected, azimuthally-averaged half-light radius to the three-dimensional stellar half-mass radius, used
      when evaluating the observational selection function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('orbital')</defaultValue>
      <description>
      The type of position to use when determining the galactocentric radius of a satellite.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>logLikelihoodZero</name>
      <source>parameters</source>
      <defaultValue>logImprobable</defaultValue>
      <description>
      The log-likelihood to assign to bins where the model expectation is zero.
      </description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=outputAnalysisLocalGroupLuminosityFunction(outputTimes_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),radiusHalfLightRatio,negativeBinomialScatterFractional,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,logLikelihoodZero)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function localGroupLuminosityFunctionConstructorParameters

  function localGroupLuminosityFunctionConstructorInternal(outputTimes_,positionType,radiusHalfLightRatio,negativeBinomialScatterFractional,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,logLikelihoodZero) result (self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisLocalGroupLuminosityFunction` output analysis class for internal use.
    !!}
    use :: Display                                 , only : displayMessage                                     , verbosityLevelStandard
    use :: Galactic_Filters                        , only : filterList                                         , galacticFilterAll                             , galacticFilterHaloIsolated             , galacticFilterHaloNotIsolated                  , &
          &                                                 galacticFilterHostMassRange                        , enumerationPositionTypeType
    use :: Galactic_Structure_Options              , only : componentTypeAll
    use :: Input_Paths                             , only : inputPath                                          , pathTypeDataStatic
    use :: Interface_Local_Group_DB                , only : comparisonEquals                                   , comparisonLessThan                            , localGroupDB                           , setOperatorIntersection                        , &
          &                                                 setOperatorRelativeComplement                      , setOperatorUnion
    use :: Node_Property_Extractors                , only : nodePropertyExtractorLuminosityStellar
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Ranges                        , only : Make_Range                                         , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorIdentity             , outputAnalysisPropertyOperatorMagnitude        , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling            , outputAnalysisWeightOperatorLocalGroupDetection, outputAnalysisWeightOperatorSequence  , weightOperatorList
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupLuminosityFunction         )                                :: self
    integer                                                              , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )                 :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum                             , &
         &                                                                                                  negativeBinomialScatterFractional                          , logLikelihoodZero
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                         , randomErrorMaximum
    double precision                                                     , intent(in   )                 :: radiusHalfLightRatio
    double precision                                                     , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                        ), intent(in   )                 :: positionType
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    type            (nodePropertyExtractorLuminosityStellar             )               , pointer        :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorMagnitude            )               , pointer        :: outputAnalysisPropertyOperatorMagnitude_
    type            (outputAnalysisPropertyOperatorSequence             )               , pointer        :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorIdentity             )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling            )               , pointer        :: outputAnalysisWeightOperatorSubsampling_
    type            (outputAnalysisWeightOperatorLocalGroupDetection    )               , pointer        :: outputAnalysisWeightOperatorDetection_
    type            (outputAnalysisWeightOperatorSequence               )               , pointer        :: outputAnalysisWeightOperatorSatellites_
    type            (weightOperatorList                                 )               , pointer        :: weightOperators_
    type            (outputAnalysisDistributionNormalizerIdentity       )               , pointer        :: outputAnalysisDistributionNormalizerCentrals_              , outputAnalysisDistributionNormalizerSatellites_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer        :: outputAnalysisDistributionOperator_
    type            (galacticFilterHaloIsolated                         )               , pointer        :: galacticFilterHaloIsolated_
    type            (galacticFilterHaloNotIsolated                      )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                        )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterAll                                  )               , pointer        :: galacticFilterSatellites_                                  , galacticFilterCentrals_
    type            (filterList                                         )               , pointer        :: filtersSatellites_                                         , filtersCentrals_
    type            (propertyOperatorList                               )               , pointer        :: operators_
    double precision                                                     , allocatable  , dimension(:  ) :: magnitudesSatellites                                       , magnitudesCentrals                                            , &
         &                                                                                                  magnitudesTarget
    logical                                                              , allocatable  , dimension(:  ) :: isPresentTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeightSatellites                                      , outputWeightCentrals
    double precision                                                     , parameter                     :: bufferWidth                                    = +3.0d+0    , errorZeroPoint                                  = -5.0d0
    integer         (c_size_t                                           ), parameter                     :: binCountSatellites                             = 20_c_size_t, binCountCentrals                                =  2_c_size_t, &
         &                                                                                                  bufferCountMinimum                             =  5_c_size_t, bufferCountCentrals                             =  0_c_size_t
    double precision                                                     , parameter                     :: magnitudeSatelliteMinimum                      =-19.0d+0    , magnitudeSatelliteMaximum                       = +0.0d0     , &
         &                                                                                                  magnitudeCentralMinimum                        =-50.0d+0    , magnitudeCentralMaximum                         =+50.0d0     , &
         &                                                                                                  radiusOuter                                    = +3.0d-1
    integer         (c_size_t                                           )                                :: i                                                           , j                                                            , &
         &                                                                                                  bufferCountSatellites
    type            (localGroupDB                                       )                                :: localGroupDB_
    !![
    <constructorAssign variables="*outputTimes_, positionType, radiusHalfLightRatio, negativeBinomialScatterFractional, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, logLikelihoodZero"/>
    !!]

    ! Initialize.
    self%finalized=.false.
    ! Compute failure count for the negative binomial distribution used in likelihood calculations.
    self%countFailures=1.0d0/negativeBinomialScatterFractional**2
    ! Construct magnitude bins.
    allocate(magnitudesSatellites(binCountSatellites))
    allocate(magnitudesCentrals  (binCountCentrals  ))
    magnitudesSatellites=Make_Range(magnitudeSatelliteMinimum,magnitudeSatelliteMaximum,int(binCountSatellites),rangeTypeLinear)
    magnitudesCentrals  =Make_Range(magnitudeCentralMinimum  ,magnitudeCentralMaximum  ,int(binCountCentrals  ),rangeTypeLinear)
    ! Construct the target distribution.
    !! Select galaxies belonging to the DELVE Milky Way Census I, within 300 kpc of the Milky Way, and excluding the Milky Way
    !! itself, then retrieve absolute V-band magnitudes for the selected galaxies. See the corresponding analysis of the stellar
    !! mass function for a discussion of this selection.
    localGroupDB_=localGroupDB()
    call localGroupDB_%select     ('census'            ,var_str('DELVE-MW-I'),comparisonEquals  ,setOperatorUnion             ,propertyRequired=.false.)
    call localGroupDB_%select     ('classification'    ,var_str('galaxy'    ),comparisonEquals  ,setOperatorIntersection                               )
    call localGroupDB_%select     ('distanceMilkyWay'  ,         radiusOuter ,comparisonLessThan,setOperatorIntersection                               )
    call localGroupDB_%select     ('name'              ,var_str('The Galaxy'),comparisonEquals  ,setOperatorRelativeComplement                         )
    call localGroupDB_%getProperty('magnitudeAbsoluteV',magnitudesTarget                                 ,isPresentTarget                              )
    allocate(self%luminosityFunctionTarget(binCountSatellites))
    self%luminosityFunctionTarget=0.0d0
    do i=1,size(magnitudesTarget)
       if (.not.isPresentTarget(i)) then
          call displayMessage('localGroupLuminosityFunction: excluding a census galaxy which has no measured absolute magnitude',verbosityLevelStandard)
          cycle
       end if
       j=int((magnitudesTarget(i)-magnitudeSatelliteMinimum)/(magnitudesSatellites(2)-magnitudesSatellites(1))+0.5d0,kind=c_size_t)+1_c_size_t
       if (j > 0 .and. j <= binCountSatellites) self%luminosityFunctionTarget(j)=self%luminosityFunctionTarget(j)+1.0d0
    end do
    ! Create a stellar luminosity property extractor. The model provides magnitudes in the AB system, while the observed
    ! absolute magnitudes are on the Vega system. Galacticus defines its Vega system such that the offset between the two
    ! systems is zero in the V-band (by which the Vega system is normalized), and the true offset in V is in any case very much
    ! smaller than the uncertainties in these measurements, so no conversion is applied.
    allocate(nodePropertyExtractor_                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorLuminosityStellar        ('Buser_V','rest',componentTypeAll,outputTimes_     )"/>
    !!]
    ! Create property operators to convert luminosity to absolute magnitude, and to apply any systematic shift. Since the
    ! binned property is itself the absolute magnitude, the unoperator is the identity.
    allocate(outputAnalysisPropertyOperatorMagnitude_        )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorMagnitude_"         constructor="outputAnalysisPropertyOperatorMagnitude        (                                                   )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(operators_     )
    allocate(operators_%next)
    operators_     %operator_ => outputAnalysisPropertyOperatorMagnitude_
    operators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (operators_                                         )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorIdentity         (                                                   )"/>
    !!]
    ! Create a weight operator which corrects for subsampling of merger tree branches, and then applies the observational
    ! selection function of the census against which this analysis compares. Only the satellites are subject to the selection
    ! function; the central galaxies, which serve only to normalize the luminosity function, are not.
    allocate(outputAnalysisWeightOperatorSubsampling_        )
    !![
    <referenceConstruct object="outputAnalysisWeightOperatorSubsampling_"         constructor="outputAnalysisWeightOperatorSubsampling        (                                                   )"/>
    !!]
    allocate(outputAnalysisWeightOperatorDetection_          )
    !![
    <referenceConstruct object="outputAnalysisWeightOperatorDetection_">
     <constructor>
      outputAnalysisWeightOperatorLocalGroupDetection(                                                                                                &amp;
       &amp;                                          inputPath(pathTypeDataStatic)//'observations/localGroup/localGroupSelectionFunctionDELVE.hdf5', &amp;
       &amp;                                          var_str('Buser_V')                                                                            , &amp;
       &amp;                                          var_str('rest'   )                                                                            , &amp;
       &amp;                                          0.0d0                                                                                         , &amp;
       &amp;                                          radiusHalfLightRatio                                                                          , &amp;
       &amp;                                          radiusOuter                                                                                   , &amp;
       &amp;                                          positionType                                                                                    &amp;
       &amp;                                         )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(weightOperators_     )
    allocate(weightOperators_%next)
    weightOperators_     %operator_ => outputAnalysisWeightOperatorSubsampling_
    weightOperators_%next%operator_ => outputAnalysisWeightOperatorDetection_
    allocate(outputAnalysisWeightOperatorSatellites_         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperatorSatellites_"          constructor="outputAnalysisWeightOperatorSequence           (weightOperators_                                   )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
    <constructor>
    outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
         &amp;                                           randomErrorMinimum              , &amp;
         &amp;                                           randomErrorMaximum              , &amp;
         &amp;                                           errorZeroPoint                  , &amp;
         &amp;                                           randomErrorPolynomialCoefficient  &amp;
         &amp;                                          )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build filters which select satellites/centrals in a specified range of host halo mass. Detectability of satellites is
    ! imposed not by a filter, but by the weight operator constructed above, which weights each satellite by the probability
    ! that it would be detected in the census.
    allocate(galacticFilterHaloIsolated_   )
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"    constructor="galacticFilterHaloIsolated   ()"/>
    !!]
    allocate(galacticFilterHaloNotIsolated_)
    !![
    <referenceConstruct object="galacticFilterHaloNotIsolated_" constructor="galacticFilterHaloNotIsolated()"/>
    !!]
    allocate(galacticFilterHostMassRange_  )
    !![
    <referenceConstruct object="galacticFilterHostMassRange_">
     <constructor>
      galacticFilterHostMassRange(                      &amp;
        &amp;                     massMinimum =1.00d12, &amp;
        &amp;                     massMaximum =2.00d12, &amp;
        &amp;                     useFinalHost=.true.   &amp;
        &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(filtersCentrals_            )
    allocate(filtersCentrals_  %next     )
    filtersCentrals_            %filter_ => galacticFilterHaloIsolated_
    filtersCentrals_       %next%filter_ => galacticFilterHostMassRange_
    allocate(galacticFilterCentrals_  )
    !![
    <referenceConstruct object="galacticFilterCentrals_"   constructor="galacticFilterAll(filtersCentrals_  )"/>
    !!]
    allocate(filtersSatellites_          )
    allocate(filtersSatellites_%next     )
    filtersSatellites_          %filter_ => galacticFilterHaloNotIsolated_
    filtersSatellites_%next     %filter_ => galacticFilterHostMassRange_
    allocate(galacticFilterSatellites_)
    !![
    <referenceConstruct object="galacticFilterSatellites_" constructor="galacticFilterAll(filtersSatellites_)"/>
    !!]
    ! Build an identity distribution normalizers for centrals and satellites.
    allocate(outputAnalysisDistributionNormalizerCentrals_  )
    allocate(outputAnalysisDistributionNormalizerSatellites_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerCentrals_"   constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerSatellites_" constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeightSatellites(binCountSatellites,outputTimes_%count()))
    allocate(outputWeightCentrals  (binCountCentrals  ,outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(i),0.0d0,absTol=1.0d-6)) then
          outputWeightSatellites(:,i)=1.0d0
          outputWeightCentrals  (:,i)=1.0d0
       else
          outputWeightSatellites(:,i)=0.0d0
          outputWeightCentrals  (:,i)=0.0d0
       end if
    end do
    if (any(sum(outputWeightSatellites,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity for satellites'//{introspection:location})
    if (any(sum(outputWeightCentrals  ,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity for centrals'  //{introspection:location})
    ! Compute the number of buffer bins to add to either side of the luminosity function - these are needed to ensure that,
    ! e.g., convolution operations on the distribution function are unaffected by edge effects.
    bufferCountSatellites=max(int(bufferWidth/abs(magnitudesSatellites(2)-magnitudesSatellites(1)),kind=c_size_t)+1_c_size_t,bufferCountMinimum)
    ! Construct the volume function 1D objects.
    allocate(self%volumeFunctionSatellites)
    allocate(self%volumeFunctionCentrals  )
    !![
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionSatellites">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                              &amp;
       &amp;                         var_str('localGroupLuminosityFunction')                     , &amp;
       &amp;                         var_str('Luminosity function of Local Group satellites')    , &amp;
       &amp;                         var_str('magnitudeAbsoluteV')                               , &amp;
       &amp;                         var_str('Absolute V-band magnitude at the bin center')      , &amp;
       &amp;                         var_str(' ')                                                , &amp;
       &amp;                         var_str(' ')                                                , &amp;
       &amp;                         .false.                                                     , &amp;
       &amp;                         1.0d0                                                       , &amp;
       &amp;                         var_str('luminosityFunction')                               , &amp;
       &amp;                         var_str('Differential satellite V-band luminosity function'), &amp;
       &amp;                         var_str(' ')                                                , &amp;
       &amp;                         var_str(' ')                                                , &amp;
       &amp;                         .false.                                                     , &amp;
       &amp;                         1.0d0                                                       , &amp;
       &amp;                         magnitudesSatellites                                        , &amp;
       &amp;                         bufferCountSatellites                                       , &amp;
       &amp;                         outputWeightSatellites                                      , &amp;
       &amp;                         nodePropertyExtractor_                                      , &amp;
       &amp;                         outputAnalysisPropertyOperator_                             , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                           , &amp;
       &amp;                         outputAnalysisWeightOperatorSatellites_                     , &amp;
       &amp;                         outputAnalysisDistributionOperator_                         , &amp;
       &amp;                         outputAnalysisDistributionNormalizerSatellites_             , &amp;
       &amp;                         galacticFilterSatellites_                                   , &amp;
       &amp;                         outputTimes_                                                , &amp;
       &amp;                         outputAnalysisCovarianceModelBinomial                       , &amp;
       &amp;                         covarianceBinomialBinsPerDecade                             , &amp;
       &amp;                         covarianceBinomialMassHaloMinimum                           , &amp;
       &amp;                         covarianceBinomialMassHaloMaximum                             &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionCentrals">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                             &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         .false.                                                    , &amp;
       &amp;                         1.0d0                                                      , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         var_str(' ')                                               , &amp;
       &amp;                         .false.                                                    , &amp;
       &amp;                         1.0d0                                                      , &amp;
       &amp;                         magnitudesCentrals                                         , &amp;
       &amp;                         bufferCountCentrals                                        , &amp;
       &amp;                         outputWeightCentrals                                       , &amp;
       &amp;                         nodePropertyExtractor_                                     , &amp;
       &amp;                         outputAnalysisPropertyOperator_                            , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                          , &amp;
       &amp;                         outputAnalysisWeightOperatorSubsampling_                   , &amp;
       &amp;                         outputAnalysisDistributionOperator_                        , &amp;
       &amp;                         outputAnalysisDistributionNormalizerCentrals_              , &amp;
       &amp;                         galacticFilterCentrals_                                    , &amp;
       &amp;                         outputTimes_                                               , &amp;
       &amp;                         outputAnalysisCovarianceModelBinomial                      , &amp;
       &amp;                         covarianceBinomialBinsPerDecade                            , &amp;
       &amp;                         covarianceBinomialMassHaloMinimum                          , &amp;
       &amp;                         covarianceBinomialMassHaloMaximum                            &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperatorMagnitude_"        />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="outputAnalysisWeightOperatorSatellites_"         />
    <objectDestructor name="outputAnalysisWeightOperatorSubsampling_"        />
    <objectDestructor name="outputAnalysisWeightOperatorDetection_"          />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="galacticFilterHaloIsolated_"                     />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                  />
    <objectDestructor name="galacticFilterHostMassRange_"                    />
    <objectDestructor name="galacticFilterCentrals_"                         />
    <objectDestructor name="galacticFilterSatellites_"                       />
    <objectDestructor name="outputAnalysisDistributionNormalizerSatellites_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerCentrals_"   />
    !!]
    nullify(weightOperators_  )
    nullify(filtersSatellites_)
    nullify(filtersCentrals_  )
    nullify(operators_        )
    return
  end function localGroupLuminosityFunctionConstructorInternal

  subroutine localGroupLuminosityFunctionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`outputAnalysisLocalGroupLuminosityFunction` output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupLuminosityFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%volumeFunctionSatellites"/>
    <objectDestructor name="self%volumeFunctionCentrals"  />
    <objectDestructor name="self%outputTimes_"            />
    !!]
    return
  end subroutine localGroupLuminosityFunctionDestructor

  subroutine localGroupLuminosityFunctionAnalyze(self,node,iOutput)
    !!{RST
    Implement a ``localGroupLuminosityFunction`` output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupLuminosityFunction), intent(inout) :: self
    type   (treeNode                                  ), intent(inout) :: node
    integer(c_size_t                                  ), intent(in   ) :: iOutput

    call self%volumeFunctionSatellites%analyze(node,iOutput)
    call self%volumeFunctionCentrals  %analyze(node,iOutput)
    return
  end subroutine localGroupLuminosityFunctionAnalyze

  subroutine localGroupLuminosityFunctionReduce(self,reduced)
    !!{RST
    Implement a ``localGroupLuminosityFunction`` output analysis reduction.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    class(outputAnalysisLocalGroupLuminosityFunction), intent(inout) :: self
    class(outputAnalysisClass                       ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupLuminosityFunction)
       call self%volumeFunctionSatellites%reduce(reduced%volumeFunctionSatellites)
       call self%volumeFunctionCentrals  %reduce(reduced%volumeFunctionCentrals  )
    class default
       call Error_Report('object is not of [outputAnalysisLocalGroupLuminosityFunction] class, but of ['//char(reduced%objectType())//'] class'//{introspection:location})
    end select
    return
  end subroutine localGroupLuminosityFunctionReduce

  subroutine localGroupLuminosityFunctionFinalizeAnalysis(self)
    !!{RST
    Finalize analysis of a ``localGroupLuminosityFunction`` output analysis.
    !!}
    implicit none
    class           (outputAnalysisLocalGroupLuminosityFunction), intent(inout)                 :: self
    double precision                                            , allocatable  , dimension(:  ) :: luminosityFunctionCentrals
    double precision                                            , allocatable  , dimension(:,:) :: covarianceCentrals
    double precision                                                                            :: weight                    , weightVariance
    integer         (c_size_t                                  )                                :: i                         , j

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! Retrieve results from our 1-D volume functions.
    call self%volumeFunctionSatellites%results(binCenter=self%magnitudes,functionValue=self%luminosityFunction  ,functionCovariance=self%covariance        )
    call self%volumeFunctionCentrals  %results(                          functionValue=luminosityFunctionCentrals,functionCovariance=    covarianceCentrals)
    ! Normalize the luminosity function to the number of host galaxies.
    weight           = sum(luminosityFunctionCentrals)
    weightVariance   = sum(covarianceCentrals        )
    if (weight > 0.0d0) then
       self%luminosityFunction=+self%luminosityFunction/weight
       self%covariance        =+self%covariance        /weight**2
       do    i=1_c_size_t,size(self%luminosityFunction)
          do j=1_c_size_t,size(self%luminosityFunction)
             self%covariance(i,j)=+self%covariance        (i,j)    &
                  &               +self%luminosityFunction(i  )    &
                  &               *self%luminosityFunction(  j)    &
                  &               *weightVariance                  &
                  &               /weight                      **2
          end do
       end do
    end if
    return
  end subroutine localGroupLuminosityFunctionFinalizeAnalysis

  subroutine localGroupLuminosityFunctionFinalize(self,groupName)
    !!{RST
    Implement a ``localGroupLuminosityFunction`` output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Group
    implicit none
    class(outputAnalysisLocalGroupLuminosityFunction), intent(inout)           :: self
    type (varying_string                            ), intent(in   ), optional :: groupName
    type (hdf5Group                                 )               , target   :: analysesGroup, subGroup
    type (hdf5Group                                 )               , pointer  :: inGroup
    type (hdf5Group                                 )                          :: analysisGroup

    ! Finalize analysis.
    call self%finalizeAnalysis()
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('localGroupLuminosityFunction','Analysis of the V-band luminosity function of Local Group satellite galaxies')
    call analysisGroup%writeAttribute('Local Group luminosity function'  ,'description'                                                                )
    call analysisGroup%writeAttribute('function1D'                       ,'type'                                                                       )
    call analysisGroup%writeAttribute('$M_\mathrm{V}$'                   ,'xAxisLabel'                                                                 )
    call analysisGroup%writeAttribute('$N$'                              ,'yAxisLabel'                                                                 )
    call analysisGroup%writeAttribute(.false.                            ,'xAxisIsLog'                                                                 )
    call analysisGroup%writeAttribute(.true.                             ,'yAxisIsLog'                                                                 )
    call analysisGroup%writeAttribute('magnitudeAbsoluteV'               ,'xDataset'                                                                   )
    call analysisGroup%writeAttribute('luminosityFunction'               ,'yDataset'                                                                   )
    call analysisGroup%writeAttribute('luminosityFunctionTarget'         ,'yDatasetTarget'                                                             )
    call analysisGroup%writeAttribute('luminosityFunctionCovariance'     ,'yCovariance'                                                                )
    call analysisGroup%writeAttribute('Observed'                         ,'targetLabel'                                                                )
    call analysisGroup%writeDataset  (self%magnitudes                    ,'magnitudeAbsoluteV'          ,'Absolute V-band magnitude at the bin center' )
    call analysisGroup%writeDataset  (self%luminosityFunction            ,'luminosityFunction'          ,'Satellite number per bin [model]'            )
    call analysisGroup%writeDataset  (self%covariance                    ,'luminosityFunctionCovariance','Satellite number per bin [model; covariance]')
    call analysisGroup%writeDataset  (self%luminosityFunctionTarget      ,'luminosityFunctionTarget'    ,'Satellite number per bin [observed]'         )
    call analysisGroup%writeAttribute(self%logLikelihood               (),'logLikelihood'                                                              )
    !$ call hdf5Access%unset()
    return
  end subroutine localGroupLuminosityFunctionFinalize

  double precision function localGroupLuminosityFunctionLogLikelihood(self)
    !!{RST
    Return the log-likelihood of a ``localGroupLuminosityFunction`` output analysis. The likelihood function assumes that the model prediction for the number of satellite galaxies in any given magnitude bin follows a negative binomial distribution as was found for dark matter subhalos :cite:p:`boylan-kolchin_theres_2010`.
    !!}
    use :: Numerical_Constants_Math         , only : Pi
    use :: Statistics_Distributions_Discrete, only : distributionFunctionDiscrete1DNegativeBinomial
    implicit none
    class           (outputAnalysisLocalGroupLuminosityFunction    ), intent(inout) :: self
    type            (distributionFunctionDiscrete1DNegativeBinomial)                :: distribution
    integer                                                                         :: i
    double precision                                                                :: negativeBinomialProbabilitySuccess, countEffective, &
         &                                                                             variance

    call self%finalizeAnalysis()
    localGroupLuminosityFunctionLogLikelihood=0.0d0
    do i=1,size(self%magnitudes)
       if (self%luminosityFunction(i) <= 0.0d0) then
          if (nint(self%luminosityFunctionTarget(i)) > 0) localGroupLuminosityFunctionLogLikelihood=+localGroupLuminosityFunctionLogLikelihood &
               &                                                                                    +self%logLikelihoodZero
       else
          negativeBinomialProbabilitySuccess =+  1.0d0                                                                &
               &                              /(                                                                      &
               &                                +1.0d0                                                                &
               &                                +self%negativeBinomialScatterFractional**2*self%luminosityFunction(i) &
               &                               )
          if (negativeBinomialProbabilitySuccess >= 1.0d0) then
             if (nint(self%luminosityFunctionTarget(i)) > 0) localGroupLuminosityFunctionLogLikelihood=+localGroupLuminosityFunctionLogLikelihood &
                  &                                                                                    +self%logLikelihoodZero
          else
             ! Compute the likelihood assuming a negative binomial distribution. Note that we "de-normalize" the likelihood by
             ! multiplying by √[2πσᵢ²] (the normalization term in the corresponding normal distribution). This is useful to allow
             ! (-logℒ) to be used as a metric for significant shifts in the model results, without changing the relative
             ! likelihood of models (as this de-normalization shift is a constant multiplicative factor).
             countEffective                           = dble(max(1.0d0,self%luminosityFunctionTarget(i)))
             variance                                 =+       countEffective                       &
                  &                                    *(                                           &
                  &                                      +     1.0d0                                &
                  &                                      +self%negativeBinomialScatterFractional**2 &
                  &                                      *     countEffective                       &
                  &                                     )
             distribution                             = distributionFunctionDiscrete1DNegativeBinomial                (negativeBinomialProbabilitySuccess,     self%countFailures               )
             localGroupLuminosityFunctionLogLikelihood=+localGroupLuminosityFunctionLogLikelihood                                                                                                 &
                  &                                    +distribution                                 %massLogarithmic(                                   nint(self%luminosityFunctionTarget(i)))  &
                  &                                    +0.5d0                                                                                                                                     &
                  &                                    *log(                                                                                                                                      &
                  &                                         +2.0d0                                                                                                                                &
                  &                                         *Pi                                                                                                                                   &
                  &                                         *variance                                                                                                                             &
                  &                                        )
          end if
       end if
    end do
    return
  end function localGroupLuminosityFunctionLogLikelihood
