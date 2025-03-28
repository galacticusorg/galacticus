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
Implements a generic two-point correlation function output analysis class.
!!}

  use               :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use               :: Dark_Matter_Halo_Biases               , only : darkMatterHaloBiasClass
  use               :: Dark_Matter_Halo_Scales               , only : darkMatterHaloScaleClass
  use               :: Dark_Matter_Profiles_DMO              , only : darkMatterProfileDMOClass
  use               :: Galactic_Filters                      , only : galacticFilterClass
  use               :: Geometry_Surveys                      , only : surveyGeometryClass
  use               :: Halo_Model_Power_Spectrum_Modifiers   , only : haloModelPowerSpectrumModifierClass
  use   , intrinsic :: ISO_C_Binding                         , only : c_size_t
  use               :: Node_Property_Extractors              , only : nodePropertyExtractorClass
  !$ use            :: OMP_Lib                               , only : omp_lock_kind
  use               :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorClass
  use               :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorClass
  use               :: Output_Times                          , only : outputTimesClass
  use               :: Power_Spectra                         , only : powerSpectrumClass

  !![
  <outputAnalysis name="outputAnalysisCorrelationFunction">
   <description>
    A generic two-point correlation function output analysis class.
  
    For constraints corresponding to (possibly, projected) correlation functions, the model expectation is computed using the
    halo model \cite{cooray_halo_2002}. For each model halo, each galaxy (satellite and central) is assessed to see if it meets
    the criteria for inclusion in the sample. Where the sample includes mass limits (either just a lower limit, or lower and
    upper limits) each galaxy is assigned a probability of inclusion in the sample based on its mass and the random error in
    mass. Thus, each halo is characterized by the probability of having a central galaxy in the sample, $p^\mathrm{(c)}$, and
    $N$ probabilities, $p_i^\mathrm{(s)}$, of each satellite galaxy being in the sample. We assume binomial statistics for each
    galaxy's probability of inclusion, and further assume that these probabilities are uncorrelated. Therefore, the
    contribution of the halo to the one- and two-halo terms of the power spectrum in the halo model are:
    \begin{equation}
     \Delta P^\mathrm{1h}(k) = {w \over n_\mathrm{gal}^2} \left[ p^\mathrm{(c)} \sum_{i=1}^N p^\mathrm{(s)} u(k|M) + \sum_{k=0}^N k(k-1) P\left(p_i^\mathrm{(s)},\ldots,p_N^{(s)}\right) u(k|M)^2 \right]
    \end{equation}
    and
    \begin{equation}
     \Delta \sqrt{P^\mathrm{2h}}(k) = {w \over n_\mathrm{gal}} \sqrt{P^\mathrm{lin}}(k) b(M) u(k|M) \left[ p^\mathrm{(c)} + \sum_{i=1}^N p_i^\mathrm{(s)} \right],
    \end{equation}
    respectively, where $w$ is the weight of the halo (i.e. the number of such model halos expected per unit volume), $b(M)$ is
    the bias of halos of mass $M$, $u(k|M)$ is the Fourier-transform of the halo density profile, and $P_\mathrm{lin}(k)$ is
    the linear theory power spectrum, and $P(p_1,\ldots,p_N)$ is the Poisson binomial distribution for $N$ events with
    probabilities $p_1,\ldots,p_N$. The contribution of the halo to the galaxy density, $n_\mathrm{gal}$, is simply $\Delta
    n_\mathrm{gal} = w \left[ p^\mathrm{(c)} + \sum_{i=1}^N p_i^\mathrm{(s)} \right]$.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisCorrelationFunction
     !!{
     A generic two-point correlation function output analysis class.
     !!}
     private
     type            (varying_string                         )                                :: label                                                , comment                                     , &
          &                                                                                      targetLabel
     class           (galacticFilterClass                    ), pointer                       :: galacticFilter_                             => null()
     class           (outputTimesClass                       ), pointer                       :: outputTimes_                                => null()
     class           (surveyGeometryClass                    ), pointer                       :: surveyGeometry_                             => null()
     class           (cosmologyFunctionsClass                ), pointer                       :: cosmologyFunctions_                         => null()
     class           (outputAnalysisDistributionOperatorClass), pointer                       :: massDistributionOperator_                   => null()
     class           (outputAnalysisPropertyOperatorClass    ), pointer                       :: massPropertyOperator_                       => null(), separationPropertyOperator_        => null()
     class           (nodePropertyExtractorClass             ), pointer                       :: massPropertyExtractor_                      => null()
     class           (darkMatterProfileDMOClass              ), pointer                       :: darkMatterProfileDMO_                       => null()
     class           (darkMatterHaloBiasClass                ), pointer                       :: darkMatterHaloBias_                         => null()
     class           (darkMatterHaloScaleClass               ), pointer                       :: darkMatterHaloScale_                        => null()
     class           (haloModelPowerSpectrumModifierClass    ), pointer                       :: haloModelPowerSpectrumModifier_             => null()
     class           (powerSpectrumClass                     ), pointer                       :: powerSpectrum_                              => null()
     double precision                                         , allocatable, dimension(:    ) :: separations                                          , wavenumber                                  , &
          &                                                                                      meanDensity                                          , massMinima                                  , &
          &                                                                                      massMaxima                                           , massMinimaLogarithmic                       , &
          &                                                                                      massMaximaLogarithmic                                , binnedProjectedCorrelationTarget1D          , &
          &                                                                                      binnedProjectedCorrelationCovarianceTarget1D
     double precision                                         , allocatable, dimension(:,:  ) :: outputWeight                                         , meanDensityMainBranch                       , &
          &                                                                                      oneHaloTerm                                          , twoHaloTerm                                 , &
          &                                                                                      termCovariance                                       , integralConstraint                          , &
          &                                                                                      binnedProjectedCorrelation                           , binnedProjectedCorrelationCovariance        , &
          &                                                                                      binnedProjectedCorrelationTarget                     , binnedProjectedCorrelationCovarianceTarget
     double precision                                         , allocatable, dimension(:,:,:) :: oneHaloTermMainBranch                                , twoHaloTermMainBranch
     integer                                                  , allocatable, dimension(:    ) :: countMainBranch
     double precision                                                                         :: massHaloLogarithmicMinimum                           , massHaloIntervalLogarithmicInverse          , &
          &                                                                                      wavenumberMinimum                                    , wavenumberMaximum                           , &
          &                                                                                      depthLineOfSight                                     , massHaloMinimum                             , &
          &                                                                                      massHaloMaximum
     integer         (c_size_t                               )                                :: binCount                                             , massCount                                   , &
          &                                                                                      countBinsMassHalo                                    , wavenumberCount
     logical                                                                                  :: finalized                                            , halfIntegral
     !$ integer      (omp_lock_kind                          )                                :: accumulateLock
     ! Workspace used while accumulating the correlation function.
     double precision                                         , allocatable, dimension(:    ) :: probabilityCentral
     double precision                                         , allocatable, dimension(:,:  ) :: probabilitySatellite
     integer                                                                                  :: countSatellites                                      , massHaloBinsPerDecade
   contains
     !![
     <methods>
       <method description="Accumulate a node to the correlation function." method="accumulateNode" />
       <method description="Accumulate a halo to the correlation function." method="accumulateHalo" />
     </methods>
     !!]
     final     ::                   correlationFunctionDestructor
     procedure :: analyze        => correlationFunctionAnalyze
     procedure :: finalize       => correlationFunctionFinalize
     procedure :: reduce         => correlationFunctionReduce
     procedure :: logLikelihood  => correlationFunctionLogLikelihood
     procedure :: accumulateNode => correlationFunctionAccumulateNode
     procedure :: accumulateHalo => correlationFunctionAccumulateHalo
  end type outputAnalysisCorrelationFunction

  interface outputAnalysisCorrelationFunction
     !!{
     Constructors for the {\normalfont \ttfamily correlationFunction} output analysis class.
     !!}
     module procedure correlationFunctionConstructorParameters
     module procedure correlationFunctionConstructorFile
     module procedure correlationFunctionConstructorInternal
  end interface outputAnalysisCorrelationFunction

contains

  function correlationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily correlationFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisCorrelationFunction      )                                :: self
    type            (inputParameters                        ), intent(inout)                 :: parameters
    class           (galacticFilterClass                    ), pointer                       :: galacticFilter_
    class           (outputTimesClass                       ), pointer                       :: outputTimes_
    class           (surveyGeometryClass                    ), pointer                       :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                       :: cosmologyFunctions_
    class           (outputAnalysisDistributionOperatorClass), pointer                       :: massDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                       :: massPropertyOperator_             , separationPropertyOperator_
    class           (nodePropertyExtractorClass             ), pointer                       :: massPropertyExtractor_
    class           (darkMatterProfileDMOClass              ), pointer                       :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                ), pointer                       :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass               ), pointer                       :: darkMatterHaloScale_
    class           (haloModelPowerSpectrumModifierClass    ), pointer                       :: haloModelPowerSpectrumModifier_
    class           (powerSpectrumClass                     ), pointer                       :: powerSpectrum_
    double precision                                         , allocatable  , dimension(:  ) :: separations                       , massMinima                                  , &
         &                                                                                      massMaxima                        , integralConstraint                          , &
         &                                                                                      binnedProjectedCorrelationTarget1D, binnedProjectedCorrelationCovarianceTarget1D
    double precision                                         , allocatable  , dimension(:,:) :: binnedProjectedCorrelationTarge  t, binnedProjectedCorrelationCovarianceTarget
    double precision                                                                         :: massHaloMinimum                   , massHaloMaximum                             , &
         &                                                                                      wavenumberMinimum                 , wavenumberMaximum                           , &
         &                                                                                      depthLineOfSight
    logical                                                                                  :: halfIntegral
    integer         (c_size_t                               )                                :: wavenumberCount
    integer                                                                                  :: massHaloBinsPerDecade
    type            (varying_string                         )                                :: label                             , comment                                     , &
         &                                                                                      targetLabel

    allocate(separations(parameters%count('separations')))
    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>A label for the mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>comment</name>
      <source>parameters</source>
      <description>A descriptive comment for the mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>separations</name>
      <source>parameters</source>
      <description>The separations corresponding to bin centers.</description>
    </inputParameter>
    <inputParameter>
      <name>massMinima</name>
      <source>parameters</source>
      <description>The minimum mass of each mass sample.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaxima</name>
      <source>parameters</source>
      <description>The maximum mass of each mass sample.</description>
    </inputParameter>
    <inputParameter>
      <name>massHaloBinsPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMinimum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMaximum</name>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberCount</name>
      <defaultValue>60_c_size_t</defaultValue>
      <description>The number of bins in wavenumber to use in computing the correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The minimum wavenumber to use when computing the correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMaximum</name>
      <defaultValue>1.0d4</defaultValue>
      <description>The maximum wavenumber to use when computing the correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>integralConstraint</name>
      <description>The integral constraint for these correlation functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>depthLineOfSight</name>
      <description>The line-of-sight depth over which the correlation function was projected.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>halfIntegral</name>
      <description>Set to true if the projection integrand should be over line-of-sight depths greater than zero.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('binnedProjectedCorrelationTarget')) then
       if (parameters%isPresent('binnedProjectedCorrelationCovarianceTarget')) then
          !![
          <inputParameter>
            <name>binnedProjectedCorrelationTarget</name>
            <source>parameters</source>
            <description>The target function for likelihood calculations.</description>
            <variable>binnedProjectedCorrelationTarget1D</variable>
          </inputParameter>
          <inputParameter>
            <name>binnedProjectedCorrelationCovarianceTarget</name>
            <source>parameters</source>
            <variable>binnedProjectedCorrelationCovarianceTarget1D</variable>
            <description>The target function covariance for likelihood calculations.</description>
          </inputParameter>
          !!]
          if (size(binnedProjectedCorrelationCovarianceTarget1D) == size(binnedProjectedCorrelationTarget1D)**2) then
             allocate(binnedProjectedCorrelationTarget          (size(separations                       ),size(massMinima                        )))
             allocate(binnedProjectedCorrelationCovarianceTarget(size(binnedProjectedCorrelationTarget1D),size(binnedProjectedCorrelationTarget1D)))
             binnedProjectedCorrelationTarget          =reshape(binnedProjectedCorrelationTarget1D          ,shape(binnedProjectedCorrelationTarget          ))
             binnedProjectedCorrelationCovarianceTarget=reshape(binnedProjectedCorrelationCovarianceTarget1D,shape(binnedProjectedCorrelationCovarianceTarget))
          else
             call Error_Report('binnedProjectedCorrelationCovarianceTarget has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('binnedProjectedCorrelationCovarianceTarget must be specified if binnedProjectedCorrelationTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('binnedProjectedCorrelationCovarianceTarget')) call Error_Report('binnedProjectedCorrelationTarget must be specified if binnedProjectedCorrelationCovarianceTarget is present'//{introspection:location})
    end if
    !![
    <inputParameter>
      <name>targetLabel</name>
      <source>parameters</source>
      <description>A label for the target dataset in a plot of this analysis.</description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    <objectBuilder class="galacticFilter"                     name="galacticFilter_"                                                            source="parameters"/>
    <objectBuilder class="outputTimes"                        name="outputTimes_"                                                               source="parameters"/>
    <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                                                            source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                                                        source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"               name="darkMatterProfileDMO_"                                                      source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"                 name="darkMatterHaloBias_"                                                        source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"                name="darkMatterHaloScale_"                                                       source="parameters"/>
    <objectBuilder class="haloModelPowerSpectrumModifier"     name="haloModelPowerSpectrumModifier_"                                            source="parameters"/>
    <objectBuilder class="powerSpectrum"                      name="powerSpectrum_"                                                             source="parameters"/>
    <objectBuilder class="outputAnalysisDistributionOperator" name="massDistributionOperator_"       parameterName="massDistributionOperator"   source="parameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator"     name="massPropertyOperator_"           parameterName="massPropertyOperator"       source="parameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator"     name="separationPropertyOperator_"     parameterName="separationPropertyOperator" source="parameters"/>
    <objectBuilder class="nodePropertyExtractor"              name="massPropertyExtractor_"          parameterName="massPropertyExtractor"      source="parameters"/>
    <conditionalCall>
     <call>
     self=outputAnalysisCorrelationFunction(                                                                  &amp;
           &amp;                            label                                                           , &amp;
           &amp;                            comment                                                         , &amp;
           &amp;                            separations                                                     , &amp;
           &amp;                            massMinima                                                      , &amp;
           &amp;                            massMaxima                                                      , &amp;
           &amp;                            massHaloBinsPerDecade                                           , &amp;
           &amp;                            massHaloMinimum                                                 , &amp;
           &amp;                            massHaloMaximum                                                 , &amp;
           &amp;                            reshape(integralConstraint,[size(separations),size(massMinima)]), &amp;
           &amp;                            wavenumberCount                                                 , &amp;
           &amp;                            wavenumberMinimum                                               , &amp;
           &amp;                            wavenumberMaximum                                               , &amp;
           &amp;                            depthLineOfSight                                                , &amp;
           &amp;                            halfIntegral                                                    , &amp;
           &amp;                            galacticFilter_                                                 , &amp;
           &amp;                            surveyGeometry_                                                 , &amp;
           &amp;                            cosmologyFunctions_                                             , &amp;
           &amp;                            outputTimes_                                                    , &amp;
           &amp;                            darkMatterProfileDMO_                                           , &amp;
           &amp;                            darkMatterHaloBias_                                             , &amp;
           &amp;                            darkMatterHaloScale_                                            , &amp;
           &amp;                            haloModelPowerSpectrumModifier_                                 , &amp;
           &amp;                            powerSpectrum_                                                  , &amp;
           &amp;                            massDistributionOperator_                                       , &amp;
           &amp;                            massPropertyOperator_                                           , &amp;
           &amp;                            separationPropertyOperator_                                     , &amp;
           &amp;                            massPropertyExtractor_                                          , &amp;
           &amp;                            targetLabel                                                       &amp;
           &amp;                            {conditions}                                                      &amp;
           &amp;                           )
     </call>
     <argument name="binnedProjectedCorrelationTarget"           value="binnedProjectedCorrelationTarget"           parameterPresent="parameters"/>
     <argument name="binnedProjectedCorrelationCovarianceTarget" value="binnedProjectedCorrelationCovarianceTarget" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"          />
    <objectDestructor name="darkMatterHaloBias_"            />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="haloModelPowerSpectrumModifier_"/>
    <objectDestructor name="powerSpectrum_"                 />
    <objectDestructor name="galacticFilter_"                />
    <objectDestructor name="outputTimes_"                   />
    <objectDestructor name="surveyGeometry_"                />
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="massDistributionOperator_"      />
    <objectDestructor name="massPropertyOperator_"          />
    <objectDestructor name="separationPropertyOperator_"    />
    <objectDestructor name="massPropertyExtractor_"         />
    !!]
    return
  end function correlationFunctionConstructorParameters

  function correlationFunctionConstructorFile(label,comment,fileName,massHaloBinsPerDecade,massHaloMinimum,massHaloMaximum,wavenumberCount,wavenumberMinimum,wavenumberMaximum,halfIntegral,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputTimes_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterHaloScale_,haloModelPowerSpectrumModifier_,powerSpectrum_,massDistributionOperator_,massPropertyOperator_,separationPropertyOperator_,massPropertyExtractor_) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily correlationFunction} output analysis class which reads bin information from a standard format file.
    !!}
    use :: Cosmology_Functions , only : cosmologyFunctionsClass  , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters, only : cosmologyParametersSimple
    use :: HDF5_Access         , only : hdf5Access
    use :: IO_HDF5             , only : hdf5Object
    implicit none
    type            (outputAnalysisCorrelationFunction      )                                :: self
    type            (varying_string                         ), intent(in   )                 :: label                                     , comment
    character       (len=*                                  ), intent(in   )                 :: fileName
    integer         (c_size_t                               ), intent(in   )                 :: wavenumberCount
    double precision                                         , intent(in   )                 :: massHaloMinimum                           , massHaloMaximum                 , &
         &                                                                                      wavenumberMinimum                         , wavenumberMaximum
    logical                                                  , intent(in   )                 :: halfIntegral
    integer                                                  , intent(in   )                 :: massHaloBinsPerDecade
    class           (galacticFilterClass                    ), intent(in   ), target         :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ), target         :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ), target         :: cosmologyFunctions_
    class           (outputTimesClass                       ), intent(in   ), target         :: outputTimes_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ), target         :: massDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), intent(in   ), target         :: massPropertyOperator_                     , separationPropertyOperator_
    class           (nodePropertyExtractorClass             ), intent(in   ), target         :: massPropertyExtractor_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target         :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                ), intent(in   ), target         :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target         :: darkMatterHaloScale_
    class           (haloModelPowerSpectrumModifierClass    ), intent(in   ), target         :: haloModelPowerSpectrumModifier_
    class           (powerSpectrumClass                     ), intent(in   ), target         :: powerSpectrum_
    double precision                                         , allocatable  , dimension(:  ) :: separations                               , massMinima                      , &
         &                                                                                      massMaxima
    double precision                                         , allocatable  , dimension(:,:) :: integralConstraint                        , binnedProjectedCorrelationTarget, &
         &                                                                                      binnedProjectedCorrelationCovarianceTarget
    type            (varying_string                         )                                :: targetLabel
    type            (hdf5Object                             )                                :: dataFile                                  , parametersGroup
    double precision                                                                         :: hubbleParameterData                       , omegaMatterData                 , &
         &                                                                                      omegaDarkEnergyData                       , depthLineOfSight

    !$ call hdf5Access%set()
    call dataFile%openFile(fileName,readOnly=.true.)
    ! Extract parameters.
    parametersGroup=dataFile%openGroup('Parameters')
    call parametersGroup%readAttribute('H_0'                                         ,hubbleParameterData)
    call parametersGroup%readAttribute('Omega_Matter'                                ,omegaMatterData    )
    call parametersGroup%readAttribute('Omega_DE'                                    ,omegaDarkEnergyData)
    call parametersGroup%readAttribute('projectedCorrelationFunctionLineOfSightDepth',depthLineOfSight   )
    call parametersGroup%close()
    ! Extract separations.
    call dataFile%readDataset('separationObserved',separations       )
    ! Extract observed datasets.
    call dataFile%readAttribute('label'                               ,targetLabel                               )
    call dataFile%readDataset  ('projectedCorrelationFunctionObserved',binnedProjectedCorrelationTarget          )
    call dataFile%readDataset  ('covariance'                          ,binnedProjectedCorrelationCovarianceTarget)
    ! Extract integral constraint.
    call dataFile%readDataset('integralConstraint',integralConstraint)
    ! Read the minimum and maximum masses.
    call dataFile%readDataset("massMinimum"       ,massMinima        )
    call dataFile%readDataset("massMaximum"       ,massMaxima        )
    ! Finish reading data file.
    call dataFile%close      (                                       )
    !$ call hdf5Access%unset()
    self=outputAnalysisCorrelationFunction(label,comment,separations,massMinima,massMaxima,massHaloBinsPerDecade,massHaloMinimum,massHaloMaximum,integralConstraint,wavenumberCount,wavenumberMinimum,wavenumberMaximum,depthLineOfSight,halfIntegral,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputTimes_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterHaloScale_,haloModelPowerSpectrumModifier_,powerSpectrum_,massDistributionOperator_,massPropertyOperator_,separationPropertyOperator_,massPropertyExtractor_,targetLabel,binnedProjectedCorrelationTarget,binnedProjectedCorrelationCovarianceTarget)
   return
  end function correlationFunctionConstructorFile

  function correlationFunctionConstructorInternal(label,comment,separations,massMinima,massMaxima,massHaloBinsPerDecade,massHaloMinimum,massHaloMaximum,integralConstraint,wavenumberCount,wavenumberMinimum,wavenumberMaximum,depthLineOfSight,halfIntegral,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputTimes_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterHaloScale_,haloModelPowerSpectrumModifier_,powerSpectrum_,massDistributionOperator_,massPropertyOperator_,separationPropertyOperator_,massPropertyExtractor_,targetLabel,binnedProjectedCorrelationTarget,binnedProjectedCorrelationCovarianceTarget) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily correlationFunction} output analysis class for internal use.
    !!}
    use, intrinsic :: ISO_C_Binding            , only : c_size_t
    use            :: Numerical_Ranges         , only : Make_Range                                 , rangeTypeLogarithmic
    use            :: Output_Analysis_Utilities, only : Output_Analysis_Output_Weight_Survey_Volume
    implicit none
    type            (outputAnalysisCorrelationFunction      )                                          :: self
    type            (varying_string                         ), intent(in   )                           :: label                          , comment
    double precision                                         , intent(in   ), dimension(:  )           :: separations                    , massMinima                                   , &
         &                                                                                                massMaxima
    double precision                                         , intent(in   ), dimension(:,:)           :: integralConstraint
    integer         (c_size_t                               ), intent(in   )                           :: wavenumberCount
    double precision                                         , intent(in   )                           :: massHaloMinimum                , massHaloMaximum                              , &
         &                                                                                                wavenumberMinimum              , wavenumberMaximum                            , &
         &                                                                                                depthLineOfSight
    logical                                                  , intent(in   )                           :: halfIntegral
    integer                                                  , intent(in   )                           :: massHaloBinsPerDecade
    class           (galacticFilterClass                    ), intent(in   ), target                   :: galacticFilter_
    class           (outputTimesClass                       ), intent(in   ), target                   :: outputTimes_
    class           (surveyGeometryClass                    ), intent(in   ), target                   :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ), target                   :: cosmologyFunctions_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ), target                   :: massDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), intent(in   ), target                   :: massPropertyOperator_          , separationPropertyOperator_
    class           (nodePropertyExtractorClass             ), intent(in   ), target                   :: massPropertyExtractor_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target                   :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                ), intent(in   ), target                   :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target                   :: darkMatterHaloScale_
    class           (haloModelPowerSpectrumModifierClass    ), intent(in   ), target                   :: haloModelPowerSpectrumModifier_
    class           (powerSpectrumClass                     ), intent(in   ), target                   :: powerSpectrum_
    type            (varying_string                         ), intent(in   )                , optional :: targetLabel
    double precision                                         , intent(in   ), dimension(:,:), optional :: binnedProjectedCorrelationTarget  , binnedProjectedCorrelationCovarianceTarget
    integer         (c_size_t                               )                                          :: i
    !![
    <constructorAssign variables="label, comment, separations, massMinima, massMaxima, massHaloBinsPerDecade, massHaloMinimum, massHaloMaximum, integralConstraint, wavenumberCount, wavenumberMinimum, wavenumberMaximum, depthLineOfSight, halfIntegral, *galacticFilter_, *outputTimes_, *darkMatterProfileDMO_, *darkMatterHaloBias_, *darkMatterHaloScale_, *haloModelPowerSpectrumModifier_, *surveyGeometry_, *cosmologyFunctions_, *powerSpectrum_, *massDistributionOperator_, *massPropertyOperator_, *separationPropertyOperator_, *massPropertyExtractor_, targetLabel, binnedProjectedCorrelationTarget, binnedProjectedCorrelationCovarianceTarget"/>
    !!]

    ! Assign 1D versions of target for use in descriptor.
    self%binnedProjectedCorrelationTarget1D          =reshape(binnedProjectedCorrelationTarget,[size(binnedProjectedCorrelationTarget          )])
    self%binnedProjectedCorrelationCovarianceTarget1D=reshape(binnedProjectedCorrelationTarget,[size(binnedProjectedCorrelationCovarianceTarget)])
    ! Compute weight that applies to each output redshift.
    self%massCount=size(self%massMinima               )
    self%binCount =size(self%separations,kind=c_size_t)
    allocate(self%outputWeight(self%massCount,self%outputTimes_%count()))
    do i=1,self%massCount
       self%outputWeight(i,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,self%outputTimes_,massMinima(i))
    end do
    ! Construct halo mass bins.
    self%massHaloLogarithmicMinimum        =                                 log10(                massHaloMinimum)
    self%countBinsMassHalo                 =                             int(log10(massHaloMaximum/massHaloMinimum)*dble(massHaloBinsPerDecade)+0.5d0)
    self%massHaloIntervalLogarithmicInverse=dble(self%countBinsMassHalo)/    log10(massHaloMaximum/massHaloMinimum)
    ! Allocate wavenumbers.
    allocate(self%wavenumber           (self%wavenumberCount                                      ))
    allocate(self%probabilityCentral   (                     self%massCount                       ))
    allocate(self%probabilitySatellite (     1              ,self%massCount                       ))
    allocate(self%meanDensity          (                     self%massCount                       ))
    allocate(self%meanDensityMainBranch(                     self%massCount,self%countBinsMassHalo))
    allocate(self%countMainBranch      (                                    self%countBinsMassHalo))
    allocate(self%oneHaloTermMainBranch(self%wavenumberCount,self%massCount,self%countBinsMassHalo))
    allocate(self%twoHaloTermMainBranch(self%wavenumberCount,self%massCount,self%countBinsMassHalo))
    allocate(self%oneHaloTerm          (self%wavenumberCount,self%massCount                       ))
    allocate(self%twoHaloTerm          (self%wavenumberCount,self%massCount                       ))
    allocate(self%termCovariance       (self%massCount*(2*self%wavenumberCount+1),self%massCount*(2*self%wavenumberCount+1)))
    self%wavenumber=Make_Range(self%wavenumberMinimum,self%wavenumberMaximum,int(self%wavenumberCount),rangeTypeLogarithmic)
    ! Compute logarithmic masses.
    allocate(self%massMinimaLogarithmic(self%massCount))
    allocate(self%massMaximaLogarithmic(self%massCount))
    self%massMinimaLogarithmic=log10(self%massMinima)
    self%massMaximaLogarithmic=log10(self%massMaxima)
    ! Initialize population statistics.
    self%countMainBranch      =0
    self%meanDensity          =0.0d0
    self%meanDensityMainBranch=0.0d0
    self%oneHaloTermMainBranch=0.0d0
    self%twoHaloTermMainBranch=0.0d0
    self%oneHaloTerm          =0.0d0
    self%twoHaloTerm          =0.0d0
    self%termCovariance       =0.0d0
    self%probabilityCentral   =0.0d0
    self%probabilitySatellite =0.0d0
    self%countSatellites      =0
    self%finalized            =.false.
    ! Initialize accumulation lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
   return
  end function correlationFunctionConstructorInternal

  subroutine correlationFunctionDestructor(self)
    !!{
    Destructor for  the {\normalfont \ttfamily correlationFunction} output analysis class.
    !!}
    implicit none
    type(outputAnalysisCorrelationFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"          />
    <objectDestructor name="self%darkMatterHaloBias_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%haloModelPowerSpectrumModifier_"/>
    <objectDestructor name="self%galacticFilter_"                />
    <objectDestructor name="self%outputTimes_"                   />
    <objectDestructor name="self%surveyGeometry_"                />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%massDistributionOperator_"      />
    <objectDestructor name="self%massPropertyOperator_"          />
    <objectDestructor name="self%separationPropertyOperator_"    />
    <objectDestructor name="self%massPropertyExtractor_"         />
    <objectDestructor name="self%powerSpectrum_"                 />
    !!]
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine correlationFunctionDestructor

  subroutine correlationFunctionAnalyze(self,node,iOutput)
    !!{
    Implement a correlationFunction output analysis.
    !!}
    use :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    implicit none
    class           (outputAnalysisCorrelationFunction        ), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    integer         (c_size_t                                 ), intent(in   ) :: iOutput
    double precision                                                           :: mass
    type            (enumerationOutputAnalysisPropertyTypeType)                :: massType

    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    ! Filter this node.
    ! Extract mass property.
    massType=self%massPropertyExtractor_%type()
    select type (extractor_ => self%massPropertyExtractor_)
       class is (nodePropertyExtractorScalar)
       mass=extractor_%extract(node)
    class default
       mass=0.0d0
    end select
    ! Apply operator to mass.
    mass=self%massPropertyOperator_%operate(mass,node,massType,iOutput)
    ! Accumulate the node.
    call self%accumulateNode(mass,massType,iOutput,node)
    ! Accumulate the halo. Since nodes are visited depth first we can be sure that all satellite of this isolated halo have
    ! already been visited.
    if (.not.node%isSatellite()) call self%accumulateHalo(iOutput,node)
    return
  end subroutine correlationFunctionAnalyze

  subroutine correlationFunctionAccumulateNode(self,mass,massType,indexOutput,node)
    !!{
    Accumulate a single galaxy to the population of the current halo. Since galaxy masses
    have random errors, each galaxy added is assigned an inclusion probability, which will be
    taken into account when evaluating the one- and two-halo terms from this halo in the halo
    model.
    !!}
    implicit none
    class           (outputAnalysisCorrelationFunction        ), intent(inout)                 :: self
    double precision                                           , intent(in   )                 :: mass
    type            (enumerationOutputAnalysisPropertyTypeType), intent(in   )                 :: massType
    integer         (c_size_t                                 ), intent(in   )                 :: indexOutput
    type            (treeNode                                 ), intent(inout)                 :: node
    double precision                                           , allocatable  , dimension(:,:) :: probabilitySatelliteTmp
    double precision                                                          , dimension(1  ) :: massDistribution
    logical                                                                                    :: satelliteIncluded
    integer                                                                                    :: j

    ! Evaluate, for each mass bin, the probability of inclusion of the galaxy in that bin. Store any such non-zero probabilities
    ! for central and satellite galaxies separately.
    satelliteIncluded=.false.
    do j=1,size(self%massMinimaLogarithmic)
       ! Find the probability that this galaxy is included in the sample.
       massDistribution=self%massDistributionOperator_%operateScalar(mass,massType,self%massMinimaLogarithmic(j:j),self%massMaximaLogarithmic(j:j),indexOutput,node)
       if (node%isSatellite()) then
          if (massDistribution(1) > 0.0d0) then
             if (.not.satelliteIncluded) then
                satelliteIncluded   =.true.
                self%countSatellites=self%countSatellites+1
                if (size(self%probabilitySatellite,dim=1) < self%countSatellites) then
                   call move_alloc(self%probabilitySatellite,probabilitySatelliteTmp)
                   allocate(self%probabilitySatellite(2*size(probabilitySatelliteTmp,dim=1),size(self%massMinimaLogarithmic)))
                   self%probabilitySatellite(1:size(probabilitySatelliteTmp,dim=1),:)=probabilitySatelliteTmp
                   deallocate(probabilitySatelliteTmp)
                end if
                self%probabilitySatellite(self%countSatellites,:)=0.0d0
             end if
             self%probabilitySatellite(self%countSatellites,j)=massDistribution(1)
          end if
       else
          self%probabilityCentral(j)=massDistribution(1)
       end if
    end do
    return
  end subroutine correlationFunctionAccumulateNode

  subroutine correlationFunctionAccumulateHalo(self,indexOutput,node)
    !!{
    Accumulate a single halo's contributions to the halo model one- and two-halo terms. For
    the one-halo term we count contributions from central-satellite pairs, and from
    satellite-satellite pairs. Contributions differ in the scalings applied to the
    Fourier-transformed dark matter halo density profile---see
    \cite[][\S6.1]{cooray_halo_2002} for a discussion of this. The number of satellites in
    the halo is assumed to follow a Poisson binomial distribution.
    !!}
    use :: Galacticus_Nodes                   , only : nodeComponentBasic                , treeNode
    use :: Halo_Model_Power_Spectrum_Modifiers, only : haloModelTermOneHalo              , haloModelTermTwoHalo
    use :: Mass_Distributions                 , only : massDistributionClass
    use :: Math_Distributions_Poisson_Binomial, only : Poisson_Binomial_Distribution_Mean, Poisson_Binomial_Distribution_Mean_Pairs , Poisson_Binomial_Distribution_Mean_Pairs_Jacobian
    use :: Output_Analyses_Options            , only : outputAnalysisPropertyTypeLinear  , enumerationOutputAnalysisPropertyTypeType
    use :: Vectors                            , only : Vector_Outer_Product
    use :: Linear_Algebra                     , only : assignment(=)                     , matrix                                   , operator(*)
    implicit none
    class           (outputAnalysisCorrelationFunction        ), intent(inout)                                                 :: self
    integer         (c_size_t                                 ), intent(in   )                                                 :: indexOutput
    type            (treeNode                                 ), intent(inout)                                                 :: node
    class           (nodeComponentBasic                       ), pointer                                                       :: basic                   , basicRoot
    class           (massDistributionClass                    ), pointer                                                       :: massDistribution_
    double precision                                                          , dimension(self%wavenumberCount,self%massCount) :: oneHaloTerm             , twoHaloTerm
    double precision                                                          , dimension(                     self%massCount) :: galaxyDensity
    logical                                                                   , dimension(                     self%massCount) :: oneHaloTermActive       , twoHaloTermActive
    double precision                                           , allocatable  , dimension(:                   ,:             ) :: termJacobian            , termCovariance            , &
         &                                                                                                                        mainBranchTermCovariance, modifierCovariance
    double precision                                           , allocatable  , dimension(:                                  ) :: satelliteJacobian       , modifierCovarianceDiagonal, &
         &                                                                                                                        fourierProfile          , wavenumber
    double precision                                                                                                           :: countSatellitePairsMean , countSatellitesMean       , &
         &                                                                                                                        haloWeightOutput        , expansionFactor           , &
         &                                                                                                                        biasHalo                , massHalo                  , &
         &                                                                                                                        radiusVirial
    integer         (c_size_t                                 )                                                                :: i                       , j                         , &
         &                                                                                                                        indexOneHalo            , indexTwoHalo              , &
         &                                                                                                                        indexDensity
    integer                                                                                                                    :: haloMassBin
    type            (enumerationOutputAnalysisPropertyTypeType)                                                                :: scaleType
    logical                                                                                                                    :: mainBranchCounted
    type            (matrix                                   )                                                                :: jacobianMatrix

    ! Return immediately if no nodes have been accumulated.
    if (all(self%probabilityCentral == 0.0d0) .and. self%countSatellites == 0) return
    ! Construct the Fourier profile of the host halo. We include the weighting by the square-root of the power spectrum here.
    massDistribution_ => self%darkMatterProfileDMO_%get            (                       node        )
    radiusVirial      =  self%darkMatterHaloScale_ %radiusVirial   (                       node        )
    expansionFactor   =  self%cosmologyFunctions_  %expansionFactor(self%outputTimes_%time(indexOutput))
    allocate(wavenumber    (self%wavenumberCount))
    allocate(fourierProfile(self%wavenumberCount))
    do i=1,self%wavenumberCount
       ! Note that wavenumbers must be converted from comoving to physical units for the dark matter profile k-space function.
       scaleType        =outputAnalysisPropertyTypeLinear
       wavenumber    (i)=+1.0d0                                                                                         &
            &            /self%separationPropertyOperator_%operate(1.0d0/self%waveNumber(i),node,scaleType,indexOutput)
       fourierProfile(i)=+      massDistribution_%fourierTransform(                                                           &
            &                                                           radiusVirial                                        , &
            &                                                           wavenumber  (i)/expansionFactor                       &
            &                                                     )                                                           &
            &            *sqrt(                                                                                               &
            &                  +self%powerSpectrum_%power         (self%wavenumber  (i),self%outputTimes_%time(indexOutput))  &
            &                 )
    end do
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Get the mass of this halo.
    basic    => node %basic()
    massHalo =  basic%mass ()
    ! Get the bias of this halo.
    biasHalo=self%darkMatterHaloBias_%bias(node)
    ! Accumulate.
    oneHaloTermActive=.false.
    twoHaloTermActive=.false.
    mainBranchCounted=.false.
    allocate(termJacobian              (self%massCount*(2*self%wavenumberCount+1),self%countSatellites+1))
    allocate(satelliteJacobian         (                                          self%countSatellites  ))
    allocate(modifierCovariance        (                  self%wavenumberCount   ,self%wavenumberCount  ))
    allocate(modifierCovarianceDiagonal(self%massCount*(2*self%wavenumberCount+1)                       ))
    termJacobian              =0.0d0
    modifierCovariance        =0.0d0
    modifierCovarianceDiagonal=0.0d0
    ! Iterate over masses.
    do i=1,self%massCount
       ! Find mean number of satellites and satellite pairs.
       if (self%countSatellites > 0) then
          countSatellitesMean    =Poisson_Binomial_Distribution_Mean      (self%probabilitySatellite(1:self%countSatellites,i))
          countSatellitePairsMean=Poisson_Binomial_Distribution_Mean_Pairs(self%probabilitySatellite(1:self%countSatellites,i))
       else
          countSatellitesMean    =0.0d0
          countSatellitePairsMean=0.0d0
       end if
       ! Skip if this halo contains no galaxies.
       if (self%probabilityCentral(i) > 0.0d0 .or. countSatellitesMean > 0.0d0) then
          ! Compute output halo weight.
          haloWeightOutput=node%hostTree%volumeWeight*self%outputWeight(i,indexOutput)
          ! Compute contribution to galaxy density.
          galaxyDensity(i)=+haloWeightOutput             &
               &           *(                            &
               &             +self%probabilityCentral(i) &
               &             +countSatellitesMean        &
               &            )
          ! For main branch galaxies, accumulate their contribution to the density as a function of halo mass, so that we can later subtract this from the variance.
          if (node%isOnMainBranch()) then
             basicRoot            => node     %hostTree%nodeBase%basic()
             haloMassBin          =  floor((log10(basicRoot%mass())-self%massHaloLogarithmicMinimum)*self%massHaloIntervalLogarithmicInverse)+1
             ! Accumulate weights to halo mass arrays.
             if (haloMassBin >= 1 .and. haloMassBin <= self%countBinsMassHalo) then
                self        %meanDensityMainBranch(  i,haloMassBin)=   &
                     & +self%meanDensityMainBranch(  i,haloMassBin)    &
                     & +     haloWeightOutput                          &
                     & *self%probabilityCentral   (  i            )
                self        %oneHaloTermMainBranch(:,i,haloMassBin)=   &
                     & +self%oneHaloTermMainBranch(:,i,haloMassBin)    &
                     & +     haloWeightOutput                          &
                     & *self%probabilityCentral   (  i            )    &
                     & *     countSatellitesMean                       &
                     & *     fourierProfile
                self        %twoHaloTermMainBranch(:,i,haloMassBin)=   &
                     & +self%twoHaloTermMainBranch(:,i,haloMassBin)    &
                     & +     haloWeightOutput                          &
                     & *self%probabilityCentral   (  i            )    &
                     & *     biasHalo                                  &
                     & *     fourierProfile
                ! If this is the first mass bin in which the central, main branch galaxy is seen, increment the number of main branch galaxies.
                if (.not.mainBranchCounted) then
                   mainBranchCounted                =.true.
                   self%countMainBranch(haloMassBin)=self%countMainBranch(haloMassBin)+1
                end if
             end if
          end if
          ! Accumulate contribution to galaxy density.
          self%meanDensity(i)=+self%meanDensity(i) &
               &              +galaxyDensity   (i)
          ! Compute and accumulate one-halo term.
          if (countSatellitesMean > 0.0d0) then
             oneHaloTermActive(  i)=.true.
             oneHaloTerm      (:,i)=+haloWeightOutput                     &
                  &                 *(                                    &
                  &                   +self%probabilityCentral     (i)    &
                  &                   *     countSatellitesMean           &
                  &                   *     fourierProfile                &
                  &                   +     countSatellitePairsMean       &
                  &                   *     fourierProfile            **2 &
                  &                  )
             call self%haloModelPowerSpectrumModifier_%modify(                           &
                  &                                                wavenumber          , &
                  &                                                haloModelTermOneHalo, &
                  &                                                oneHaloTerm(:,i)    , &
                  &                                                modifierCovariance  , &
                  &                                           mass=massHalo              &
                  &                                          )
             call correlationFunctionTermIndices(i,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
             forall(j=1:self%wavenumberCount)
                modifierCovarianceDiagonal(indexOneHalo+j-1)=modifierCovariance(j,j)
             end forall
             modifierCovarianceDiagonal(indexDensity)=0.0d0
             self%oneHaloTerm(:,i)=+self%oneHaloTerm(:,i) &
                  &                +     oneHaloTerm(:,i)
          end if
          ! Compute and accumulate two-halo term.
          twoHaloTermActive(  i)=.true.
          twoHaloTerm      (:,i)=+galaxyDensity (i) &
               &                 *biasHalo          &
               &                 *fourierProfile
          call self%haloModelPowerSpectrumModifier_%modify(                           &
               &                                                wavenumber          , &
               &                                                haloModelTermTwoHalo, &
               &                                                twoHaloTerm(:,i)    , &
               &                                                modifierCovariance  , &
               &                                           mass=massHalo              &
               &                                          )
          call correlationFunctionTermIndices(i,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
          forall(j=1:self%wavenumberCount)
             modifierCovarianceDiagonal(indexTwoHalo+j-1)=modifierCovariance(j,j)
          end forall
          modifierCovarianceDiagonal(indexDensity)=0.0d0
          self%twoHaloTerm(:,i)=+self%twoHaloTerm(:,i) &
               &                +     twoHaloTerm(:,i)
          ! Construct Jacobian of the terms being accumulated. The Jacobian here is an MxN matrix, where M=massCount*(2*wavenumberCount+1)
          ! (the number of terms in halo model quantities being accumulated {wavenumberCount for 1- and 2-halo terms, plus a density, for
          ! each mass bin}), and N is the total number of galaxies in the halo (number of satellites plus 1 central).
          ! Compute indices.
          call correlationFunctionTermIndices(i,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
          ! One halo terms.
          if (self%countSatellites > 0) then
             satelliteJacobian=Poisson_Binomial_Distribution_Mean_Pairs_Jacobian(self%probabilitySatellite(1:self%countSatellites,i))*self%probabilitySatellite(1:self%countSatellites,i)
             do j=1,self%wavenumberCount
                termJacobian(indexOneHalo             +j                   -1,1:self%countSatellites  )=haloWeightOutput         *fourierProfile(j)**2*     satelliteJacobian
             end do
          end if
          termJacobian      (indexOneHalo:indexOneHalo+self%wavenumberCount-1,  self%countSatellites+1)=haloWeightOutput         *fourierProfile      *self%probabilityCentral  (                       i)*countSatellitesMean
          ! Two halo terms.
          do j=1,self%wavenumberCount
             termJacobian   (indexTwoHalo             +j                   -1,1:self%countSatellites  )=haloWeightOutput*biasHalo*fourierProfile(j)   *self%probabilitySatellite(1:self%countSatellites,i)
          end do
          termJacobian      (indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,  self%countSatellites+1)=haloWeightOutput*biasHalo*fourierProfile      *self%probabilityCentral  (                       i)
          ! Compute density terms.
          termJacobian      (indexDensity                                    ,1:self%countSatellites  )=haloWeightOutput                              *self%probabilitySatellite(1:self%countSatellites,i)
          termJacobian      (indexDensity                                    ,  self%countSatellites+1)=haloWeightOutput                              *self%probabilityCentral  (                       i)
       end if
    end do
    ! Construct and accumulate term covariance.
    allocate(termCovariance(self%massCount*(2*self%wavenumberCount+1),self%massCount*(2*self%wavenumberCount+1)))
    jacobianMatrix=termJacobian
    termCovariance=jacobianMatrix*jacobianMatrix%transpose()
    ! Add modifier covariance.
    termCovariance=termCovariance+Vector_Outer_Product(modifierCovarianceDiagonal)
    ! For main branch galaxies, zero all off-diagonal contributions.
    if (node%isOnMainBranch()) then
       termJacobian(:,1:self%countSatellites)=0.0d0
       jacobianMatrix=termJacobian
       allocate(mainBranchTermCovariance(self%massCount*(2*self%wavenumberCount+1),self%massCount*(2*self%wavenumberCount+1)))
       mainBranchTermCovariance=jacobianMatrix*jacobianMatrix%transpose()
       do i=1,self%massCount
          mainBranchTermCovariance(                                                                 &
               &                   (i-1)*(2*self%wavenumberCount+1)+1:i*(2*self%wavenumberCount+1), &
               &                   (i-1)*(2*self%wavenumberCount+1)+1:i*(2*self%wavenumberCount+1)  &
               &                  )                                                                 &
               &                  =0.0d0
       end do
       termCovariance=termCovariance-mainBranchTermCovariance
       deallocate(mainBranchTermCovariance)
    end if
    self%termCovariance=self%termCovariance+termCovariance
    deallocate(termJacobian     )
    deallocate(satelliteJacobian)
    ! Reset counts.
    self%probabilityCentral=0.0d0
    self%countSatellites   =0
    return
  end subroutine correlationFunctionAccumulateHalo

  subroutine correlationFunctionTermIndices(iMass,wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
    !!{
    Return the indices in the term covariances array at which one-halo, two-halo, and density terms are stored for the given
    mass.
    !!}
    implicit none
    integer(c_size_t), intent(in   ) :: iMass       , wavenumberCount
    integer(c_size_t), intent(  out) :: indexOneHalo, indexTwoHalo   , &
         &                              indexDensity

    indexOneHalo=(iMass-1)*(2*wavenumberCount+1)                  +1
    indexTwoHalo=(iMass-1)*(2*wavenumberCount+1)+  wavenumberCount+1
    indexDensity=(iMass-1)*(2*wavenumberCount+1)+2*wavenumberCount+1
    return
  end subroutine correlationFunctionTermIndices

  subroutine correlationFunctionReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily correlationFunction} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisCorrelationFunction), intent(inout) :: self
    class(outputAnalysisClass              ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisCorrelationFunction)
       !$ call OMP_Set_Lock(reduced%accumulateLock)
       reduced%meanDensity          =reduced%meanDensity          +self%meanDensity
       reduced%oneHaloTerm          =reduced%oneHaloTerm          +self%oneHaloTerm
       reduced%twoHaloTerm          =reduced%twoHaloTerm          +self%twoHaloTerm
       reduced%countMainBranch      =reduced%countMainBranch      +self%countMainBranch
       reduced%meanDensityMainBranch=reduced%meanDensityMainBranch+self%meanDensityMainBranch
       reduced%oneHaloTermMainBranch=reduced%oneHaloTermMainBranch+self%oneHaloTermMainBranch
       reduced%twoHaloTermMainBranch=reduced%twoHaloTermMainBranch+self%twoHaloTermMainBranch
       reduced%termCovariance       =reduced%termCovariance       +self%termCovariance
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine correlationFunctionReduce

  subroutine correlationFunctionFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily correlationFunction} output analysis finalization.
    !!}
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(outputAnalysisCorrelationFunction), intent(inout)           :: self
    type (varying_string                   ), intent(in   ), optional :: groupName
    type (hdf5Object                       )               , target   :: analysesGroup, subGroup
    type (hdf5Object                       )               , pointer  :: inGroup
    type (hdf5Object                       )                          :: analysisGroup, dataset

    ! Finalize analysis.
    call correlationFunctionFinalizeAnalysis(self)
    ! Output the correlation function.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('correlationFunction'//char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call    analysisGroup%writeAttribute(char(self%comment)                                  ,'description'                                                                                          )
    call    analysisGroup%writeAttribute("function1DSequence"                                ,'type'                                                                                                 )
    call    analysisGroup%writeAttribute('$s\, [_\chi\mathrm{Mpc}]$'                         ,'xAxisLabel'                                                                                           )
    call    analysisGroup%writeAttribute('$w(s)\, [_\chi\mathrm{Mpc}]$'                      ,'yAxisLabel'                                                                                           )
    call    analysisGroup%writeAttribute(.true.                                              ,'xAxisIsLog'                                                                                           )
    call    analysisGroup%writeAttribute(.true.                                              ,'yAxisIsLog'                                                                                           )
    call    analysisGroup%writeAttribute('separation'                                        ,'xDataset'                                                                                             )
    call    analysisGroup%writeAttribute('correlationFunction'                               ,'yDataset'                                                                                             )
    call    analysisGroup%writeAttribute('correlationFunctionTarget'                         ,'yDatasetTarget'                                                                                       )
    call    analysisGroup%writeAttribute('correlationFunctionCovariance'                     ,'yCovariance'                                                                                          )
    call    analysisGroup%writeAttribute('correlationFunctionCovarianceTarget'               ,'yCovarianceTarget'                                                                                    )
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (self%separations                                    ,'separation'                         ,'Separation'                             ,datasetReturned=dataset)
    call    dataset      %writeAttribute('ᵪMpc'                                              ,'units'                                                                                                )
    call    dataset      %writeAttribute(megaParsec                                          ,'unitsInSI'                                                                                            )
    call    dataset      %close         (                                                                                                                                                            )
    call    analysisGroup%writeDataset  (self%binnedProjectedCorrelation                     ,'correlationFunction'                ,'Projected correlation'                  ,datasetReturned=dataset)
    call    dataset      %writeAttribute('ᵪMpc'                                              ,'units'                                                                                                )
    call    dataset      %writeAttribute(megaParsec                                          ,'unitsInSI'                                                                                            )
    call    dataset      %close         (                                                                                                                                                            )
    call    analysisGroup%writeDataset  (self%binnedProjectedCorrelationCovariance           ,'correlationFunctionCovariance'      ,'Projected correlation covariance'       ,datasetReturned=dataset)
    call    dataset      %writeAttribute('ᵪMpc²'                                             ,'units'                                                                                                )
    call    dataset      %writeAttribute(megaParsec**2                                       ,'unitsInSI'                                                                                            )
    call    dataset      %close         (                                                                                                                                                            )
    ! If available, include the log-likelihood and target dataset.
    if (allocated(self%binnedProjectedCorrelationTarget)) then
       call analysisGroup%writeAttribute(     self%logLikelihood()                           ,'logLikelihood'                                                                                        )
       call analysisGroup%writeAttribute(char(self%targetLabel)                              ,'targetLabel'                                                                                          )
       call analysisGroup%writeDataset  (     self%binnedProjectedCorrelationTarget          ,"correlationFunctionTarget"          ,'Projected correlation target'           ,datasetReturned=dataset)
       call dataset      %writeAttribute(     "ᵪMpc"                                         ,'units'                                                                                                )
       call dataset      %writeAttribute(     megaParsec                                     ,'unitsInSI'                                                                                            )
       call dataset      %close         (                                                                                                                                                            )
       call analysisGroup%writeDataset  (     self%binnedProjectedCorrelationCovarianceTarget,"correlationFunctionCovarianceTarget",'Projected correlation covariance target',datasetReturned=dataset)
       call dataset      %writeAttribute(     "ᵪMpc²"                                        ,'units'                                                                                                )
       call dataset      %writeAttribute(      megaParsec**2                                 ,'unitsInSI'                                                                                            )
       call dataset      %close         (                                                                                                                                                            )
    end if
    call    analysisGroup%close         (                                                                                                                                                            )
    if (present(groupName)) &
         & call subGroup %close         (                                                                                                                                                            )
    call    analysesGroup%close         (                                                                                                                                                            )
    !$ call hdf5Access%unset()
    return
  end subroutine correlationFunctionFinalize

  subroutine correlationFunctionFinalizeAnalysis(self)
    !!{
    Compute final covariances and normalize.
    !!}
    use :: FFTLogs                 , only : FFTLogSineTransform                , fftLogForward
#ifdef USEMPI
    use :: MPI_Utilities           , only : mpiSelf
#endif
    use :: Numerical_Constants_Math, only : Pi
    use :: Table_Labels            , only : extrapolationTypeExtrapolate
    use :: Tables                  , only : table1DLogarithmicLinear           , tablesIntegrationWeightFunction
    use :: Vectors                 , only : Matrix_Copy_Upper_To_Lower_Triangle, Vector_Outer_Product
    use :: Linear_Algebra          , only : assignment(=)                      , matrix                         , operator(*)
    implicit none
    class           (outputAnalysisCorrelationFunction), intent(inout)                 :: self
    double precision                                   , allocatable  , dimension(:  ) :: separation
    double precision                                   , allocatable  , dimension(:,:) :: powerSpectrumCovariance       , jacobian            , &
         &                                                                                correlationCovariance         , covarianceTmp       , &
         &                                                                                projectedCorrelationCovariance, oneTwoHaloCovariance, &
         &                                                                                powerSpectrumValue            , correlation         , &
         &                                                                                projectedCorrelation
    integer         (c_size_t                         )                                :: i                             , j                   , &
         &                                                                                m                             , n                   , &
         &                                                                                indexOneHalo                  , indexTwoHalo        , &
         &                                                                                indexDensity
    type            (table1DLogarithmicLinear         )                                :: correlationTable
    double precision                                                                   :: projectedSeparation           , binSeparationMinimum, &
         &                                                                                binSeparationMaximum          , binWidthLogarithmic
    type            (matrix                           )                                :: jacobianMatrix                , covarianceMatrix
    procedure       (tablesIntegrationWeightFunction  ), pointer                       :: integrandWeightFunction

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
#ifdef USEMPI
    ! If running under MPI, perform a summation reduction across all processes.
    self%meanDensity          =mpiSelf%sum(self%meanDensity          )
    self%oneHaloTerm          =mpiSelf%sum(self%oneHaloTerm          )
    self%twoHaloTerm          =mpiSelf%sum(self%twoHaloTerm          )
    self%countMainBranch      =mpiSelf%sum(self%countMainBranch      )
    self%meanDensityMainBranch=mpiSelf%sum(self%meanDensityMainBranch)
    self%oneHaloTermMainBranch=mpiSelf%sum(self%oneHaloTermMainBranch)
    self%twoHaloTermMainBranch=mpiSelf%sum(self%twoHaloTermMainBranch)
    self%termCovariance       =mpiSelf%sum(self%termCovariance       )
#endif
   ! Copy upper to lower triangle of covariance matrix (we've accumulated only the upper triangle).
    self%termCovariance=Matrix_Copy_Upper_To_Lower_Triangle(self%termCovariance)
    ! Find average density contribution of main branch galaxies in each halo mass bin.
    do n=1,self%massCount
       where    (self%countMainBranch(:) > 0)
          self                 %meanDensityMainBranch(  n,:)  &
               &    =+     self%meanDensityMainBranch(  n,:)  &
               &     /dble(self%countMainBranch      (    :))
       end where
       do i=1,self%wavenumberCount
          where (self%countMainBranch(:) > 0)
             self              %oneHaloTermMainBranch(i,n,:)  &
                  & =+     self%oneHaloTermMainBranch(i,n,:)  &
                  &  /dble(self%countMainBranch      (    :))
             self              %twoHaloTermMainBranch(i,n,:)  &
                  & =+     self%twoHaloTermMainBranch(i,n,:)  &
                  &  /dble(self%countMainBranch      (    :))
          end where
       end do
    end do
    ! Subtract out Poisson component of main branch galaxy variance (since these galaxies are not Poisson distributed).
    do m=1,self%massCount
       call correlationFunctionTermIndices(m,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
       do i=1,self%countBinsMassHalo
          ! Density-density.
          self                             %termCovariance       (  indexDensity                                    ,indexDensity                                    )=    &
               & +                     self%termCovariance       (  indexDensity                                    ,indexDensity                                    )     &
               & -                     self%meanDensityMainBranch(  m                                               ,i                                               ) **2 &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          ! One-halo-one-halo.
          self                             %termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexOneHalo:indexOneHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexOneHalo:indexOneHalo+self%wavenumberCount-1)     &
               & -Vector_Outer_Product(                                                                                                                                    &
               &                       self%oneHaloTermMainBranch(:,m                                               ,i                                               ),    &
               &                       self%oneHaloTermMainBranch(:,m                                               ,i                                               )     &
               &                      )                                                                                                                                    &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          ! Two-halo-two-halo.
          self                             %termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)     &
               & -Vector_Outer_Product(                                                                                                                                    &
               &                       self%twoHaloTermMainBranch(:,m                                               ,i                                               ),    &
               &                       self%twoHaloTermMainBranch(:,m                                               ,i                                               )     &
               &                      )                                                                                                                                    &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          ! Density-one-halo.
          self                             %termCovariance       (  indexDensity                                    ,indexOneHalo:indexOneHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexDensity                                    ,indexOneHalo:indexOneHalo+self%wavenumberCount-1)     &
               & -                     self%meanDensityMainBranch(  m                                               ,i                                               )     &
               & *                     self%oneHaloTermMainBranch(:,m                                               ,i                                               )     &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          self                             %termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexDensity                                    )=    &
               & +                     self%termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexDensity                                    )     &
               & -                     self%oneHaloTermMainBranch(:,m                                               ,i                                               )     &
               & *                     self%meanDensityMainBranch(  m                                               ,i                                               )     &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          ! Density-two-halo.
          self                             %termCovariance       (  indexDensity                                    ,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexDensity                                    ,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)     &
               & -                     self%meanDensityMainBranch(  m                                               ,i                                               )     &
               & *                     self%twoHaloTermMainBranch(:,m                                               ,i                                               )     &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          self                             %termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexDensity                                    )=    &
               & +                     self%termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexDensity                                    )     &
               & -                     self%twoHaloTermMainBranch(:,m                                               ,i                                               )     &
               & *                     self%meanDensityMainBranch(  m                                               ,i                                               )     &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          ! One-halo-two-halo
          self                             %termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexOneHalo:indexOneHalo+self%wavenumberCount-1,indexTwoHalo:indexTwoHalo+self%wavenumberCount-1)     &
               & -Vector_Outer_Product(                                                                                                                                    &
               &                       self%oneHaloTermMainBranch(:,m                                               ,i                                               ),    &
               &                       self%twoHaloTermMainBranch(:,m                                               ,i                                               )     &
               &                      )                                                                                                                                    &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
          self                             %termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexOneHalo:indexOneHalo+self%wavenumberCount-1)=    &
               & +                     self%termCovariance       (  indexTwoHalo:indexTwoHalo+self%wavenumberCount-1,indexOneHalo:indexOneHalo+self%wavenumberCount-1)     &
               & -Vector_Outer_Product(                                                                                                                                    &
               &                       self%twoHaloTermMainBranch(:,m                                               ,i                                               ),    &
               &                       self%oneHaloTermMainBranch(:,m                                               ,i                                               )     &
               &                      )                                                                                                                                    &
               & *                dble(self%countMainBranch      (                                                   i                                               ))
       end do
    end do
    ! Normalize one- and two-halo terms.
    allocate(jacobian            (self%massCount*(2*self%wavenumberCount),self%massCount*(2*self%wavenumberCount+1)))
    allocate(oneTwoHaloCovariance(self%massCount*(2*self%wavenumberCount),self%massCount*(2*self%wavenumberCount  )))
    ! One-halo term.
    jacobian=0.0d0
    do n=1,self%massCount
       call correlationFunctionTermIndices(n,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
       if (self%meanDensity(n) > 0.0d0) then
          do i=1,self%wavenumberCount
             jacobian((n-1)*(2*self%wavenumberCount)+i,indexOneHalo+i-1)=1.0d0/self%meanDensity(n)**2
          end do
          jacobian((n-1)*(2*self%wavenumberCount)+1:(n-1)*(2*self%wavenumberCount)+self%wavenumberCount,indexDensity)=-2.0d0*self%oneHaloTerm(:,n)/self%meanDensity(n)**3
       end if
    end do
    ! Two-halo term.
    do n=1,self%massCount
       call correlationFunctionTermIndices(n,self%wavenumberCount,indexOneHalo,indexTwoHalo,indexDensity)
       if (self%meanDensity(n) > 0.0d0) then
          do i=1,self%wavenumberCount
             jacobian((n-1)*(2*self%wavenumberCount)+self%wavenumberCount+i,indexTwoHalo+i-1)=1.0d0/self%meanDensity(n)
          end do
          jacobian((n-1)*(2*self%wavenumberCount)+self%wavenumberCount+1:(n-1)*(2*self%wavenumberCount)+2*self%wavenumberCount,indexDensity)=-self%twoHaloTerm(:,n)/self%meanDensity(n)**2
       end if
    end do
    jacobianMatrix                     =jacobian
    covarianceMatrix                   =self%termCovariance
    oneTwoHaloCovariance               =jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
    do n=1,self%massCount
       if (self%meanDensity(n) > 0.0d0) then
          self%oneHaloTerm(:,n)=self%oneHaloTerm(:,n)/self%meanDensity(n)**2
          self%twoHaloTerm(:,n)=self%twoHaloTerm(:,n)/self%meanDensity(n)
       end if
    end do
    deallocate(jacobian)
    ! Square the two halo term, and multiply by the linear theory power spectrum.
    allocate(jacobian            (self%massCount*(2*self%wavenumberCount),self%massCount*(2*self%wavenumberCount)))
    jacobian=0.0d0
    do n=1,self%massCount
       do i=1,self%wavenumberCount
          jacobian((n-1)*(2*self%wavenumberCount)                     +i,(n-1)*(2*self%wavenumberCount)                     +i)=+1.0d0
          jacobian((n-1)*(2*self%wavenumberCount)+self%wavenumberCount+i,(n-1)*(2*self%wavenumberCount)+self%wavenumberCount+i)=+2.0d0                 &
               &                                                                                                                *self%twoHaloTerm(i,n)
       end do
    end do
    jacobianMatrix      =jacobian
    covarianceMatrix    =oneTwoHaloCovariance
    oneTwoHaloCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
    do n=1,self%massCount
       do i=1,self%wavenumberCount
          self%twoHaloTerm(i,n)=+self%twoHaloTerm(i,n)**2
       end do
    end do
    deallocate(jacobian)
    ! Construct the final power spectra.
    allocate(powerSpectrumValue     (               self%wavenumberCount,self%massCount                         ))
    allocate(powerSpectrumCovariance(self%massCount*self%wavenumberCount,self%massCount*   self%wavenumberCount ))
    allocate(jacobian               (self%massCount*self%wavenumberCount,self%massCount*(2*self%wavenumberCount)))
    jacobian=0.0d0
    do n=1,self%massCount
       do i=1,self%wavenumberCount
          jacobian((n-1)*self%wavenumberCount+i,(n-1)*(2*self%wavenumberCount)                     +i)=1.0d0
          jacobian((n-1)*self%wavenumberCount+i,(n-1)*(2*self%wavenumberCount)+self%wavenumberCount+i)=1.0d0
       end do
    end do
    jacobianMatrix         =jacobian
    covarianceMatrix       =oneTwoHaloCovariance
    powerSpectrumCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
    do n=1,self%massCount
       powerSpectrumValue(:,n)=self%oneHaloTerm(:,n)+self%twoHaloTerm(:,n)
    end do
    deallocate(jacobian            )
    deallocate(oneTwoHaloCovariance)
    ! Allocate correlation function and separation arrays.
    allocate(correlation,mold=powerSpectrumValue)
    allocate(separation (self%wavenumberCount))
    ! Fourier transform the power spectrum to get the correlation function.
    do n=1,self%massCount
       call FFTLogSineTransform(                          &
            &                    self%wavenumber        , &
            &                    separation             , &
            &                   +powerSpectrumValue(:,n)  &
            &                   *self%wavenumber          &
            &                   * 4.0d0*Pi                &
            &                   /(2.0d0*Pi)**3          , &
            &                    correlation       (:,n), &
            &                    fftLogForward            &
            &                  )
       correlation(:,n)=correlation(:,n)/separation
    end do
    ! Compute the covariance of the correlation function.
    allocate(covarianceTmp        (self%massCount*self%wavenumberCount,self%massCount*self%wavenumberCount))
    allocate(correlationCovariance(self%massCount*self%wavenumberCount,self%massCount*self%wavenumberCount))
    ! Apply wavenumber weighting to the power spectrum covariance.
    do n=1,self%massCount
       do m=1,self%massCount
          do i=1,self%wavenumberCount
             do j=1,self%wavenumberCount
                powerSpectrumCovariance         ((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+j) &
                     & =+powerSpectrumCovariance((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+j) &
                     &  *self%wavenumber        (                           i                             ) &
                     &  *self%wavenumber        (                                                        j) &
                     &  *(                                                                                  &
                     &    + 4.0d0*Pi                                                                        &
                     &    /(2.0d0*Pi)**3                                                                    &
                     &   )**2
             end do
          end do
       end do
    end do
    ! Derive the covariance of the correlation function by first Fourier transforming each row of the power spectrum covariance
    ! matrix, and then Fourier transforming each column.
    do n=1,self%massCount
       do m=1,self%massCount
          do i=1,self%wavenumberCount
             call FFTLogSineTransform(                                                                                                           &
                  &                   self%wavenumber                                                                                          , &
                  &                   separation                                                                                               , &
                  &                   powerSpectrumCovariance((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+1:m*self%wavenumberCount), &
                  &                   covarianceTmp          ((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+1:m*self%wavenumberCount), &
                  &                   fftLogForward                                                                                              &
                  &                  )
          end do
       end do
    end do
    do n=1,self%massCount
       do m=1,self%massCount
          do i=1,self%wavenumberCount
             call FFTLogSineTransform(                                                                                                           &
                  &                   self%wavenumber                                                                                          , &
                  &                   separation                                                                                               , &
                  &                   covarianceTmp          ((n-1)*self%wavenumberCount+1:n*self%wavenumberCount,(m-1)*self%wavenumberCount+i), &
                  &                   correlationCovariance  ((n-1)*self%wavenumberCount+1:n*self%wavenumberCount,(m-1)*self%wavenumberCount+i), &
                  &                   fftLogForward                                                                                              &
                  &                  )
          end do
       end do
    end do
    do n=1,self%massCount
       do m=1,self%massCount
          do i=1,self%wavenumberCount
             do j=1,self%wavenumberCount
                correlationCovariance        ((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+j) &
                     & =correlationCovariance((n-1)*self%wavenumberCount+i,(m-1)*self%wavenumberCount+j) &
                     & /separation           (                           i                             ) &
                     & /separation           (                                                        j)
             end do
          end do
       end do
    end do
    deallocate(covarianceTmp)
    ! Construct correlation table.
    call correlationTable%create(separation(1),separation(self%wavenumberCount),size(separation),extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
    ! Project the correlation function.
    allocate(jacobian                      (self%massCount*self%wavenumberCount,self%massCount*self%wavenumberCount))
    allocate(projectedCorrelationCovariance(self%massCount*self%wavenumberCount,self%massCount*self%wavenumberCount))
    allocate(projectedCorrelation          (               self%wavenumberCount,self%massCount                     ))
    jacobian=0.0d0
    integrandWeightFunction => projectionIntegrandWeight
    do i=1,self%wavenumberCount
       projectedSeparation=correlationTable%x(int(i))
       jacobian(i,1:self%wavenumberCount)=correlationTable%integrationWeights(                                &
            &                                                                 projectedSeparation           , &
            &                                                                 sqrt(                           &
            &                                                                      +projectedSeparation  **2  &
            &                                                                      +self%depthLineOfSight**2  &
            &                                                                     )                         , &
            &                                                                 integrandWeightFunction         &
            &                                                                )
       do n=1,self%massCount
          if (n > 1) jacobian((n-1)*self%wavenumberCount+i,(n-1)*self%wavenumberCount+1:n*self%wavenumberCount)=jacobian(i,1:self%wavenumberCount)
          projectedCorrelation(i,n)=sum(jacobian(i,1:self%wavenumberCount)*correlation(:,n))
       end do
    end do
    jacobianMatrix                =jacobian
    covarianceMatrix              =correlationCovariance
    projectedCorrelationCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
    deallocate(jacobian)
    ! If the integral was taken over the half range, 0 < π < πₘₐₓ, rather than the full range, -πₘₐₓ < π < πₘₐₓ, then divide
    ! the projected correlation function by two.
    if (self%halfIntegral) then
       projectedCorrelation          =projectedCorrelation          /2.0d0
       projectedCorrelationCovariance=projectedCorrelationCovariance/2.0d0**2
    end if
    ! Integrate the projected correlation function over bins.
    allocate(self%binnedProjectedCorrelation          (               size(self%separations),self%massCount                           ))
    allocate(self%binnedProjectedCorrelationCovariance(self%massCount*size(self%separations),self%massCount*size(self%separations    )))
    allocate(     jacobian                            (self%massCount*size(self%separations),self%massCount*     self%wavenumberCount ))
    jacobian=0.0d0
    integrandWeightFunction => binningIntegrandWeight
    binWidthLogarithmic=log(self%separations(2)/self%separations(1))
    do i=1,size(self%separations)
       binSeparationMinimum         =self%separations(i)*exp(-0.5d0*binWidthLogarithmic)
       binSeparationMaximum         =self%separations(i)*exp(+0.5d0*binWidthLogarithmic)
       jacobian(i,1:self%wavenumberCount)=correlationTable%integrationWeights(                         &
            &                                                                 binSeparationMinimum   , &
            &                                                                 binSeparationMaximum   , &
            &                                                                 integrandWeightFunction  &
            &                                                                )                         &
            &                             /Pi                                                          &
            &                             /(                                                           &
            &                               +binSeparationMaximum**2                                   &
            &                               -binSeparationMinimum**2                                   &
            &                              )
       do n=1,self%massCount
          if (n > 1) jacobian((n-1)*size(self%separations)+i,(n-1)*self%wavenumberCount+1:n*self%wavenumberCount)=jacobian(i,1:self%wavenumberCount)
          self%binnedProjectedCorrelation(i,n)=sum(jacobian(i,1:self%wavenumberCount)*projectedCorrelation(:,n))
       end do
    end do
    jacobianMatrix                           =jacobian
    covarianceMatrix                         =projectedCorrelationCovariance
    self%binnedProjectedCorrelationCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
    deallocate(jacobian)
    call correlationTable%destroy()
    ! Apply the integral constraint.
    self%binnedProjectedCorrelation=self%binnedProjectedCorrelation/self%integralConstraint
    return

  contains

    double precision function projectionIntegrandWeight(separation)
      !!{
      The weight function applied to the correlation function when integrating to get the projected correlation function.
      !!}
      implicit none
      double precision, intent(in   ) :: separation

      if (separation > projectedSeparation) then
         projectionIntegrandWeight=2.0d0*separation/sqrt(separation**2-projectedSeparation**2)
      else
         projectionIntegrandWeight=0.0d0
      end if
      return
    end function projectionIntegrandWeight

    double precision function binningIntegrandWeight(separation)
      !!{
      The weight function applied to the projected correlation function when integrating into bins.
      !!}
      implicit none
      double precision, intent(in   ) :: separation

      binningIntegrandWeight=2.0d0*Pi*separation
      return
    end function binningIntegrandWeight

  end subroutine correlationFunctionFinalizeAnalysis

  double precision function correlationFunctionLogLikelihood(self)
    !!{
    Return the log-likelihood of a correlationFunction output analysis.
    !!}
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisCorrelationFunction), intent(inout)                 :: self
    double precision                                   , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                   , allocatable  , dimension(:  ) :: functionValueDifference
    type            (vector                           )                                :: residual
    type            (matrix                           )                                :: covariance
    integer                                                                            :: status

    ! Check for existence of a target distribution.
    if (allocated(self%binnedProjectedCorrelationTarget)) then
       ! Finalize analysis.
       call correlationFunctionFinalizeAnalysis(self)
       ! Allocate workspaces.
       allocate(functionCovarianceCombined(self%binCount*self%massCount,self%binCount*self%massCount))
       allocate(functionValueDifference   (self%binCount*self%massCount                             ))
       ! Find combined covariance and difference between model and target.
       functionValueDifference   =reshape(                                                   &
            &                              +self%binnedProjectedCorrelation                  &
            &                              -self%binnedProjectedCorrelationTarget          , &
            &                             [                                                  &
            &                              +self%binCount                                    &
            &                              *self%massCount                                   &
            &                             ]                                                  &
            &                            )
       functionCovarianceCombined=         +self%binnedProjectedCorrelationCovariance        &
            &                              +self%binnedProjectedCorrelationCovarianceTarget
       residual                  =vector(functionValueDifference   )
       covariance                =matrix(functionCovarianceCombined)
       ! Compute the log-likelihood.
       correlationFunctionLogLikelihood=-0.5d0*covariance%covarianceProduct     (residual,status) &
            &                           -0.5d0*covariance%logarithmicDeterminant(               ) &
            &                           -0.5d0*dble(self%binCount)*log(2.0d0*Pi)
       if (status /= GSL_Success) correlationFunctionLogLikelihood=logImprobable
    else
       correlationFunctionLogLikelihood=0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function correlationFunctionLogLikelihood
