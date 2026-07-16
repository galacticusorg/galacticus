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

!!{RST
Implements a generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.
!!}

  use               :: Galactic_Filters                        , only : galacticFilterClass
  use   , intrinsic :: ISO_C_Binding                           , only : c_size_t
  use               :: ISO_Varying_String                      , only : varying_string
  use               :: Node_Property_Extractors                , only : nodePropertyExtractorClass
  !$ use            :: Locks                                   , only : ompLock
  use               :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerClass
  use               :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorClass
  use               :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorClass
  use               :: Output_Analysis_Target_Data             , only : outputAnalysisTargetDataClass             , outputAnalysisTargetDataStandard
  use               :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorClass
  use               :: Output_Times                            , only : outputTimesClass
  use               :: Output_Analyses_Options                 , only : enumerationOutputAnalysisCovarianceModelType, enumerationOutputAnalysisStateType
  !![
  <outputAnalysis name="outputAnalysisVolumeFunction1D" docformat="rst">
   <description>
   A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.

   In addition to the volume function itself, the covariance matrix, :math:`\mathbf{C}_\mathrm{model}`, of the mass function is also computed. The assumptions used when constructing the covariance matrix are controlled by the parameter ``[covarianceModel]``. If set to ``binomial``, then to construct :math:`\mathbf{C}_\mathrm{model}` we make use of the fact that Galacticus works by sampling a set of tree "root masses" from the :math:`z=0` dark matter halo mass function. From each root, a tree is grown, within which the physics of galaxy formation is then solved. Root masses are sampled uniformly from the halo mass function. That is, the cumulative halo mass function, :math:`N(M)`, is constructed between the maximum and minimum halo masses to be simulated. The number of root masses, :math:`N_\mathrm{r}`, to be used in a model evaluation is then determined. Root masses are then chosen such that

   .. math::

      N(M_i) = N(M_\mathrm{min}) {i-1 \over N_\mathrm{r}-1}

   for :math:`i=1\ldots N_\mathrm{r}` (noting that :math:`N(M_\mathrm{max})=0` by construction).

   Consider first those galaxies which form in the main branch of each tree (i.e. those galaxies which are destined to become the central galaxy of the :math:`z=0` halo). Suppose that we simulate :math:`N_k` halos of root mass :math:`M_k` at :math:`z=0`. In such halos the main branch galaxies will, at any time, have property values drawn from some distribution :math:`p_k(M_\star|t)`. The number of such galaxies contributing to bin :math:`i` of the mass function is therefore binomially distributed with success probability :math:`p_{ik} = \int_{M_{i,\mathrm min}}^{M_{i,\mathrm max}} p_k(M_\star|t) \d M_\star` and a sample size of :math:`N_k`.

   Generalizing to consider all bins in our volume function, the number of galaxies in each bin will jointly follow a `multinomial distribution &lt;https://en.wikipedia.org/wiki/Multinomial_distribution&gt;`_. The contribution to the covariance matrix from these main branch galaxies is therefore:

   .. math::

      \mathcal{C}_{ij} = \left\{ \begin{array}{ll} p_{ik}(1-p_{ik}) N_k w_k^2 &amp; \hbox{ if } i = j \\ -p_{ik} p_{jk} N_k w_k^2 &amp; \hbox{ otherwise,} \end{array} \right.

   where :math:`w_k` is the weight to be assigned to each tree. To compute this covariance requires knowledge of the probabilities, :math:`p_{ik}`. We estimate these directly from the model. To do this, we bin trees into narrow bins of root mass and assume that :math:`p_{ik}` does not vary significantly across the mass range of each bin. Using all realizations of trees that fall within a given bin, :math:`k`, we can directly estimate :math:`p_{ik}`. Similarly, :math:`N_k w_k^2` is found by accumulating squared weights in bins of root mass. In computing :math:`p_{ik}` and :math:`N_k`, the range of halo masses considered and the fineness of binning in halo mass are determined by the parameters ``[covarianceBinomialMassHaloMinimum]``, ``[covarianceBinomialMassHaloMaximum]``, and ``[covarianceBinomialBinsPerDecade]``.

   If instead, ``[covarianceModel]``\ :math:`=`\ ``Poisson``, the main branch galaxies are modeled as being sampled from a Poisson distribution (and so off-diagonal terms in the covariance matrix will be zero).

   In addition to the main branch galaxies, each tree will contain a number of other galaxies (these will be "satellite" galaxies at :math:`z=0`, but at higher redshifts may still be central galaxies in their own halos). Tests have established that the number of satellites in halos is well described by a Poisson process. Note that, as described above, each galaxy contributes a Gaussian distribution to the mass function due to modeling of random errors in property value determinations. For main branch galaxies this is simply accounted for when accumulating the probabilities, :math:`p_{ik}`. For satellite galaxies, off-diagonal contributions to the covariance matrix arise as a result, :math:`C_{ij} = w_k f_i f_j`, where :math:`f_i` is the fraction of the galaxy contributing to bin :math:`i` of the mass function.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisVolumeFunction1D
     !!{RST
     A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.
     !!}
     private
     type            (varying_string                              )                              :: label                                          , comment                                          , &
          &                                                                                         propertyLabel                                  , propertyComment                                  , &
          &                                                                                         distributionLabel                              , distributionComment                              , &
          &                                                                                         propertyUnits                                  , distributionUnits                                , &
          &                                                                                         propertyQuantity                               , distributionQuantity                             , &
          &                                                                                         xAxisLabel                                     , yAxisLabel                                       , &
          &                                                                                         targetLabel                                    , reportLabel
     ! Axis labels, log-scale flags, and the (optional) target dataset are also bundled into a
     ! single `outputAnalysisTargetDataStandard` instance so the wrapper-pipeline doesn't have
     ! to enumerate 2^N presence combinations for the otherwise individually-optional fields.
     ! The shadow `xAxisLabel` / `yAxisLabel` / `targetLabel` / `xAxisIsLog` / `yAxisIsLog` /
     ! `functionValueTarget` / `functionCovarianceTarget1D` fields are kept on the outer type so
     ! the auto-built descriptor (which walks the type definition, not contained types) can
     ! reconstruct a parameter block that recreates this object.
     type            (outputAnalysisTargetDataStandard            )                              :: targetData_
     double precision                                                                            :: propertyUnitsInSI                              , distributionUnitsInSI
     logical                                                                                     :: propertyIsComoving                             , distributionIsComoving
     class           (nodePropertyExtractorClass                  ), pointer                     :: nodePropertyExtractor_                => null()
     class           (outputAnalysisPropertyOperatorClass         ), pointer                     :: outputAnalysisPropertyOperator_       => null()                                                   , &
          &                                                                                         outputAnalysisPropertyUnoperator_     => null()
     class           (outputAnalysisWeightOperatorClass           ), pointer                     :: outputAnalysisWeightOperator_         => null()
     class           (outputAnalysisDistributionOperatorClass     ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (outputAnalysisDistributionNormalizerClass   ), pointer                     :: outputAnalysisDistributionNormalizer_ => null()
     class           (galacticFilterClass                         ), pointer                     :: galacticFilter_                       => null()
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                          => null()
     double precision                                              , dimension(:,:), allocatable :: outputWeight                                   , functionCovariance                               , &
          &                                                                                         weightMainBranch
     double precision                                              , dimension(:  ), allocatable :: binCenter                                      , functionValue                                    , &
          &                                                                                         functionValueTarget                            , weightMainBranchSquared                          , &
          &                                                                                         functionCovarianceTarget1D
     double precision                                              , dimension(:  ), allocatable :: binMinimum                                     , binMaximum
     integer         (c_size_t                                    )                              :: binCount                                       , bufferCount                                      , &
          &                                                                                         binCountTotal                                  , covarianceModelBinomialBinCount                  , &
          &                                                                                         reportCount                                    , reportCountInRange
     type            (enumerationOutputAnalysisCovarianceModelType)                              :: covarianceModel
     type            (enumerationOutputAnalysisStateType          )                              :: state
     integer                                                                                     :: covarianceBinomialBinsPerDecade
     double precision                                                                            :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum                , &
          &                                                                                         covarianceModelHaloMassMinimumLogarithmic      , covarianceModelHaloMassIntervalLogarithmicInverse, &
          &                                                                                         binWidth
     logical                                                                                     :: finalized                                      , xAxisIsLog                                       , &
          &                                                                                         yAxisIsLog                                     , likelihoodNormalize                              , &
          &                                                                                         report
     !$ type         (ompLock                                     )                              :: accumulateLock
   contains
     !![
     <methods docformat="rst">
       <method description="Return the results of the volume function operator." method="results"           />
       <method description="Finalize the analysis of this function."             method="finalizeAnalysis"  />
       <method description="Activate/deactivate reporting."                      method="setReporting"      />
       <method description="Write the log-likelihood of this analysis to the output group. Child classes that compute their own log-likelihood should override this to avoid evaluating the parent-class logLikelihood method." method="logLikelihoodWrite"/>
       <method description="Write class-specific metadata to the analysis output group. The default implementation does nothing; child classes may override it to add further attributes or datasets." method="metadataWrite"/>
     </methods>
     !!]
     final     ::                     volumeFunction1DDestructor
     procedure :: analyze          => volumeFunction1DAnalyze
     procedure :: finalize         => volumeFunction1DFinalize
     procedure :: results          => volumeFunction1DResults
     procedure :: reduce           => volumeFunction1DReduce
     procedure :: logLikelihood    => volumeFunction1DLogLikelihood
     procedure :: finalizeAnalysis => volumeFunction1DFinalizeAnalysis
     procedure :: setReporting     => volumeFunction1DSetReporting
     procedure :: logLikelihoodWrite => volumeFunction1DLogLikelihoodWrite
     procedure :: metadataWrite      => volumeFunction1DMetadataWrite
  end type outputAnalysisVolumeFunction1D

  interface outputAnalysisVolumeFunction1D
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisVolumeFunction1D` output analysis class.
     !!}
     module procedure volumeFunction1DConstructorParameters
     module procedure volumeFunction1DConstructorInternal
  end interface outputAnalysisVolumeFunction1D

contains

  function volumeFunction1DConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisVolumeFunction1D` output analysis class which takes a parameter set as input.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter                                , inputParameters
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisCovarianceModelEncode
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass               ), pointer                     :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass        ), pointer                     :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    type            (outputAnalysisTargetDataStandard         )                              :: targetData_
    double precision                                           , dimension(:  ), allocatable :: binCenter                            , outputWeight                     , &
         &                                                                                      functionValueTarget                  , functionCovarianceTarget1D
    double precision                                           , dimension(:,:), allocatable :: functionCovarianceTarget
    integer         (c_size_t                                 )                              :: bufferCount
    type            (varying_string                           )                              :: label                                , comment                          , &
         &                                                                                      propertyLabel                        , propertyComment                  , &
         &                                                                                      distributionLabel                    , distributionComment              , &
         &                                                                                      propertyUnits                        , distributionUnits                , &
         &                                                                                      propertyQuantity                     , distributionQuantity             , &
         &                                                                                      covarianceModel                      , targetLabel                      , &
         &                                                                                      xAxisLabel                           , yAxisLabel
    integer                                                                                  :: covarianceBinomialBinsPerDecade
    type            (inputParameters                          )                              :: unoperatorParameters
    double precision                                                                         :: propertyUnitsInSI                    , distributionUnitsInSI            , &
         &                                                                                      covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum, &
         &                                                                                      binWidth
    logical                                                                                  :: xAxisIsLog                           , yAxisIsLog                       , &
         &                                                                                      propertyIsComoving                   , distributionIsComoving           , &
         &                                                                                      likelihoodNormalize

    ! Check and read parameters.
    unoperatorParameters=parameters%subParameters('unoperatorParameters',requireValue=.false.)
    !![
    <objectBuilder class="nodePropertyExtractor"                name="nodePropertyExtractor_"                source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"     source="unoperatorParameters"/>
    <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"         source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionNormalizer" name="outputAnalysisDistributionNormalizer_" source="parameters"          />
    <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"          />
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"          />
    !!]
    allocate(binCenter   (int(parameters%count('binCenter'))                          ))
    allocate(outputWeight(int(parameters%count('binCenter'))*self%outputTimes_%count()))
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*self%outputTimes_%count()) &
         & call Error_Report('incorrect number of output weights provided'//{introspection:location})
    !![
    <inputParameter docformat="rst">
      <name>label</name>
      <source>parameters</source>
      <variable>label</variable>
      <description>
      A label for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>xAxisLabel</name>
      <source>parameters</source>
      <description>
      A label for the :math:`x`-axis in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisLabel</name>
      <source>parameters</source>
      <description>
      A label for the :math:`y`-axis in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>xAxisIsLog</name>
      <source>parameters</source>
      <description>
      If true, indicates that the :math:`x`-axis should be logarithmic in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisIsLog</name>
      <source>parameters</source>
      <description>
      If true, indicates that the :math:`y`-axis should be logarithmic in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>comment</name>
      <source>parameters</source>
      <variable>comment</variable>
      <description>
      A descriptive comment for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyLabel</name>
      <source>parameters</source>
      <variable>propertyLabel</variable>
      <description>
      A label for the property variable.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyComment</name>
      <source>parameters</source>
      <variable>propertyComment</variable>
      <description>
      A descriptive comment for the property variable.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyUnits</name>
      <source>parameters</source>
      <variable>propertyUnits</variable>
      <description>
      A human-readable description of the units for the property.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyQuantity</name>
      <source>parameters</source>
      <variable>propertyQuantity</variable>
      <description>
      An ``astropy.units``-parseable units string for the property.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyIsComoving</name>
      <source>parameters</source>
      <variable>propertyIsComoving</variable>
      <description>
      If true, the property is in comoving units.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyUnitsInSI</name>
      <source>parameters</source>
      <variable>propertyUnitsInSI</variable>
      <description>
      A units for the property in the SI system.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionLabel</name>
      <source>parameters</source>
      <variable>distributionLabel</variable>
      <description>
      A label for the distribution.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionComment</name>
      <source>parameters</source>
      <variable>distributionComment</variable>
      <description>
      A descriptive comment for the distribution.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionUnits</name>
      <source>parameters</source>
      <variable>distributionUnits</variable>
      <description>
      A human-readable description of the units for the distribution.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionQuantity</name>
      <source>parameters</source>
      <variable>distributionQuantity</variable>
      <description>
      An ``astropy.units``-parseable units string for the distribution.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionIsComoving</name>
      <source>parameters</source>
      <variable>distributionIsComoving</variable>
      <description>
      If true, the distribution is in comoving units.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>distributionUnitsInSI</name>
      <source>parameters</source>
      <variable>distributionUnitsInSI</variable>
      <description>
      A units for the distribution in the SI system.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>
      The value of the property at the center of each bin.
      </description>
    </inputParameter>
    !!]
    if (size(binCenter) == 1) then
       !![
       <inputParameter docformat="rst">
	 <name>binWidth</name>
	 <source>parameters</source>
	 <variable>binWidth</variable>
	 <description>
	 The width of the bins.
	 </description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter docformat="rst">
      <name>bufferCount</name>
      <source>parameters</source>
      <variable>bufferCount</variable>
      <description>
      The number of buffer bins to include below and above the range of actual bins.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputWeight</name>
      <source>parameters</source>
      <variable>outputWeight</variable>
      <description>
      The weight to assign to each bin at each output.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceModel</name>
      <source>parameters</source>
      <variable>covarianceModel</variable>
      <description>
      The model to use for computing covariances.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>
      The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>
      The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true then normalize the likelihood to make it a probability density.
      </description>
    </inputParameter>
    !!]
    if (parameters%isPresent('functionValueTarget')) then
       if (parameters%isPresent('functionCovarianceTarget')) then
          !![
          <inputParameter docformat="rst">
            <name>functionValueTarget</name>
            <source>parameters</source>
            <description>
            The target function for likelihood calculations.
            </description>
          </inputParameter>
          <inputParameter docformat="rst">
            <name>functionCovarianceTarget</name>
            <source>parameters</source>
            <variable>functionCovarianceTarget1D</variable>
            <description>
            The target function covariance for likelihood calculations.
            </description>
          </inputParameter>
          !!]
          if (size(functionCovarianceTarget1D) == size(functionValueTarget)**2) then
             allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
             functionCovarianceTarget=reshape(functionCovarianceTarget1D,shape(functionCovarianceTarget))
          else
             call Error_Report('functionCovariance has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('functionCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('functionCovariance')) call Error_Report('functionTarget must be specified if functionCovariance is present'//{introspection:location})
    end if
    !![
    <inputParameter docformat="rst">
      <name>targetLabel</name>
      <source>parameters</source>
      <description>
      A label for the target dataset in a plot of this analysis.
      </description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    !!]
    ! Bundle the (potentially partial) target data into a single object for the internal constructor.
    targetData_=outputAnalysisTargetDataStandard(                                          &
         &                                       xAxisLabel      =xAxisLabel             , &
         &                                       yAxisLabel      =yAxisLabel             , &
         &                                       targetLabel     =targetLabel            , &
         &                                       xAxisIsLog      =xAxisIsLog             , &
         &                                       yAxisIsLog      =yAxisIsLog             , &
         &                                       valueTarget     =functionValueTarget    , &
         &                                       covarianceTarget=functionCovarianceTarget &
         &                                      )
    ! Build the object.
    !![
    <conditionalCall>
     <call>
      self=outputAnalysisVolumeFunction1D(                                                                                                    &amp;
           &amp;                          label                                                                                             , &amp;
           &amp;                          comment                                                                                           , &amp;
           &amp;                          propertyLabel                                                                                     , &amp;
           &amp;                          propertyComment                                                                                   , &amp;
           &amp;                          propertyUnits                                                                                     , &amp;
           &amp;                          propertyQuantity                                                                                  , &amp;
           &amp;                          propertyIsComoving                                                                                , &amp;
           &amp;                          propertyUnitsInSI                                                                                 , &amp;
           &amp;                          distributionLabel                                                                                 , &amp;
           &amp;                          distributionComment                                                                               , &amp;
           &amp;                          distributionUnits                                                                                 , &amp;
           &amp;                          distributionQuantity                                                                              , &amp;
           &amp;                          distributionIsComoving                                                                            , &amp;
           &amp;                          distributionUnitsInSI                                                                             , &amp;
           &amp;                          binCenter                                                                                         , &amp;
           &amp;                          bufferCount                                                                                       , &amp;
           &amp;                          reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),self%outputTimes_%count()]), &amp;
           &amp;                          nodePropertyExtractor_                                                                            , &amp;
           &amp;                          outputAnalysisPropertyOperator_                                                                   , &amp;
           &amp;                          outputAnalysisPropertyUnoperator_                                                                 , &amp;
           &amp;                          outputAnalysisWeightOperator_                                                                     , &amp;
           &amp;                          outputAnalysisDistributionOperator_                                                               , &amp;
           &amp;                          outputAnalysisDistributionNormalizer_                                                             , &amp;
           &amp;                          galacticFilter_                                                                                   , &amp;
           &amp;                          outputTimes_                                                                                      , &amp;
           &amp;                          enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.)      , &amp;
           &amp;                          covarianceBinomialBinsPerDecade                                                                   , &amp;
           &amp;                          covarianceBinomialMassHaloMinimum                                                                 , &amp;
           &amp;                          covarianceBinomialMassHaloMaximum                                                                 , &amp;
           &amp;                          likelihoodNormalize                                                                               , &amp;
           &amp;                          targetData_                                                                                         &amp;
           &amp;                          {conditions}                                                                                        &amp;
           &amp;                         )
     </call>
     <argument name="binWidth"                 value="binWidth"                        condition="size(binCenter) == 1"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"               />
    <objectDestructor name="outputAnalysisPropertyOperator_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="outputAnalysisWeightOperator_"        />
    <objectDestructor name="outputAnalysisDistributionOperator_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="galacticFilter_"                      />
    <objectDestructor name="outputTimes_"                         />
    !!]
    return
  end function volumeFunction1DConstructorParameters

  function volumeFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyQuantity,propertyIsComoving,propertyUnitsInSI,distributionLabel,distributionComment,distributionUnits,distributionQuantity,distributionIsComoving,distributionUnitsInSI,binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,outputAnalysisDistributionNormalizer_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,targetData_,binWidth) result (self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisVolumeFunction1D` output analysis class for internal use.
    !!}
    use :: Error                   , only : Error_Report
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass           , nodePropertyExtractorScalar
    use :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial
    implicit none
    type            (outputAnalysisVolumeFunction1D              )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                , comment                          , &
         &                                                                                                     propertyLabel                        , propertyComment                  , &
         &                                                                                                     distributionLabel                    , distributionComment              , &
         &                                                                                                     propertyUnits                        , distributionUnits                , &
         &                                                                                                     propertyQuantity                     , distributionQuantity
    logical                                                       , intent(in   ), optional                 :: likelihoodNormalize
    logical                                                       , intent(in   )                           :: propertyIsComoving                   , distributionIsComoving
    double precision                                              , intent(in   )                           :: propertyUnitsInSI                    , distributionUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    class           (nodePropertyExtractorClass                  ), intent(in   ), target                   :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(in   ), target                   :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(in   ), target                   :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass   ), intent(in   ), target                   :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                         ), intent(in   ), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(in   ), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum, &
         &                                                                                                     binWidth
    class           (outputAnalysisTargetDataClass               ), intent(in   ), optional                 :: targetData_
    integer         (c_size_t                                    )                                          :: i
    !![
    <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyQuantity, propertyIsComoving, propertyUnitsInSI, distributionLabel, distributionComment, distributionUnits, distributionQuantity, distributionIsComoving, distributionUnitsInSI, binCenter, bufferCount, outputWeight, *nodePropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *outputAnalysisDistributionNormalizer_, *galacticFilter_, *outputTimes_, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, binWidth"/>
    !!]

    ! Initialise the bundled target-data fields.  An explicit `targetData_` must be of the
    ! concrete `outputAnalysisTargetDataStandard` type (the only impl in the project today);
    ! without one we default-construct, matching the per-arg defaults the previous signature
    ! exposed.
    if (present(targetData_)) then
       select type (targetData_)
       type is (outputAnalysisTargetDataStandard)
          self%targetData_=targetData_
       class default
          call Error_Report('targetData_ must be of type outputAnalysisTargetDataStandard'//{introspection:location})
       end select
    else
       self%targetData_=outputAnalysisTargetDataStandard()
    end if
    ! Mirror the bundled target-data fields onto the outer object so the auto-built descriptor
    ! (which walks the type definition, not contained types) can reconstruct a parameter block
    ! that recreates this object.  Reshape the 2D covariance into the 1D form the parameter
    ! reader produces.
    self%xAxisLabel =self%targetData_%xAxisLabel
    self%yAxisLabel =self%targetData_%yAxisLabel
    self%targetLabel=self%targetData_%targetLabel
    self%xAxisIsLog =self%targetData_%xAxisIsLog
    self%yAxisIsLog =self%targetData_%yAxisIsLog
    if (self%targetData_%hasTarget()) then
       self%functionValueTarget       =self%targetData_%valueTarget
       self%functionCovarianceTarget1D=reshape(self%targetData_%covarianceTarget,[size(self%targetData_%covarianceTarget)])
    end if
    ! Validate.
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    ! Set normalization state for likelihood.
    self%likelihoodNormalize=.true.
    if (present(likelihoodNormalize)) self%likelihoodNormalize=likelihoodNormalize
    ! Count bins.
    self%binCount     =size(binCenter,kind=c_size_t)
    self%binCountTotal=self%binCount+2*bufferCount
    ! Determine bin minima and maxima. Allocate arrays for bins that include buffer regions.
    allocate(self%binMinimum(-bufferCount+1:-bufferCount+self%binCountTotal))
    allocate(self%binMaximum(-bufferCount+1:-bufferCount+self%binCountTotal))
    if (present(binWidth)) then
       self%binMinimum(1:self%binCount)=+binCenter-0.5d0*binWidth
       self%binMaximum(1:self%binCount)=+binCenter+0.5d0*binWidth
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*binWidth
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*binWidth
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*binWidth
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*binWidth
       end do
    else
       do i=1,self%binCount
          if (i == 1) then
             self%binMinimum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i+1))
          else
             self%binMinimum(i)=             +0.5d0*(self%binCenter(i-1)+self%binCenter(i  ))
          end if
          if (i == self%binCount) then
             self%binMaximum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i-1))
          else
             self%binMaximum(i)=             +0.5d0*(self%binCenter(i  )+self%binCenter(i+1))
          end if
       end do
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
       end do
    end if
    ! Allocate and initialize function values.
    allocate(self%functionValue     (self%binCount              ))
    allocate(self%functionCovariance(self%binCount,self%binCount))
    self%functionValue     =0.0d0
    self%functionCovariance=0.0d0
    ! Allocate and initialize binomial covariance model halo arrays if necessary.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       self%covarianceModelBinomialBinCount                  =int(log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)*dble(self%covarianceBinomialBinsPerDecade)+0.5d0)
       self%covarianceModelHaloMassMinimumLogarithmic        =log10(self%covarianceBinomialMassHaloMinimum)
       self%covarianceModelHaloMassIntervalLogarithmicInverse=dble(self%covarianceModelBinomialBinCount)/log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)
       allocate(self%weightMainBranch       (self%binCount,self%covarianceModelBinomialBinCount))
       allocate(self%weightMainBranchSquared(self%covarianceModelBinomialBinCount))
       self%weightMainBranch       =0.0d0
       self%weightMainBranchSquared=0.0d0
    end if
    ! Initialize finalization status.
    self%finalized=.false.
    ! Initialize OpenMP accumulation lock.
    !$ self%accumulateLock=ompLock()
    ! Initialize reporting state.
    self%report            =.false.
    self%reportLabel       ="unknown"
    self%reportCount       =0_c_size_t
    self%reportCountInRange=0_c_size_t
   return
  end function volumeFunction1DConstructorInternal

  subroutine volumeFunction1DDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`outputAnalysisVolumeFunction1D` output analysis class.
    !!}
    implicit none
    type(outputAnalysisVolumeFunction1D), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"               />
    <objectDestructor name="self%outputAnalysisPropertyOperator_"      />
    <objectDestructor name="self%outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="self%outputAnalysisWeightOperator_"        />
    <objectDestructor name="self%outputAnalysisDistributionOperator_"  />
    <objectDestructor name="self%outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="self%galacticFilter_"                      />
    <objectDestructor name="self%outputTimes_"                         />
    !!]
    return
  end subroutine volumeFunction1DDestructor

  subroutine volumeFunction1DAnalyze(self,node,iOutput)
    !!{RST
    Implement a volumeFunction1D output analysis.
    !!}
    use :: Display                 , only : displayMessage
    use :: Galacticus_Nodes        , only : nodeComponentBasic                   , treeNode
    use :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    use :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial, enumerationOutputAnalysisPropertyTypeType, enumerationOutputAnalysisPropertyQuantityType, outputAnalysisState, &
         &                                  enumerationOutputAnalysisStateDecode , enumerationOutputAnalysisStateType
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (outputAnalysisVolumeFunction1D               ), intent(inout)                 :: self
    type            (treeNode                                     ), intent(inout)                 :: node
    integer         (c_size_t                                     ), intent(in   )                 :: iOutput
    double precision                                               , allocatable  , dimension(:  ) :: distribution
    double precision                                               , allocatable  , dimension(:,:) :: covariance
    class           (nodeComponentBasic                           ), pointer                       :: basic
    double precision                                                                               :: propertyValue         , weightValue  , &
         &                                                                                            propertyValueIntrinsic
    type            (enumerationOutputAnalysisPropertyTypeType    )                                :: propertyType 
    type            (enumerationOutputAnalysisPropertyQuantityType)                                :: propertyQuantity
    type            (enumerationOutputAnalysisStateType           )                                :: state
    integer         (c_size_t                                     )                                :: j                     , k            , &
         &                                                                                            indexHaloMass

    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    ! Filter this node.
    if (.not.self%galacticFilter_%passes(node)) return
    ! Allocate work arrays.
    allocate(distribution(-self%bufferCount+1:self%binCount+self%bufferCount))
    ! Extract the property from the node.
    propertyType             =self%nodePropertyExtractor_%type    (    )
    propertyQuantity         =self%nodePropertyExtractor_%quantity(    )
    select type (extractor_ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       propertyValue         =                           extractor_%extract (node)
    class default
       propertyValue         =0.0d0
    end select
    propertyValueIntrinsic   =propertyValue
    ! Apply property operators.
    propertyValue=self%outputAnalysisPropertyOperator_%operate(propertyValue,node,propertyType,iOutput)
    ! Apply distribution operators.
    distribution=self%outputAnalysisDistributionOperator_%operateScalar(propertyValue,propertyType,self%binMinimum,self%binMaximum,iOutput,node)
    ! Compute the weight.
    weightValue=node%hostTree%volumeWeight
    ! Apply weight operators.
    weightValue=self%outputAnalysisWeightOperator_%operate(weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,iOutput)
    ! Apply weights.
    distribution(1:self%binCount)=+distribution              (1:self%binCount        ) &
         &                        *self        %outputWeight ( :             ,iOutput) &
         &                        *weightValue
    ! Accumulate the property, including weights from both the host tree and the output. Note that we accumulate only the
    ! non-buffer bins of the distribution.
    self%functionValue=+self%functionValue+distribution(1:self%binCount)
    ! Report on state changes and count in-range nodes.
    if (self%report) then
       self         %reportCount       =self%reportCount       +1_c_size_t
       if     (                                                            &
            &   propertyValue > self%binMinimum(     1       )             &
            &  .and.                                                       &
            &   propertyValue < self%binMaximum(self%binCount)             &
            & ) self%reportCountInRange=self%reportCountInRange+1_c_size_t
       ! Check for state changes.
       state=outputAnalysisState(self%functionValue)
       if (state /= self%state) then
          block
            type(varying_string) :: message
            message="report: "//self%reportLabel//": state change: "//enumerationOutputAnalysisStateDecode(self%state,includePrefix=.false.)//" → "//enumerationOutputAnalysisStateDecode(state,includePrefix=.false.)//" [node: "//node%uniqueID()//"]"
            call displayMessage(message)
          end block
          self%state=state
       end if
    end if
    ! Accumulate covariance. If using the binomial model for main branch galaxies, handle them separately.
    if (node%isOnMainBranch() .and. self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       ! Find the bin to which this halo mass belongs.
       basic         => node%basic()
       indexHaloMass =  floor((log10(basic%mass())-self%covarianceModelHaloMassMinimumLogarithmic)*self%covarianceModelHaloMassIntervalLogarithmicInverse)+1
       ! Accumulate weights to halo mass arrays.
       if (indexHaloMass >= 1 .and. indexHaloMass <= self%covarianceModelBinomialBinCount) then
          self                     %weightMainBranch       ( :             ,indexHaloMass)=    &
               &  +    self        %weightMainBranch       ( :             ,indexHaloMass)     &
               &  +    distribution                        (1:self%binCount)
          self                     %weightMainBranchSquared(                indexHaloMass)=    &
               &  +    self        %weightMainBranchSquared(                indexHaloMass)     &
               &  +sum(distribution                        (1:self%binCount              ))**2
       end if
    else
       ! Construct contribution to the covariance matrix assuming Poisson statistics.
       allocate(covariance(self%binCount,self%binCount))
       forall(j=1:self%binCount)
          forall(k=j:self%binCount)
             covariance(j,k)=+distribution(j  ) &
                  &          *distribution(  k)
             covariance(k,j)=+covariance  (j,k)
          end forall
       end forall
       ! Accumulate covariance.
       self        %functionCovariance= &
            & +self%functionCovariance  &
            & +             covariance
       deallocate(covariance)
    end if
    ! Deallocate workspace.
    deallocate(distribution)
    return
  end subroutine volumeFunction1DAnalyze

  subroutine volumeFunction1DReduce(self,reduced)
    !!{RST
    Implement a volumeFunction1D output analysis reduction.
    !!}
    use :: Display                , only : displayMessage                       , displayIndent , displayUnindent
    use :: Error                  , only : Error_Report
    use :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    use :: String_Handling        , only : operator(//)
    implicit none
    class    (outputAnalysisVolumeFunction1D), intent(inout) :: self
    class    (outputAnalysisClass           ), intent(inout) :: reduced
    type     (varying_string                )                :: message
    character(len=30                        )                :: label
    integer                                                  :: i
    
    select type (reduced)
    class is (outputAnalysisVolumeFunction1D)
       !$ call reduced%accumulateLock%set()
       if (self%report) call displayIndent('begin reduction lock: '//char(self%reportLabel))
       if (self%report) then
          message="report: "//self%reportLabel//": reduce: pre [nodes total/in-range = "//self%reportCount//" / "//self%reportCountInRange//"]"
          call displayIndent(message)
          call displayMessage("i    value        reduced     ")
          call displayMessage("---- ------------ ------------")
          do i=1,size(self%functionValue)
             write (label,'(i4,1x,e12.6,1x,e12.6)') i,self%functionValue(i),reduced%functionValue(i)
             call displayMessage(label)
          end do
          call displayUnindent("done")
       end if
       if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
          reduced%weightMainBranch       =reduced%weightMainBranch       +self%weightMainBranch
          reduced%weightMainBranchSquared=reduced%weightMainBranchSquared+self%weightMainBranchSquared
       end if
       reduced%functionCovariance        =reduced%functionCovariance     +self%functionCovariance
       reduced%functionValue             =reduced%functionValue          +self%functionValue
       !$ call reduced%accumulateLock%unset()
       if (self%report) then
          message="report: "//self%reportLabel//": reduce: post"
          call displayIndent(message)
          call displayMessage("i    value        reduced     ")
          call displayMessage("---- ------------ ------------")
          do i=1,size(self%functionValue)
             write (label,'(i4,1x,e12.6,1x,e12.6)') i,self%functionValue(i),reduced%functionValue(i)
             call displayMessage(label)
          end do
          call displayUnindent("done")
       end if
       if (self%report) call displayUnindent('end reduction lock: '//char(self%reportLabel))
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine volumeFunction1DReduce

  subroutine volumeFunction1DFinalize(self,groupName)
    !!{RST
    Implement a volumeFunction1D output analysis finalization.
    !!}
    use :: Output_HDF5   , only : outputFile
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5File, hdf5Group, hdf5Dataset
    use :: Units_MetaData, only : unitType
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout)           :: self
    type (varying_string                ), intent(in   ), optional :: groupName
    type (hdf5Group                     )               , target   :: analysesGroup, subGroup
    type (hdf5Group                     )               , pointer  :: inGroup
    type (hdf5Group                     )                          :: analysisGroup
    type (hdf5Dataset                   )                          :: dataset

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Output.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup =  inGroup      %openGroup(char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call    analysisGroup%writeAttribute(     char(self%             comment    )                ,'description'                                                                                                          )
    call    analysisGroup%writeAttribute("function1D"                                            ,'type'                                                                                                                 )
    call    analysisGroup%writeAttribute(     char(self%targetData_%xAxisLabel  )                ,'xAxisLabel'                                                                                                           )
    call    analysisGroup%writeAttribute(     char(self%targetData_%yAxisLabel  )                ,'yAxisLabel'                                                                                                           )
    call    analysisGroup%writeAttribute(          self%targetData_%xAxisIsLog                   ,'xAxisIsLog'                                                                                                           )
    call    analysisGroup%writeAttribute(          self%targetData_%yAxisIsLog                   ,'yAxisIsLog'                                                                                                           )
    call    analysisGroup%writeAttribute(     char(self%       propertyLabel    )                ,'xDataset'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%   distributionLabel    )                ,'yDataset'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%   distributionLabel    )//"Covariance"  ,'yCovariance'                                                                                                          )
    if (self%targetData_%hasTarget()) then
       call analysisGroup%writeAttribute(     char(self%   distributionLabel    )//"Target"            ,'yDatasetTarget'                                                                                                 )
       call analysisGroup%writeAttribute(     char(self%   distributionLabel    )//"CovarianceTarget"  ,'yCovarianceTarget'                                                                                              )
    end if
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (self%binCenter    (1:self%binCount                     ),char(self%    propertyLabel)                   ,char(self%   propertyComment)                  ,datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%propertyUnitsInSI       ,description=     char(self%propertyUnits    ),      quantity=     char(self%propertyQuantity    )       ,isComoving=self%propertyIsComoving    ),'units')
    call    analysisGroup%writeDataset  (self%functionValue(1:self%binCount                     ),char(self%distributionLabel)                   ,char(self%distributionComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%distributionUnitsInSI   ,description=     char(self%distributionUnits)      ,quantity=     char(self%distributionQuantity)       ,isComoving=self%distributionIsComoving),'units')
    call    analysisGroup%writeDataset  (self%functionCovariance(1:self%binCount,1:self%binCount),char(self%distributionLabel)//"Covariance"     ,char(self%distributionComment)//" [covariance]",datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%distributionUnitsInSI**2,description="["//char(self%distributionUnits)//"]²",quantity="("//char(self%distributionQuantity)//")^2",isComoving=self%distributionIsComoving),'units')
    ! Write the log-likelihood. This is delegated to the "logLikelihoodWrite" method so that child classes which compute their
    ! own log-likelihood can override it, and thereby avoid evaluating the (possibly invalid) parent-class "logLikelihood"
    ! method.
    call self%logLikelihoodWrite(analysisGroup)
    ! If available, include the target dataset.
    if (self%targetData_%hasTarget()) then
       call analysisGroup%writeAttribute(     char(self%targetData_%targetLabel)                     ,'targetLabel'                                                                                                      )
       call analysisGroup%writeDataset  (          self%targetData_%valueTarget                      ,char(self%distributionLabel)//"Target"          ,char(self%distributionComment)                 ,datasetReturned=dataset)
       call dataset      %writeAttribute(unitType(self%distributionUnitsInSI   ,description=     char(self%distributionUnits)      ,quantity=     char(self%distributionQuantity)       ,isComoving=self%distributionIsComoving),'units')
       call analysisGroup%writeDataset  (          self%targetData_%covarianceTarget                 ,char(self%distributionLabel)//"CovarianceTarget",char(self%distributionComment)//" [covariance]",datasetReturned=dataset)
       call dataset      %writeAttribute(unitType(self%distributionUnitsInSI**2,description="["//char(self%distributionUnits)//"]²",quantity="("//char(self%distributionQuantity)//")^2",isComoving=self%distributionIsComoving),'units')
    end if
    ! Write any class-specific metadata.
    call self%metadataWrite(analysisGroup)
    !$ call hdf5Access%unset()
    return
  end subroutine volumeFunction1DFinalize

  subroutine volumeFunction1DLogLikelihoodWrite(self,analysisGroup)
    !!{RST
    Write the log-likelihood of this analysis to the output group. This default implementation writes the log-likelihood returned by the :galacticus-class:`outputAnalysisVolumeFunction1D`  logLikelihood method whenever a target dataset is available. Child classes that compute their own log-likelihood (and which may not, e.g., initialize the covariance matrix used by the default  logLikelihood method) should override this method so that the parent-class  logLikelihood is never evaluated.
    !!}
    use :: IO_HDF5, only : hdf5File, hdf5Group
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout) :: self
    type (hdf5Group                     ), intent(inout) :: analysisGroup

    if (self%targetData_%hasTarget()) &
         & call analysisGroup%writeAttribute(self%logLikelihood(),'logLikelihood')
    return
  end subroutine volumeFunction1DLogLikelihoodWrite

  subroutine volumeFunction1DMetadataWrite(self,analysisGroup)
    !!{RST
    Write class-specific metadata to the analysis output group. This default implementation does nothing; child classes may override it to add further attributes or datasets to the analysis group.
    !!}
    use :: IO_HDF5, only : hdf5File, hdf5Group
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout) :: self
    type (hdf5Group                     ), intent(inout) :: analysisGroup
    !$GLC attributes unused :: self, analysisGroup

    return
  end subroutine volumeFunction1DMetadataWrite

  subroutine volumeFunction1DFinalizeAnalysis(self)
    !!{RST
    Compute final covariances and normalize.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities          , only : mpiSelf
#endif
    use :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout) :: self
    integer         (c_size_t                      )                :: i                    , j, &
         &                                                             m
    double precision                                                :: weightMainBranchTotal

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! If using the binomial covariance model, add on the contribution to covariance from main branch galaxies.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       do m=1,self%covarianceModelBinomialBinCount
          weightMainBranchTotal=sum(self%weightMainBranch(:,m))
          if (weightMainBranchTotal > 0.0d0) then
             do i=1,self%binCount
                ! For on-diagonal terms the variance is just n pᵢ (1-pᵢ). Here, pᵢ is the probability of contributing to bin i in
                ! our function. We find this by dividing the weight of main branch galaxies in bin i by the sum of weights over
                ! all bins. Finally, n is just the sum of squared weights for this halo mass bin.
                self               %functionCovariance     (i,i)=                       &
                     &         self%functionCovariance     (i,i)                        &
                     & +       self%weightMainBranchSquared(  m)                        & ! +   n
                     & *       self%weightMainBranch       (i,m)/weightMainBranchTotal  & !     pᵢ
                     & *(1.0d0-self%weightMainBranch       (i,m)/weightMainBranchTotal)   !  (1-pᵢ)
                do j=1,self%binCount
                   if (i == j) cycle
                   ! For off-diagonal terms the covariance is - n pᵢ pⱼ.
                   self            %functionCovariance     (i,j)=                       &
                        &  +   self%functionCovariance     (i,j)                        &
                        &  -   self%weightMainBranchSquared(  m)                        & ! -   n
                        &  *   self%weightMainBranch       (i,m)/weightMainBranchTotal  & !     pᵢ
                        &  *   self%weightMainBranch       (j,m)/weightMainBranchTotal    !     pⱼ
                end do
             end do
          end if
       end do
    end if
#ifdef USEMPI
    ! If running under MPI, perform a summation reduction across all processes.
    self%functionValue     =mpiSelf%sum(self%functionValue     )
    self%functionCovariance=mpiSelf%sum(self%functionCovariance)
#endif
    ! Apply final distribution operators - pass only the non-buffer bin values here.
    call self%outputAnalysisDistributionNormalizer_%normalize(self%functionValue,self%functionCovariance,self%binMinimum(1:self%binCount),self%binMaximum(1:self%binCount))
    ! Apply any "unoperator" to output bin values. This can be used to reverse transformations (e.g. if masses were converted to
    ! log10 for analysis, that can be undone here).
    do i=1,self%binCount
       self%binCenter(i)=self%outputAnalysisPropertyUnoperator_%operate(self%binCenter(i))
    end do
    return
  end subroutine volumeFunction1DFinalizeAnalysis

  subroutine volumeFunction1DResults(self,binCenter,functionValue,functionCovariance)
    !!{RST
    Implement a volumeFunction1D output analysis finalization.
    !!}
    implicit none
    class           (outputAnalysisVolumeFunction1D)                             , intent(inout)           :: self
    double precision                                , allocatable, dimension(:  ), intent(inout), optional :: binCenter         , functionValue
    double precision                                , allocatable, dimension(:,:), intent(inout), optional :: functionCovariance

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Return results.
    if (present(binCenter         )) then
       if (allocated(binCenter         )) deallocate(binCenter         )
       allocate(binCenter(size(self%binCenter)))
       binCenter         =self%binCenter
    end if
    if (present(functionValue     )) then
       if (allocated(functionValue     )) deallocate(functionValue     )
       allocate(functionValue(size(self%functionValue)))
       functionValue     =self%functionValue
    end if
    if (present(functionCovariance)) then
       if (allocated(functionCovariance)) deallocate(functionCovariance)
        allocate(functionCovariance(size(self%functionCovariance,dim=1),size(self%functionCovariance,dim=2)))
       !!]
       functionCovariance=self%functionCovariance
    end if
    return
  end subroutine volumeFunction1DResults

  double precision function volumeFunction1DLogLikelihood(self)
    !!{RST
    Return the log-likelihood of a volumeFunction1D output analysis.
    !!}
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Models_Likelihoods_Constants, only : logImprobable
    use :: Error                       , only : Error_Report
    use :: Interface_GSL               , only : GSL_Success
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout)                 :: self
    double precision                                , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                , allocatable  , dimension(:  ) :: functionValueDifference
    type            (vector                        )                                :: residual
    type            (matrix                        )                                :: covariance
    integer                                                                         :: status

    ! Check for existence of a target distribution.
    if (self%targetData_%hasTarget()) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! If model has everywhere zero return an improbable likelihood.
       if (all(self%functionValue == 0.0d0) .and. all(self%functionCovariance == 0.0d0)) then
          volumeFunction1DLogLikelihood=logImprobable
       else
          ! Allocate workspaces.
          allocate(functionCovarianceCombined(self%binCount,self%binCount))
          allocate(functionValueDifference   (self%binCount              ))
          ! Find combined covariance and difference between model and target.
          functionValueDifference   =+self%functionValue                  &
               &                     -self%targetData_%valueTarget
          functionCovarianceCombined=+self%functionCovariance              &
               &                     +self%targetData_%covarianceTarget
          residual                  = vector(functionValueDifference   )
          covariance                = matrix(functionCovarianceCombined)
          ! Compute the log-likelihood.
          volumeFunction1DLogLikelihood          =-0.5d0*covariance%covarianceProduct(residual,status)
          if (status == GSL_Success) then
             if (self%likelihoodNormalize)                                                   &
                  & volumeFunction1DLogLikelihood=+volumeFunction1DLogLikelihood             &
                  &                               -0.5d0*covariance%logarithmicDeterminant() &
                  &                               -0.5d0*dble(self%binCount)                 &
                  &                               *log(2.0d0*Pi)
          else
             volumeFunction1DLogLikelihood       =+logImprobable
          end if
       end if
    else
       volumeFunction1DLogLikelihood=0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function volumeFunction1DLogLikelihood

  subroutine volumeFunction1DSetReporting(self,report,reportLabel)
    !!{RST
    Activate/deactivate reporting.
    !!}
    use :: ISO_Varying_String     , only : assignment(=)
    use :: Output_Analyses_Options, only : outputAnalysisStateUnknown
    implicit none
    class    (outputAnalysisVolumeFunction1D), intent(inout)           :: self
    logical                                  , intent(in   )           :: report
    character(len=*                         ), intent(in   ), optional :: reportLabel

    self%report=report
    if (present(reportLabel)) self%reportLabel=reportLabel
    if (        report      ) self%state      =outputAnalysisStateUnknown
    return
  end subroutine volumeFunction1DSetReporting
