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
  Implements an output analysis class that computes subhalo mean maximum velocity as a function of mass.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <outputAnalysis name="outputAnalysisSubhaloVMaxVsMass">
   <description>An output analysis class that computes subhalo mean maximum velocity as a function of mass.</description>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisSubhaloVMaxVsMass
     !!{
     An output analysis class that computes subhalo mean maximum velocity as a function of mass.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null(), virialDensityContrastDefinition_ => null()
     class           (darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
     double precision                                      :: massMinimum                     , massMaximum                               , &
          &                                                   redshift
     integer         (c_size_t                  )          :: countMasses
     type            (varying_string            )          :: fileName
   contains
     final     ::                  subhaloVMaxVsMassDestructor
     procedure :: logLikelihood => subhaloVMaxVsMassLogLikelihood
  end type outputAnalysisSubhaloVMaxVsMass

  interface outputAnalysisSubhaloVMaxVsMass
     !!{
     Constructors for the \refClass{outputAnalysisSubhaloVMaxVsMass} output analysis class.
     !!}
     module procedure subhaloVMaxVsMassConstructorParameters
     module procedure subhaloVMaxVsMassConstructorFile
     module procedure subhaloVMaxVsMassConstructorInternal
  end interface outputAnalysisSubhaloVMaxVsMass

contains

  function subhaloVMaxVsMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisSubhaloVMaxVsMass} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters        , only : inputParameter            , inputParameters
    use :: Output_Times            , only : outputTimesClass
    use :: Cosmology_Functions     , only : cosmologyFunctionsClass
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
    implicit none
    type            (outputAnalysisSubhaloVMaxVsMass)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass               ), pointer       :: outputTimes_
    class           (virialDensityContrastClass     ), pointer       :: virialDensityContrast_ , virialDensityContrastDefinition_
    class           (darkMatterProfileDMOClass      ), pointer       :: darkMatterProfileDMO_
    double precision                                                 :: massMinimum            , massMaximum                     , &
         &                                                              redshift
    integer         (c_size_t                       )                :: countMasses
    type            (varying_string                 )                :: fileName

    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read the target dataset.</description>
       </inputParameter>
       !!]
    else
       !![
       <inputParameter>
         <name>massMinimum</name>
         <source>parameters</source>
         <defaultValue>1.0d6</defaultValue>
         <description>The minimum mass to consider.</description>
       </inputParameter>
       <inputParameter>
         <name>massMaximum</name>
         <source>parameters</source>
         <defaultValue>1.0d12</defaultValue>
         <description>The maximum mass to consider.</description>
       </inputParameter>
       <inputParameter>
         <name>countMasses</name>
         <source>parameters</source>
         <defaultValue>12_c_size_t</defaultValue>
         <description>The number of bins in mass to use.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which to compute the subhalo $V_\mathrm{max}$--$M$ relation.</description>
    </inputParameter>
    <objectBuilder class="outputTimes"           name="outputTimes_"                     source="parameters"                                                />
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <conditionalCall>
        <call>self=outputAnalysisSubhaloVMaxVsMass(outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,fileName{conditions})</call>
         <argument name="redshift" value="redshift" parameterPresent="parameters"/>
       </conditionalCall>
       !!]
    else
       self=outputAnalysisSubhaloVMaxVsMass(outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),massMinimum,massMaximum,countMasses)
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function subhaloVMaxVsMassConstructorParameters
  
  function subhaloVMaxVsMassConstructorFile(outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,fileName,redshift) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSubhaloVMaxVsMass} output analysis class for internal use.
    !!}
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: Output_Times            , only : outputTimesClass
    use :: Cosmology_Functions     , only : cosmologyFunctionsClass
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
    implicit none
    type            (outputAnalysisSubhaloVMaxVsMass)                                :: self
    type            (varying_string                 ), intent(in   )                 :: fileName
    class           (outputTimesClass               ), intent(inout)                 :: outputTimes_
    class           (virialDensityContrastClass     ), intent(in   )                 :: virialDensityContrast_  , virialDensityContrastDefinition_
    class           (cosmologyParametersClass       ), intent(inout)                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), intent(inout), target         :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass      ), intent(in   )                 :: darkMatterProfileDMO_
    double precision                                 , intent(in   ), optional       :: redshift
    double precision                                 , allocatable  , dimension(:  ) :: massesTarget            , functionTarget                  , &
         &                                                                              functionErrorTarget
    double precision                                 , allocatable  , dimension(:,:) :: functionCovarianceTarget
    double precision                                                                 :: massMinimum             , massMaximum                     , &
         &                                                                              time                    , redshift_
    integer         (c_size_t                         )                              :: countMasses             , i
    type            (varying_string                   )                              :: labelTarget
    type            (hdf5Object                       )                              :: file                    , velocityMaximumVsMassGroup
    !![
    <constructorAssign variables="redshift"/>
    !!]

    ! Read properties from the file.
    !$ call hdf5Access%set()
    call file                      %openFile     (fileName                  ,readOnly=.true.             )
    call file                      %readAttribute('label'                   ,         labelTarget        )
    call file                      %readAttribute('redshift'                ,         redshift_          )
    velocityMaximumVsMassGroup=file%openGroup('velocityMaximum')
    call velocityMaximumVsMassGroup%readDataset  ('mass'                    ,         massesTarget       )
    call velocityMaximumVsMassGroup%readDataset  ('velocityMaximumMean'     ,         functionTarget     )
    call velocityMaximumVsMassGroup%readDataset  ('velocityMaximumMeanError',         functionErrorTarget)
    call velocityMaximumVsMassGroup%close        (                                                       )
    call file                      %close        (                                                       )
    !$ call hdf5Access%unset()
    ! Override the redshift if one is provided.
    if (present(redshift)) redshift_=redshift
    ! Construct the velocity maximum function.
    time       =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift_))
    countMasses=size(massesTarget)
    massMinimum=massesTarget(          1)
    massMaximum=massesTarget(countMasses)
    allocate(functionCovarianceTarget(countMasses,countMasses))
    functionCovarianceTarget=0.0d0
    do i=1_c_size_t,countMasses
       functionCovarianceTarget(i,i)=functionErrorTarget(i)**2
    end do
    self=outputAnalysisSubhaloVMaxVsMass(outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,time,massMinimum,massMaximum,countMasses,functionTarget,functionCovarianceTarget,labelTarget)
    return
  end function subhaloVMaxVsMassConstructorFile

  function subhaloVMaxVsMassConstructorInternal(outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,time,massMinimum,massMaximum,countMasses,functionTarget,functionCovarianceTarget,labelTarget) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSubhaloVMaxVsMass} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                      , only : galacticFilterHaloIsolated                , galacticFilterHaloNotIsolated      , galacticFilterLowPass                 , galacticFilterAll                , &
         &                                                filterList
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassBound            , nodePropertyExtractorRadiusOrbital , nodePropertyExtractorRatio            , nodePropertyExtractorRadiusVirial, &
         &                                                nodePropertyExtractorVelocityMaximum      , nodePropertyExtractorHostNode
    use :: Numerical_Comparison                  , only : Values_Agree
    use :: Numerical_Constants_Prefixes          , only : kilo
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Numerical_Ranges                      , only : Make_Range                                , rangeTypeLinear
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10   , outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorIdentity
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                          , only : outputTimesClass
    use :: Virial_Density_Contrast               , only : virialDensityContrastClass
    use :: Dark_Matter_Profiles_DMO              , only : darkMatterProfileDMOClass
    implicit none
    type            (outputAnalysisSubhaloVMaxVsMass           )                                          :: self
    double precision                                            , intent(in   )                           :: massMinimum                          , massMaximum                     , &
         &                                                                                                   time
    integer         (c_size_t                                  ), intent(in   )                           :: countMasses
    class           (outputTimesClass                          ), intent(inout), target                   :: outputTimes_
    class           (virialDensityContrastClass                ), intent(in   ), target                   :: virialDensityContrast_               , virialDensityContrastDefinition_
    class           (cosmologyParametersClass                  ), intent(in   ), target                   :: cosmologyParameters_
    class           (cosmologyFunctionsClass                   ), intent(in   ), target                   :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                 ), intent(in   ), target                   :: darkMatterProfileDMO_
    double precision                                            , intent(in   ), dimension(:  ), optional :: functionTarget
    double precision                                            , intent(in   ), dimension(:,:), optional :: functionCovarianceTarget
    type            (varying_string                            ), intent(in   )                , optional :: labelTarget
    type            (nodePropertyExtractorMassBound            )               , pointer                  :: nodePropertyExtractorMassBound_
    type            (nodePropertyExtractorRadiusOrbital        )               , pointer                  :: nodePropertyExtractorRadiusOrbital_
    type            (nodePropertyExtractorRadiusVirial         )               , pointer                  :: nodePropertyExtractorRadiusVirial_
    type            (nodePropertyExtractorRatio                )               , pointer                  :: nodePropertyExtractorRadiusFractional_
    type            (nodePropertyExtractorVelocityMaximum      )               , pointer                  :: nodeWeightPropertyExtractor_
    type            (nodePropertyExtractorHostNode             )               , pointer                  :: nodePropertyExtractorRadiusVirialHost_
    type            (outputAnalysisPropertyOperatorIdentity    )               , pointer                  :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10       )               , pointer                  :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10   )               , pointer                  :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling   )               , pointer                  :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionOperatorIdentity)               , pointer                  :: outputAnalysisDistributionOperator_
    type            (galacticFilterHaloNotIsolated             )               , pointer                  :: galacticFilterIsSubhalo_
    type            (galacticFilterLowPass                     ), pointer                                 :: galacticFilterVirialRadius_
    type            (galacticFilterAll                         ), pointer                                 :: galacticFilterSubhalos_
    type            (filterList                                ), pointer                                 :: filters_
    double precision                                            , allocatable  , dimension(:  )           :: masses
    double precision                                            , allocatable  , dimension(:,:)           :: outputWeight
    integer         (c_size_t                                  )                                          :: i
    !![
    <constructorAssign variables="*cosmologyFunctions_, *virialDensityContrastDefinition_, *cosmologyParameters_, *virialDensityContrast_, *darkMatterProfileDMO_, massMinimum, massMaximum, countMasses"/>
    !!]

    ! Initialize.
    self%redshift =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Construct mass bins.
    allocate(masses(countMasses))
    masses=Make_Range(log10(massMinimum),log10(massMaximum),int(countMasses),rangeTypeLinear)
    ! Create a mass ratio property extractor.
    allocate(nodePropertyExtractorMassBound_       )
    allocate(nodePropertyExtractorRadiusOrbital_   )
    allocate(nodePropertyExtractorRadiusVirial_    )
    allocate(nodePropertyExtractorRadiusVirialHost_)
    allocate(nodePropertyExtractorRadiusFractional_)
    allocate(nodeWeightPropertyExtractor_          )
    !![
    <referenceConstruct object="nodePropertyExtractorMassBound_"        constructor="nodePropertyExtractorMassBound      (                                                                                                                                                   )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusOrbital_"    constructor="nodePropertyExtractorRadiusOrbital  (                                                                                                                                                   )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusVirial_"     constructor="nodePropertyExtractorRadiusVirial   (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                     )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusVirialHost_" constructor="nodePropertyExtractorHostNode       (nodePropertyExtractorRadiusVirial_                                                                                                                 )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusFractional_" constructor="nodePropertyExtractorRatio          ('radiusFraction','Ratio of subhalo orbital radius to host virial radius',nodePropertyExtractorRadiusOrbital_,nodePropertyExtractorRadiusVirialHost_)"/>
    <referenceConstruct object="nodeWeightPropertyExtractor_"           constructor="nodePropertyExtractorVelocityMaximum(var_str('VelocityMaximum'),darkMatterProfileDMO_                                                                                                   )"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass.
    allocate(outputAnalysisPropertyOperator_  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"       constructor="outputAnalysisPropertyOperatorLog10       (                                            )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_)
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"     constructor="outputAnalysisPropertyOperatorAntiLog10   (                                            )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_" constructor="outputAnalysisPropertyOperatorIdentity    (                                            )"/>
    !!]
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"         constructor="outputAnalysisWeightOperatorSubsampling   (                                            )"/>
    !!]
    ! Build filters which select subhalos/hosts.
    allocate(galacticFilterIsSubhalo_)
    !![
    <referenceConstruct object="galacticFilterIsSubhalo_"              constructor="galacticFilterHaloNotIsolated             (                                            )"/>
    !!]
    allocate(galacticFilterVirialRadius_)
    !![
    <referenceConstruct object="galacticFilterVirialRadius_"           constructor="galacticFilterLowPass                     (1.0d0,nodePropertyExtractorRadiusFractional_)"/>
    !!]
    allocate(galacticFilterSubhalos_     )
    allocate(filters_                    )
    allocate(filters_               %next)
    filters_     %filter_ => galacticFilterIsSubhalo_
    filters_%next%filter_ => galacticFilterVirialRadius_
    !![
    <referenceConstruct object="galacticFilterSubhalos_"               constructor="galacticFilterAll                         (filters_                                    )"/>
    !!]
    ! Build an identity distribution operator.
    allocate(outputAnalysisDistributionOperator_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"   constructor="outputAnalysisDistributionOperatorIdentity(                                            )"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countMasses,outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(i),time,absTol=1.0d-6)) then
          outputWeight(:,i)=1.0d0
       else
          outputWeight(:,i)=0.0d0
       end if
    end do
    ! Construct the mean function analyzer.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                                       &
         &                                                                              var_str('subhaloVelocityMaximumMean'                           ), &
         &                                                                              var_str('Subhalo mean maximum velocity vs. bound mass relation'), &
         &                                                                              var_str('massBound'                                            ), &
         &                                                                              var_str('Halo bound mass'                                      ), &
         &                                                                              var_str('M☉'                                                  ), &
         &                                                                              massSolar                                                       , &
         &                                                                              var_str('velocityMaximumMean'                                  ), &
         &                                                                              var_str('Mean velocity maximum'                                ), &
         &                                                                              var_str('km/s'                                                 ), &
         &                                                                              kilo                                                            , &
         &                                                                              masses                                                          , &
         &                                                                              0_c_size_t                                                      , &
         &                                                                              outputWeight                                                    , &
         &                                                                              nodePropertyExtractorMassBound_                                 , &
         &                                                                              nodeWeightPropertyExtractor_                                    , &
         &                                                                              outputAnalysisPropertyOperator_                                 , &
         &                                                                              outputAnalysisWeightPropertyOperator_                           , &
         &                                                                              outputAnalysisPropertyUnoperator_                               , &
         &                                                                              outputAnalysisWeightOperator_                                   , &
         &                                                                              outputAnalysisDistributionOperator_                             , &
         &                                                                              galacticFilterSubhalos_                                         , &
         &                                                                              outputTimes_                                                    , &
         &                                                                              outputAnalysisCovarianceModelPoisson                            , &
         &                                                         likelihoodNormalize =.false.                                                         , &
         &                                                         xAxisLabel          =var_str('$M_\mathrm{bound}/\mathrm{M}_\odot$'                  ), &
         &                                                         yAxisLabel          =var_str('$\langle V_\mathrm{max} \rangle / \hbox{km s}^{-1}$'  ), &
         &                                                         xAxisIsLog          =.true.                                                          , &
         &                                                         yAxisIsLog          =.true.                                                          , &
         &                                                         targetLabel         =labelTarget                                                     , &
         &                                                         meanValueTarget     =functionTarget                                                  , &
         &                                                         meanCovarianceTarget=functionCovarianceTarget                                          &
         &                                                        )
    !![
    <objectDestructor name="nodePropertyExtractorMassBound_"       />
    <objectDestructor name="nodePropertyExtractorRadiusOrbital_"   />
    <objectDestructor name="nodePropertyExtractorRadiusVirial_"    />
    <objectDestructor name="nodePropertyExtractorRadiusVirialHost_"/>
    <objectDestructor name="nodePropertyExtractorRadiusFractional_"/>
    <objectDestructor name="nodeWeightPropertyExtractor_"          />
    <objectDestructor name="outputAnalysisPropertyOperator_"       />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_" />
    <objectDestructor name="outputAnalysisWeightOperator_"         />
    <objectDestructor name="outputAnalysisDistributionOperator_"   />
    <objectDestructor name="galacticFilterSubhalos_"               />
    <objectDestructor name="galacticFilterIsSubhalo_"              />
    <objectDestructor name="galacticFilterVirialRadius_"           />
    !!]
    nullify(filters_)
    return
  end function subhaloVMaxVsMassConstructorInternal

  subroutine subhaloVMaxVsMassDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSubhaloVMaxVsMass} output analysis class.
    !!}
    implicit none
    type(outputAnalysisSubhaloVMaxVsMass), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"                    />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine subhaloVMaxVsMassDestructor

  double precision function subhaloVMaxVsMassLogLikelihood(self) result(logLikelihood)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily outputAnalysisSubhaloVMaxVsMass} output analysis.
    !!}
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisSubhaloVMaxVsMass), intent(inout)                 :: self
    double precision                                 , allocatable  , dimension(:,:) :: velocityMeanCovariance
    double precision                                 , allocatable  , dimension(:  ) :: velocityMeanDifference
    type            (vector                         )                                :: residual
    type            (matrix                         )                                :: covariance
    integer                                                                          :: i                     , j     , &
         &                                                                              countNonZero          , status

    ! Count the number of non-zero bins.
    countNonZero=0
    do i=1,size(self%meanValueTarget)
       if     (                                         &
            &   self%meanValueTarget     (i  ) <= 0.0d0 &
            &  .or.                                     &
            &   self%meanCovarianceTarget(i,i) <= 0.0d0 &
            & ) cycle
       countNonZero=countNonZero+1
    end do
    allocate(velocityMeanDifference(countNonZero             ))
    allocate(velocityMeanCovariance(countNonZero,countNonZero))
    ! Populate reduced bins.
    velocityMeanDifference=0.0d0
    velocityMeanCovariance=0.0d0
    j=0
    do i=1,size(self%meanValueTarget)
       if     (                                         &
            &   self%meanValueTarget     (i  ) <= 0.0d0 &
            &  .or.                                     &
            &   self%meanCovarianceTarget(i,i) <= 0.0d0 &
            & ) cycle
       j=j+1
       velocityMeanDifference(j  )=self%meanValue     (i  )-self%meanValueTarget     (i  )
       velocityMeanCovariance(j,j)=self%meanCovariance(i,i)+self%meanCovarianceTarget(i,i)
    end do
    ! Construct residual vector and covariance matrix.
    residual  =vector(velocityMeanDifference)
    covariance=matrix(velocityMeanCovariance)
    ! Compute the log-likelihood.
    logLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
    if (status /= GSL_Success) logLikelihood=logImprobable
    return
  end function subhaloVMaxVsMassLogLikelihood
