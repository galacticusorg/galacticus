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
  Implementation of a posterior sampling likelihood class which implements a likelihood for SED fitting.
  !!}

  use :: Cosmology_Functions                       , only : cosmologyFunctionsClass
  use :: Stellar_Population_Selectors              , only : stellarPopulationSelectorClass
  use :: Stellar_Population_Spectra_Postprocess    , only : stellarPopulationSpectraPostprocessorBuilderClass, stellarPopulationSpectraPostprocessorList
  use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass

  !![
  <enumeration>
   <name>sedFitDustType</name>
   <description>Used to specify the type of dust model to use in SED fitting likelihoods.</description>
   <visibility>private</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="null"           />
   <entry label="charlotFall2000"/>
   <entry label="cardelli1989"   />
   <entry label="gordon2003"     />
   <entry label="calzetti2000"   />
   <entry label="wittGordon2000" />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>sedFitStartTime</name>
   <description>Used to specify the type of start time to use in SED fitting likelihoods.</description>
   <visibility>private</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="time"/>
   <entry label="age" />
  </enumeration>
  !!]

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodSEDFit">
   <description>A posterior sampling likelihood class which implements a likelihood for SED fitting.</description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodSEDFit
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for SED fitting.
     !!}
     private
     class           (cosmologyFunctionsClass                          ), pointer                   :: cosmologyFunctions_                           => null()
     class           (stellarPopulationSelectorClass                   ), pointer                   :: stellarPopulationSelector_                    => null()
     class           (stellarPopulationSpectraPostprocessorBuilderClass), pointer                   :: stellarPopulationSpectraPostprocessorBuilder_ => null()
     class           (stellarPopulationBroadBandLuminositiesClass      ), pointer                   :: stellarPopulationBroadBandLuminosities_       => null()
     type            (enumerationSEDFitDustTypeType                    )                            :: dustType
     type            (enumerationSEDFitStartTimeType                   )                            :: startTimeType
     integer                                                                                        :: photometryCount                                        , burstCount
     integer                                                            , allocatable, dimension(:) :: filterIndex                                            , luminosityIndex
     type            (stellarPopulationSpectraPostprocessorList        ), allocatable, dimension(:) :: postprocessor
     double precision                                                   , allocatable, dimension(:) :: magnitude                                              , error              , &
          &                                                                                            redshift                                               , burstFraction      , &
          &                                                                                            age                                                    , wavelengthEffective, &
          &                                                                                            burstTimeStart                                         , burstTimescale
     type            (varying_string                                   ), allocatable, dimension(:) :: filter                                                 , system
     type            (varying_string                                   )                            :: startTime
     double precision                                                                , dimension(1) :: massToLightRatio
   contains
     final     ::                    sedFitDestructor
     procedure :: evaluate        => sedFitEvaluate
     procedure :: functionChanged => sedFitFunctionChanged
  end type posteriorSampleLikelihoodSEDFit

  interface posteriorSampleLikelihoodSEDFit
     !!{
     Constructors for the \refClass{posteriorSampleLikelihoodSEDFit} posterior sampling convergence class.
     !!}
     module procedure sedFitConstructorParameters
     module procedure sedFitConstructorInternal
  end interface posteriorSampleLikelihoodSEDFit

contains

  function sedFitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodSEDFit} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodSEDFit                  )                              :: self
    type            (inputParameters                                  ), intent(inout)               :: parameters
    double precision                                                   , allocatable  , dimension(:) :: magnitude                                    , error
    type            (varying_string                                   ), allocatable  , dimension(:) :: filter                                       , system
    integer                                                                                          :: burstCount                                   , dustType, &
         &                                                                                              startTime
    class           (cosmologyFunctionsClass                          ), pointer                     :: cosmologyFunctions_
    class           (stellarPopulationSelectorClass                   ), pointer                     :: stellarPopulationSelector_
    class           (stellarPopulationSpectraPostprocessorBuilderClass), pointer                     :: stellarPopulationSpectraPostprocessorBuilder_
    class           (stellarPopulationBroadBandLuminositiesClass      ), pointer                     :: stellarPopulationBroadBandLuminosities_

    !![
    <inputParameter>
      <name>magnitude</name>
      <description>The magnitudes of the broad-band SED.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>error</name>
      <description>The errors on the magnitudes of the broad-band SED.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>filter</name>
      <description>The names of the filters in the broad-band SED.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>system</name>
      <description>The photometric system (AB or Vega) of the broad-band SED.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>burstCount</name>
      <description>The number of bursts events to include in the star formation history.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dustType</name>
      <description>The type of dust model to apply to the SED.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>startTime</name>
      <description>The definition of start time (absolute {\normalfont \ttfamily time} or {\normalfont \ttfamily age}).</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                           name="cosmologyFunctions_"                           source="parameters"/>
    <objectBuilder class="stellarPopulationSelector"                    name="stellarPopulationSelector_"                    source="parameters"/>
    <objectBuilder class="stellarPopulationSpectraPostprocessorBuilder" name="stellarPopulationSpectraPostprocessorBuilder_" source="parameters"/>
    <objectBuilder class="stellarPopulationBroadBandLuminosities"       name="stellarPopulationBroadBandLuminosities_"       source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodSEDFit(                                                                                                              &
         &                                                                     magnitude                                                             , &
         &                                                                     error                                                                 , &
         &                                                                     filter                                                                , &
         &                                                                     system                                                                , &
         &                                                                     burstCount                                                            , &
         &                               enumerationSedFitDustTypeEncode (char(dustType                                     ),includesPrefix=.false.), &
         &                               enumerationSedFitStartTimeEncode(char(startTime                                    ),includesPrefix=.false.), &
         &                                                                     cosmologyFunctions_                                                   , &
         &                                                                     stellarPopulationSelector_                                            , &
         &                                                                     stellarPopulationSpectraPostprocessorBuilder_                         , &
         &                                                                     stellarPopulationBroadBandLuminosities_                                 &
         &                              )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                          />
    <objectDestructor name="stellarPopulationSelector_"                   />
    <objectDestructor name="stellarPopulationSpectraPostprocessorBuilder_"/>
    <objectDestructor name="stellarPopulationBroadBandLuminosities_"      />
    !!]
    return
  end function sedFitConstructorParameters

  function sedFitConstructorInternal(magnitude,error,filter,system,burstCount,dustType,startTimeType,cosmologyFunctions_,stellarPopulationSelector_,stellarPopulationSpectraPostprocessorBuilder_,stellarPopulationBroadBandLuminosities_) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodSEDFit} posterior sampling likelihood class.
    !!}
    use :: ISO_Varying_String , only : var_str         , varying_string
    use :: Instruments_Filters, only : Filter_Get_Index, Filter_Vega_Offset, Filter_Wavelength_Effective
    implicit none
    type            (posteriorSampleLikelihoodSEDFit                  )                              :: self
    double precision                                                   , intent(in   ), dimension(:) :: magnitude                                    , error
    type            (varying_string                                   ), intent(in   ), dimension(:) :: filter                                       , system
    type            (enumerationSEDFitDustTypeType                    ), intent(in   )               :: dustType
    type            (enumerationSEDFitStartTimeType                   ), intent(in   )               :: startTimeType
    integer                                                            , intent(in   )               :: burstCount
    class           (cosmologyFunctionsClass                          ), intent(in   ), target       :: cosmologyFunctions_
    class           (stellarPopulationSelectorClass                   ), intent(in   ), target       :: stellarPopulationSelector_
    class           (stellarPopulationSpectraPostprocessorBuilderClass), intent(in   ), target       :: stellarPopulationSpectraPostprocessorBuilder_
    class           (stellarPopulationBroadBandLuminositiesClass      ), intent(in   ), target       :: stellarPopulationBroadBandLuminosities_
    integer                                                                                          :: i
    !![
    <constructorAssign variables="magnitude, error, filter, system, burstCount, dustType, startTimeType, *cosmologyFunctions_, *stellarPopulationSelector_, *stellarPopulationSpectraPostprocessorBuilder_, *stellarPopulationBroadBandLuminosities_"/>
    !!]

    self%photometryCount=size(magnitude)
    allocate          (self%postprocessor            (self%photometryCount))
    allocate(self%filterIndex             (self%photometryCount))
    allocate(self%luminosityIndex         (self%photometryCount))
    allocate(self%redshift                (self%photometryCount))
    allocate(self%age                     (self%photometryCount))
    allocate(self%wavelengthEffective     (self%photometryCount))
    ! Determine indices.
    do i=1,self%photometryCount
       ! Set a luminosity index.
       self%luminosityIndex(i)=i
       ! Find the index for this filter.
       self%filterIndex    (i)=Filter_Get_Index(filter(i))
       ! Convert to Vega system if required.
       if (system(i) == "vega")                         &
            & self                     %magnitude  (i)  &
            &  =                   self%magnitude  (i)  &
            &  -Filter_Vega_Offset(self%filterIndex(i))
       ! Get the effective wavelength.
       self%wavelengthEffective(i)=Filter_Wavelength_Effective(self%filterIndex(i))
    end do
    ! Create burst arrays.
    allocate(self%burstTimeStart(self%burstCount))
    allocate(self%burstTimescale(self%burstCount))
    allocate(self%burstFraction (self%burstCount))
    ! Find stellar spectra postprocessing chain to use.
    do i=1,self%photometryCount
       self%postprocessor(i)%stellarPopulationSpectraPostprocessor_ => self%stellarPopulationSpectraPostprocessorBuilder_%build(var_str('default'))
    end do
    return
  end function sedFitConstructorInternal

  subroutine sedFitDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleLikelihoodSEDFit} posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodSEDFit), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                          />
    <objectDestructor name="self%stellarPopulationSelector_"                   />
    <objectDestructor name="self%stellarPopulationSpectraPostprocessorBuilder_"/>
    <objectDestructor name="self%stellarPopulationBroadBandLuminosities_"      />
    !!]
    return
  end subroutine sedFitDestructor

  double precision function sedFitEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the SED fitting likelihood function.
    !!}
    use            :: Abundances_Structure             , only : abundances                          , max                                      , metallicityTypeLinearByMassSolar
    use            :: Error                            , only : Error_Report
    use            :: Galacticus_Nodes                 , only : nodeComponentDisk
    use, intrinsic :: ISO_C_Binding                    , only : c_size_t
    use            :: Models_Likelihoods_Constants     , only : logImpossible
    use            :: Numerical_Integration            , only : integrator                          , GSL_Integ_Gauss61
    use            :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use            :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    use            :: Stellar_Populations              , only : stellarPopulationClass
    use            :: Stellar_Spectra_Dust_Attenuations, only : gordon2003SampleLMC                 , stellarSpectraDustAttenuationCalzetti2000, stellarSpectraDustAttenuationCardelli1989  , stellarSpectraDustAttenuationCharlotFall2000, &
          &                                                     stellarSpectraDustAttenuationClass  , stellarSpectraDustAttenuationGordon2003  , stellarSpectraDustAttenuationWittGordon2000, stellarSpectraDustAttenuationZero           , &
          &                                                     wittGordon2000ModelMilkyWayShellTau3
    implicit none
    class           (posteriorSampleLikelihoodSEDFit   ), intent(inout), target         :: self
    class           (posteriorSampleStateClass         ), intent(inout)                 :: simulationState
    type            (modelParameterList                ), intent(inout), dimension(:  ) :: modelParametersActive_                          , modelParametersInactive_
    class           (posteriorSampleConvergenceClass   ), intent(inout)                 :: simulationConvergence
    double precision                                    , intent(in   )                 :: temperature                                     , logLikelihoodCurrent         , &
         &                                                                                 logPriorCurrent                                 , logPriorProposed
    real                                                , intent(inout)                 :: timeEvaluate
    double precision                                    , intent(  out), optional       :: logLikelihoodVariance
    logical                                             , intent(inout), optional       :: forceAcceptance
    double precision                                    , allocatable  , dimension(:  ) :: stateVector                                     , ages                         , &
         &                                                                                 weights
    double precision                                    , allocatable  , dimension(:,:) :: massToLightRatios
    double precision                                    , parameter                     :: logUnlikely                =-7.00000000000000d+0
    double precision                                    , parameter                     :: toleranceRelativeFractional=+1.00000000000000d-2
    double precision                                    , parameter                     :: stellarAgeArbitrary        =+1.00000000000000d+0
    double precision                                    , parameter                     :: vBandWavelength            =+5.50461227375652d+3
    class           (stellarSpectraDustAttenuationClass), allocatable                   :: dust
    class           (stellarPopulationClass            ), pointer                       :: stellarPopulation_
    type            (nodeComponentDisk                 )                                :: disk
    integer         (c_size_t                          )                                :: i
    integer                                                                             :: iMagnitude                                      , burstIndexOffset
    double precision                                                                    :: mass                                            , timeScale                    , &
         &                                                                                 starFormationRateNormalization                  , magnitude                    , &
         &                                                                                 luminosity                                      , metallicity                  , &
         &                                                                                 redshift                                        , timeObserved                 , &
         &                                                                                 vBandAttenuation                                , opticalDepthBirthClouds      , &
         &                                                                                 timeStart                                       , Rv                           , &
         &                                                                                 toleranceRelative                               , termLinear                   , &
         &                                                                                 termConstant                                    , recycledFractionInstantaneous, &
         &                                                                                 agePrevious                                     , ageNow                       , &
         &                                                                                 ageNext
    type            (abundances                )                                        :: abundancesStars
    type            (integrator                )                                        :: integrator_
    logical                                                                             :: useRapidEvaluation
    !$GLC attributes unused :: simulationConvergence, timeEvaluate, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       sedFitEvaluate=0.0d0
       return
    end if
    ! Get the simulation state.
    stateVector=simulationState%get()
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    mass        =      stateVector(1)
    timeScale   =1.0d0/stateVector(3)
    metallicity =      stateVector(4)
    redshift    =      stateVector(5)
    ! Determine time observed.
    timeObserved=self%cosmologyFunctions_%cosmicTime                 (          &
         &       self%cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                             redshift &
         &                                                            )         &
         &                                                           )
    ! Determine start time.
    select case (self%startTimeType%ID)
    case (sedFitStartTimeTime%ID)
       timeStart=             stateVector(2)
    case (sedFitStartTimeAge %ID)
       timeStart=timeObserved-stateVector(2)
    end select
    ! Return impossibility if start time is after the observed time or before the Big Bang.
    if (timeStart >= timeObserved*(1.0d0-toleranceRelativeFractional) .or. timeStart <= 0.0d0) then
       sedFitEvaluate=logImpossible
       return
    end if
    ! Construct dust attenuation object.
    select case (self%dustType%ID)
    case (sedFitDustTypeNull           %ID)
       allocate(stellarSpectraDustAttenuationZero :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationZero)
          vBandAttenuation=0.0d0
          dust=stellarSpectraDustAttenuationZero()
       end select
       burstIndexOffset=5
    case (sedFitDustTypeCharlotFall2000%ID)
       allocate(stellarSpectraDustAttenuationCharlotFall2000 :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationCharlotFall2000)
          vBandAttenuation       =     stateVector(6)
          opticalDepthBirthClouds=-log(stateVector(7))
          dust=stellarSpectraDustAttenuationCharlotFall2000(                                                              &
               &                                            opacityExponent        =0.7d+0                              , &
               &                                            birthCloudLifetime     =1.0d-2                              , &
               &                                            opticalDepthISM        =1.0d+0                              , &
               &                                            opticalDepthBirthClouds=opticalDepthBirthClouds               &
               &                                           )
       end select
       burstIndexOffset=7
    case (sedFitDustTypeCardelli1989%ID)
       allocate(stellarSpectraDustAttenuationCardelli1989    :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationCardelli1989)
          vBandAttenuation=stateVector(6)
          Rv              =stateVector(7)
          dust=stellarSpectraDustAttenuationCardelli1989   (                                                              &
               &                                            Rv                     =Rv                                    &
               &                                           )
       end select
       burstIndexOffset=7
    case (sedFitDustTypeGordon2003%ID)
       allocate(stellarSpectraDustAttenuationGordon2003      :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationGordon2003)
          vBandAttenuation=stateVector(6)
          dust=stellarSpectraDustAttenuationGordon2003     (                                                              &
               &                                            sample                 =gordon2003SampleLMC                   &
               &                                           )
       end select
       burstIndexOffset=6
    case (sedFitDustTypeCalzetti2000%ID)
       allocate(stellarSpectraDustAttenuationCalzetti2000    :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationCalzetti2000)
          vBandAttenuation=stateVector(6)
          dust=stellarSpectraDustAttenuationCalzetti2000   (                                                              &
               &                                           )
       end select
       burstIndexOffset=6
    case (sedFitDustTypeWittGordon2000%ID)
       allocate(stellarSpectraDustAttenuationWittGordon2000  :: dust)
       select type (dust)
       type is (stellarSpectraDustAttenuationWittGordon2000)
          vBandAttenuation=stateVector(6)
          dust=stellarSpectraDustAttenuationWittGordon2000 (                                                             &
               &                                            model                 =wittGordon2000ModelMilkyWayShellTau3  &
               &                                           )
       end select
       burstIndexOffset=6
    case default
       burstIndexOffset=-1
       call Error_Report('unknown dust type'//{introspection:location})
    end select
    ! Extract bursts.
    do i=1,self%burstCount
       self%burstTimeStart(i)=stateVector(burstIndexOffset+(i-1)*3+1)
       self%burstTimescale(i)=stateVector(burstIndexOffset+(i-1)*3+2)
       self%burstFraction (i)=stateVector(burstIndexOffset+(i-1)*3+3)
       ! Return impossibility if the burst begins before the galaxy forms.
       if (self%burstTimeStart(i) < timeStart) then
          sedFitEvaluate=logImpossible
          return
       end if
    end do
    ! Set redshift.
    self%redshift=redshift
    ! Construct the abundances object.
    call abundancesStars%metallicitySet(metallicity,metallicityType=metallicityTypeLinearByMassSolar)
    ! Determine if we can use a rapid summation over tabulated mass-to-light ratios.
    stellarPopulation_ =>  self%stellarPopulationSelector_%select(1.0d0,abundancesStars,disk)
    useRapidEvaluation =   self%burstCount == 0                                                  &
         &                .and.                                                                  &
         &                 (                                                                     &
         &                        dust                           %isSeparable                 () &
         &                  .or.                                                                 &
         &                   .not.dust                           %isAgeDependent              () &
         &                 )                                                                     &
         &                .and.                                                                  &
         &                   .not.self%stellarPopulationSelector_%isStarFormationRateDependent()
    if (useRapidEvaluation) then
       ! Get tables of ages and luminosities.
       call self%stellarPopulationBroadBandLuminosities_%luminosityTracks(                      &
            &                                                             self%luminosityIndex, &
            &                                                             self%filterIndex    , &
            &                                                             self%postprocessor  , &
            &                                                             stellarPopulation_  , &
            &                                                             abundancesStars     , &
            &                                                             self%redshift       , &
            &                                                             ages                , &
            &                                                             massToLightRatios     &
            &                                                            )
       ! Evaluate the star formation rate normalization.
       starFormationRateNormalization=+mass                  &
            &                         /timeScale             &
            &                         /(                     &
            &                           +1.0d0               &
            &                           -exp(                &
            &                                -(              &
            &                                  +timeObserved &
            &                                  -timeStart    &
            &                                 )              &
            &                                /  timeScale    &
            &                               )                &
            &                          )
       ! Evaluate weight for each tabulated age.
       allocate(weights(size(ages)))
       weights(size(ages))=0.0d0
       do i=1,size(ages,kind=c_size_t)-1
          if (i == 1) then
             agePrevious=0.0d0
          else
             agePrevious=min(ages(i-1),timeObserved-timeStart)
          end if
          ageNow =min(ages(i  ),timeObserved-timeStart)
          ageNext=min(ages(i+1),timeObserved-timeStart)
          termConstant=0.0d0
          if (ageNext     > ageNow) termConstant=                                                                                  &
               &                                 +termConstant                                                                     &
               &                                 +(1.0d0+ageNow/(ageNext-ageNow))                                                  &
               &                                 *(                                                                                &
               &                                   -exp(-((timeObserved-ageNow     )-timeStart)/timeScale)                         &
               &                                   +exp(-((timeObserved-ageNext    )-timeStart)/timeScale)                         &
               &                                  )
          if (agePrevious < ageNow) termConstant=                                                                                  &
               &                                 +termConstant                                                                     &
               &                                 -(agePrevious/(ageNow-agePrevious))                                               &
               &                                 *(                                                                                &
               &                                   -exp(-((timeObserved-agePrevious)-timeStart)/timeScale)                         &
               &                                   +exp(-((timeObserved-ageNow     )-timeStart)/timeScale)                         &
               &                                 )
          termLinear=0.0d0
          if (ageNext     > ageNow) termLinear  =                                                                                  &
               &                                 +termLinear                                                                       &
               &                                 -(                                                                                &
               &                                   -exp(-((timeObserved-ageNow     )-timeStart)/timeScale)*(ageNow     -timeScale) &
               &                                   +exp(-((timeObserved-ageNext    )-timeStart)/timeScale)*(ageNext    -timeScale) &
               &                                 )                                                                                 &
               &                                 /(ageNext-ageNow)
          if (agePrevious < ageNow) termLinear  =                                                                                  &
               &                                 +termLinear                                                                       &
               &                                 +(                                                                                &
               &                                   -exp(-((timeObserved-agePrevious)-timeStart)/timeScale)*(agePrevious-timeScale) &
               &                                   +exp(-((timeObserved-ageNow     )-timeStart)/timeScale)*(ageNow     -timeScale) &
               &                                 )                                                                                 &
               &                                 /(ageNow-agePrevious)
          ! Find the recycled fraction.
          stellarPopulation_            => self              %stellarPopulationSelector_%select                       (1.0d0,abundancesStars,disk)
          recycledFractionInstantaneous =  stellarPopulation_                           %recycledFractionInstantaneous(                          )
          weights(i)=(termConstant+termLinear)*starFormationRateNormalization*timeScale/(1.0d0-recycledFractionInstantaneous)
          if (.not.dust%isSeparable())                                       &
               & weights(i)=                                                 &
               &            +weights(i)                                      &
               &            *10.0d0**(                                       &
               &                      -0.4d0                                 &
               &                      *dust%attenuation(                     &
               &                                        vBandWavelength    , &
               &                                        ageNow             , &
               &                                        vBandAttenuation     &
               &                                       )                     &
               &                      +0.4d0                                 &
               &                      *dust%attenuation(                     &
               &                                        vBandWavelength    , &
               &                                        stellarAgeArbitrary, &
               &                                        vBandAttenuation     &
               &                                       )                     &
               &                     )
       end do
    else
       allocate(weights(0))
    end if
    ! Iterate over bands.
    sedFitEvaluate                =0.0d0
    starFormationRateNormalization=1.0d0
    do iMagnitude=0,size(self%magnitude)
       ! Compute luminosity.
       if (useRapidEvaluation) then
          if (iMagnitude > 0) then
             luminosity=+sum(weights*massToLightRatios(:,iMagnitude))                     &
                  &     *10.0d0**(                                                        &
                  &               -0.4d0                                                  &
                  &               *dust%attenuation(                                      &
                  &                                 self%wavelengthEffective(iMagnitude), &
                  &                                 stellarAgeArbitrary                 , &
                  &                                 vBandAttenuation                      &
                  &                                )                                      &
                  &              )
          else
             luminosity=0.0d0
          end if
       else
          ! Integrate star formation history to get luminosity.
          integrator_=integrator(luminosityIntegrand,toleranceRelative=toleranceRelative,integrationRule=GSL_Integ_Gauss61)
          if (iMagnitude == 0) then
             toleranceRelative=    1.0d-3
          else
             toleranceRelative=max(1.0d-3,toleranceRelativeFractional*self%error(imagnitude)*log(10.0d0)/2.5d0)
          end if
          luminosity=integrator_%integrate(timeStart,timeObserved)
       end if
       ! Determine if we're computing mass or luminosity.
       if (iMagnitude > 0) then
          ! Accumulate to likelihood.
          magnitude=-2.5d0*log10(luminosity)
          sedFitEvaluate=                       &
               & +sedFitEvaluate                &
               & -0.5d0                         &
               & *(                             &
               &   (                            &
               &    +     magnitude             &
               &    -self%magnitude(iMagnitude) &
               &   )                            &
               &   / self%error    (iMagnitude) &
               &  )**2
          ! Exit as soon as the proposed likelihood is found to be sufficiently below the current likelihood.
          if     (                                          &
               &    (sedFitEvaluate      +logPriorProposed) &
               &   -                                        &
               &    (logLikelihoodCurrent+logPriorCurrent ) &
               &  <                                         &
               &    logUnlikely                             &
               &   *temperature                             &
               & ) exit
       else if (.not.useRapidEvaluation) then
          ! Normalize the star formation rate.
          starFormationRateNormalization=mass/luminosity
       end if
    end do
    ! Destroy abundances object.
    call abundancesStars%destroy()
    return

  contains

    double precision function luminosityIntegrand(time)
      !!{
      Star formation rate integrand.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: time
      double precision                :: recycledFractionInstantaneous, starFormationRate

      ! Find stellar population age.
      self%age         =+timeObserved&
           &            -time
      ! Find continuous star formation rate.
      starFormationRate=+starFormationRateNormalization &
           &            *exp(                           &
           &                 -(                         &
           &                   +time                    &
           &                   -timeStart               &
           &                  )                         &
           &                 /  timeScale               &
           &                )
      ! Add burst star formation rate.
      do i=1,self%burstCount
         starFormationRate=+starFormationRate                &
              &            +starFormationRateNormalization   &
              &            *timeScale                        &
              &            *(                                &
              &              +1.0d0                          &
              &              -exp(                           &
              &                   -self%burstTimeStart(i)    &
              &                   /timeScale                 &
              &                  )                           &
              &             )                                &
              &            *       self%burstFraction (i)    &
              &            *exp(                             &
              &                 -0.5d0                       &
              &                 *(                           &
              &                   +time                      &
              &                   -self%burstTimeStart(i)    &
              &                  )**2                        &
              &                 /  self%burstTimescale(i)**2 &
              &                )                             &
              &            /sqrt(                            &
              &                  +2.0d0                      &
              &                  *Pi                         &
              &                 )
      end do
      ! Find the stellar population.
      stellarPopulation_ => self%stellarPopulationSelector_%select(1.0d0,abundancesStars,disk)
      ! Determine if we're computing luminosity or mass.
      if (iMagnitude > 0) then
         ! Find the mass-to-light ratio.
         self%massToLightRatio=self%stellarPopulationBroadBandLuminosities_%luminosities(                                             &
              &                                                                          self%luminosityIndex(iMagnitude:iMagnitude), &
              &                                                                          self%filterIndex    (iMagnitude:iMagnitude), &
              &                                                                          self%postprocessor  (iMagnitude:iMagnitude), &
              &                                                                          stellarPopulation_                         , &
              &                                                                          abundancesStars                            , &
              &                                                                          self%age            (iMagnitude:iMagnitude), &
              &                                                                          self%redshift       (iMagnitude:iMagnitude)  &
              &                                                                         )
         luminosityIntegrand=                                                                  &
              &              +starFormationRate                                                &
              &              *self%massToLightRatio(1)                                         &
              &              *10.0d0**(                                                        &
              &                        -0.4d0                                                  &
              &                        *dust%attenuation(                                      &
              &                                          self%wavelengthEffective(iMagnitude), &
              &                                          self%age                (iMagnitude), &
              &                                          vBandAttenuation                      &
              &                                         )                                      &
              &                       )
      else
         ! Find the recycled fraction.
          stellarPopulation_            => self              %stellarPopulationSelector_%select                       (starFormationRate,abundancesStars,disk)
          recycledFractionInstantaneous =  stellarPopulation_                           %recycledFractionInstantaneous(                                      )
         ! Find the contribution to final mass.
         luminosityIntegrand          =starFormationRate*(1.0d0-recycledFractionInstantaneous)
      end if
      return
    end function luminosityIntegrand

  end function sedFitEvaluate

  subroutine sedFitFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodSEDFit), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine sedFitFunctionChanged
