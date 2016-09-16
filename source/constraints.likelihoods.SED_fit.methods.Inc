!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a SED fitting likelihood function.

function likelihoodSEDFitConstructor(definition)
  !% Constructor for the SED fitting likelihood class.
  use FoX_DOM
  use IO_XML
  use Instruments_Filters
  use ISO_Varying_String
  use Memory_Management
  use Stellar_Population_Spectra_Postprocess
  use Galacticus_Error
  implicit none
  type     (likelihoodSEDFit       )                         :: likelihoodSEDFitConstructor
  type     (node                   ), pointer, intent(in   ) :: definition
  type     (node                   ), pointer                :: thisPhotometry             , thisElement
  type     (nodeList               ), pointer                :: photometry
  character(len=128                )                         :: filterName                 , systemName
  integer                                                    :: i

  ! Find the photometry for the galaxy to fit.
  photometry => getElementsByTagName(definition,"photometry")
  ! Allocate storage for photometry.
  likelihoodSEDFitConstructor%photometryCount=getLength(photometry)
  call allocateArray(likelihoodSEDFitConstructor%magnitude               ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%error                   ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%filterIndex             ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%luminosityIndex         ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%redshift                ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%age                     ,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%postprocessingChainIndex,[likelihoodSEDFitConstructor%photometryCount])
  call allocateArray(likelihoodSEDFitConstructor%wavelengthEffective     ,[likelihoodSEDFitConstructor%photometryCount])
  ! Extract photometry information.
  do i=1,likelihoodSEDFitConstructor%photometryCount
     thisPhotometry => item(photometry,i-1)
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"magnitude")
     call extractDataContent(thisElement,likelihoodSEDFitConstructor%magnitude(i))
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"error"    )
     call extractDataContent(thisElement,likelihoodSEDFitConstructor%error    (i))
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"filter"   )
     call extractDataContent(thisElement,filterName                              )
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"system"   )
     call extractDataContent(thisElement,systemName                              )
     ! Set a luminosity index.
     likelihoodSEDFitConstructor%luminosityIndex(i)=i
     ! Find the index for this filter.
     likelihoodSEDFitConstructor%filterIndex    (i)=Filter_Get_Index(var_str(trim(filterName)))
     ! Convert to AB system.
     if (trim(systemName) == "vega")                                         &
          & likelihoodSEDFitConstructor                     %magnitude  (i)  &
          &  =                   likelihoodSEDFitConstructor%magnitude  (i)  &
          &  -Filter_Vega_Offset(likelihoodSEDFitConstructor%filterIndex(i))
     ! Get the effective wavelength.
     likelihoodSEDFitConstructor%wavelengthEffective(i)=Filter_Wavelength_Effective(likelihoodSEDFitConstructor%filterIndex(i))
  end do
  ! Extract number of bursts.
  thisElement => XML_Get_First_Element_By_Tag_Name(definition,"burstCount")
  call extractDataContent(thisElement,likelihoodSEDFitConstructor%burstCount)
  call allocateArray(likelihoodSEDFitConstructor%burstTimeStart,[likelihoodSEDFitConstructor%burstCount])
  call allocateArray(likelihoodSEDFitConstructor%burstTimescale,[likelihoodSEDFitConstructor%burstCount])
  call allocateArray(likelihoodSEDFitConstructor%burstFraction ,[likelihoodSEDFitConstructor%burstCount])
  ! Extract dust model information.
  if (XML_Path_Exists(definition,"dustType")) then
     select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"dustType"))))
     case ("null"           )
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeNull
     case ("CharlotFall2000")
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeCharlotFall2000
     case ("Cardelli1989"   )
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeCardelli1989
     case ("Gordon2003"     )
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeGordon2003
     case ("Calzetti2000"   )
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeCalzetti2000
     case ("WittGordon2000" )
        likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeWittGordon2000
     case default
        call Galacticus_Error_Report('likelihoodSEDFitConstructor','unrecognized dust type')
     end select
  else
     likelihoodSEDFitConstructor%dustType=likelihoodSEDFitDustTypeNull
  end if
  ! Extract age definition.
  select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"startTimeType"))))
  case ("time")
     likelihoodSEDFitConstructor%startTimeType=likelihoodSEDFitStartTimeTypeTime
  case ("age" )
     likelihoodSEDFitConstructor%startTimeType=likelihoodSEDFitStartTimeTypeAge
  end select
  ! Find stellar spectra postprocessing chain to use.
  likelihoodSEDFitConstructor%postprocessingChainIndex=Stellar_Population_Spectrum_Postprocess_Index(var_str('default'))
  return
end function likelihoodSEDFitConstructor

double precision function likelihoodSEDFitEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
  !% Return the log-likelihood for the SED fitting likelihood function.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Error
  use FGSL
  use Constraints_Constants
  use Constraints_State
  use Constraints_Convergence
  use Numerical_Integration
  use Galactic_Structure_Options
  use Abundances_Structure
  use Cosmology_Functions
  use Stellar_Spectra_Dust_Attenuations
  use Star_Formation_IMF
  use Stellar_Population_Luminosities
  implicit none
  class           (likelihoodSEDFit                  ), intent(inout)                 :: self
  class           (state                             ), intent(in   )                 :: simulationState
  type            (mappingList                       ), intent(in   ), dimension(:  ) :: parameterMappings
  class           (convergence                       ), intent(inout)                 :: simulationConvergence
  double precision                                    , intent(in   )                 :: temperature                         , logLikelihoodCurrent   , &
       &                                                                                 logPriorCurrent                     , logPriorProposed
  real                                                , intent(inout)                 :: timeEvaluate
  double precision                                    , intent(  out), optional       :: logLikelihoodVariance
  double precision                                    , allocatable  , dimension(:  ) :: stateVector                         , ages                   , &
       &                                                                                 weights
  double precision                                    , allocatable  , dimension(:,:) :: massToLightRatios
  double precision                                    , parameter                     :: logImprobable              =-7.0d+0
  double precision                                    , parameter                     :: toleranceRelativeFractional=+1.0d-2
  double precision                                    , parameter                     :: stellarAgeArbitrary        =+1.0d+0
  double precision                                    , parameter                     :: vBandWavelength            =5504.61227375652d0
  class           (cosmologyFunctionsClass           ), pointer                       :: cosmologyFunctions_
  class           (stellarSpectraDustAttenuationClass), allocatable                   :: dust
  integer         (c_size_t                          )                                :: i
  integer                                                                             :: iMagnitude                          , imfIndex                     , &
       &                                                                                 burstIndexOffset
  double precision                                                                    :: mass                                , timeScale                    , &
       &                                                                                 starFormationRateNormalization      , magnitude                    , &
       &                                                                                 luminosity                          , metallicity                  , &
       &                                                                                 redshift                            , timeObserved                 , &
       &                                                                                 vBandAttenuation                    , opticalDepthBirthClouds      , &
       &                                                                                 timeStart                           , Rv                           , &
       &                                                                                 toleranceRelative                   , termLinear                   , &
       &                                                                                 termConstant                        , recycledFractionInstantaneous, &
       &                                                                                 agePrevious                         , ageNow                       , &
       &                                                                                 ageNext
  type            (abundances                )                                        :: abundancesStars
  type            (fgsl_function             )                                        :: integrandFunction
  type            (fgsl_integration_workspace)                                        :: integrationWorkspace
  logical                                                                             :: integrationReset                    , useRapidEvaluation
  !GCC$ attributes unused :: simulationConvergence, timeEvaluate
  
  ! There is no variance in our likelihood estimate.
  if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodSEDFitEvaluate=0.0d0
     return
  end if
  ! Get the simulation state.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  mass       =      stateVector(1)
  timeScale  =1.0d0/stateVector(3)
  metallicity=      stateVector(4)
  redshift   =      stateVector(5)
  ! Determine time observed.
  cosmologyFunctions_ => cosmologyFunctions()
  timeObserved        =  cosmologyFunctions_ %cosmicTime                 (          &
       &                  cosmologyFunctions_%expansionFactorFromRedshift (         &
       &                                                                   redshift &
       &                                                                  )         &
       &                                                                 )
  ! Determine start time.
  select case (self%startTimeType)
  case (likelihoodSEDFitStartTimeTypeTime)
     timeStart=             stateVector(2)
  case (likelihoodSEDFitStartTimeTypeAge )
     timeStart=timeObserved-stateVector(2)
  end select
  ! Return impossibility if start time is after the observed time or before the Big Bang.
  if (timeStart >= timeObserved*(1.0d0-toleranceRelativeFractional) .or. timeStart <= 0.0d0) then
     likelihoodSEDFitEvaluate=logImpossible
     return
  end if
  ! Construct dust attenuation object.
  select case (self%dustType)
  case (likelihoodSEDFitDustTypeNull           )
     allocate(stellarSpectraDustAttenuationNull :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationNull)
        vBandAttenuation=0.0d0
        dust=stellarSpectraDustAttenuationNull()
     end select
     burstIndexOffset=5
  case (likelihoodSEDFitDustTypeCharlotFall2000)
     allocate(stellarSpectraDustAttenuationCharlotFall2000 :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationCharlotFall2000)
        vBandAttenuation       =     stateVector(6)
        opticalDepthBirthClouds=-log(stateVector(7))
        dust=stellarSpectraDustAttenuationCharlotFall2000(                                                 &
             &                                            opacityExponent        =0.7d+0                 , &
             &                                            birthCloudLifetime     =1.0d-2                 , &
             &                                            opticalDepthISM        =1.0d+0                 , &
             &                                            opticalDepthBirthClouds=opticalDepthBirthClouds  &
             &                                           )
     end select
     burstIndexOffset=7
  case (likelihoodSEDFitDustTypeCardelli1989)
     allocate(stellarSpectraDustAttenuationCardelli1989    :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationCardelli1989)
        vBandAttenuation=stateVector(6)
        Rv              =stateVector(7)
        dust=stellarSpectraDustAttenuationCardelli1989   (                                                 &
             &                                            Rv                     =Rv                       &
             &                                           )
     end select
     burstIndexOffset=7
  case (likelihoodSEDFitDustTypeGordon2003)
     allocate(stellarSpectraDustAttenuationGordon2003      :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationGordon2003)
        vBandAttenuation=stateVector(6)
        dust=stellarSpectraDustAttenuationGordon2003     (                                                 &
             &                                            sample                 ='LMC'                    &
             &                                           )
     end select
     burstIndexOffset=6
  case (likelihoodSEDFitDustTypeCalzetti2000)
     allocate(stellarSpectraDustAttenuationCalzetti2000    :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationCalzetti2000)
        vBandAttenuation=stateVector(6)
        dust=stellarSpectraDustAttenuationCalzetti2000   (                                                 &
             &                                           )
     end select
     burstIndexOffset=6
  case (likelihoodSEDFitDustTypeWittGordon2000)
     allocate(stellarSpectraDustAttenuationWittGordon2000  :: dust)
     select type (dust)
     type is (stellarSpectraDustAttenuationWittGordon2000)
        vBandAttenuation=stateVector(6)
        dust=stellarSpectraDustAttenuationWittGordon2000 (                                                &
             &                                            model                 ='MilkyWayShellTau3.0'    &
             &                                           )
     end select
     burstIndexOffset=6
  case default
     burstIndexOffset=-1
     call Galacticus_Error_Report('likelihoodSEDFitEvaluate','unknown dust type')
  end select
  ! Extract bursts.
  do i=1,self%burstCount
     self%burstTimeStart(i)=stateVector(burstIndexOffset+(i-1)*3+1)
     self%burstTimescale(i)=stateVector(burstIndexOffset+(i-1)*3+2)
     self%burstFraction (i)=stateVector(burstIndexOffset+(i-1)*3+3)
     ! Return impossibility if the burst begins before the galaxy forms.
     if (self%burstTimeStart(i) < timeStart) then
        likelihoodSEDFitEvaluate=logImpossible
        return
     end if
  end do
  ! Set redshift.
  self%redshift=redshift
  ! Construct the abundances object.
  call abundancesStars%metallicitySet(metallicity,metallicityType=metallicityTypeLinearByMassSolar)
  ! Determine if we can use a rapid summation over tabulated mass-to-light ratios.
  useRapidEvaluation= self%burstCount == 0                        &
       &             .and.                                        &
       &              (                                           &
       &                     dust%isSeparable   ()                &
       &               .or.                                       &
       &                .not.dust%isAgeDependent()                &
       &              )                                           &
       &             .and.                                        &
       &              .not.IMF_Is_Star_Formation_Rate_Dependent()
  if (useRapidEvaluation) then
     ! Find the IMF.
     imfIndex         =IMF_Select(1.0d0,abundancesStars,componentTypeDisk)
     ! Get tables of ages and luminosities.
     call Stellar_Population_Luminosity_Track(                               &
          &                                   self%luminosityIndex         , &
          &                                   self%filterIndex             , &
          &                                   self%postprocessingChainIndex, &
          &                                   imfIndex                     , &
          &                                   abundancesStars              , &
          &                                   self%redshift                , &
          &                                   ages                         , &
          &                                   massToLightRatios              &
          &                                  )
     ! Evaluate the star formation rate normalization.
     starFormationRateNormalization=mass/timeScale/(1.0d0-exp(-(timeObserved-timeStart)/timeScale))
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
        recycledFractionInstantaneous=IMF_Recycled_Fraction_Instantaneous(1.0d0,abundancesStars,componentTypeDisk)
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
  likelihoodSEDFitEvaluate      =0.0d0
  starFormationRateNormalization=1.0d0
  do iMagnitude=0,size(self%magnitude)
     ! Compute luminosity.
     if (useRapidEvaluation) then
        if (iMagnitude > 0.0d0)                                                             &
             & luminosity=+sum(weights*massToLightRatios(:,iMagnitude))                     &
             &            *10.0d0**(                                                        &
             &                      -0.4d0                                                  &
             &                      *dust%attenuation(                                      &
             &                                        self%wavelengthEffective(iMagnitude), &
             &                                        stellarAgeArbitrary                 , &
             &                                        vBandAttenuation                      &
             &                                       )                                      &
             &                     )
     else
        ! Integrate star formation history to get luminosity.
        integrationReset=.true.
        if (iMagnitude == 0) then
           toleranceRelative=    1.0d-3
        else
           toleranceRelative=max(1.0d-3,toleranceRelativeFractional*self%error(imagnitude)*log(10.0d0)/2.5d0)
        end if
        luminosity=Integrate(                                         &
             &               timeStart                              , &
             &               timeObserved                           , &
             &               luminosityIntegrand                    , &
             &               integrandFunction                      , &
             &               integrationWorkspace                   , &
             &               toleranceRelative   =toleranceRelative , &
             &               integrationRule     =FGSL_Integ_Gauss61, &
             &               reset               =integrationReset    &
             &              )
        call Integrate_Done(integrandFunction,integrationWorkspace)
     end if
     ! Determine if we're computing mass or luminosity.
     if (iMagnitude > 0) then
        ! Accumulate to likelihood.
        magnitude=-2.5d0*log10(luminosity)
        likelihoodSEDFitEvaluate=             &
             & +likelihoodSEDFitEvaluate      &
             & -0.5d0                         &
             & *(                             &
             &   (                            &
             &    +     magnitude             &
             &    -self%magnitude(iMagnitude) &
             &   )                            &
             &   / self%error    (iMagnitude) &
             &  )**2
        ! Exit as soon as the proposed likelihood is found to be sufficiently below the current likelihood.
        if     (                                              &
             &    (likelihoodSEDFitEvaluate+logPriorProposed) &
             &   -                                            &
             &    (logLikelihoodCurrent    +logPriorCurrent ) &
             &  <                                             &
             &    logImprobable                               &
             &   *temperature                                 &
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
    !% Star formation rate integrand.
    use Star_Formation_IMF
    use Stellar_Population_Luminosities
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: time
    integer                         :: imfIndex
    double precision                :: recycledFractionInstantaneous, starFormationRate
    
    ! Find stellar population age.
    self%age         =timeObserved-time
    ! Find continuous star formation rate.
    starFormationRate=starFormationRateNormalization*exp(-(time-timeStart)/timeScale)
    ! Add burst star formation rate.
    do i=1,self%burstCount
       starFormationRate=starFormationRate+starFormationRateNormalization*timeScale*(1.0d0-exp(-self%burstTimeStart(i)/timeScale))*self%burstFraction(i)*exp(-0.5d0*(time-self%burstTimeStart(i))**2/self%burstTimescale(i)**2)/sqrt(2.0d0*Pi)
    end do
    ! Find the IMF.
    imfIndex         =IMF_Select(starFormationRate,abundancesStars,componentTypeDisk)
    ! Determine if we're computing luminosity or mass.
    if (iMagnitude > 0) then
       ! Find the mass-to-light ratio.
       self%massToLightRatio        =Stellar_Population_Luminosity      (                                                      &
            &                                                            self%luminosityIndex         (iMagnitude:iMagnitude), &
            &                                                            self%filterIndex             (iMagnitude:iMagnitude), &
            &                                                            self%postprocessingChainIndex(iMagnitude:iMagnitude), &
            &                                                            imfIndex                                            , &
            &                                                            abundancesStars                                     , &
            &                                                            self%age                     (iMagnitude:iMagnitude), &
            &                                                            self%redshift                (iMagnitude:iMagnitude)  &
            &                                                           )
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
       recycledFractionInstantaneous=IMF_Recycled_Fraction_Instantaneous(starFormationRate,abundancesStars,componentTypeDisk)
       ! Find the contribution to final mass.
       luminosityIntegrand          =starFormationRate*(1.0d0-recycledFractionInstantaneous)
    end if
    return
  end function luminosityIntegrand
   
end function likelihoodSEDFitEvaluate

subroutine likelihoodSEDFitFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodSEDFit), intent(inout) :: self
  !GCC$ attributes unused :: self
  
  return
end subroutine likelihoodSEDFitFunctionChanged
