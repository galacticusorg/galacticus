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
  An implementation of cosmological density field mass variance which modifies another member of the class by offsetting for
  the peak-background split.
  !!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <cosmologicalMassVariance name="cosmologicalMassVariancePeakBackgroundSplit">
   <description>
    The cosmological mass variance is computed by taking the variance from some other mass variance class, $\sigma^2(M)$, and
    offsetting it by the variance of the background in the peak-background split model, $\sigma^2(M_\mathrm{e})$, where
    $M_\mathrm{e}$ is the mass contained within the region defined as the background.
   </description>
  </cosmologicalMassVariance>
  !!]

  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVariancePeakBackgroundSplit
     !!{
     A cosmological mass variance class computing variance from a filtered power spectrum.
     !!}
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (haloEnvironmentClass         ), pointer :: haloEnvironment_          => null()
     double precision                                         :: massBackground                     , factorMassEnvironment     , &
          &                                                      timePrevious                       , varianceBackgroundPrevious
   contains
     !![
     <methods>
       <method description="Compute the variance of the background at this time." method="varianceBackground" />
     </methods>
     !!]
     final     ::                                       variancePeakBackgroundSplitDestructor
     procedure :: sigma8                             => variancePeakBackgroundSplitSigma8
     procedure :: powerNormalization                 => variancePeakBackgroundSplitPowerNormalization
     procedure :: rootVariance                       => variancePeakBackgroundSplitRootVariance
     procedure :: rootVarianceLogarithmicGradient    => variancePeakBackgroundSplitRootVarianceLogarithmicGradient
     procedure :: rootVarianceAndLogarithmicGradient => variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient
     procedure :: mass                               => variancePeakBackgroundSplitMass
     procedure :: growthIsMassDependent              => variancePeakBackgroundSplitGrowthIsMassDependent
     procedure :: varianceBackground                 => variancePeakBackgroundSplitVarianceBackground
  end type cosmologicalMassVariancePeakBackgroundSplit

  interface cosmologicalMassVariancePeakBackgroundSplit
     !!{
     Constructors for the \refClass{cosmologicalMassVariancePeakBackgroundSplit} cosmological mass variance class.
     !!}
     module procedure variancePeakBackgroundSplitConstructorParameters
     module procedure variancePeakBackgroundSplitConstructorInternal
  end interface cosmologicalMassVariancePeakBackgroundSplit

contains

  function variancePeakBackgroundSplitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{cosmologicalMassVariancePeakBackgroundSplit} cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologicalMassVariancePeakBackgroundSplit)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (haloEnvironmentClass                       ), pointer       :: haloEnvironment_
    class           (cosmologicalMassVarianceClass              ), pointer       :: cosmologicalMassVariance_
    class           (cosmologyParametersClass                   ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                    ), pointer       :: cosmologyFunctions_
    double precision                                                             :: factorMassEnvironment

    !![
    <inputParameter>
      <name>factorMassEnvironment</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The background variance is computed as $\sigma(\alpha M_\mathrm{env})^2$ where $M_\mathrm{env}$ is the environment mass and $\alpha=${\normalfont \ttfamily [factorMassEnvironment]}.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !!]
    self=cosmologicalMassVariancePeakBackgroundSplit(factorMassEnvironment,haloEnvironment_,cosmologicalMassVariance_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloEnvironment_"         />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    !!]
    return
  end function variancePeakBackgroundSplitConstructorParameters

  function variancePeakBackgroundSplitConstructorInternal(factorMassEnvironment,haloEnvironment_,cosmologicalMassVariance_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{cosmologicalMassVariancePeakBackgroundSplit} linear growth class.
    !!}
    implicit none
    type            (cosmologicalMassVariancePeakBackgroundSplit)                        :: self
    class           (haloEnvironmentClass                       ), intent(in   ), target :: haloEnvironment_
    class           (cosmologicalMassVarianceClass              ), intent(in   ), target :: cosmologicalMassVariance_
    class           (cosmologyParametersClass                   ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                    ), intent(in   ), target :: cosmologyFunctions_
    double precision                                             , intent(in   )         :: factorMassEnvironment
    !![
    <constructorAssign variables="factorMassEnvironment, *haloEnvironment_, *cosmologicalMassVariance_, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Evaluate and store the mass of the background.
    self%massBackground=+                      factorMassEnvironment   &
         &              *self%haloEnvironment_%environmentMass      ()
    ! Initialize memoized results.
    self%timePrevious              =-huge(0.0d0)
    self%varianceBackgroundPrevious=-huge(0.0d0)
    return
  end function variancePeakBackgroundSplitConstructorInternal

  subroutine variancePeakBackgroundSplitDestructor(self)
    !!{
    Destructor for the \refClass{cosmologicalMassVariancePeakBackgroundSplit} linear growth class.
    !!}
    implicit none
    type(cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self

    !![
    <objectDestructor name="self%haloEnvironment_"         />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    !!]
    return
  end subroutine variancePeakBackgroundSplitDestructor

  double precision function variancePeakBackgroundSplitPowerNormalization(self)
    !!{
    Return the normalization of the power spectrum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    !$GLC attributes unused :: self

    variancePeakBackgroundSplitPowerNormalization=0.0d0
    call Error_Report('power spectrum normalization is not well-defined in peak-background split model'//{introspection:location})
    return
  end function variancePeakBackgroundSplitPowerNormalization

  double precision function variancePeakBackgroundSplitSigma8(self)
    !!{
    Return the value of $\sigma_8$.
    !!}
    use :: Cosmology_Parameters    , only : hubbleUnitsLittleH
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    ! Radius for σ(M) normalization in Mpc/h.
    double precision                                             , parameter     :: radiusNormalization=8.0d0

    variancePeakBackgroundSplitSigma8=self%rootVariance(                                                                  &
         &                                              +(                                                                &
         &                                                +4.0d0                                                          &
         &                                                /3.0d0                                                          &
         &                                                *Pi                                                             &
         &                                               )                                                                &
         &                                              *  self%cosmologyParameters_%OmegaMatter    (                  )  &
         &                                              *  self%cosmologyParameters_%densityCritical(                  )  &
         &                                              *(                                                                &
         &                                                +radiusNormalization                                            &
         &                                                /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)  &
         &                                               )**3                                                           , &
         &                                                 self%cosmologyFunctions_ %cosmicTime     (1.0d0             )  &
         &                                             )
    return
  end function variancePeakBackgroundSplitSigma8

  double precision function variancePeakBackgroundSplitRootVariance(self,mass,time)
    !!{
    Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass         , time
    double precision                                                             :: varianceTotal, varianceBackground

    varianceTotal     =self%cosmologicalMassVariance_%rootVariance      (mass,time)**2
    varianceBackground=self                          %varianceBackground(     time)
    if (varianceTotal > varianceBackground) then
       variancePeakBackgroundSplitRootVariance=+sqrt(                    &
            &                                        +varianceTotal      &
            &                                        -varianceBackground &
            &                                       )
    else
       variancePeakBackgroundSplitRootVariance=+0.0d0
    end if
    return
  end function variancePeakBackgroundSplitRootVariance

  double precision function variancePeakBackgroundSplitRootVarianceLogarithmicGradient(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass         , time
    double precision                                                             :: varianceTotal, varianceBackground

    varianceTotal     =self%cosmologicalMassVariance_%rootVariance      (mass,time)**2
    varianceBackground=self                          %varianceBackground(     time)
    if (varianceTotal > varianceBackground) then
       variancePeakBackgroundSplitRootVarianceLogarithmicGradient=+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time) &
            &                                                     *                               varianceTotal                              &
            &                                                     /(                                                                         &
            &                                                       +                             varianceTotal                              &
            &                                                       -                             varianceBackground                         &
            &                                                      )
    else
       variancePeakBackgroundSplitRootVarianceLogarithmicGradient=+0.0d0
    end if
    return
  end function variancePeakBackgroundSplitRootVarianceLogarithmicGradient

  subroutine variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient)
    !!{
    Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    spherical region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass             , time
    double precision                                             , intent(  out) :: rootVariance     , rootVarianceLogarithmicGradient
    double precision                                                             :: varianceTotal    , varianceBackground                  , &
         &                                                                          variance         , rootVarianceTotalLogarithmicGradient, &
         &                                                                          rootVarianceTotal

    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,time,rootVarianceTotal,rootVarianceTotalLogarithmicGradient)
    varianceTotal     =     rootVarianceTotal       **2
    varianceBackground=self%varianceBackground(time)
    if (varianceTotal > varianceBackground) then
       variance                       =+(                    &
            &                            +varianceTotal      &
            &                            -varianceBackground &
            &                           )
       rootVariance                   =+sqrt(variance)
       rootVarianceLogarithmicGradient=+rootVarianceTotalLogarithmicGradient &
            &                          *varianceTotal                        &
            &                          /variance
    else
       rootVariance                   =+0.0d0
       rootVarianceLogarithmicGradient=+0.0d0
    end if
    return
  end subroutine variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient

  double precision function variancePeakBackgroundSplitMass(self,rootVariance,time)
    !!{
    Return the mass corresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    !!}
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: rootVariance      , time
    double precision                                                             :: varianceBackground

    varianceBackground             =self                          %varianceBackground(                              &
         &                                                                                   time                   &
         &                                                                           )
    variancePeakBackgroundSplitMass=self%cosmologicalMassVariance_%mass              (                              &
         &                                                                            +sqrt(                        &
         &                                                                                  +rootVariance      **2  &
         &                                                                                  +varianceBackground     &
         &                                                                                 )                      , &
         &                                                                                   time                   &
         &                                                                           )
    return
  end function variancePeakBackgroundSplitMass

  double precision function variancePeakBackgroundSplitVarianceBackground(self,time)
    !!{
    Return the variance of the background at this time.
    !!}
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: time

    if (time /= self%timePrevious) then
       self%timePrevious=time
       if (self%massBackground < huge(0.0d0)) then
          self%varianceBackgroundPrevious=self%cosmologicalMassVariance_%rootVariance(self%massBackground,time)**2
       else
          self%varianceBackgroundPrevious=0.0d0
       end if
    end if
    variancePeakBackgroundSplitVarianceBackground=self%varianceBackgroundPrevious
    return
  end function variancePeakBackgroundSplitVarianceBackground

  logical function variancePeakBackgroundSplitGrowthIsMassDependent(self)
    !!{
    Return true if the growth rate of the variance is mass-dependent.
    !!}
    implicit none
    class(cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self

    variancePeakBackgroundSplitGrowthIsMassDependent=self%cosmologicalMassVariance_%growthIsMassDependent()
    return
  end function variancePeakBackgroundSplitGrowthIsMassDependent
