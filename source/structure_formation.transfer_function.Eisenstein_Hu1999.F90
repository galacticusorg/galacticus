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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implements a transfer function class using the fitting function of \cite{eisenstein_power_1999}.
  !!}

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionEisensteinHu1999">
   <description>
    Provides the \cite{eisenstein_power_1999} fitting function to the transfer function. The effective number of neutrino
    species and the summed mass (in electron volts) of all neutrino species are specified via the {\normalfont \ttfamily
    neutrinoNumberEffective} and {\normalfont \ttfamily neutrinoMassSummed} parameters respectively.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionEisensteinHu1999
     !!{
     The {\normalfont \ttfamily eisensteinHu1999} transfer function class.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
     double precision                                    :: temperatureCMB27              , distanceSoundWave      , &
          &                                                 neutrinoMassFraction          , neutrinoNumberEffective, &
          &                                                 neutrinoFactor                , betaDarkMatter         , &
          &                                                 neutrinoMassSummed            , shapeParameterEffective
     double precision                                    :: C                             , L                      , &
          &                                                 wavenumberEffective           , wavenumberNeutrino     , &
          &                                                 wavenumberEffectivePow        , wavenumberPrevious     , &
          &                                                 time
   contains
     !![
     <methods>
       <method description="Compute common factors needed by \cite{eisenstein_power_1999} transfer function calculations." method="computeFactors" />
     </methods>
     !!]
     final     ::                          eisensteinHu1999Destructor
     procedure :: value                 => eisensteinHu1999Value
     procedure :: logarithmicDerivative => eisensteinHu1999LogarithmicDerivative
     procedure :: computeFactors        => eisensteinHu1999ComputeFactors
     procedure :: halfModeMass          => eisensteinHu1999HalfModeMass
     procedure :: quarterModeMass       => eisensteinHu1999QuarterModeMass
     procedure :: epochTime             => eisensteinHu1999EpochTime
  end type transferFunctionEisensteinHu1999

  interface transferFunctionEisensteinHu1999
     !!{
     Constructors for the {\normalfont \ttfamily eisensteinHu1999} transfer function class.
     !!}
     module procedure eisensteinHu1999ConstructorParameters
     module procedure eisensteinHu1999ConstructorInternal
  end interface transferFunctionEisensteinHu1999

contains

  function eisensteinHu1999ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily eisensteinHu1999} transfer function class
    which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (transferFunctionEisensteinHu1999)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (darkMatterParticleClass         ), pointer       :: darkMatterParticle_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    double precision                                                  :: neutrinoNumberEffective             , neutrinoMassSummed

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>neutrinoNumberEffective</name>
      <source>parameters</source>
      <defaultValue>3.046d0</defaultValue>
      <defaultSource>\citep{mangano_relic_2005}</defaultSource>
      <description>The effective number of neutrino species.</description>
    </inputParameter>
    <inputParameter>
      <name>neutrinoMassSummed</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The summed mass (in electron volts) of all neutrino species.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    ! Call the internal constructor.
    self=transferFunctionEisensteinHu1999(neutrinoNumberEffective,neutrinoMassSummed,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function eisensteinHu1999ConstructorParameters

  function eisensteinHu1999ConstructorInternal(neutrinoNumberEffective,neutrinoMassSummed,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily eisensteinHu1999} transfer function class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM
    use :: Error                , only : Error_Report
    implicit none
    type            (transferFunctionEisensteinHu1999)                        :: self
    double precision                                  , intent(in   )         :: neutrinoNumberEffective     , neutrinoMassSummed
    class           (darkMatterParticleClass         ), intent(in   ), target :: darkMatterParticle_
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                          :: redshiftEquality            , redshiftComptonDrag   , &
         &                                                                       b1                          , b2                    , &
         &                                                                       massFractionBaryonic        , massFractionDarkMatter, &
         &                                                                       expansionFactorRatio        , massFractionMatter    , &
         &                                                                       massFractionBaryonsNeutrinos, suppressionDarkMatter , &
         &                                                                       suppressionMatter
    !![
    <constructorAssign variables="*darkMatterParticle_, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
       class is (darkMatterParticleCDM)
          ! Cold dark matter particle - this is as expected.
       class default
       call Error_Report('transfer function expects a cold dark matter particle'//{introspection:location})
    end select
    ! Compute the epoch - the transfer function is assumed to be for z=0.
    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))
    ! Present day CMB temperature [in units of 2.7K].
    self%temperatureCMB27       =+self%cosmologyParameters_%temperatureCMB(                  )    &
         &                                                       /2.7d0
    ! Redshift of matter-radiation equality.
    redshiftEquality            =+2.50d4                                                          &
         &                       *self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                       *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                       /self%temperatureCMB27                                       **4
    ! Compute redshift at which baryons are released from Compton drag of photons (eqn. 2)
    b1                          =+0.313d0                                                             &
         &                       /(                                                                   &
         &                         +  self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                         *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                        )**0.419d0                                                          &
         &                       *(                                                                   &
         &                         +1.000d0                                                           &
         &                         +0.607d0                                                           &
         &                         *(                                                                 &
         &                           +self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                           *self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                          )**0.674d0                                                        &
         &                        )
    b2                          =+0.238d0                                                             &
         &                       *(                                                                   &
         &                         +  self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                         *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                        )**0.223d0
    redshiftComptonDrag         =+1291.0d0                                                            &
         &                       *(                                                                   &
         &                         +  self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                         *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                        )**0.251d0                                                          &
         &                       *(                                                                   &
         &                         +1.0d0                                                             &
         &                         +b1                                                                &
         &                         *(                                                                 &
         &                           +self%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                           *self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                          )**b2                                                             &
         &                        )                                                                   &
         &                       /(                                                                   &
         &                         +1.000d0                                                           &
         &                         +0.659d0                                                           &
         &                         *(                                                                 &
         &                           +self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                           *self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                          )**0.828d0                                                        &
         &                        )
    ! Relative expansion factor between previous two computed redshifts.
    expansionFactorRatio        =+(1.0d0+redshiftEquality   ) &
         &                       /(1.0d0+redshiftComptonDrag)
    ! Compute the comoving distance that a sound wave can propagate prior to redshiftComptonDrag (i.e. sound horizon; eq. 4)
    self%distanceSoundWave      =+44.5d0                                                                  &
         &                       *log (                                                                   &
         &                             +9.83d0                                                            &
         &                             /  self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                             /  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                            )                                                                   &
         &                       /sqrt(                                                                   &
         &                             + 1.0d0                                                            &
         &                             +10.0d0                                                            &
         &                             *(                                                                 &
         &                               +self%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                               *self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                              )**0.75d0                                                         &
         &                            )
    ! Specify properties of neutrinos. Mass fraction formula is from Komatsu et al. (2007; http://adsabs.harvard.edu/abs/2010arXiv1001.4538K).
    self%neutrinoMassSummed     =+neutrinoMassSummed
    self%neutrinoMassFraction   =+neutrinoMassSummed                                              &
         &                       /94.0d0                                                          &
         &                       /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                       /self%cosmologyParameters_%OmegaMatter   (                  )
    self%neutrinoNumberEffective=+neutrinoNumberEffective
    ! Compute baryonic and cold dark matter fractions.
    massFractionBaryonic        =+  self%cosmologyParameters_%OmegaBaryon() &
         &                       /  self%cosmologyParameters_%OmegaMatter()
    massFractionDarkMatter      =+(                                         &
         &                         +self%cosmologyParameters_%OmegaMatter() &
         &                         -self%cosmologyParameters_%OmegaBaryon() &
         &                        )                                         &
         &                       /  self%cosmologyParameters_%OmegaMatter()
    ! Total matter fraction.
    massFractionMatter          =+massFractionBaryonic   &
         &                       +massFractionDarkMatter
    ! Baryonic + neutrino fraction.
    massFractionBaryonsNeutrinos=+self%neutrinoMassFraction &
         &                       +massFractionBaryonic
    ! Compute small scale suppression factor (eqn. 15).
    suppressionDarkMatter  =+0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*massFractionDarkMatter))
    suppressionMatter      =+0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*massFractionMatter    ))
    self%neutrinoFactor    =+(                                                         &
         &                    +massFractionDarkMatter                                  &
         &                    /massFractionMatter                                      &
         &                   )                                                         &
         &                  *(                                                         &
         &                    +(5.0d0-2.0d0*(suppressionDarkMatter+suppressionMatter)) &
         &                    /(5.0d0-4.0d0*                       suppressionMatter ) &
         &                   )                                                         &
         &                  *(                                                         &
         &                    +(                                                       &
         &                      +1.000d0                                               &
         &                      -0.533d0                                               &
         &                      *massFractionBaryonsNeutrinos                          &
         &                      +0.126d0                                               &
         &                      *massFractionBaryonsNeutrinos**3                       &
         &                     )                                                       &
         &                    *(                                                       &
         &                      +1.0d0                                                 &
         &                      +expansionFactorRatio                                  &
         &                     )**(                                                    &
         &                         +suppressionMatter                                  &
         &                         -suppressionDarkMatter                              &
         &                        )                                                    &
         &                    /(                                                       &
         &                      +1.000d0                                               &
         &                      -0.193d0                                               &
         &                      *sqrt(                                                 &
         &                            +self%neutrinoMassFraction                       &
         &                            *self%neutrinoNumberEffective                    &
         &                           )                                                 &
         &                      +0.169d0                                               &
         &                      *self%neutrinoMassFraction                             &
         &                      *self%neutrinoNumberEffective**0.2d0                   &
         &                     )                                                       &
         &                   )                                                         &
         &                  *(                                                         &
         &                    +1.0d0                                                   &
         &                    +0.5d0                                                   &
         &                    *(                                                       &
         &                      +suppressionDarkMatter                                 &
         &                      -suppressionMatter                                     &
         &                     )                                                       &
         &                    *(                                                       &
         &                      +1.0d0                                                 &
         &                      +1.0d0                                                 &
         &                      /(3.0d0-4.0d0*suppressionDarkMatter)                   &
         &                      /(7.0d0-4.0d0*suppressionMatter    )                   &
         &                     )                                                       &
         &                    /(                                                       &
         &                      +1.0d0                                                 &
         &                      +expansionFactorRatio                                  &
         &                     )                                                       &
         &                   )
    self%betaDarkMatter    =+1.0d0                                                     & ! Eqn. 21.
         &                  /(                                                         &
         &                    +1.0d0                                                   &
         &                    -0.949d0                                                 &
         &                    *massFractionBaryonsNeutrinos                            &
         &                   )
    ! Initialize wavenumber for which results were computed to a non-physical value.
    self%wavenumberPrevious=-1.0d0
    return
  end function eisensteinHu1999ConstructorInternal

  subroutine eisensteinHu1999Destructor(self)
    !!{
    Destructor for the eisensteinHu1999 transfer function class.
    !!}
    implicit none
    type(transferFunctionEisensteinHu1999), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine eisensteinHu1999Destructor

  subroutine eisensteinHu1999ComputeFactors(self,wavenumber)
    !!{
    Compute common factors required by {\normalfont \ttfamily eisensteinHu1999} transfer function class.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: wavenumberScaleFree


    ! If called again with the same wavenumber, return without recomputing result.
    if (wavenumber == self%wavenumberPrevious) return
    self%wavenumberPrevious=wavenumber
    ! Compute rescaled shape parameter (eqn. 16)
    self%shapeParameterEffective=+self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                       *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                       *(                                                               &
         &                         +       sqrt(self%neutrinoFactor)                              &
         &                         +(1.0d0-sqrt(self%neutrinoFactor))                             &
         &                         /(                                                             &
         &                           +1.0d0                                                       &
         &                           +(                                                           &
         &                             +0.43d0                                                    &
         &                             *wavenumber                                                &
         &                             *self%distanceSoundWave                                    &
         &                            )**4                                                        &
         &                          )                                                             &
         &                       )
    self%wavenumberEffective    =+wavenumber                                                 &
         &                       *self%temperatureCMB27**2                                   &
         &                       /self%shapeParameterEffective
    self%wavenumberEffectivePow=+self%wavenumberEffective**1.11d0
    self%L                      =+log(                                                       & ! Eqn. 19.
         &                       +exp(1.0d0)                                                 &
         &                       +1.84d0                                                     &
         &                       *self%betaDarkMatter                                        &
         &                       *sqrt(self%neutrinoFactor)                                  &
         &                       *self%wavenumberEffective                                   &
         &                      )
    self%C                      =+ 14.4d0                                                    & ! Eqn. 20.
         &                       +325.0d0                                                    &
         &                       /(                                                          &
         &                         + 1.0d0                                                   &
         &                         +60.5d0                                                   &
         &                         *self%wavenumberEffectivePow                              &
         &                        )
    ! Compute wavenumbers needed for horizon scale modes.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) then
       ! Compute effective q.
       wavenumberScaleFree      =+wavenumber                                                      &
            &                    *self%temperatureCMB27                                       **2 &
            &                    /self%cosmologyParameters_%OmegaMatter   (                  )    &
            &                    /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
       self%wavenumberNeutrino =+3.92d0                             &
            &                   *wavenumberScaleFree                &
            &                   *sqrt(self%neutrinoNumberEffective) &
            &                   /     self%neutrinoMassFraction
    else
       self%wavenumberNeutrino =+0.0d0
    end if
    return
  end subroutine eisensteinHu1999ComputeFactors

  double precision function eisensteinHu1999Value(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: suppressionNeutrino

    ! Compute common factors.
    call self%computeFactors(wavenumber)
    ! Evaluate the transfer function.
    eisensteinHu1999Value  =+  self%L                      & ! Zero baryon form of the transfer function (eqn. 18).
         &                  /(                             &
         &                    +self%L                      &
         &                    +self%C                      &
         &                    *self%wavenumberEffective**2 &
         &                   )
    ! Apply correction for scales close to horizon.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) then
       suppressionNeutrino=+1.0d0                                                                    &
            &              +(                                                                        &
            &                +1.2d0                                                                  &
            &                *self%neutrinoMassFraction   ** 0.64d0                                  &
            &                *self%neutrinoNumberEffective**(0.30d0+0.6d0*self%neutrinoMassFraction) &
            &               )                                                                        &
            &              /(                                                                        &
            &                +self%wavenumberNeutrino**(-1.6d0)                                      &
            &                +self%wavenumberNeutrino**(+0.8d0)                                      &
            &               )
       eisensteinHu1999Value=eisensteinHu1999Value*suppressionNeutrino
    end if
    return
  end function eisensteinHu1999Value

  double precision function eisensteinHu1999LogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: transferFunction             , dLdwavenumberEffective       , &
         &                                                               dCdwavenumberEffective       , suppressionNeutrino          , &
         &                                                               wavenumberEffectiveDerivative, suppressionNeutrinoDerivative, &
         &                                                               wavenumberNeutrinoDerivative

    ! Compute common factors.
    call self%computeFactors(wavenumber)
    ! Get the transfer function itself.
    transferFunction=self%value(wavenumber)
    ! Compute derivatives of transfer function terms.
    dCdwavenumberEffective               =-325.00d0                      &
         &                                * 60.50d0                      &
         &                                *  1.11d0                      &
         &                                /(                             &
         &                                  + 1.0d0                      &
         &                                  +60.5d0                      &
         &                                  *self%wavenumberEffectivePow &
         &                                 )**2                          &
         &                                *  self%wavenumberEffectivePow &
         &                                /  self%waveNumberEffective
    dLdwavenumberEffective               =+1.0d0                         &
         &                                /(                             &
         &                                  +exp(1.0d0)                  &
         &                                  /1.84d0                      &
         &                                  /self%betaDarkMatter         &
         &                                  /sqrt(self%neutrinoFactor)   &
         &                                  +self%wavenumberEffective    &
         &                                 )
    wavenumberEffectiveDerivative=+self%wavenumberEffective                                        &
         &                        /     wavenumber                                                 &
         &                        +self%wavenumberEffective                                        &
         &                        /self%shapeParameterEffective                                    &
         &                        *self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                        *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                        *(1.0d0-sqrt(self%neutrinoFactor))                               &
         &                        /(                                                               &
         &                          +1.0d0                                                         &
         &                          +(                                                             &
         &                            +0.43d0                                                      &
         &                            *wavenumber                                                  &
         &                            *self%distanceSoundWave                                      &
         &                           )**4                                                          &
         &                         )                                                           **2 &
         &                        *4.0d0                                                           &
         &                        *(                                                               &
         &                          +0.43d0                                                        &
         &                          *wavenumber                                                    &
         &                          *self%distanceSoundWave                                        &
         &                         )**3                                                            &
         &                        *0.43d0                                                          &
         &                        *self%distanceSoundWave
    ! Compute logarithmic derivative of transfer function.
    eisensteinHu1999LogarithmicDerivative=+(                                                                                                                       &
         &                                  +dLdwavenumberEffective                                                                                                &
         &                                  /(                                      + self%L               + self%C               *self%wavenumberEffective**2)    &
         &                                                                          - self%L                                                                       &
         &                                  *(+2.0d0*self%C*self%wavenumberEffective+dLdwavenumberEffective+dCdwavenumberEffective*self%wavenumberEffective**2)    &
         &                                  /(                                      + self%L               + self%C               *self%wavenumberEffective**2)**2 &
         &                                 )                                                                                                                       &
         &                                *wavenumberEffectiveDerivative                                                                                           &
         &                                *wavenumber                                                                                                              &
         &                                /transferFunction
    ! Apply correction for scales close to horizon.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) then
       suppressionNeutrino                  =+1.0d0                                                                    &
            &                                +(                                                                        &
            &                                  +1.2d0                                                                  &
            &                                  *self%neutrinoMassFraction   ** 0.64d0                                  &
            &                                  *self%neutrinoNumberEffective**(0.30d0+0.6d0*self%neutrinoMassFraction) &
            &                                 )                                                                        &
            &                                /(                                                                        &
            &                                  +self%wavenumberNeutrino**(-1.6d0)                                      &
            &                                  +self%wavenumberNeutrino**(+0.8d0)                                      &
            &                                 )
       suppressionNeutrinoDerivative        =-(                                                                        &
            &                                  +1.2d0                                                                  &
            &                                  *self%neutrinoMassFraction   ** 0.64d0                                  &
            &                                  *self%neutrinoNumberEffective**(0.30d0+0.6d0*self%neutrinoMassFraction) &
            &                                 )                                                                        &
            &                                /(                                                                        &
            &                                  +      self%wavenumberNeutrino**(-1.6d0)                                &
            &                                  +      self%wavenumberNeutrino**(+0.8d0)                                &
            &                                 )**2                                                                     &
            &                                *(                                                                        &
            &                                  -1.6d0*self%wavenumberNeutrino**(-2.6d0)                                &
            &                                  +0.8d0*self%wavenumberNeutrino**(-0.2d0)                                &
            &                                 )
       wavenumberNeutrinoDerivative         =+3.92d0&
            &                                *sqrt(self%neutrinoNumberEffective)                              &
            &                                /     self%neutrinoMassFraction                                  &
            &                                *self%temperatureCMB27                                       **2 &
            &                                /self%cosmologyParameters_%OmegaMatter   (                  )    &
            &                                /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
       eisensteinHu1999LogarithmicDerivative=+eisensteinHu1999LogarithmicDerivative   &
            &                                *suppressionNeutrino                     &
            &                                +suppressionNeutrinoDerivative           &
            &                                *wavenumberNeutrinoDerivative            &
            &                                *wavenumber
    end if
    return
  end function eisensteinHu1999LogarithmicDerivative

  double precision function eisensteinHu1999HalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionEisensteinHu1999), intent(inout), target   :: self
    integer                                  , intent(  out), optional :: status
    !$GLC attributes unused :: self

    eisensteinHu1999HalfModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function eisensteinHu1999HalfModeMass

  double precision function eisensteinHu1999QuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionEisensteinHu1999), intent(inout), target   :: self
    integer                                  , intent(  out), optional :: status
    !$GLC attributes unused :: self

    eisensteinHu1999QuarterModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function eisensteinHu1999QuarterModeMass

  double precision function eisensteinHu1999EpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionEisensteinHu1999), intent(inout) :: self

    eisensteinHu1999EpochTime=self%time
    return
  end function eisensteinHu1999EpochTime
