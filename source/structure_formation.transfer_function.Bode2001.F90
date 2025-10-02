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
Implements a transfer function class based on the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.
!!}

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <enumeration>
   <name>scaleCutOffModel</name>
   <description>The model to use for the cut-off scale in the \cite{bode_halo_2001} transfer function.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="bode2001"              />
   <entry label="barkana2001"           />
   <entry label="viel05"                />
   <entry label="vogel23SpinHalf"       />
   <entry label="vogel23SpinThreeHalves"/>
  </enumeration>
  !!]

  !![
  <transferFunction name="transferFunctionBode2001">
   <description>
    A transfer function class which applies the thermal \gls{wdm} modifier of \cite{bode_halo_2001} to the provided \gls{cdm}
    transfer function.  The modifier is given by:
    \begin{equation}
    T(k) \rightarrow T(k) (1+[\epsilon k R_\mathrm{c}^0]^{2\nu})^{-\eta/\nu},
    \end{equation}
    where $\epsilon=${\normalfont \ttfamily [epsilon]}, $\eta=${\normalfont \ttfamily [eta]}, $\nu=${\normalfont \ttfamily
    [nu]}. The cut-off scale is computed from the dark matter particle (which must be of the {\normalfont \ttfamily
    darkMatterParticleWDMThermal} class) properties.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionBode2001
     !!{
     A transfer function class which modifies another transfer function using the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.
     !!}
     private
     type            (enumerationScaleCutOffModelType)          :: scaleCutOffModel
     double precision                                           :: epsilon                       , eta        , &
          &                                                        nu                            , scaleCutOff, &
          &                                                        time                          , redshift
     class           (transferFunctionClass          ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass        ), pointer :: darkMatterParticle_  => null()
   contains
     !![
     <methods>
      <method description="Compute the wavenumber at which the transfer function is suppressed by the given factor relative to the large-scale value." method="wavenumberAtSuppression"/>
     </methods>
     !!]
     final     ::                            bode2001Destructor
     procedure :: value                   => bode2001Value
     procedure :: logarithmicDerivative   => bode2001LogarithmicDerivative
     procedure :: halfModeMass            => bode2001HalfModeMass
     procedure :: quarterModeMass         => bode2001QuarterModeMass
     procedure :: fractionModeMass        => bode2001FractionModeMass
     procedure :: epochTime               => bode2001EpochTime
     procedure :: wavenumberAtSuppression => bode2001WavenumberAtSupression
  end type transferFunctionBode2001

  interface transferFunctionBode2001
     !!{
     Constructors for the \refClass{transferFunctionBode2001} transfer function class.
     !!}
     module procedure bode2001ConstructorParameters
     module procedure bode2001ConstructorInternal
  end interface transferFunctionBode2001

contains

  function bode2001ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionBode2001} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Error                         , only : Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionBode2001)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_
    type            (varying_string          )                :: scaleCutOffModel
    double precision                                          :: epsilon             , eta     , &
         &                                                       nu                  , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
    <inputParameter>
      <name>scaleCutOffModel</name>
      <source>parameters</source>
      <defaultValue>var_str('barkana2001')</defaultValue>
      <description>The model to use to compute the cut-off scale, either ``{\normalfont \ttfamily bode2001}'' to use the fitting function given by equation~A9 of \cite{bode_halo_2001}, ``{\normalfont \ttfamily barkana2001}'' to use the fitting function given by equation~(4) of \cite{barkana_constraints_2001}, ``{\normalfont \ttfamily viel05}'' to use the fitting function given by equation~(7) of \cite{viel_constraining_2005}, or ``{\normalfont \ttfamily vogel23spinHalf}'' or ``{\normalfont \ttfamily vogel23spinThreeHalves}'' to use the fitting function given by equation~(9) of \cite{vogel_entering_2023} for spin-1/2 or spin-3/2 particles respectively.</description>
    </inputParameter>
    <inputParameter>
      <name>epsilon</name>
      <source>parameters</source>
      <defaultValue>0.359d0</defaultValue>
      <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
      <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    </inputParameter>
    <inputParameter>
      <name>eta</name>
      <source>parameters</source>
      <defaultValue>3.81d0</defaultValue>
      <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
      <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    </inputParameter>
    <inputParameter>
      <name>nu</name>
      <source>parameters</source>
      <defaultValue>1.1d0</defaultValue>
      <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
      <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
      <description>The redshift of the epoch at which the transfer function is defined.</description>
    </inputParameter>
    !!]
    self=transferFunctionBode2001(transferFunctionCDM,enumerationScaleCutOffModelEncode(char(scaleCutOffModel),includesPrefix=.false.),epsilon,eta,nu,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="transferFunctionCDM" />
    !!]
    return
  end function bode2001ConstructorParameters

  function bode2001ConstructorInternal(transferFunctionCDM,scaleCutOffModel,epsilon,eta,nu,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionBode2001} transfer function class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Error                , only : Error_Report
    implicit none
    type            (transferFunctionBode2001       )                           :: self
    class           (transferFunctionClass          ), target   , intent(in   ) :: transferFunctionCDM
    type            (enumerationScaleCutOffModelType)           , intent(in   ) :: scaleCutOffModel
    double precision                                            , intent(in   ) :: epsilon                   , eta                            , &
         &                                                                         nu                        , time
    class           (cosmologyParametersClass       ), target   , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), target   , intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass        ), target   , intent(in   ) :: darkMatterParticle_
    double precision                                 , parameter                :: massReference       =1.0d0, degreesOfFreedomReference=1.5d0
    !![
    <constructorAssign variables="*transferFunctionCDM, scaleCutOffModel, epsilon, eta, nu, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>
    !!]

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Compute the comoving cut-off scale.
    select type (particle => self%darkMatterParticle_)
    class is (darkMatterParticleWDMThermal)
       select case (scaleCutOffModel%ID)
       case (scaleCutOffModelBarkana2001           %ID)
          ! This uses equation (4) from Barkana et al. (2001; http://adsabs.harvard.edu/abs/2001ApJ...558..482B), with the
          ! prefactor of 0.932 to give the cut-off scale at the epoch of matter-radiation equality as discussed in the paragraph
          ! following their equation (4).
          self%scaleCutOff=+0.932d0                                                                  &
               &           *0.201d0                                                                  &
               &           *(                                                                        &
               &             +(                                                                      &
               &               +self%cosmologyParameters_%OmegaMatter   (                  )         &
               &               -self%cosmologyParameters_%OmegaBaryon   (                  )         &
               &              )                                                                      &
               &             *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2      &
               &             /0.15d0                                                                 &
               &            )                                                               **0.15d0 &
               &           /(particle%degreesOfFreedomEffective()/degreesOfFreedomReference)**0.29d0 &
               &           /(particle%mass                     ()/            massReference)**1.15d0
       case (scaleCutOffModelBode2001              %ID)
          ! This uses equation (A9) from Bode et al. (2001; https://ui.adsabs.harvard.edu/abs/2001ApJ...556...93B).
          self%scaleCutOff=+0.048d0                                                                  &
               &           *(                                                                        &
               &             +(                                                                      &
               &               +self%cosmologyParameters_%OmegaMatter   (                  )         &
               &               -self%cosmologyParameters_%OmegaBaryon   (                  )         &
               &              )                                                                      &
               &             /0.4d0                                                                  &
               &            )**0.15d0                                                                &
               &           *(                                                                        &
               &             +self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)           &
               &             /0.65d0                                                                 &
               &            )**1.3d0                                                                 &
               &           /  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)           &
               &           /(particle%degreesOfFreedomEffective()/degreesOfFreedomReference)**0.29d0 &
               &           /(particle%mass                     ()/            massReference)**1.15d0
       case (scaleCutOffModelViel05                %ID)
          ! This uses equation (7) from Viel et al. (2005; https://ui.adsabs.harvard.edu/abs/2005PhRvD..71f3534V).
          self%scaleCutOff=+0.049d0                                                                  &
               &           *(                                                                        &
               &             +(                                                                      &
               &               +self%cosmologyParameters_%OmegaMatter   (                  )         &
               &               -self%cosmologyParameters_%OmegaBaryon   (                  )         &
               &              )                                                                      &
               &             /0.25d0                                                                 &
               &            )**0.11d0                                                                &
               &           *(                                                                        &
               &             +self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)           &
               &             /0.7d0                                                                  &
               &            )**1.22d0                                                                &
               &           /  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)           &
               &           /(particle%mass                     ()/            massReference)**1.11d0

       case (scaleCutOffModelVogel23SpinHalf       %ID)
          ! This uses equation (9) from Vogel & Azabajian (2023; https://ui.adsabs.harvard.edu/abs/2023PhRvD.108d3520V) with parameters for spin-1/2 particles.
          self%scaleCutOff=+0.0437d0                                                                  &
               &           *(                                                                         &
               &             +(                                                                       &
               &               +self%cosmologyParameters_%OmegaMatter   (                  )          &
               &               -self%cosmologyParameters_%OmegaBaryon   (                  )          &
               &              )                                                                       &
               &             *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2       &
               &             /0.12d0                                                                  &
               &            )**0.2463d0                                                               &
               &           *(                                                                         &
               &             +  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)          &
               &             /0.6736d0                                                                &
               &            )**2.012d0                                                                &
               &           /    self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)          &
               &           /(particle%mass                     ()/            massReference)**1.188d0
       case (scaleCutOffModelVogel23SpinThreeHalves%ID)
          ! This uses equation (9) from Vogel & Azabajian (2023; https://ui.adsabs.harvard.edu/abs/2023PhRvD.108d3520V) with parameters for spin-3/2 particles.
          self%scaleCutOff=+0.0345d0                                                                  &
               &           *(                                                                         &
               &             +(                                                                       &
               &               +self%cosmologyParameters_%OmegaMatter   (                  )          &
               &               -self%cosmologyParameters_%OmegaBaryon   (                  )          &
               &              )                                                                       &
               &             *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2       &
               &             /0.12d0                                                                  &
               &            )**0.2463d0                                                               &
               &           *(                                                                         &
               &             +  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)          &
               &             /0.6736d0                                                                &
               &            )**2.012d0                                                                &
               &           /    self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)          &
               &           /(particle%mass                     ()/            massReference)**1.195d0
       case default
          call Error_Report('invalid cut-off scale model'//{introspection:location})
       end select
    class default
       call Error_Report('transfer function expects a thermal warm dark matter particle'//{introspection:location})
    end select
    return
  end function bode2001ConstructorInternal

  subroutine bode2001Destructor(self)
    !!{
    Destructor for the \refClass{transferFunctionBode2001} transfer function class.
    !!}
    implicit none
    type(transferFunctionBode2001), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterParticle_" />
    <objectDestructor name="self%transferFunctionCDM" />
    !!]
    return
  end subroutine bode2001Destructor

  double precision function bode2001Value(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: wavenumberDimensionless

    bode2001Value=+self%transferFunctionCDM%value(wavenumber)
    if (self%scaleCutOff > 0.0d0) then
       wavenumberDimensionless=+     wavenumber  &
            &                  *self%scaleCutOff &
            &                  *self%epsilon
       if (exponent(wavenumberDimensionless)*2.0d0*self%eta < maxExponent(0.0d0)) then
          bode2001Value=+bode2001Value                                 &
               &        /(                                             &
               &          +1.0d0                                       &
               &          +wavenumberDimensionless**(2.0d0   *self%nu) &
               &         )                        **(self%eta/self%nu)
       else
          bode2001Value=+0.0d0
       end if
    end if
    return
  end function bode2001Value

  double precision function bode2001LogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: wavenumberDimensionless

    bode2001LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%scaleCutOff > 0.0d0) then
       wavenumberDimensionless=+     wavenumber  &
            &                  *self%scaleCutOff &
            &                  *self%epsilon
       if (exponent(wavenumberDimensionless)*2.0d0*self%eta < maxExponent(0.0d0)) then          
          bode2001LogarithmicDerivative=+bode2001LogarithmicDerivative                &
               &                        +2.0d0                                        &
               &                        *self%eta                                     &
               &                        *(                                            &
               &                          +1.0d0                                      &
               &                          /(                                          &
               &                            +1.0d0                                    &
               &                            +wavenumberDimensionless**(2.0d0*self%nu) &
               &                           )                                          &
               &                          -1.0d0                                      &
               &                         )
       else
          bode2001LogarithmicDerivative=+0.0d0
       end if
    end if
    return
  end function bode2001LogarithmicDerivative

  double precision function bode2001WavenumberAtSupression(self,factorSuppression)
    !!{
    Compute the wavenumber at which the transfer function is suppressed by the given factor relative to the large-scale value.
    !!}
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: factorSuppression

    bode2001WavenumberAtSupression=+(                                        &
         &                           +factorSuppression**(+self%nu/self%eta) &
         &                           -1.0d0                                  &
         &                          )                  **(+0.5d0  /self%nu ) &
         &                          /self%epsilon                            &
         &                          /self%scaleCutOff
    return
  end function bode2001WavenumberAtSupression
  
  double precision function bode2001HalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionBode2001), intent(inout), target   :: self
    integer                                   , intent(  out), optional :: status
    double precision                                                    :: matterDensity

    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of \cite{schneider_non-linear_2012}].
    bode2001HalfModeMass=+4.0d0                                 &
         &               *Pi                                    &
         &               /3.0d0                                 &
         &               *matterDensity                         &
         &               *(                                     &
         &                 +Pi                                  &
         &                 /self%wavenumberAtSuppression(2.0d0) &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function bode2001HalfModeMass

  double precision function bode2001QuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionBode2001), intent(inout), target   :: self
    integer                                   , intent(  out), optional :: status
    double precision                                                    :: matterDensity

    matterDensity          =+self%cosmologyParameters_%OmegaMatter    () &
         &                  *self%cosmologyParameters_%densityCritical()
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of \cite{schneider_non-linear_2012}].
    bode2001QuarterModeMass=+4.0d0                                 &
         &                  *Pi                                    &
         &                  /3.0d0                                 &
         &                  *matterDensity                         &
         &                  *(                                     &
         &                    +Pi                                  &
         &                    /self%wavenumberAtSuppression(4.0d0) &
         &                  )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function bode2001QuarterModeMass

  double precision function bode2001FractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is reduced by {\normalfont \ttfamily fraction} relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionBode2001), intent(inout), target   :: self
    double precision                          , intent(in   )           :: fraction
    integer                                   , intent(  out), optional :: status
    double precision                                                    :: matterDensity

    matterDensity           =+self%cosmologyParameters_%OmegaMatter    ()    &
         &                   *self%cosmologyParameters_%densityCritical()
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of \cite{schneider_non-linear_2012}].
    bode2001FractionModeMass=+4.0d0                                          &
         &                   *Pi                                             &
         &                   /3.0d0                                          &
         &                   *matterDensity                                  &
         &                   *(                                              &
         &                     +Pi                                           &
         &                     /self%wavenumberAtSuppression(1.0d0/fraction) &
         &                   )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function bode2001FractionModeMass

  double precision function bode2001EpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionBode2001), intent(inout) :: self

    bode2001EpochTime=self%time
    return
  end function bode2001EpochTime
