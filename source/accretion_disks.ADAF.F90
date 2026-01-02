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
  Implementation of an ADAF accretion disk.
  !!}

  use :: Tables, only : table1DLogarithmicLinear

  !![
  <enumeration>
   <name>adafRadiativeEfficiencyType</name>
   <description>Type of radiative efficiency model to use for ADAFs.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="fixed"   />
   <entry label="thinDisk"/>
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>adafViscosity</name>
   <description>Type of viscosity model to use for ADAFs.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="fit"  />
   <entry label="fixed"/>
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>adafFieldEnhancement</name>
   <description>Type of field enhancement model to use for ADAFs.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="exponential"/>
   <entry label="linear"     />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>adafEnergy</name>
   <description>Type of energy model to use for ADAFs.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="pureADAF"/>
   <entry label="ISCO"    />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>adafTable</name>
   <description>Enumeration of ADAF look-up tables.</description>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="powerJet"  />
   <entry label="rateSpinUp"/>
  </enumeration>
  !!]

  !![
  <accretionDisks name="accretionDisksADAF">
   <description>
    A circumnuclear accretion disk class, in which accretion is via an \gls{adaf} \citep{narayan_advection-dominated_1994}
    which is radiatively inefficient and geometrically thick. The radiative efficiency of the flow, which will be zero for a
    pure \gls{adaf}, is controlled by {\normalfont \ttfamily [efficiencyRadiationType]}. If set to {\normalfont \ttfamily
    fixed}, then the radiative efficiency is set to the value of the input parameter {\normalfont \ttfamily
    [efficiencyRadiation]}. Alternatively, if set to {\normalfont \ttfamily thinDisk} the radiative efficiency will be set to
    that of a Shakura-Sunyaev thin disk. The spin up rate of the black hole and the jet power produced as material accretes
    into the black hole are computed using the method of \cite{benson_maximum_2009}. The maximum efficiency of the jet (in
    units of the accretion power $\dot{M} \mathrm{c}^2$) is set by {\normalfont \ttfamily [efficiencyJetMaximum]}---in the
    model of \cite{benson_maximum_2009} the jet efficiency diverges as $j\rightarrow 1$, setting a maximum is important to
    avoid numerical instabilities. The energy of the accreted material can be set equal to the energy at infinity (as expected
    for a pure \gls{adaf}) or the energy at the \gls{isco} by use of the {\normalfont \ttfamily [energyOption]} parameter (set
    to {\normalfont \ttfamily pureADAF} or {\normalfont \ttfamily ISCO} respectively). The \gls{adaf} structure is controlled
    by the adiabatic index, $\gamma$, and viscosity parameter, $\alpha$, which are specified via the {\normalfont \ttfamily
    [adiabaticIndex]} and {\normalfont \ttfamily [viscosityOption]} input parameters respectively. The field-enhancing shear,
    $g$, is computed using $g=\exp(\omega \tau)$ if {\normalfont \ttfamily [fieldEnhancementOption]} is set to ``exponential''
    where $\omega$ is the frame-dragging frequency and $\tau$ is the smaller of the radial inflow and azimuthal velocity
    timescales. If {\normalfont \ttfamily [fieldEnhancementOption]} is set to ``linear'' then the alternative version,
    $g=1+\omega \tau$ is used instead. {\normalfont \ttfamily [viscosityOption]} may be set to ``{\normalfont \ttfamily fit}'',
    in which case the fitting function for $\alpha$ as a function of black hole spin is used:
    \begin{eqnarray}
     \alpha(j)=0.015+0.02 j^4 &amp; \hbox{ if  }&amp; g=\exp(\omega\tau) \hbox{ and } E=E_\mathrm{ISCO}, \\
     \alpha(j)=0.025+0.08 j^4 &amp; \hbox{ if } &amp; g=1+\omega\tau \hbox{ and } E=E_\mathrm{ISCO}, \\
     \alpha(j)=0.010+0.00 j^4 &amp; \hbox{ if } &amp; g=\exp(\omega\tau) \hbox{ and } E=1, \\
     \alpha(j)=0.025+0.02 j^4 &amp; \hbox{ if } &amp; g=1+\omega\tau \hbox{ and } E=1.  
    \end{eqnarray}
   </description>
  </accretionDisks>
  !!]
  type, extends(accretionDisksClass) :: accretionDisksADAF
     !!{
     Implementation of an ADAF accretion disk class.
     !!}
     private
     type            (enumerationADAFRadiativeEfficiencyTypeType) :: efficiencyRadiationType
     type            (enumerationADAFViscosityType              ) :: viscosityOption
     type            (enumerationADAFFieldEnhancementType       ) :: fieldEnhancementOption
     type            (enumerationADAFEnergyType                 ) :: energyOption     
     double precision                                             :: efficiencyRadiation                    , adiabaticIndex                       , &
          &                                                          pressureThermalFractional              , efficiencyJetMaximum                 , &
          &                                                          viscosityAlpha
     type            (table1DLogarithmicLinear                  ) :: tabulations
     type            (accretionDisksShakuraSunyaev              ) :: thinDisk
     ! Stored solutions.
     !! Height.
     double precision                                             :: heightStored                           , heightSpinPrevious                   , &
          &                                                          heightRadiusPrevious
     !! Velocity.
     double precision                                             :: velocityRadiusPrevious                 , velocitySpinPrevious                 , &
          &                                                          velocityStored
     !! Temperature.
     double precision                                             :: temperatureRadiusPrevious              , temperatureSpinPrevious              , &
          &                                                          temperatureStored
     !! Enthalpy.
     double precision                                             :: enthalpyRadiusPrevious                 , enthalpySpinPrevious                 , &
          &                                                          enthalpyStored
     !! Enthalpy-angular momentum product.
     double precision                                             :: enthlpyAngMPrdctRadiusPrevious         , enthlpyAngMPrdctSpinPrevious         , &
          &                                                          enthlpyAngMPrdctStored
     !! Angular momentum product.
     double precision                                             :: angularMomentumRadiusPrevious         , angularMomentumSpinPrevious          , &
          &                                                          angularMomentumStored
     !! Radial gamma factor.
     double precision                                             :: gammaRadialRadiusPrevious              , gammaRadialSpinPrevious              , &
          &                                                          gammaRadialStored
     !! Azimuthal gamma factor.
     double precision                                             :: gammaAzimuthalRadiusPrevious           , gammaAzimuthalSpinPrevious           , &
          &                                                          gammaAzimuthalStored
     !! Total gamma factor.
     double precision                                             :: gammaRadiusPrevious                    , gammaSpinPrevious                    , &
          &                                                          gammaStored
     !! Viscosity parameter.
     double precision                                             :: alphaStored                            , alphaSpinPrevious
     !! Fluid angular velocity.
     double precision                                             :: fluidAngularVelocityRadiusPrevious     , fluidAngularVelocitySpinPrevious     , &
          &                                                          fluidAngularVelocityStored
     !! Field enhancement.
     double precision                                             :: fieldEnhancementRadiusPrevious         , fieldEnhancementSpinPrevious         , &
          &                                                          fieldEnhancementStored
     !! Black hole-launched jet power.
     double precision                                             :: jetPowerBlackHoleRadiusPrevious        , jetPowerBlackHoleSpinPrevious        , &
          &                                                          jetPowerBlackHoleStored
     !! Disk-launched jet power.
     double precision                                             :: jetPowerDiskRadiusPrevious             , jetPowerDiskSpinPrevious             , &
          &                                                          jetPowerDiskStored
     !! Disk-launched black hole jet power.
     double precision                                             :: jetPowerDiskFromBlackHoleRadiusPrevious, jetPowerDiskFromBlackHoleSpinPrevious, &
          &                                                          jetPowerDiskFromBlackHoleStored
   contains
     !![
     <methods>
       <method description="Return the dimensionless height of the ADAF at a given radius." method="height" />
       <method description="Return the dimensionless velocity of the ADAF at a given radius." method="velocity" />
       <method description="Return the dimensionless temperature of the ADAF at a given radius." method="temperature" />
       <method description="Return the dimensionless enthalpy of the ADAF at a given radius." method="enthalpy" />
       <method description="Return the product of dimensionless enthalpy and angular momentum of the ADAF at a given radius." method="enthalpyAngularMomentumProduct" />
       <method description="Return the dimensionless angular momentum of the ADAF at a given radius." method="angularMomentum" />
       <method description="Return the radial part of the relativistic boost factor in the ADAF at a given radius." method="gammaRadial" />
       <method description="Return the azimuthal part of the relativistic boost factor in the ADAF at a given radius." method="gammaAzimuthal" />
       <method description="Return the relativistic boost factor in the ADAF at a given radius." method="gamma" />
       <method description="Return the viscosity parameter, $\alpha$, in the ADAF." method="viscosityParameter" />
       <method description="Return the dimensionless angular velocity of the ADAF fluid at the given radius." method="fluidAngularVelocity" />
       <method description="Return the magnetic field enhancement factor in the ADAF at the given radius." method="fieldEnhancement" />
       <method description="Return the power of the jet launched by the black hole for the ADAF." method="jetPowerBlackHole" />
       <method description="Return the power of the jet launched from the disk for the ADAF." method="jetPowerDisk" />
       <method description="Return the power of the jet launched from the disk which is derived frm the black hole." method="jetPowerDiskFromBlackHole" />
     </methods>
     !!]
     procedure :: efficiencyRadiative            => adafEfficiencyRadiative
     procedure :: powerJet                       => adafPowerJet
     procedure :: rateSpinUp                     => adafRateSpinUp
     procedure :: height                         => adafHeight
     procedure :: velocity                       => adafVelocity
     procedure :: temperature                    => adafTemperature
     procedure :: enthalpy                       => adafEnthalpy
     procedure :: enthalpyAngularMomentumProduct => adafEnthalpyAngularMomentumProduct
     procedure :: angularMomentum                => adafAngularMomentum
     procedure :: gammaRadial                    => adafGammaRadial
     procedure :: gammaAzimuthal                 => adafGammaAzimuthal
     procedure :: gamma                          => adafGamma
     procedure :: viscosityParameter             => adafViscosityParameter
     procedure :: fluidAngularVelocity           => adafFluidAngularVelocity
     procedure :: fieldEnhancement               => adafFieldEnhancement
     procedure :: jetPowerBlackHole              => adafJetPowerBlackHole
     procedure :: jetPowerDisk                   => adafJetPowerDisk
     procedure :: jetPowerDiskFromBlackHole      => adafJetPowerDiskFromBlackHole
  end type accretionDisksADAF

  interface accretionDisksADAF
     !!{
     Constructors for the ADAF accretion disk class.
     !!}
     module procedure adafConstructorParameters
     module procedure adafConstructorInternal
  end interface accretionDisksADAF

  ! Number of points to use in ADAF look-up tables.
  integer, parameter :: countTable=10000

contains

  function adafConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ADAF accretion disk class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionDisksADAF)                :: self
    type            (inputParameters   ), intent(inout) :: parameters
    type            (varying_string    )                :: efficiencyRadiationType, energyOption        , &
         &                                                 fieldEnhancementOption , viscosityOption
    double precision                                    :: efficiencyRadiation    , adiabaticIndex      , &
         &                                                 viscosityAlpha         , efficiencyJetMaximum

    !![
    <inputParameter>
      <name>efficiencyRadiationType</name>
      <source>parameters</source>
      <defaultValue>var_str('thinDisk')</defaultValue>
      <description>Specifies the specific energy of material at the inner edge of an ADAF. {\normalfont \ttfamily pureADAF} makes the specific energy equal
        to 1 (i.e. all energy is advected with the flow); {\normalfont \ttfamily ISCO} makes the specific energy equal to that for the innermost
        stable circular orbit.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiencyRadiation</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>Specifies the radiative efficiency of an ADAF (i.e. the fraction of $\dot{M}\clight^2$ that is emitted in radiation).</description>
    </inputParameter>
    <inputParameter>
      <name>energyOption</name>
      <source>parameters</source>
      <defaultValue>var_str('pureADAF')</defaultValue>
      <description>Specifies the specific energy of material at the inner edge of an ADAF. {\normalfont \ttfamily pureADAF} makes the specific energy equal
        to 1 (i.e. all energy is advected with the flow); {\normalfont \ttfamily ISCO} makes the specific energy equal to that for the innermost
        stable circular orbit.</description>
    </inputParameter>
    <inputParameter>
      <name>fieldEnhancementOption</name>
      <source>parameters</source>
      <defaultValue>var_str('exponential')</defaultValue>
      <description>Controls how the field enhancing shear is determined. {\normalfont \ttfamily exponential} will cause the form $g=\exp(\omega t)$ \citep{benson_maximum_2009}
       to be used, while {\normalfont \ttfamily linear} will cause $g=1+\omega t$ to be used instead. The functional form of $\alpha(j)$ (if used) will be adjusted
       to achieve a sensible spin-up function in each case.</description>
    </inputParameter>
    <inputParameter>
      <name>adiabaticIndex</name>
      <source>parameters</source>
      <defaultValue>adafAdiabaticIndexDefault(enumerationAdafFieldEnhancementEncode(char(fieldEnhancementOption),includesPrefix=.false.))</defaultValue>
      <description>Specifies the effective adiabatic index of gas in an ADAF.</description>
    </inputParameter>
    <inputParameter>
      <name>viscosityOption</name>
      <source>parameters</source>
      <defaultValue>var_str('fit')</defaultValue>
      <description>Controls how the viscosity parameter $\alpha$ in an ADAF is determined. {\normalfont \ttfamily fit} will cause $\alpha$ to be computed
       using the fitting function of \cite{benson_maximum_2009}; {\normalfont \ttfamily fixed} will cause $\alpha=${\normalfont \ttfamily [adafViscosityFixedAlpha]}
       to be used.</description>
    </inputParameter>
    <inputParameter>
      <name>viscosityAlpha</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The value for the viscosity parameter $\alpha$ in an ADAF to be used if {\normalfont \ttfamily [adafViscosityOption]}$=${\normalfont \ttfamily fixed}.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiencyJetMaximum</name>
      <source>parameters</source>
      <defaultValue>2.0d0</defaultValue>
      <description>The maximum efficiency allowed for ADAF-driven jets (in units of the accretion power).</description>
    </inputParameter>
    !!]
    self=accretionDisksADAF(                                                                                                    &
         &                  enumerationAdafEnergyEncode                 (char(energyOption           ),includesPrefix=.false.), &
         &                  enumerationAdafFieldEnhancementEncode       (char(fieldEnhancementOption ),includesPrefix=.false.), &
         &                  enumerationAdafRadiativeEfficiencyTypeEncode(char(efficiencyRadiationType),includesPrefix=.false.), &
         &                  enumerationAdafViscosityEncode              (char(viscosityOption        ),includesPrefix=.false.), &
         &                  efficiencyJetMaximum                                                                              , &
         &                  efficiencyRadiation                                                                               , &
         &                  adiabaticIndex                                                                                    , &
         &                  viscosityAlpha                                                                                      &
         &                 )
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adafConstructorParameters

  function adafConstructorInternal(energyOption,fieldEnhancementOption,efficiencyRadiationType,viscosityOption,efficiencyJetMaximum,efficiencyRadiation,adiabaticIndex,viscosityAlpha) result(self)
    !!{
    Internal constructor for the ADAF accretion disk class.
    !!}
    use :: Black_Hole_Fundamentals     , only : Black_Hole_ISCO_Radius, Black_Hole_ISCO_Specific_Energy, Black_Hole_Rotational_Energy_Spin_Down, Black_Hole_Static_Radius, &
          &                                     orbitPrograde
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Table_Labels                , only : extrapolationTypeFix
    implicit none
    type            (accretionDisksADAF                        )                          :: self
    type            (enumerationADAFRadiativeEfficiencyTypeType), intent(in   )           :: efficiencyRadiationType
    type            (enumerationADAFViscosityType              ), intent(in   )           :: viscosityOption
    type            (enumerationADAFFieldEnhancementType       ), intent(in   )           :: fieldEnhancementOption
    type            (enumerationADAFEnergyType                 ), intent(in   )           :: energyOption     
    double precision                                            , intent(in   )           :: efficiencyJetMaximum
    double precision                                            , intent(in   ), optional :: efficiencyRadiation              , adiabaticIndex                    , &
         &                                                                                   viscosityAlpha
    double precision                                            , parameter               :: spinBlackHoleInverseMaximum=1.0d0, spinBlackHoleInverseMinimum=1.0d-6
    integer                                                                               :: iSpin
    double precision                                                                      :: adafEnergy                       , radiusStatic                      , &
         &                                                                                   spinBlackHole                    , radiusISCO

    ! Validate arguments.
    if (efficiencyRadiationType == adafRadiativeEfficiencyTypeFixed .and. .not.present(efficiencyRadiation)) &
         & call Error_Report('radiation efficiency must be provided'//{introspection:location})
    if (viscosityOption         == adafViscosityFixed               .and. .not.present(viscosityAlpha     )) &
         & call Error_Report('viscosity parameter must be provided'//{introspection:location})
    ! Make assignments.
    !![
    <constructorAssign variables="energyOption, fieldEnhancementOption, efficiencyRadiationType, viscosityOption, efficiencyJetMaximum, efficiencyRadiation, adiabaticIndex, viscosityAlpha"/>
    !!]
    ! Set the default adiabatic index if none was provided.
    if (.not.present(adiabaticIndex)) self%adiabaticIndex=adafAdiabaticIndexDefault(fieldEnhancementOption)
    ! Set the thermal pressure fraction.
    self%pressureThermalFractional=(8.0d0-6.0d0*self%adiabaticIndex)/3.0d0/(1.0d0-self%adiabaticIndex)
    ! Initialize stored solutions.
    self%heightStored                           =-     1.0d0
    self%heightSpinPrevious                     =-huge(1.0d0)
    self%heightRadiusPrevious                   =-     1.0d0
    self%velocityRadiusPrevious                 =-     1.0d0
    self%velocitySpinPrevious                   =-huge(1.0d0)
    self%velocityStored                         =-     1.0d0
    self%temperatureRadiusPrevious              =-     1.0d0
    self%temperatureSpinPrevious                =-huge(1.0d0)
    self%temperatureStored                      =-     1.0d0
    self%enthalpyRadiusPrevious                 =-     1.0d0
    self%enthalpySpinPrevious                   =-huge(1.0d0)
    self%enthalpyStored                         =-     1.0d0
    self%enthlpyAngMPrdctRadiusPrevious         =-     1.0d0
    self%enthlpyAngMPrdctSpinPrevious           =-huge(1.0d0)
    self%enthlpyAngMPrdctStored                 =-     1.0d0
    self%angularMomentumRadiusPrevious          =-     1.0d0
    self%angularMomentumSpinPrevious            =-huge(1.0d0)
    self%angularMomentumStored                  =-     1.0d0
    self%gammaRadialRadiusPrevious              =-     1.0d0
    self%gammaRadialSpinPrevious                =-huge(1.0d0)
    self%gammaRadialStored                      =-     1.0d0
    self%gammaAzimuthalRadiusPrevious           =-     1.0d0
    self%gammaAzimuthalSpinPrevious             =-huge(1.0d0)
    self%gammaAzimuthalStored                   =-     1.0d0
    self%gammaRadiusPrevious                    =-     1.0d0
    self%gammaSpinPrevious                      =-huge(1.0d0)
    self%gammaStored                            =-     1.0d0
    self%alphaSpinPrevious                      =-huge(1.0d0)
    self%alphaStored                            =-     1.0d0
    self%fluidAngularVelocityRadiusPrevious     =-     1.0d0
    self%fluidAngularVelocitySpinPrevious       =-huge(1.0d0)
    self%fluidAngularVelocityStored             =-     1.0d0
    self%fieldEnhancementRadiusPrevious         =-     1.0d0
    self%fieldEnhancementSpinPrevious           =-huge(1.0d0)
    self%fieldEnhancementStored                 =-     1.0d0
    self%jetPowerBlackHoleRadiusPrevious        =-     1.0d0
    self%jetPowerBlackHoleSpinPrevious          =-huge(1.0d0)
    self%jetPowerBlackHoleStored                =-     1.0d0
    self%jetPowerDiskRadiusPrevious             =-     1.0d0
    self%jetPowerDiskSpinPrevious               =-huge(1.0d0)
    self%jetPowerDiskStored                     =-     1.0d0
    self%jetPowerDiskFromBlackHoleRadiusPrevious=-     1.0d0
    self%jetPowerDiskFromBlackHoleSpinPrevious  =-huge(1.0d0)
    self%jetPowerDiskFromBlackHoleStored        =-     1.0d0
    ! Build tabulations of jet power and spin.
    call self%tabulations%destroy()
    call self%tabulations%create (                                                                           &
         &                        spinBlackHoleInverseMinimum                                              , &
         &                        spinBlackHoleInverseMaximum                                              , &
         &                        countTable                                                               , &
         &                        tableCount                   =2                                          , &
         &                        extrapolationType            =[extrapolationTypeFix,extrapolationTypeFix]  &
         &                       )
    do iSpin=1,countTable
       ! Get the black hole spin. The "inverse spin parameter" that we tabulate in is 1-j so that we can easily pack
       ! many points close to j=1.
       spinBlackHole=1.0d0-self%tabulations%x(iSpin)
       ! Determine the ADAF energy.
       select case (self%energyOption%ID)
       case (adafEnergyPureADAF%ID)
          adafEnergy=1.0d0
       case (adafEnergyISCO    %ID)
          adafEnergy=Black_Hole_ISCO_Specific_Energy(spinBlackHole,orbitPrograde)
       case default
          adafEnergy=0.0d0
          call Error_Report('unknown energy type'//{introspection:location})
       end select
       ! Compute jet launch radii.
       radiusISCO  =Black_Hole_ISCO_Radius  (spinBlackHole)
       radiusStatic=Black_Hole_Static_Radius(spinBlackHole)
       ! Compute the jet power.
       call self%tabulations%populate(                                                                           &
            &                               +min(                                                                &
            &                                    (                                                               &
            &                                     +self%jetPowerBlackHole          (spinBlackHole,radiusStatic)  &
            &                                     +self%jetPowerDisk               (spinBlackHole,radiusISCO  )  &
            &                                    ),                                                              &
            &                                    self%efficiencyJetMaximum                                       &
            &                                   )                                                                &
            &                               *(speedLight/kilo)**2                                              , &
            &                               iSpin                                                              , &
            &                         table=adafTablePowerJet%ID                                                 &
            &                 )
       ! Compute the rate of spin up to mass rate of change ratio.
       call self%tabulations%populate(                                                                           &
            &                               +self%angularMomentum                  (spinBlackHole,radiusISCO  )  &
            &                               -2.0d0                                                               &
            &                               *                                       spinBlackHole                &
            &                               *adafEnergy                                                          &
            &                               -Black_Hole_Rotational_Energy_Spin_Down(spinBlackHole             )  &
            &                               *(                                                                   &
            &                                 +self%jetPowerBlackHole              (spinBlackHole,radiusStatic)  &
            &                                 +self%jetPowerDiskFromBlackHole      (spinBlackHole,radiusISCO  )  &
            &                                )                                                                 , &
            &                               iSpin                                                              , &
            &                         table=adafTableRateSpinUp%ID                                               &
            &                        )
    end do
    ! Initialize the thin disk component used for radiative efficiency calculations.
    self%thinDisk=accretionDisksShakuraSunyaev()
    return
  end function adafConstructorInternal

  double precision function adafEfficiencyRadiative(self,blackHole,accretionRateMass)
    !!{
    Computes the radiative efficiency for an ADAF.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (accretionDisksADAF    ), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass

    select case (self%efficiencyRadiationType%ID)
    case (adafRadiativeEfficiencyTypeFixed   %ID)
       adafEfficiencyRadiative=self%efficiencyRadiation
    case (adafRadiativeEfficiencyTypeThinDisk%ID)
       adafEfficiencyRadiative=self%thinDisk%efficiencyRadiative(blackHole,accretionRateMass)
    case default
       adafEfficiencyRadiative=0.0d0
       call Error_Report('unknown radiative efficiency type'//{introspection:location})
    end select
    return
  end function adafEfficiencyRadiative

  double precision function adafPowerJet(self,blackHole,accretionRateMass)
    !!{
    Computes the jet power of the given black hole in due to accretion from an ADAF disk.
    !!}
    implicit none
    class           (accretionDisksADAF    ), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                                        :: spinBlackHole    , spinInverseBlackHole

    ! Get the black hole spin.
    spinBlackHole       =blackHole%spin()
    ! Get the "inverse spin parameter".
    spinInverseBlackHole=+1.0d0-spinBlackHole
    ! Compute the jet power.
    adafPowerJet        =+accretionRateMass                                                       &
         &               *self%tabulations%interpolate(spinInverseBlackHole,adafTablePowerJet%ID)
    return
  end function adafPowerJet

  double precision function adafRateSpinUp(self,blackHole,accretionRateMass)
    !!{
    Computes the spin up rate of the black hole in {\normalfont \ttfamily blackHole} due to accretion from an ADAF.
    disk.
    !!}
    implicit none
    class           (accretionDisksADAF    ), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                                        :: spinBlackHole              , spinInverseBlackHole, &
         &                                                     spinToMassRateOfChangeRatio

    ! Get the black hole spin.
    spinBlackHole              =blackHole%spin()
    ! Get the "inverse spin parameter".
    spinInverseBlackHole       =+1.0d0-spinBlackHole
    ! Compute the ratio of spin and mass rates of change.
    spinToMassRateOfChangeRatio=self%tabulations%interpolate(spinInverseBlackHole,adafTableRateSpinUp%ID)
    ! Scale to the mass rate of change.
    adafRateSpinUp             =+spinToMassRateOfChangeRatio &
         &                      *accretionRateMass           &
         &                      /blackHole%mass()
    return
  end function adafRateSpinUp

  double precision function adafJetPowerDiskFromBlackHole(self,spinBlackHole,radius)
    !!{
    Returns the power extracted from the black hole by the disk-launched jet from an ADAF.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%jetPowerDiskFromBlackHoleRadiusPrevious .or. spinBlackHole /= self%jetPowerDiskFromBlackHoleSpinPrevious) then
       self%jetPowerDiskFromBlackHoleStored        =+  self%jetPowerDisk    (spinBlackHole,radius)    &
            &                                       *(                                                &
            &                                         +1.0d0                                          &
            &                                         -1.0d0                                          &
            &                                         /self%fieldEnhancement(spinBlackHole,radius)**2 &
            &                                        )
       self%jetPowerDiskFromBlackHoleRadiusPrevious=+                                      radius
       self%jetPowerDiskFromBlackHoleSpinPrevious  =+                        spinBlackHole
    end if
    adafJetPowerDiskFromBlackHole=self%jetPowerDiskFromBlackHoleStored
    return
  end function adafJetPowerDiskFromBlackHole

  double precision function adafJetPowerDisk(self,spinBlackHole,radius)
    !!{
    Returns the power of the disk-launched jet from an ADAF.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Frame_Dragging_Frequency, Black_Hole_Metric_A_Factor, Black_Hole_Metric_D_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius       , spinBlackHole
    double precision                                    :: betaAzimuthal

    ! Check if arguments are the same as on the previous call.
    if (radius /= self%jetPowerDiskRadiusPrevious .or. spinBlackHole /= self%jetPowerDiskSpinPrevious) then
       ! They are not, so compute (and store) a new value.
       betaAzimuthal                  =+sqrt(                                                                     &
            &                                +1.0d0                                                               &
            &                                -1.0d0                                                               &
            &                                /  self%gammaAzimuthal                     (spinBlackHole,radius)**2 &
            &                               )
       self%jetPowerDiskStored        =+(                                                                         &
            &                            + 3.0d0                                                                  &
            &                            /80.0d0                                                                  &
            &                           )                                                                         &
            &                          *                                                               radius **2 &
            &                          *(                                                                         &
            &                            +2.0d0                                                                   &
            &                            *                                               spinBlackHole            &
            &                            *betaAzimuthal                                                           &
            &                            /                                                             radius **2 &
            &                            +sqrt(                                                                   &
            &                                  +     Black_Hole_Metric_D_Factor         (spinBlackHole,radius)    &
            &                                 )                                                                   &
            &                           )**2                                                                      &
            &                          *(                                                                         &
            &                            +1.0d0                                                                   &
            &                            -      self%pressureThermalFractional                                    &
            &                           )                                                                         &
            &                          *(                                                                         &
            &                            +      self%fieldEnhancement                   (spinBlackHole,radius)    &
            &                            *      self%gamma                              (spinBlackHole,radius)    &
            &                           )**2                                                                      &
            &                          *sqrt(                                                                     &
            &                                +(                                                                   &
            &                                  +1.0d0                                                             &
            &                                  -self%velocity                           (spinBlackHole,radius)**2 &
            &                                 )                                                                   &
            &                                /       Black_Hole_Metric_D_Factor         (spinBlackHole,radius)    &
            &                               )                                                                     &
            &                          *(                                                                         &
            &                            +      self%fluidAngularVelocity               (spinBlackHole,radius)    &
            &                            +           Black_Hole_Frame_Dragging_Frequency(spinBlackHole,radius)    &
            &                           )**2                                                                      &
            &                          *        self%temperature                        (spinBlackHole,radius)    &
            &                          /             Black_Hole_Metric_A_Factor         (spinBlackHole,radius)    &
            &                          /        self%velocity                           (spinBlackHole,radius)    &
            &                          /        self%height                             (spinBlackHole,radius)
       self%jetPowerDiskRadiusPrevious=+                                                               radius
       self%jetPowerDiskSpinPrevious  =+                                                 spinBlackHole
    end if
    adafJetPowerDisk=self%jetPowerDiskStored
    return
  end function adafJetPowerDisk

  double precision function adafJetPowerBlackHole(self,spinBlackHole,radius)
    !!{
    Returns the power of the black hole-launched jet from an ADAF.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Frame_Dragging_Frequency, Black_Hole_Metric_A_Factor, Black_Hole_Metric_D_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius                     , spinBlackHole
    double precision                    , parameter     :: spinBlackHoleMinimum=5.0d-8
    double precision                                    :: betaAzimuthal

    if (spinBlackHole > spinBlackHoleMinimum) then
       ! Check if arguments are the same as on the previous call.
       if (radius /= self%jetPowerBlackHoleRadiusPrevious .or. spinBlackHole /= self%jetPowerBlackHoleSpinPrevious) then
          ! They are not, so compute (and store) a new value.
          betaAzimuthal                       =+sqrt(                                                                     &
               &                                     +1.0d0                                                               &
               &                                     -1.0d0                                                               &
               &                                     /  self%gammaAzimuthal                     (spinBlackHole,radius)**2 &
               &                                    )
          self%jetPowerBlackHoleStored        =+(                                                                         &
               &                                 + 3.0d0                                                                  &
               &                                 /80.0d0                                                                  &
               &                                )                                                                         &
               &                               *                                                               radius **2 &
               &                               *(                                                                         &
               &                                 +2.0d0                                                                   &
               &                                 *                                               spinBlackHole            &
               &                                 *betaAzimuthal                                                           &
               &                                 /                                                             radius **2 &
               &                                 +sqrt(                                                                   &
               &                                       +     Black_Hole_Metric_D_Factor         (spinBlackHole,radius)    &
               &                                      )                                                                   &
               &                                )**2                                                                      &
               &                               *(                                                                         &
               &                                 +1.0d0                                                                   &
               &                                 -      self%pressureThermalFractional                                    &
               &                                )                                                                         &
               &                               *(                                                                         &
               &                                 +      self%fieldEnhancement                   (spinBlackHole,radius)    &
               &                                 *      self%gamma                              (spinBlackHole,radius)    &
               &                                )**2                                                                      &
               &                               *sqrt(                                                                     &
               &                                     +(                                                                   &
               &                                       +1.0d0                                                             &
               &                                       -self%velocity                           (spinBlackHole,radius)**2 &
               &                                      )                                                                   &
               &                                     /       Black_Hole_Metric_D_Factor         (spinBlackHole,radius)    &
               &                                    )                                                                     &
               &                               *             Black_Hole_Frame_Dragging_Frequency(spinBlackHole,radius)**2 &
               &                               *        self%temperature                        (spinBlackHole,radius)    &
               &                               /             Black_Hole_Metric_A_Factor         (spinBlackHole,radius)    &
               &                               /        self%velocity                           (spinBlackHole,radius)    &
               &                               /        self%height                             (spinBlackHole,radius)
          self%jetPowerBlackHoleRadiusPrevious=+                                                               radius
          self%jetPowerBlackHoleSpinPrevious  =+                                                 spinBlackHole
       end if
       adafJetPowerBlackHole=self%jetPowerBlackHoleStored
    else
       adafJetPowerBlackHole=0.0d0
    end if
    return
  end function adafJetPowerBlackHole

  double precision function adafFieldEnhancement(self,spinBlackHole,radius)
    !!{
    Returns the field enhancement factor, $g$, in the ADAF.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Frame_Dragging_Frequency, Black_Hole_Metric_D_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole
    double precision                                    :: timescale   , timescaleAzimuthal       , &
         &                                                 timescaleRadial

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%fieldEnhancementRadiusPrevious .or. spinBlackHole /= self%fieldEnhancementSpinPrevious) then
       timescaleAzimuthal=+1.0d0                                                       &
            &             /     self%fluidAngularVelocity      (spinBlackHole,radius)
       timescaleRadial   =+                                                   radius   &
            &             *     self%gammaAzimuthal            (spinBlackHole,radius)  &
            &             /     self%velocity                  (spinBlackHole,radius)  &
            &             /sqrt(     Black_Hole_Metric_D_Factor(spinBlackHole,radius))
       timescale         =min(timescaleAzimuthal,timescaleRadial)
       select case (self%fieldEnhancementOption%ID)
       case (adafFieldEnhancementExponential%ID)
          self%fieldEnhancementStored=      +exp(Black_Hole_Frame_Dragging_Frequency(spinBlackHole,radius)*timescale)
       case (adafFieldEnhancementLinear     %ID)
          self%fieldEnhancementStored=+1.0d0+    Black_Hole_Frame_Dragging_Frequency(spinBlackHole,radius)*timescale
       end select
       self%fieldEnhancementRadiusPrevious            =radius
       self%fieldEnhancementSpinPrevious     =spinBlackHole
    end if
    adafFieldEnhancement=self%fieldEnhancementStored
    return
  end function adafFieldEnhancement

  double precision function adafFluidAngularVelocity(self,spinBlackHole,radius)
    !!{
    Returns the angular velocity of the rotating fluid with respect to the local inertial observer (ZAMO).
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Metric_A_Factor, Black_Hole_Metric_D_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%fluidAngularVelocityRadiusPrevious .or. spinBlackHole /= self%fluidAngularVelocitySpinPrevious) then
       self%fluidAngularVelocityStored        =+     self%angularMomentum           (spinBlackHole,radius)    &
            &                                  *sqrt(                                                         &
            &                                             Black_Hole_Metric_D_Factor(spinBlackHole,radius)    &
            &                                        /    Black_Hole_Metric_A_Factor(spinBlackHole,radius)**3 &
            &                                       )                                                         &
            &                                  /                                                   radius **2 &
            &                                  /     self%gammaAzimuthal            (spinBlackHole,radius)    &
            &                                  /     self%gammaRadial               (spinBlackHole,radius)
       self%fluidAngularVelocityRadiusPrevious=+radius
       self%fluidAngularVelocitySpinPrevious  =+spinBlackHole
    end if
    adafFluidAngularVelocity=self%fluidAngularVelocityStored
    return
  end function adafFluidAngularVelocity

  double precision function adafViscosityParameter(self,spinBlackHole)
    !!{
    Returns the effective value of the $\alpha$ viscosity parameter for an ADAF.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (spinBlackHole /= self%alphaSpinPrevious) then
       select case (self%energyOption%ID)
       case (adafEnergyISCO%ID)
          select case (self%fieldEnhancementOption%ID)
          case (adafFieldEnhancementExponential%ID)
             self%alphaStored=0.015d0+0.02d0*spinBlackHole**4
          case (adafFieldEnhancementLinear     %ID)
             self%alphaStored=0.025d0+0.08d0*spinBlackHole**4
          end select
       case (adafEnergyPureADAF   %ID)
          select case (self%fieldEnhancementOption%ID)
          case (adafFieldEnhancementExponential%ID)
             self%alphaStored=0.010d0
          case (adafFieldEnhancementLinear     %ID)
             self%alphaStored=0.025d0+0.02d0*spinBlackHole**4
          end select
       end select
       self%alphaSpinPrevious=spinBlackHole
    end if
    adafViscosityParameter=self%alphaStored
    return
  end function adafViscosityParameter

  double precision function adafGamma(self,spinBlackHole,radius)
    !!{
    Returns the net relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    The input quantities are in natural units.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%gammaRadiusPrevious .or. spinBlackHole /= self%gammaSpinPrevious) then
       self%gammaStored        =+self%gammaRadial   (spinBlackHole,radius) &
            &                   *self%gammaAzimuthal(spinBlackHole,radius)
       self%gammaRadiusPrevious=+radius
       self%gammaSpinPrevious  =+spinBlackHole
    end if
    adafGamma=self%gammaStored
    return
  end function adafGamma

  double precision function adafGammaAzimuthal(self,spinBlackHole,radius)
    !!{
    Returns the $\phi$ component relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    The input quantities are in natural units.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Metric_A_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%gammaAzimuthalRadiusPrevious .or. spinBlackHole /= self%gammaAzimuthalSpinPrevious) then
       self%gammaAzimuthalStored        =+sqrt(                                                  &
            &                                  +1.0d0                                            &
            &                                  +(                                                &
            &                                    +(                                              &
            &                                      +self%angularMomentum  (spinBlackHole,radius) &
            &                                      /self%gammaRadial      (spinBlackHole,radius) &
            &                                      /                                     radius  &
            &                                     )**2                                           &
            &                                   )                                                &
            &                                  /Black_Hole_Metric_A_Factor(spinBlackHole,radius) &
            &                                 )
       self%gammaAzimuthalRadiusPrevious=+radius
       self%gammaAzimuthalSpinPrevious  =+spinBlackHole
    end if
    adafGammaAzimuthal=self%gammaAzimuthalStored
    return
  end function adafGammaAzimuthal

  double precision function adafGammaRadial(self,spinBlackHole,radius)
    !!{
    Returns the $r$ component relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    The input quantities are in natural units.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius, spinBlackHole

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%gammaRadialRadiusPrevious .or. spinBlackHole /= self%gammaRadialSpinPrevious) then
       self%gammaRadialStored        =+sqrt(                                          &
            &                               +1.0d0                                    &
            &                               /(                                        &
            &                                 +1.0d0                                  &
            &                                 -self%velocity(spinBlackHole,radius)**2 &
            &                                )                                        &
            &                              )
       self%gammaRadialRadiusPrevious=+radius
       self%gammaRadialSpinPrevious  =+spinBlackHole
    end if
    adafGammaRadial=self%gammaRadialStored
    return
  end function adafGammaRadial

  double precision function adafAngularMomentum(self,spinBlackHole,radius)
    !!{
    Returns the specific angular momentum of accreted material in the ADAF.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: spinBlackHole, radius

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%angularMomentumRadiusPrevious .or. spinBlackHole /= self%angularMomentumSpinPrevious) then
       self%angularMomentumStored        =+self%enthalpyAngularMomentumProduct(spinBlackHole,radius) &
            &                             /self%enthalpy                      (spinBlackHole,radius)
       self%angularMomentumRadiusPrevious=+radius
       self%angularMomentumSpinPrevious  =+spinBlackHole
    end if
    adafAngularMomentum=self%angularMomentumStored
    return
  end function adafAngularMomentum

  double precision function adafEnthalpyAngularMomentumProduct(self,spinBlackHole,radius)
    !!{
    Return the product of enthalpy and angular momentum for the ADAF.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_ISCO_Radius, unitsGravitational
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius                       , spinBlackHole
    double precision                                    :: enthalpyAngularMomentum1     , enthalpyAngularMomentum2, &
         &                                                 enthalpyAngularMomentum3     , enthalpyAngularMomentum4, &
         &                                                 enthalpyAngularMomentum5     , enthalpyAngularMomentum6, &
         &                                                 viscosityParameterLogarithmic, radiusISCO

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%enthlpyAngMPrdctRadiusPrevious .or. spinBlackHole /= self%enthlpyAngMPrdctSpinPrevious) then
       ! We are not, so compute (and store) the value.
       viscosityParameterLogarithmic     =+log10(                                                   &
            &                                    +self%viscosityParameter(spinBlackHole)            &
            &                                   )
       radiusISCO                        =+Black_Hole_ISCO_Radius(spinBlackHole,unitsGravitational)
       enthalpyAngularMomentum1          =+0.0871d0                                                 &
            &                             *radiusISCO                                               &
            &                             -0.10282d0
       enthalpyAngularMomentum2          =+0.5000d0                                                 &
            &                             -7.7983d0                                                 &
            &                             *(                                                        &
            &                               +self%adiabaticIndex                                    &
            &                               -1.333d0                                                &
            &                              )**1.26d0
       enthalpyAngularMomentum3          =+0.153d0                                                  &
            &                             *(                                                        &
            &                               +radiusISCO                                             &
            &                               -0.6d0                                                  &
            &                              )**0.30d0                                                &
            &                             +0.105d0
       enthalpyAngularMomentum4          =+enthalpyAngularMomentum3                                 &
            &                             *(                                                        &
            &                               +0.9000d0                                               &
            &                               *self%adiabaticIndex                                    &
            &                               -0.2996d0                                               &
            &                              )                                                        &
            &                             *(                                                        &
            &                               +1.202d0                                                &
            &                               -0.080d0                                                &
            &                               *(                                                      &
            &                                 +viscosityParameterLogarithmic                        &
            &                                 +2.5d0                                                &
            &                                )**2.6d0                                               &
            &                              )
       enthalpyAngularMomentum5          =-1.800d0                                                  &
            &                             *self%adiabaticIndex                                      &
            &                             +4.299d0                                                  &
            &                             -0.018d0                                                  &
            &                             +0.018d0                                                  &
            &                             *(                                                        &
            &                               +viscosityParameterLogarithmic                          &
            &                               +2.0d0                                                  &
            &                              )**3.571d0
       enthalpyAngularMomentum6          =+enthalpyAngularMomentum4                                 &
            &                             *(                                                        &
            &                               +(                                                      &
            &                                 +(                                                    &
            &                                   +0.14d0                                             &
            &                                   *log10(radius)**enthalpyAngularMomentum5            &
            &                                   +0.23d0                                             &
            &                                  )                                                    &
            &                                 /enthalpyAngularMomentum4                             &
            &                                )**10.0d0                                              &
            &                               +1.0d0                                                  &
            &                              )**0.1d0
       self%enthlpyAngMPrdctStored       =+enthalpyAngularMomentum2                                 &
            &                             +(                                                        &
            &                               +enthalpyAngularMomentum1                               &
            &                               +10.0d0**enthalpyAngularMomentum6                       &
            &                              )                                                        &
            &                             *(                                                        &
            &                               +1.15d0                                                 &
            &                               -0.03d0                                                 &
            &                               *(                                                      &
            &                                 +3.0d0                                                &
            &                                 +viscosityParameterLogarithmic                        &
            &                                )**2.37d0                                              &
            &                              )
       self%enthlpyAngMPrdctRadiusPrevious=+radius
       self%enthlpyAngMPrdctSpinPrevious  =+spinBlackHole
    end if
    adafEnthalpyAngularMomentumProduct=self%enthlpyAngMPrdctStored
    return
  end function adafEnthalpyAngularMomentumProduct

  double precision function adafEnthalpy(self,spinBlackHole,radius)
    !!{
    Returns the relativistic enthalpy of the ADAF.
    !!}
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: spinBlackHole, radius

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%enthalpyRadiusPrevious .or. spinBlackHole /= self%enthalpySpinPrevious) then
       self%enthalpyStored        =+1.0d0                                  &
            &                      +(                                      &
            &                        +self%adiabaticIndex                  &
            &                        /(                                    &
            &                          +self%adiabaticIndex                &
            &                          -1.0d0                              &
            &                         )                                    &
            &                       )                                      &
            &                      *self%temperature(spinBlackHole,radius)
       self%enthalpyRadiusPrevious=+radius
       self%enthalpySpinPrevious  =+spinBlackHole
    end if
    adafEnthalpy=self%enthalpyStored
    return
  end function adafEnthalpy

  double precision function adafTemperature(self,spinBlackHole,radius)
    !!{
    Return the dimensionless temperature of the ADAF
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_ISCO_Radius, unitsGravitational
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius                       , spinBlackHole
    double precision                                    :: viscosityParameterLogarithmic, radiusISCO        , &
         &                                                 temperatureFactor1           , temperatureFactor2, &
         &                                                 temperatureFactor3           , temperatureFactor4, &
         &                                                 temperatureFactor5

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= self%temperatureRadiusPrevious .or. spinBlackHole /= self%temperatureSpinPrevious) then
       viscosityParameterLogarithmic =+log10(                                        &
            &                                +self%viscosityParameter(spinBlackHole) &
            &                               )
       radiusISCO                    =+Black_Hole_ISCO_Radius(spinBlackHole,unitsGravitational)
       temperatureFactor1            =-0.270278d0                                               &
            &                         *self%adiabaticIndex                                      &
            &                         +1.360270d0
       temperatureFactor2            =-0.9400d0                                                 &
            &                         +4.4744d0                                                 &
            &                         *(                                                        &
            &                           +self%adiabaticIndex                                    &
            &                           -1.444d0                                                &
            &                          )                                                        &
            &                         -5.1402d0                                                 &
            &                         *(                                                        &
            &                           +self%adiabaticIndex                                    &
            &                           -1.444d0                                                &
            &                          )**2
       temperatureFactor3            =+0.840d0                                                  &
            &                         *viscosityParameterLogarithmic                            &
            &                         +0.919d0                                                  &
            &                         -0.643d0                                                  &
            &                         *exp(                                                     &
            &                              -0.209d0                                             &
            &                              /self%viscosityParameter(spinBlackHole)              &
            &                             )
       temperatureFactor4            =+(                                                        &
            &                           + 0.6365d0                                              &
            &                           *radiusISCO                                             &
            &                           - 0.4828d0                                              &
            &                          )                                                        &
            &                         *(                                                        &
            &                           + 1.0000d0                                              &
            &                           +11.9000d0                                              &
            &                           *exp(                                                   &
            &                                -0.838d0                                           &
            &                                *radiusISCO**4                                     &
            &                               )                                                   &
            &                          )
       temperatureFactor5            =+1.444d0                                                  &
            &                         *exp(                                                     &
            &                              -1.01d0                                              &
            &                              *radiusISCO**0.86d0                                  &
            &                             )                                                     &
            &                         +0.1d0
       self%temperatureStored        =+0.31d0                                                   &
            &                         *(                                                        &
            &                           +(                                                      &
            &                             +1.0d0                                                &
            &                             +(                                                    &
            &                               +temperatureFactor4                                 &
            &                               /radius                                             &
            &                              )**0.9d0                                             &
            &                            )**(temperatureFactor2+temperatureFactor3)             &
            &                          )                                                        &
            &                         /(                                                        &
            &                           +radius                                                 &
            &                           -temperatureFactor5                                     &
            &                          )**temperatureFactor1
       self%temperatureRadiusPrevious=+radius
       self%temperatureSpinPrevious  =+spinBlackHole
    end if
    adafTemperature=self%temperatureStored
    return
  end function adafTemperature

  double precision function adafVelocity(self,spinBlackHole,radius)
    !!{
    Return the (dimensionless) velocity in an ADAF at given {\normalfont \ttfamily radius}, for a black hole of given
    {\normalfont \ttfamily spinBlackHole}.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Horizon_Radius, Black_Hole_ISCO_Radius
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius          , spinBlackHole
    double precision                                    :: Phi             , viscosityParameterEffective, &
         &                                                 radiusISCO      , radiusEffective            , &
         &                                                 radiusHorizon   , PhiFactor1                 , &
         &                                                 PhiFactor2      , PhiFactor3                 , &
         &                                                 PhiFactor4      , PhiFactor5                 , &
         &                                                 radiusFractional, radiusHorizonFractional

    ! Check if we are being called with the same arguments as the previous call.
    if     (                                              &
         &   radius        /= self%velocityRadiusPrevious &
         &  .or.                                          &
         &   spinBlackHole /= self%velocitySpinPrevious   &
         & ) then
       radiusHorizon              =+Black_Hole_Horizon_Radius(spinBlackHole)
       radiusISCO                 =+Black_Hole_ISCO_Radius   (spinBlackHole)
       radiusFractional           =+radius                                        &
            &                      /radiusISCO
       radiusHorizonFractional    =+radiusHorizon                                 &
            &                      /radiusISCO
       viscosityParameterEffective=+self%viscosityParameter(spinBlackHole)        &
            &                      *(                                             &
            &                        +1.00d0                                      &
            &                        +6.450d0*(+self%adiabaticIndex-1.444d0)      &
            &                        +1.355d0*(+self%adiabaticIndex-1.444d0)**2   &
            &                       )
       PhiFactor1                 =+9.0d0                                         &
            &                      *log(                                          &
            &                           +9.0d0                                    &
            &                           *radiusFractional                         &
            &                          )
       PhiFactor2                 =+exp(                                          &
            &                           -0.66d0                                   &
            &                           *(                                        &
            &                             +1.0d0                                  &
            &                             -2.0d0                                  &
            &                             *viscosityParameterEffective            &
            &                            )                                        &
            &                           *log(                                     &
            &                                +viscosityParameterEffective         &
            &                                /0.1d0                               &
            &                               )                                     &
            &                           *log(                                     &
            &                                +radiusFractional                    &
            &                                /radiusHorizonFractional             &
            &                               )                                     &
            &                          )
       PhiFactor3                 =+1.0d0                                         &
            &                      -exp(                                          &
            &                           -radiusFractional                         &
            &                           *(                                        &
            &                             +0.16d0                                 &
            &                             *(                                      &
            &                               +spinBlackHole                        &
            &                               -1.0d0                                &
            &                              )                                      &
            &                             +0.76d0                                 &
            &                            )                                        &
            &                          )
       PhiFactor4                 =+1.40000d0                                     &
            &                      +0.29065d0*(+spinBlackHole-0.5d0)**4           &
            &                      -0.87560d0*(+spinBlackHole-0.5d0)**2           &
            &                      +(                                             &
            &                        -0.33000d0                                   &
            &                        *spinBlackHole                               &
            &                        +0.45035d0                                   &
            &                       )                                             &
            &                      *(                                             &
            &                        +1.0d0                                       &
            &                        -exp(                                        &
            &                             -(                                      &
            &                               +radiusFractional                     &
            &                               -radiusHorizonFractional              &
            &                              )                                      &
            &                            )                                        &
            &                       )
       PhiFactor5                 =+2.3d0                                         &
            &                      *exp(                                          &
            &                           +40.0d0                                   &
            &                           *(                                        &
            &                             +spinBlackHole                          &
            &                             -1.0d0                                  &
            &                            )                                        &
            &                          )                                          &
            &                      *exp(                                          &
            &                           -15.0d0                                   &
            &                           *radiusISCO                               &
            &                           *(                                        &
            &                             +radiusFractional                       &
            &                             -radiusHorizonFractional                &
            &                            )                                        &
            &                          )                                          &
            &                      +1.0d0
       Phi                        =+PhiFactor1                                    &
            &                      *PhiFactor2                                    &
            &                      *PhiFactor3                                    &
            &                      *PhiFactor4                                    &
            &                      *PhiFactor5
       radiusEffective            =+radiusHorizon                                 &
            &                      +Phi                                           &
            &                      *(                                             &
            &                        +radius                                      &
            &                        -radiusHorizon                               &
            &                       )
       self%velocityStored        =+sqrt(                                         &
            &                            +1.0d0                                   &
            &                            -(                                       &
            &                              +1.0d0                                 &
            &                              -2.0d0                                 &
            &                              /radiusEffective                       &
            &                              +(                                     &
            &                                +spinBlackHole                       &
            &                                /radiusEffective                     &
            &                               )**2                                  &
            &                             )                                       &
            &                           )
       self%velocityRadiusPrevious=+radius
       self%velocitySpinPrevious  =+spinBlackHole
    end if
    adafVelocity=self%velocityStored
    return
  end function adafVelocity

  double precision function adafHeight(self,spinBlackHole,radius)
    !!{
    Return the (dimensionless) height in an ADAF at given {\normalfont \ttfamily radius}, for a black hole of given {\normalfont \ttfamily spinBlackHole}.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Frame_Dragging_Frequency, Black_Hole_Metric_A_Factor, Black_Hole_Metric_D_Factor
    implicit none
    class           (accretionDisksADAF), intent(inout) :: self
    double precision                    , intent(in   ) :: radius                           , spinBlackHole
    double precision                                    :: verticalFrequencyEffectiveSquared

    ! Check if we are being called with the same arguments as the previous call.
    if     (                                            &
         &   radius        /= self%heightRadiusPrevious &
         &  .or.                                        &
         &   spinBlackHole /= self%heightSpinPrevious   &
         & ) then
       verticalFrequencyEffectiveSquared=+(                                                                  &
            &                              +                                        spinBlackHole        **2 &
            &                              +(                                                                &
            &                                +1.0d0                                                          &
            &                                -(                                                              &
            &                                  +                                    spinBlackHole            &
            &                                  *Black_Hole_Frame_Dragging_Frequency(spinBlackHole,radius)    &
            &                                 )                                                          **2 &
            &                               )                                                                &
            &                              *self%angularMomentum                   (spinBlackHole,radius)**2 &
            &                              -(                                                                &
            &                                +(                                                              &
            &                                  +                                    spinBlackHole            &
            &                                  *self%gammaAzimuthal                (spinBlackHole,radius)    &
            &                                 )                                                          **2 &
            &                                /Black_Hole_Metric_A_Factor           (spinBlackHole,radius)    &
            &                               )                                                                &
            &                              *self%gammaRadial                       (spinBlackHole,radius)**2 &
            &                              *Black_Hole_Metric_D_Factor             (spinBlackHole,radius)    &
            &                              -self%gammaRadial                       (spinBlackHole,radius)    &
            &                              *sqrt(                                                            &
            &                                    +Black_Hole_Metric_D_Factor       (spinBlackHole,radius)    &
            &                                    /Black_Hole_Metric_A_Factor       (spinBlackHole,radius)    &
            &                                   )                                                            &
            &                              *2.0d0                                                            &
            &                              *self%angularMomentum                   (spinBlackHole,radius)    &
            &                              *Black_Hole_Frame_Dragging_Frequency    (spinBlackHole,radius)    &
            &                              *self%gammaAzimuthal                    (spinBlackHole,radius)    &
            &                              *                                        spinBlackHole        **2 &
            &                             )                                                                  &
            &                            /                                                        radius **4
       self%heightStored                =+sqrt(                                                              &
            &                                  +self%temperature                   (spinBlackHole,radius)    &
            &                                  /self%enthalpy                      (spinBlackHole,radius)    &
            &                                  /                                                  radius **2 &
            &                                  /verticalFrequencyEffectiveSquared&
            &                                 )
       self%heightRadiusPrevious        =+radius
       self%heightSpinPrevious          =+spinBlackHole
    end if
    adafHeight=self%heightStored
    return
  end function adafHeight

  double precision function adafAdiabaticIndexDefault(fieldEnhancementOption)
    !!{
    Returns the default adiabatic index in an ADAF give the field enhancement option.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(enumerationADAFFieldEnhancementType), intent(in   ) :: fieldEnhancementOption

    select case (fieldEnhancementOption%ID)
    case (adafFieldEnhancementExponential%ID)
       adafAdiabaticIndexDefault=1.444d0
    case (adafFieldEnhancementLinear     %ID)
       adafAdiabaticIndexDefault=1.333d0
    case default
       adafAdiabaticIndexDefault=0.000d0
      call Error_Report('unknown field enhancement option'//{introspection:location})
    end select
    return
  end function adafAdiabaticIndexDefault

