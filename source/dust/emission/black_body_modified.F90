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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a modified blackbody spectrum of thermal emission from dust.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Dust_Properties    , only : dustPropertiesClass

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumBlackBodyModified" docformat="rst">
   <description>
   Thermal emission from optically thin dust at a single temperature---a modified blackbody, or greybody,

   .. math::

      L_\nu = 4\pi M_\mathrm{dust} \kappa_\mathrm{abs}(\nu) B_\nu(T), \qquad \kappa_\mathrm{abs}(\nu) = \kappa_\mathrm{ref} \left(\frac{\nu}{\nu_\mathrm{ref}}\right)^\beta,

   where the absorption opacity of the dust---its value :math:`\kappa_\mathrm{ref}` at a reference frequency
   :math:`\nu_\mathrm{ref}`, and emissivity index :math:`\beta`---is taken from a
   :galacticus-class:`dustPropertiesClass` object, the same one which sets the dust content for attenuation.

   By default the temperature follows from energy balance. Dust re-emits the luminosity :math:`L_\mathrm{abs}` it
   absorbs, and integrating the spectrum over frequency gives

   .. math::

      L_\mathrm{abs} = 8\pi M_\mathrm{dust} \frac{\kappa_\mathrm{ref}}{\nu_\mathrm{ref}^\beta} \frac{h}{c^2} \left(\frac{k T_\star}{h}\right)^{4+\beta} \Gamma(4+\beta)\, \zeta(4+\beta),

   which is solved for the temperature :math:`T_\star` in closed form. Setting ``temperature`` instead fixes
   :math:`T_\star`, and the mass of dust is then not used.

   The cosmic microwave background also heats dust, which matters at high redshift. With ``heatingCMB`` true, the
   default, the dust is taken to be in equilibrium with both the radiation it absorbs and the CMB, following
   :cite:t:`da_cunha_effect_2013`,

   .. math::

      T^{4+\beta} = T_\star^{4+\beta} + T_\mathrm{CMB}(z)^{4+\beta},

   and so emits :math:`L_\mathrm{abs} (T/T_\star)^{4+\beta}`: the absorbed luminosity, plus the energy absorbed from the
   CMB and re-radiated. The spectrum returned is that total emission. An observer measures it against the CMB, which
   reduces its contrast; that reduction is not applied here. For dust at typical temperatures the CMB contributes
   negligibly below a redshift of a few.

   The spectrum is normalized analytically, so that it integrates to exactly the emitted luminosity whatever wavelengths
   it is evaluated at. The dust is assumed optically thin to its own emission, which fails shortward of about
   :math:`100\,\mu\hbox{m}` in dense, compact systems. A single temperature describes dust in equilibrium with the
   radiation field, not the stochastically heated small grains which dominate the mid-infrared; several of these spectra
   may be combined through :galacticus-class:`dustEmissionSpectrumSum` to represent a range of temperatures.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumBlackBodyModified
     !!{RST
     A modified blackbody spectrum of thermal emission from dust.
     !!}
     private
     class           (dustPropertiesClass    ), pointer :: dustProperties_     => null()
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: temperature_
     logical                                            :: heatingCMB
   contains
     !![
     <methods docformat="rst">
       <method method="temperature" description="Return the temperature of the dust, and optionally the temperature it would have if heated by absorbed radiation alone."/>
     </methods>
     !!]
     final     ::                blackBodyModifiedDestructor
     procedure :: luminosity  => blackBodyModifiedLuminosity
     procedure :: temperature => blackBodyModifiedTemperature
  end type dustEmissionSpectrumBlackBodyModified

  interface dustEmissionSpectrumBlackBodyModified
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumBlackBodyModified` dust emission spectrum class.
     !!}
     module procedure blackBodyModifiedConstructorParameters
     module procedure blackBodyModifiedConstructorInternal
  end interface dustEmissionSpectrumBlackBodyModified

contains

  function blackBodyModifiedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumBlackBodyModified` dust emission spectrum class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustEmissionSpectrumBlackBodyModified)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (dustPropertiesClass                  ), pointer       :: dustProperties_
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    double precision                                                       :: temperature_
    logical                                                                :: heatingCMB

    !![
    <inputParameter docformat="rst">
      <name>temperature</name>
      <variable>temperature_</variable>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      If positive, the temperature, in K, of the dust when heated by the radiation it absorbs, overriding the
      temperature which energy balance would give. The mass of dust is then not used. Heating by the cosmic microwave
      background, if enabled, is applied on top of this.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>heatingCMB</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, the dust is also heated by the cosmic microwave background, following :cite:t:`da_cunha_effect_2013`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustProperties"     name="dustProperties_"     source="parameters"/>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=dustEmissionSpectrumBlackBodyModified(temperature_,heatingCMB,dustProperties_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustProperties_"    />
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function blackBodyModifiedConstructorParameters

  function blackBodyModifiedConstructorInternal(temperature_,heatingCMB,dustProperties_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumBlackBodyModified` dust emission spectrum class.
    !!}
    implicit none
    type            (dustEmissionSpectrumBlackBodyModified)                        :: self
    double precision                                       , intent(in   )         :: temperature_
    logical                                                , intent(in   )         :: heatingCMB
    class           (dustPropertiesClass                  ), intent(in   ), target :: dustProperties_
    class           (cosmologyFunctionsClass              ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="temperature_, heatingCMB, *dustProperties_, *cosmologyFunctions_"/>
    !!]

    return
  end function blackBodyModifiedConstructorInternal

  subroutine blackBodyModifiedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustEmissionSpectrumBlackBodyModified` dust emission spectrum class.
    !!}
    implicit none
    type(dustEmissionSpectrumBlackBodyModified), intent(inout) :: self

    !![
    <objectDestructor name="self%dustProperties_"    />
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine blackBodyModifiedDestructor

  double precision function blackBodyModifiedTemperature(self,luminosityAbsorbed,massDust,time,temperatureAbsorbed) result(temperature)
    !!{RST
    Return the temperature, in K, of dust of mass ``massDust`` (in :math:`M_\odot`) absorbing a luminosity
    ``luminosityAbsorbed`` (in :math:`L_\odot`) at cosmic ``time`` (in Gyr), including heating by the cosmic microwave
    background if enabled. If present, ``temperatureAbsorbed`` is set to the temperature the dust would have if heated
    by the absorbed luminosity alone.
    !!}
    use :: Error                           , only : Error_Report
    use :: Gamma_Functions                 , only : Gamma_Function
    use :: Numerical_Constants_Astronomical, only : luminositySolar   , massSolar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant, plancksConstant  , speedLight
    use :: Numerical_Constants_Prefixes    , only : centi             , kilo
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    use :: Zeta_Functions                  , only : Zeta_Function
    implicit none
    class           (dustEmissionSpectrumBlackBodyModified), intent(inout)           :: self
    double precision                                       , intent(in   )           :: luminosityAbsorbed , massDust           , &
         &                                                                              time
    double precision                                       , intent(  out), optional :: temperatureAbsorbed
    double precision                                                                 :: opacityReference   , wavelengthReference, &
         &                                                                              exponent           , frequencyReference , &
         &                                                                              temperatureStar    , temperatureCMB

    call self%dustProperties_%opacityAbsorptionPowerLaw(opacityReference,wavelengthReference,exponent)
    if (self%temperature_ > 0.0d0) then
       temperatureStar=self%temperature_
    else if (luminosityAbsorbed <= 0.0d0) then
       temperatureStar=0.0d0
    else
       if (massDust <= 0.0d0 .or. opacityReference <= 0.0d0)                                                               &
            & call Error_Report('dust absorbs luminosity but has no mass or no opacity, so no temperature can be found'//  &
            &                   {introspection:location}                                                                 )
       frequencyReference=+speedLight*metersToAngstroms/wavelengthReference
       ! Solve the energy balance for the temperature, working in SI units. The opacity is converted from cm² g⁻¹ to
       ! m² kg⁻¹.
       temperatureStar   =+plancksConstant                                     &
            &             /boltzmannsConstant                                  &
            &             *(                                                   &
            &               +luminosityAbsorbed                                &
            &               *luminositySolar                                   &
            &               *speedLight        **2                             &
            &               *frequencyReference**exponent                      &
            &               /8.0d0                                             &
            &               /Pi                                                &
            &               /plancksConstant                                   &
            &               /massDust                                          &
            &               /massSolar                                         &
            &               /opacityReference                                  &
            &               /centi**2                                          &
            &               /kilo                                              &
            &               /Gamma_Function(4.0d0+exponent)                    &
            &               /Zeta_Function (4.0d0+exponent)                    &
            &              )**(1.0d0/(4.0d0+exponent))
    end if
    if (present(temperatureAbsorbed)) temperatureAbsorbed=temperatureStar
    temperature=temperatureStar
    if (self%heatingCMB) then
       temperatureCMB=self%cosmologyFunctions_%temperatureCMBEpochal(time=time)
       temperature   =(temperatureStar**(4.0d0+exponent)+temperatureCMB**(4.0d0+exponent))**(1.0d0/(4.0d0+exponent))
    end if
    return
  end function blackBodyModifiedTemperature

  function blackBodyModifiedLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the dust, normalized analytically to the total it emits.

    Writing :math:`x = h\nu/kT`, the normalized spectrum is
    :math:`L_\nu = L_\mathrm{emit} (h/kT) x^{3+\beta} / [(\mathrm{e}^x-1) \Gamma(4+\beta) \zeta(4+\beta)]`, which
    integrates to :math:`L_\mathrm{emit}` over all frequencies.
    !!}
    use :: Gamma_Functions             , only : Gamma_Function
    use :: Numerical_Constants_Physical, only : boltzmannsConstant, plancksConstant, speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    use :: Zeta_Functions              , only : Zeta_Function
    implicit none
    class           (dustEmissionSpectrumBlackBodyModified), intent(inout)                               :: self
    double precision                                       , intent(in   ), dimension(:                ) :: wavelengths
    double precision                                       , intent(in   )                               :: luminosityAbsorbed         , massDust          , &
         &                                                                                                  time
    double precision                                                      , dimension(size(wavelengths)) :: luminosity
    ! Beyond this argument the exponential overflows, and the emission is in any case negligible.
    double precision                                       , parameter                                   :: argumentMaximum     =700.0d0
    ! Below this argument the denominator is evaluated by its series, avoiding cancellation in exp(x)-1.
    double precision                                       , parameter                                   :: argumentSmall       =1.0d-4
    double precision                                                                                     :: opacityReference           , wavelengthReference, &
         &                                                                                                  exponent                   , temperature        , &
         &                                                                                                  temperatureAbsorbed        , luminosityEmitted  , &
         &                                                                                                  normalization              , argument           , &
         &                                                                                                  denominator
    integer                                                                                              :: i

    luminosity=0.0d0
    if (luminosityAbsorbed <= 0.0d0) return
    call self%dustProperties_%opacityAbsorptionPowerLaw(opacityReference,wavelengthReference,exponent)
    temperature      =self%temperature(luminosityAbsorbed,massDust,time,temperatureAbsorbed)
    ! Dust heated by the CMB as well re-radiates that energy too.
    luminosityEmitted=luminosityAbsorbed*(temperature/temperatureAbsorbed)**(4.0d0+exponent)
    normalization    =+luminosityEmitted              &
         &            *plancksConstant                &
         &            /boltzmannsConstant             &
         &            /temperature                    &
         &            /Gamma_Function(4.0d0+exponent) &
         &            /Zeta_Function (4.0d0+exponent)
    do i=1,size(wavelengths)
       argument=plancksConstant*speedLight*metersToAngstroms/wavelengths(i)/boltzmannsConstant/temperature
       if (argument > argumentMaximum) cycle
       if (argument < argumentSmall) then
          denominator=argument*(1.0d0+argument/2.0d0+argument**2/6.0d0)
       else
          denominator=exp(argument)-1.0d0
       end if
       luminosity(i)=normalization*argument**(3.0d0+exponent)/denominator
    end do
    return
  end function blackBodyModifiedLuminosity
