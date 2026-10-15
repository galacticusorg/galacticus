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
  Implements a spectrum of thermal emission from the dusty torus around an AGN.
  !!}

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumAGNTorus" docformat="rst">
   <description>
   Thermal emission from the dusty torus around an AGN, approximated as a superposition of modified blackbodies, of
   emissivity index :math:`\beta` (``[exponentEmissivity]``), spanning temperatures from ``[temperatureMinimum]`` up to
   the sublimation temperature of the dust, ``[temperatureSublimation]``. The luminosity emitted per unit logarithmic
   interval of temperature is taken to be a power law, :math:`\mathrm{d}L/\mathrm{d}\ln T \propto T^\alpha`, with
   :math:`\alpha` given by ``[exponentTemperature]``, so that

   .. math::

      L_\nu = L_\mathrm{abs} \int_{T_\mathrm{min}}^{T_\mathrm{sub}} p(\ln T) \frac{h}{kT} \frac{x^{3+\beta}}{(\mathrm{e}^x-1)\Gamma(4+\beta)\zeta(4+\beta)}\,\mathrm{d}\ln T, \qquad x = \frac{h\nu}{kT},

   where :math:`p(\ln T) \propto T^\alpha` is normalized to unit integral. Each modified blackbody integrates to one
   over frequency, so :math:`L_\nu` integrates to exactly the luminosity absorbed, :math:`L_\mathrm{abs}`. The integral over
   temperature is evaluated by a Gauss-Legendre rule in :math:`\ln T` whose weights are normalized to sum to one, which
   keeps that property exact.

   The temperatures of the dust in a torus are set by its distance from the accretion disk, not by energy balance with
   the luminosity it absorbs, so the shape of this spectrum is fixed, and the mass of dust is not used. This is
   appropriate for the torus phase of :galacticus-class:`dustAttenuationAGNTorus`, whose dust is not part of the dust
   mass of the galaxy: give it a zero ``fractionsMassPhase`` in :galacticus-class:`nodePropertyExtractorSEDDustEmission`.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumAGNTorus
     !!{RST
     A spectrum of thermal emission from the dusty torus around an AGN.
     !!}
     private
     double precision                            :: temperatureSublimation, temperatureMinimum, &
          &                                         exponentTemperature   , exponentEmissivity, &
          &                                         normalizationShape
     ! Temperatures, and the fraction of the luminosity emitted at each, of the quadrature over temperature.
     double precision, allocatable, dimension(:) :: temperatures          , weights
   contains
     procedure :: luminosity => agnTorusLuminosity
  end type dustEmissionSpectrumAGNTorus

  interface dustEmissionSpectrumAGNTorus
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumAGNTorus` dust emission spectrum class.
     !!}
     module procedure agnTorusConstructorParameters
     module procedure agnTorusConstructorInternal
  end interface dustEmissionSpectrumAGNTorus

contains

  function agnTorusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorus` dust emission spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustEmissionSpectrumAGNTorus)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: temperatureSublimation, temperatureMinimum, &
         &                                                           exponentTemperature   , exponentEmissivity

    !![
    <inputParameter docformat="rst">
      <name>temperatureSublimation</name>
      <defaultValue>1500.0d0</defaultValue>
      <description>
      The sublimation temperature of the dust, in K: the temperature of the hottest dust in the torus.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>temperatureMinimum</name>
      <defaultValue>100.0d0</defaultValue>
      <description>
      The temperature of the coolest dust in the torus, in K.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentTemperature</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The exponent, :math:`\alpha`, of the distribution of emitted luminosity over temperature,
      :math:`\mathrm{d}L/\mathrm{d}\ln T \propto T^\alpha`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentEmissivity</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The emissivity index, :math:`\beta`, of the dust.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustEmissionSpectrumAGNTorus(temperatureSublimation,temperatureMinimum,exponentTemperature,exponentEmissivity)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function agnTorusConstructorParameters

  function agnTorusConstructorInternal(temperatureSublimation,temperatureMinimum,exponentTemperature,exponentEmissivity) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorus` dust emission spectrum class. The
    quadrature over temperature is built here, once.
    !!}
    use :: Dust_Attenuations, only : gaussLegendreRule
    use :: Error            , only : Error_Report
    use :: Gamma_Functions  , only : Gamma_Function
    use :: Zeta_Functions   , only : Zeta_Function
    implicit none
    type            (dustEmissionSpectrumAGNTorus)                            :: self
    double precision                              , intent(in   )             :: temperatureSublimation  , temperatureMinimum, &
         &                                                                       exponentTemperature     , exponentEmissivity
    ! Order of the quadrature over the logarithm of temperature.
    integer                                       , parameter                 :: orderQuadrature       =32
    double precision                              , allocatable, dimension(:) :: abscissae
    !![
    <constructorAssign variables="temperatureSublimation, temperatureMinimum, exponentTemperature, exponentEmissivity"/>
    !!]

    if (temperatureMinimum <= 0.0d0 .or. temperatureSublimation <= temperatureMinimum)                                      &
         & call Error_Report('require 0 < `temperatureMinimum` < `temperatureSublimation`'//{introspection:location})
    call gaussLegendreRule(orderQuadrature,abscissae,self%weights)
    self%temperatures      =exp(log(temperatureMinimum)+log(temperatureSublimation/temperatureMinimum)*abscissae)
    self%weights           =self%weights*self%temperatures**exponentTemperature
    self%weights           =self%weights/sum(self%weights)
    self%normalizationShape=1.0d0/Gamma_Function(4.0d0+exponentEmissivity)/Zeta_Function(4.0d0+exponentEmissivity)
    return
  end function agnTorusConstructorInternal

  function agnTorusLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the torus, normalized analytically to the luminosity it absorbs.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant, plancksConstant, speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    implicit none
    class           (dustEmissionSpectrumAGNTorus), intent(inout)                               :: self
    double precision                              , intent(in   ), dimension(:                ) :: wavelengths
    double precision                              , intent(in   )                               :: luminosityAbsorbed        , massDust   , &
         &                                                                                         time
    double precision                                             , dimension(size(wavelengths)) :: luminosity
    ! Beyond this argument the exponential overflows, and the emission is in any case negligible.
    double precision                              , parameter                                   :: argumentMaximum   =700.0d0
    ! Below this argument the denominator is evaluated by its series, avoiding cancellation in exp(x)-1.
    double precision                              , parameter                                   :: argumentSmall     =1.0d-4
    double precision                                                                            :: argument                  , denominator
    integer                                                                                     :: i                         , j
    !$GLC attributes unused :: massDust, time

    luminosity=0.0d0
    if (luminosityAbsorbed <= 0.0d0) return
    do i=1,size(wavelengths)
       do j=1,size(self%temperatures)
          argument=plancksConstant*speedLight*metersToAngstroms/wavelengths(i)/boltzmannsConstant/self%temperatures(j)
          if (argument > argumentMaximum) cycle
          if (argument < argumentSmall) then
             denominator=argument*(1.0d0+argument/2.0d0+argument**2/6.0d0)
          else
             denominator=exp(argument)-1.0d0
          end if
          luminosity(i)=+luminosity(i)                                           &
               &        +self%weights(j)                                         &
               &        *plancksConstant/boltzmannsConstant/self%temperatures(j) &
               &        *argument**(3.0d0+self%exponentEmissivity)               &
               &        /denominator
       end do
    end do
    luminosity=luminosity*luminosityAbsorbed*self%normalizationShape
    return
  end function agnTorusLuminosity
