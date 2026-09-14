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
  Implements a dust screen whose optical depth scales with the surface density of metals.
  !!}

  !![
  <dustAttenuation name="dustAttenuationScreenSurfaceDensityMetals" docformat="rst">
   <description>
   A uniform dust screen whose :math:`V`-band optical depth is proportional to the surface density of dust in the
   component being attenuated,

   .. math::

      \tau_\mathrm{V} = C \, \kappa_\mathrm{V} \, f_\mathrm{dust:metals} \, \Sigma_\mathrm{Z},

   where :math:`\Sigma_\mathrm{Z}` is the surface density of metals in the gas, and the :math:`V`-band extinction opacity
   per unit mass of dust, :math:`\kappa_\mathrm{V}`, and the dust-to-metals ratio, :math:`f_\mathrm{dust:metals}`, are
   supplied by a :galacticus-class:`dustPropertiesClass` object. Taking both from that object, rather than fixing them
   here, keeps the dust which attenuates a galaxy's light consistent with the dust which re-emits it. With the default
   :galacticus-class:`dustPropertiesSimple` the product :math:`\kappa_\mathrm{V} f_\mathrm{dust:metals}` reproduces the
   Milky Way relation between column density and reddening of :cite:t:`savage_observed_1979`.

   The dimensionless coefficient :math:`C` (``coefficient``) scales the optical depth *without* changing the mass of
   dust. It therefore stands for the effects of geometry---clumping of the dust, say, or a screen which covers only part
   of the emission---and should not be used to change how much dust a galaxy has: set the dust-to-metals ratio of the
   :galacticus-class:`dustPropertiesClass` object for that.

   The surface density is that of an exponential disk or of a spheroid of the same scale radius,
   :math:`\Sigma_\mathrm{Z} = Z M_\mathrm{gas} / 2\pi r^2`, and is taken to be zero for a component with no gas or no
   size.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationScreen) :: dustAttenuationScreenSurfaceDensityMetals
     !!{RST
     A dust screen whose optical depth scales with the surface density of metals.
     !!}
     private
     class           (dustPropertiesClass), pointer :: dustProperties_ => null()
     double precision                               :: coefficient
   contains
     final     ::                  screenSurfaceDensityMetalsDestructor
     procedure :: depthOpticalV => screenSurfaceDensityMetalsDepthOpticalV
  end type dustAttenuationScreenSurfaceDensityMetals

  interface dustAttenuationScreenSurfaceDensityMetals
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` dust attenuation class.
     !!}
     module procedure screenSurfaceDensityMetalsConstructorParameters
     module procedure screenSurfaceDensityMetalsConstructorInternal
  end interface dustAttenuationScreenSurfaceDensityMetals

contains

  function screenSurfaceDensityMetalsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` dust attenuation class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationScreenSurfaceDensityMetals)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (dustExtinctionCurveClass                 ), pointer       :: dustExtinctionCurve_
    class           (dustPropertiesClass                      ), pointer       :: dustProperties_
    double precision                                                           :: coefficient

    !![
    <inputParameter docformat="rst">
      <name>coefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A dimensionless multiplicative coefficient applied to the :math:`V`-band optical depth, representing the effects
      of geometry. It does not change the mass of dust, which is set by the ``dustProperties`` object.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
    <objectBuilder class="dustProperties"      name="dustProperties_"      source="parameters"/>
    !!]
    self=dustAttenuationScreenSurfaceDensityMetals(coefficient,dustExtinctionCurve_,dustProperties_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"/>
    <objectDestructor name="dustProperties_"     />
    !!]
    return
  end function screenSurfaceDensityMetalsConstructorParameters

  function screenSurfaceDensityMetalsConstructorInternal(coefficient,dustExtinctionCurve_,dustProperties_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` dust attenuation class.
    !!}
    implicit none
    type            (dustAttenuationScreenSurfaceDensityMetals)                        :: self
    double precision                                           , intent(in   )         :: coefficient
    class           (dustExtinctionCurveClass                 ), intent(in   ), target :: dustExtinctionCurve_
    class           (dustPropertiesClass                      ), intent(in   ), target :: dustProperties_
    !![
    <constructorAssign variables="coefficient, *dustExtinctionCurve_, *dustProperties_"/>
    !!]

    return
  end function screenSurfaceDensityMetalsConstructorInternal

  subroutine screenSurfaceDensityMetalsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationScreenSurfaceDensityMetals), intent(inout) :: self

    !![
    <objectDestructor name="self%dustProperties_"/>
    !!]
    return
  end subroutine screenSurfaceDensityMetalsDestructor

  double precision function screenSurfaceDensityMetalsDepthOpticalV(self,node,componentType) result(depthOpticalV)
    !!{RST
    Return the :math:`V`-band optical depth of a screen scaling with the surface density of metals.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : hecto    , kilo
    implicit none
    class           (dustAttenuationScreenSurfaceDensityMetals), intent(inout)         :: self
    type            (treeNode                                 ), intent(inout), target :: node
    type            (enumerationComponentTypeType             ), intent(in   )         :: componentType
    double precision                                                                   :: massGas             , radius, &
         &                                                                                metallicity         ,         &
         &                                                                                densitySurfaceMetals

    call componentGasProperties(node,componentType,massGas,radius,metallicity)
    ! A component with no gas, or no size, has no dust.
    if (massGas <= 0.0d0 .or. radius <= 0.0d0) then
       depthOpticalV=0.0d0
       return
    end if
    ! Surface density of metals, in g/cm².
    densitySurfaceMetals=+metallicity  &
         &               *massGas      &
         &               *massSolar    &
         &               *kilo         &
         &               /2.0d0        &
         &               /Pi           &
         &               /(            &
         &                 +radius     &
         &                 *megaParsec &
         &                 *hecto      &
         &                )**2
    depthOpticalV       =+self                %coefficient                            &
         &               *self%dustProperties_%opacityExtinctionV(                  ) &
         &               *self%dustProperties_%dustToMetalsRatio (node,componentType) &
         &               *densitySurfaceMetals
    return
  end function screenSurfaceDensityMetalsDepthOpticalV
