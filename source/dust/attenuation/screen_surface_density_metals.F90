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
   A uniform dust screen whose :math:`V`-band optical depth is proportional to the surface density of metals in the
   gas of the component being attenuated,

   .. math::

      \tau_\mathrm{V} = C \, \frac{A_\mathrm{V}/E(B-V)}{N_\mathrm{H}/E(B-V)} \, \frac{X_\mathrm{H}}{m_\mathrm{u}} \, \frac{\Sigma_\mathrm{Z}}{Z_\odot} \, \frac{1}{2.5 \log_{10} \mathrm{e}},

   with the ratios :math:`A_\mathrm{V}/E(B-V)=3.1` and :math:`N_\mathrm{H}/E(B-V)=5.8\times10^{21}\,\hbox{atoms cm}^{-2}\,\hbox{mag}^{-1}` from
   :cite:t:`savage_observed_1979`, and :math:`\Sigma_\mathrm{Z}` the surface density of gas-phase metals. The
   dimensionless coefficient :math:`C` (``coefficient``) allows the overall normalization to be adjusted.

   The scaling assumes a universal dust-to-metals ratio, so that a galaxy with the metallicity and gas surface density
   of the local interstellar medium reproduces the Milky Way relation between column density and reddening.

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
     double precision :: coefficient
   contains
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
    double precision                                                           :: coefficient

    !![
    <inputParameter docformat="rst">
      <name>coefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A dimensionless multiplicative coefficient applied to the :math:`V`-band optical depth, allowing the overall
      normalization of the dust content to be adjusted away from the Milky Way calibration.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
    !!]
    self=dustAttenuationScreenSurfaceDensityMetals(coefficient,dustExtinctionCurve_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"/>
    !!]
    return
  end function screenSurfaceDensityMetalsConstructorParameters

  function screenSurfaceDensityMetalsConstructorInternal(coefficient,dustExtinctionCurve_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` dust attenuation class.
    !!}
    implicit none
    type            (dustAttenuationScreenSurfaceDensityMetals)                        :: self
    double precision                                           , intent(in   )         :: coefficient
    class           (dustExtinctionCurveClass                 ), intent(in   ), target :: dustExtinctionCurve_
    !![
    <constructorAssign variables="coefficient, *dustExtinctionCurve_"/>
    !!]

    return
  end function screenSurfaceDensityMetalsConstructorInternal

  double precision function screenSurfaceDensityMetalsDepthOpticalV(self,node,componentType) result(depthOpticalV)
    !!{RST
    Return the :math:`V`-band optical depth of a screen scaling with the surface density of metals.
    !!}
    use :: Numerical_Constants_Astronomical, only : hydrogenByMassSolar, massSolar, opticalDepthToMagnitudes, parsec
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : hecto              , mega
    implicit none
    class           (dustAttenuationScreenSurfaceDensityMetals), intent(inout)         :: self
    type            (treeNode                                 ), intent(inout), target :: node
    type            (enumerationComponentTypeType             ), intent(in   )         :: componentType
    ! A_V/E(B-V) (Savage & Mathis 1979).
    double precision                                           , parameter             :: AVToEBV                  =+3.10d+00
    ! N_H/E(B-V), in atoms/cm²/mag (Savage & Mathis 1979).
    double precision                                           , parameter             :: NHToEBV                  =+5.80d+21
    ! Optical depth per unit surface density of metals, in units of (M☉/pc²)⁻¹. Note that the classes replaced
    ! here derived the optical-depth-to-magnitudes conversion locally as 2.5 log10(e); `opticalDepthToMagnitudes` is
    ! the algebraically identical 2.5/ln(10), and differs from it by one unit in the last place (2e-16 relative), so
    ! attenuations computed here differ from those of the originals only at that level.
    double precision                                           , parameter             :: depthOpticalNormalization=+AVToEBV                  &
         &                                                                                                          /NHToEBV                  &
         &                                                                                                          *hydrogenByMassSolar      &
         &                                                                                                          /atomicMassUnit*massSolar &
         &                                                                                                          /(                        &
         &                                                                                                            +parsec                 &
         &                                                                                                            *hecto                  &
         &                                                                                                          )**2                      &
         &                                                                                                          /metallicityISMLocal      &
         &                                                                                                          /opticalDepthToMagnitudes
    double precision                                                                   :: massGas                  , radius              , &
         &                                                                                metallicity              , densitySurfaceMetals

    call componentGasProperties(node,componentType,massGas,radius,metallicity)
    ! A component with no gas, or no size, has no dust.
    if (massGas <= 0.0d0 .or. radius <= 0.0d0) then
       depthOpticalV=0.0d0
       return
    end if
    ! Surface density of metals, in M☉/pc².
    densitySurfaceMetals=+metallicity          &
         &               *massGas              &
         &               /2.0d0                &
         &               /Pi                   &
         &               /(                    &
         &                 +mega               &
         &                 *radius             &
         &                )**2
    depthOpticalV       =+self%coefficient          &
         &               *depthOpticalNormalization &
         &               *densitySurfaceMetals
    return
  end function screenSurfaceDensityMetalsDepthOpticalV
