!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a class to describe radiation fields.

module Radiation_Fields
  !% Implements a class to describe radiation fields.
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !# <functionClass>
  !#  <name>radiationField</name>
  !#  <descriptiveName>Radiation Fields</descriptiveName>
  !#  <description>Class providing radiation fields.</description>
  !#  <default>null</default>
  !#  <method name="flux">
  !#   <description>Return the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision          , intent(in   ) :: wavelength</argument>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="integrateOverCrossSection">
  !#   <description>Integrates the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure between the wavelengths given in {\normalfont \ttfamily wavelengthRange} over a cross section specified by the function {\normalfont \ttfamily crossSectionFunction}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision          , dimension(2), intent(in   ) :: wavelengthRange</argument>
  !#   <argument>double precision          , external                    :: crossSectionFunction</argument>
  !#   <argument>type            (treeNode)              , intent(inout) :: node</argument>
  !#   <code>
  !#    radiationFieldIntegrateOverCrossSection=radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
  !#   </code>
  !#  </method>
  !# </functionClass>

  ! Module global variables for use in integrand routines.
  class    (radiationFieldClass), pointer :: selfGlobal
  type     (treeNode           ), pointer :: nodeGlobal
  procedure(double precision   ), pointer :: crossSectionFunctionGlobal
  !$omp threadprivate(selfGlobal,nodeGlobal,crossSectionFunctionGlobal)

contains

  double precision function radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
    !% Integrate the photon number of the radiation field over a given cross-section function (which should return the cross
    !% section in units of cm$^2$), i.e.:
    !% \begin{equation}
    !% {4 \pi \over \mathrm{h}} \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) j_{\nu}(\lambda) {\mathrm{d}\lambda \over \lambda},
    !% \end{equation}
    !% where $j_{\nu}$ is the flux of energy per unit area per unit solid angle and per unit frequency.
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : plancksConstant
    use :: Numerical_Constants_Units   , only : ergs
    use :: Numerical_Integration2      , only : integratorAdaptiveCompositeTrapezoidal1D
    use :: Timers
    implicit none
    class           (radiationFieldClass                     ), target      , intent(inout) :: self
    double precision                                          , dimension(2), intent(in   ) :: wavelengthRange
    double precision                                          , external                    :: crossSectionFunction
    type            (treeNode                                ), target      , intent(inout) :: node
    type            (integratorAdaptiveCompositeTrapezoidal1D)                              :: integrator_
    
    ! Set module-scope pointers to self and the cross-section function for use in the integrand routine.
    selfGlobal                 => self
    nodeGlobal                 => node
    crossSectionFunctionGlobal => crossSectionFunction
    ! Perform the integration. An adapative, composite trapezoidal integrator is used here - this is efficient since typically
    ! these integrands can involve cross-sections with sharp edges, and the radiation field can also have sharp edges (and is
    ! often approximated using a tabulated function).
    call integrator_%integrandSet(integrand        =crossSectionIntegrand)
    call integrator_%toleranceSet(toleranceRelative=1.0d-3               )
    radiationFieldIntegrateOverCrossSection_=integrator_%evaluate(wavelengthRange(1),wavelengthRange(2))    
    ! Scale result by multiplicative prefactors to give answer in units of inverse seconds.
    radiationFieldIntegrateOverCrossSection_=+radiationFieldIntegrateOverCrossSection_ &
         &                                   *4.0d0                                    &
         &                                   *Pi                                       &
         &                                   *ergs                                     &
         &                                   /plancksConstant
    return
  end function radiationFieldIntegrateOverCrossSection_

  double precision function crossSectionIntegrand(wavelength)
    !% Integrand function use in integrating a radiation field over a cross section function.
    double precision, intent(in   ) :: wavelength

    if (wavelength > 0.0d0) then
       crossSectionIntegrand=+crossSectionFunctionGlobal     (wavelength           ) &
            &                *selfGlobal                %flux(wavelength,nodeGlobal) &
            &                /                                wavelength
    else
       crossSectionIntegrand=0.0d0
    end if
    return
  end function crossSectionIntegrand

end module Radiation_Fields
