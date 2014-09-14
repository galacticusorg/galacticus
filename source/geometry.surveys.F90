!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements geometries of galaxy surveys.

module Geometry_Surveys
  !% Implements geometries of galaxy surveys.
  use ISO_Varying_String
  use, intrinsic :: ISO_C_Binding
  !# <include directive="surveyGeometry" type="functionModules" >
  include 'surveyGeometry.functionModules.inc'
  !# </include>
  implicit none
  private

  !# <include directive="surveyGeometry" type="function" >
  !#  <descriptiveName>Survey Geometry</descriptiveName>
  !#  <description>Object providing galaxy surveys geometries and related functions.</description>
  !#  <default>liWhite2009SDSS</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="fieldCount" >
  !#   <description>Returns the number of distinct fields included in the survey.</description>
  !#   <type>integer</type>
  !#   <pass>yes</pass>
  !#   <code>surveyGeometryFieldCount=1</code>
  !#  </method>
  !#  <method name="distanceMinimum" >
  !#   <description>Returns the minimum distance (in Mpc) at which a galaxy of the specified {\tt mass} (in $M_\odot$) would be included in the survey.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: mass</argument>
  !#   <argument>integer         , intent(in   ), optional :: field</argument>
  !#   <code>surveyGeometryDistanceMinimum=0.0d0</code>
  !#  </method>
  !#  <method name="distanceMaximum" >
  !#   <description>Returns the maximum distance (in Mpc) at which a galaxy of the specified {\tt mass} (in $M_\odot$) could be detected.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: mass</argument>
  !#   <argument>integer         , intent(in   ), optional :: field</argument>
  !#  </method>
  !#  <method name="solidAngle" >
  !#   <description>Return the solid angle (in steradians) of the survey.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>integer, intent(in   ), optional :: field</argument>
  !#  </method>
  !#  <method name="volumeMaximum" >
  !#   <description>Returns the maximum volume (in Mpc$^3$) at which a galaxy of the specified {\tt mass} (in $M_\odot$) could be detected.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: mass</argument>
  !#   <argument>integer         , intent(in   ), optional :: field</argument>
  !#   <code>surveyGeometryVolumeMaximum=self%solidAngle(field)*self%distanceMaximum(mass,field)**3/3.0d0</code>
  !#  </method>
  !#  <method name="windowFunctionAvailable" >
  !#   <description>Returns true if survey 3-D window functions are available.</description>
  !#   <type>logical</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="angularPowerAvailable" >
  !#   <description>Returns true if angular power spectrum of survey window function is available.</description>
  !#   <type>logical</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="windowFunctions" >
  !#   <description>Returns the window functions on a grid of the specified size ({\tt gridCount} cells in each dimension) for galaxies of the specified {\tt mass1} and {\tt mass2} (in $M_\odot$). The {\tt boxLength} should be set to an appropriate value to fully enclose (with sufficient buffering to allow for Fourier transformation) the two window functions.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision                  , intent(in   )                                           :: mass1          , mass2</argument>
  !#   <argument>integer                           , intent(in   )                                           :: gridCount</argument>
  !#   <argument>double precision                  , intent(  out)                                           :: boxLength</argument>
  !#   <argument>complex         (c_double_complex), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1, windowFunction2</argument>
  !#  </method>
  !#  <method name="angularPower" >
  !#   <description>Return $C^{ij}_\ell$, where $(2\ell+1) C^{ij}_\ell = \sum_{m=-\ell}^{+\ell} \Psi^i_{\ell m} \Psi^{j*}_{\ell m}$, and $\Psi^i_{\ell m}$ are the cofficients of the spherical harmonic expansion of the $i^{\rm th}$ field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>integer, intent(in   ) :: i,j,l</argument>
  !#  </method>
  !#  <method name="angularPowerMaximumDegree" >
  !#   <description>Return the maximum degree, $\ell_{\rm max}$, for which the angular power is available.</description>
  !#   <type>integer</type>
  !#   <pass>yes</pass>
  !#   <code>surveyGeometryAngularPowerMaximumDegree=-1</code>
  !#  </method>
  include 'surveyGeometry.type.inc'
  !# </include>

end module Geometry_Surveys
