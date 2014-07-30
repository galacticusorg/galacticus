!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides an object that implements dark matter halo profiles.

module Dark_Matter_Profiles
  !% Provides an object that implements dark matter halo profiles.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Galacticus_Nodes
  use               FGSL
  !# <include directive="darkMatterProfile" type="functionModules" >
  include 'darkMatterProfile.functionModules.inc'
  !# </include>
  private

  !# <include directive="darkMatterProfile" type="function" >
  !#  <descriptiveName>Dark Matter Halo Profiles</descriptiveName>
  !#  <description>Object providing dark matter halo profiles.</description>
  !#  <default>NFW</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <calculationReset>yes</calculationReset>
  !#  <method name="density" >
  !#   <description> Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="energy" >
  !#   <description>Return the total energy for the given {\tt node} in units of $M_\odot$ km$^2$ s$^{-1}$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  !#  <method name="energyGrowthRate" >
  !#   <description> Returns the rate of chance of the total energy of {\tt node} in units of $M_\odot$ km$^2$ s$^{-1}$ Gyr$^{-1}$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  !#  <method name="rotationNormalization" >
  !#   <description> Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for the given {\tt node}. Specifically, the normalization, $A$, returned is such that $V_{\rm rot} = A J/M$</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  !#  <method name="radiusFromSpecificAngularMomentum" >
  !#   <description> Returns the radius (in Mpc) in the dark matter profile of {\tt node} at which the specific angular momentum of a circular orbit equals {\tt specificAngularMomentum} (specified in units of km s$^{-1}$ Mpc.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: specificAngularMomentum</argument>
  !#  </method>
  !#  <method name="circularVelocity" >
  !#   <description>Returns the circular velocity (in km/s) in the dark matter profile of {\tt node} at the given {\tt radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="circularVelocityMaximum" >
  !#   <description>Returns the maximum circular velocity (in km/s) in the dark matter profile of {\tt node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  !#  <method name="potential" >
  !#   <description>Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer  :: node</argument>
  !#   <argument>double precision          , intent(in   )           :: radius</argument>
  !#   <argument>integer                   , intent(  out), optional :: status</argument>
  !#  </method>
  !#  <method name="enclosedMass" >
  !#   <description>Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in units of Mpc). for the given {\tt node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="kSpace" >
  !#   <description>Returns the normalized Fourier space density profile of the dark matter profile of {\tt node} at the given {\tt waveNumber} (given in units of Mpc$^{-1}$).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: wavenumber</argument>
  !#  </method>
  !#  <method name="freefallRadius" >
  !#   <description>Returns the freefall radius (in Mpc) corresponding to the given {\tt time} (in Gyr) in {\tt node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: time</argument>
  !#  </method>
  !#  <method name="freeFallRadiusIncreaseRate" >
  !#   <description>Returns the rate of increase of the freefall radius (in Mpc/Gyr) corresponding to the given {\tt time} (in Gyr) in {\tt node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: time</argument>
  !#  </method>
  include 'darkMatterProfile.type.inc'
  !# </include>

end module Dark_Matter_Profiles
