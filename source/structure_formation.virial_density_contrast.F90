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

!% Contains a module which provides a class implementing the virial density contrast for halos.

module Virial_Density_Contrast
  !% Provides a class implementing the virial density contrast for halos.
  use ISO_Varying_String
  use Galacticus_Nodes
  use FGSL
  !# <include directive="virialDensityContrast" type="functionModules" >
  include 'virialDensityContrast.functionModules.inc'
  !# </include>
  private

  !# <include directive="virialDensityContrast" type="function" >
  !#  <descriptiveName>Virial Density Contrasts</descriptiveName>
  !#  <description>Class providing dark matter halo virial density contrasts.</description>
  !#  <default>sphericalCollapseMatterLambda</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <method name="densityContrast" >
  !#   <description>Returns the virial density contrast at the given epoch.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#  </method>
  !#  <method name="densityContrastRateOfChange" >
  !#   <description>Returns the rate of change of virial density contrast at the given epoch.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#  </method>
  !#  <method name="turnAroundOverVirialRadii" >
  !#   <description>Returns the ratio of turnaround and virial radii at the given epoch.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Galacticus_Error</modules>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#   <code>call Galacticus_Error_Report('turnAroundOverVirialRadii','ratio is undefined for this density contrast class')</code>
  !#  </method>
  include 'virialDensityContrast.type.inc'
  !# </include>

end module Virial_Density_Contrast
