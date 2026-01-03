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
Contains a module which provides a class implementing the virial density contrast for halos.
!!}

module Virial_Density_Contrast
  !!{
  Provides a class implementing the virial density contrast for halos.
  !!}
  private

  !![
  <functionClass>
   <name>virialDensityContrast</name>
   <descriptiveName>Virial Density Contrasts</descriptiveName>
   <description>
    Class providing dark matter halo virial mean density contrasts.
   </description>
   <default>sphericalCollapseClsnlssMttrCsmlgclCnstnt</default>
   <method name="densityContrast" >
    <description>Returns the virial density contrast at the given epoch.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: mass</argument>
    <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsing</argument>
   </method>
   <method name="densityContrastRateOfChange" >
    <description>Returns the rate of change of virial density contrast at the given epoch.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: mass</argument>
    <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsing</argument>
   </method>
   <method name="turnAroundOverVirialRadii" >
    <description>Returns the ratio of turnaround and virial radii at the given epoch.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>double precision, intent(in   )           :: mass</argument>
    <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsing</argument>
    <code>
     !$GLC attributes unused :: self, mass, time, expansionFactor, collapsing
     virialDensityContrastTurnaroundOverVirialRadii=0.0d0
     call Error_Report('ratio is undefined for the "'//char(self%objectType())//'" density contrast class'//{introspection:location})
    </code>
   </method>
   <method name="isMassDependent" >
    <description>Returns true if the virial density contrast is mass-dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     virialDensityContrastIsMassDependent=.false.
    </code>
   </method>
  </functionClass>
  !!]

end module Virial_Density_Contrast
