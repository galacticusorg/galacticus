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
Contains a module which provides a class that implements linear growth of cosmological structure.
!!}

module Linear_Growth
  !!{
  Provides a class that implements linear growth of cosmological structure.
  !!}
  private

  ! Enumeration for normalization options.
  !![
  <enumeration>
   <name>normalize</name>
   <description>Specifies normalization options for linear growth factor.</description>
   <entry label="matterDominated" />
   <entry label="presentDay"      />
  </enumeration>
  !!]

  ! Enumeration for components.
  !![
  <enumeration>
   <name>component</name>
   <description>Specifies components for linear growth factor.</description>
   <entry label="darkMatter" />
   <entry label="baryons"    />
   <entry label="radiation"  />
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>linearGrowth</name>
   <descriptiveName>Linear Growth of Cosmological Structure</descriptiveName>
   <description>Object providing linear growth of cosmological structure.</description>
   <default>collisionlessMatter</default>
   <method name="value" >
    <description>Return the linear growth factor at the given time and mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>type            (enumerationNormalizeType), intent(in   ), optional :: normalize                  </argument>
    <argument>type            (enumerationComponentType), intent(in   ), optional :: component                  </argument>
    <argument>double precision                          , intent(in   ), optional :: wavenumber                 </argument>
   </method>
   <method name="logarithmicDerivativeExpansionFactor" >
    <description>Return the logarithmic derivative of linear growth factor with respect to expansion factor.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>type            (enumerationComponentType), intent(in   ), optional :: component                  </argument>
    <argument>double precision                          , intent(in   ), optional :: wavenumber                 </argument>
   </method>
   <method name="logarithmicDerivativeWavenumber" >
    <description>Return the logarithmic derivative of linear growth factor with respect to wavenumber.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>type            (enumerationComponentType), intent(in   ), optional :: component                  </argument>
    <argument>double precision                          , intent(in   ), optional :: wavenumber                 </argument>
   </method>
   <method name="isWavenumberDependent" >
    <description>Return true if the growth function is wavenumber-dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(enumerationComponentType), intent(in   ), optional :: component</argument>
   </method>
  </functionClass>
  !!]

end module Linear_Growth
