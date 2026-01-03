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
Contains a module that provides a class implementing the chemical state.
!!}

module Chemical_States
  !!{
  Provides a class implementing the chemical state.
  !!}
  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Radiation_Fields             , only : radiationFieldClass
  private

  !![
  <functionClass>
   <name>chemicalState</name>
   <descriptiveName>Chemical State</descriptiveName>
   <description>
    Class providing chemical state of gas.
   </description>
   <default>atomicCIECloudy</default>
   <method name="electronDensity" >
    <description>Return the electron density at the given temperature and hydrogen density for the specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
   </method>
   <method name="electronDensityTemperatureLogSlope" >
    <description>Return the logarithmic gradient of electron density with temperature at the given temperature and hydrogen density for the specified set of abundances and radiation field.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
   </method>
   <method name="electronDensityDensityLogSlope" >
    <description>Return the logarithmic gradient of electron density with respect to density at the given temperature and hydrogen density for the specified set of abundances and radiation field.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
   </method>
   <method name="chemicalDensities" >
    <description>Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (chemicalAbundances ), intent(inout) :: chemicalDensities</argument>
    <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
   </method>
  </functionClass>
  !!]

end module Chemical_States
