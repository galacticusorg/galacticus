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
Contains a module which provides an object that implements concentrations of dark matter halo profiles.
!!}

module Dark_Matter_Profiles_Concentration
  !!{
  Provides a class that implements concentrations of dark matter halo profiles.
  !!}
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO , darkMatterProfileDMOClass
  use :: Galacticus_Nodes        , only : treeNode
  use :: Virial_Density_Contrast , only : virialDensityContrast, virialDensityContrastClass
  implicit none
  private

  !![
  <functionClass>
   <name>darkMatterProfileConcentration</name>
   <descriptiveName>Dark Matter Profile Concentrations</descriptiveName>
   <description>
    Class providing dark matter profile concentrations.
   </description>
   <default>gao2008</default>
   <method name="concentration" >
    <description>Returns the concentration parameter for the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="concentrationMean" >
    <description>Returns the mean concentration parameter for a {\normalfont \ttfamily node} of the given mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
    <code>darkMatterProfileConcentrationConcentrationMean=self%concentration(node)</code>
   </method>
   <method name="densityContrastDefinition" >
    <description>Returns a {\normalfont \ttfamily virialDensityContrast} object describing the virial density contrast used to define this concentration.</description>
    <type>class(virialDensityContrastClass)</type>
    <pass>yes</pass>
   </method>
   <method name="darkMatterProfileDMODefinition" >
    <description>Returns a {\normalfont \ttfamily darkMatterProfileDMO} object describing the dark matter density profile used to define this concentration.</description>
    <type>class(darkMatterProfileDMOClass)</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Profiles_Concentration
