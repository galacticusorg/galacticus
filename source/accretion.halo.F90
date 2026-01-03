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
Contains a module which implements a class implementing accretion of gas from the \gls{igm} onto halos.
!!}

module Accretion_Halos
  !!{
  Implements a class implementing accretion of gas from the \gls{igm} onto halos.
  !!}
  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Galacticus_Nodes             , only : treeNode
  private

  ! Enumeration of accretion types.
  !![
  <enumeration>
   <name>accretionMode</name>
   <description>Enumeration of accretion modes for the {\normalfont \ttfamily accretionHalo} class.</description>
   <visibility>public</visibility>
   <entry label="total"/>
   <entry label="hot"  />
   <entry label="cold" />
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>accretionHalo</name>
   <descriptiveName>Accretion Onto Halos</descriptiveName>
   <description>
    Class providing rates of accretion of gas from the \gls{igm} onto halos. This is expected to depend on (at least) the rate
    at which that halo mass is growing, the depth of its potential well and the thermodynamical properties of the accreting
    gas.
   </description>
   <default>simple</default>
   <method name="branchHasBaryons" >
    <description>Returns {\normalfont \ttfamily true} if this tree branch may accrete baryons, and {\normalfont \ttfamily false} otherwise.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="accretionRate" >
    <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of mass from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="accretedMass" >
    <description>Returns the mass (in units of $M_\odot$) of accreted from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}. Used to initialize nodes.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="failedAccretionRate" >
    <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of failed accretion of mass from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="failedAccretedMass" >
    <description>Returns the mass (in units of $M_\odot$) that failed to accrete from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}. Used to initialize nodes.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="accretionRateMetals" >
    <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of metals from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}.</description>
    <type>type(abundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="accretedMassMetals" >
    <description>Returns the mass of metals (in units of $M_\odot$) of accreted from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}. Used to initialize nodes.</description>
    <type>type(abundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="failedAccretionRateMetals" >
    <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of failed accretion of metals from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}.</description>
    <type>type(abundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="failedAccretedMassMetals" >
    <description>Returns the mass of metals (in units of $M_\odot$) that failed to accrete from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}. Used to initialize nodes.</description>
    <type>type(abundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="accretionRateChemicals" >
    <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of chemicals from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}.</description>
    <type>type(chemicalAbundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
   <method name="accretedMassChemicals" >
    <description>Returns the mass of chemicals (in units of $M_\odot$) of accreted from the \gls{igm} onto {\normalfont \ttfamily node} in the given {\normalfont \ttfamily accretionMode}. Used to initialize nodes.</description>
    <type>type(chemicalAbundances)</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout) :: node</argument>
    <argument>type(enumerationAccretionModeType), intent(in   ) :: accretionMode</argument>
   </method>
  </functionClass>
  !!]

end module Accretion_Halos
