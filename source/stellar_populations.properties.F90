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
Contains a module which implements a class for computing properties of stellar populations.
!!}

module Stellar_Population_Properties
  !!{
  Implements a class for computing properties of stellar populations.
  !!}
  use :: Abundances_Structure, only : abundances
  use :: Galacticus_Nodes    , only : nodeComponent, treeNode
  use :: Histories           , only : history
  implicit none
  private

  !![
  <functionClass>
   <name>stellarPopulationProperties</name>
   <descriptiveName>Stellar Population Properties</descriptiveName>
   <description>
    Class providing stellar population properties---essentially the rates of change of stellar and gas mass and abundances given
    a star formation rate and fuel abundances (and perhaps a historical record of star formation in the component).
   </description>
   <default>instantaneous</default>
   <method name="rates" >
    <description>Returns rates of change of stellar population properties.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                     , intent(in   ) :: rateStarFormation                                                   </argument>
    <argument>type            (abundances         ), intent(in   ) :: abundancesFuel                                                      </argument>
    <argument>class           (nodeComponent      ), intent(in   ) :: component                                                           </argument>
    <argument>type            (treeNode           ), intent(inout) :: node                                                                </argument>
    <argument>type            (history            ), intent(inout) :: history_                                                            </argument>
    <argument>double precision                     , intent(  out) :: rateMassStellar             , rateMassFuel         , rateEnergyInput</argument>
    <argument>type            (abundances         ), intent(inout) :: rateAbundancesFuel          , rateAbundancesStellar                 </argument>
    <argument>type            (stellarLuminosities), intent(inout) :: rateLuminosityStellar                                               </argument>
    <argument>logical                              , intent(in   ) :: computeRateLuminosityStellar                                        </argument>
   </method>
   <method name="scales">
    <description>Return scaling factors of stellar population properties for an \gls{ode} solver.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision            , intent(in   ) :: massStellar      </argument>
    <argument>type            (abundances), intent(in   ) :: abundancesStellar</argument>
    <argument>type            (history   ), intent(inout) :: history_         </argument>
   </method>
   <method name="historyCount">
    <description>Return the number of stellar population property histories which must be stored.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="historyCreate">
    <description>Create histories needed to store stellar population properties.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node    </argument>
    <argument>type(history ), intent(inout) :: history_</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Population_Properties
