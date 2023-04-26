!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which provides a hot halo mass distribution class.
!!}

module Hot_Halo_Mass_Distributions
  !!{
  Provides an object which provides a hot halo mass distribution class.
  !!}
  use :: Galacticus_Nodes          , only : treeNode
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Galactic_Structure_Options, only : enumerationWeightByType
  private

  !![
  <functionClass>
   <name>hotHaloMassDistribution</name>
   <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
   <description>
    Object implementing hot halo mass distributions.
   </description>
   <default>betaProfile</default>
   <method name="get" >
    <description>Return the mass distribution of the hot halo.</description>
    <type>class(massDistributionClass)</type>
    <pass>yes</pass>
    <argument>type   (treeNode               ), intent(inout)           :: node       </argument>
    <argument>type   (enumerationWeightByType), intent(in   ), optional :: weightBy   </argument>
    <argument>integer                         , intent(in   ), optional :: weightIndex</argument>
   </method>
   <method name="density" >
    <description>Return the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="densityLogSlope" >
    <description>Return the logarithmic slope of the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   )          :: radius</argument>
   </method>
   <method name="enclosedMass" >
    <description>Return the mass enclosed in the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: radius</argument>
   </method>
   <method name="radialMoment" >
    <description>Return the specified radial{\normalfont \ttfamily moment} of the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: moment, radius</argument>
   </method>
   <method name="densitySquaredIntegral" >
    <description>Return the integral of the square of the density of the hot halo from zero to the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="rotationNormalization" >
    <description>Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Mass_Distributions
