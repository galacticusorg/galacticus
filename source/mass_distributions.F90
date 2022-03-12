!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a class that provides mass distributions.
!!}

module Mass_Distributions
  !!{
  Implements a class that provides mass distributions.
  !!}
  use :: Coordinates             , only : coordinate
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use :: Tensors                 , only : tensorRank2Dimension3Symmetric
  private

  !![
  <functionClass>
   <name>massDistribution</name>
   <descriptiveName>Mass Distributions</descriptiveName>
   <description>Class providing mass distributions.</description>
   <method name="symmetry" >
    <description>Return the symmetry of the distribution.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="isDimensionless" >
    <description>Return true if the distribution is dimensionless.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     massDistributionIsDimensionless=self%dimensionless
    </code>
   </method>
   <method name="massTotal" >
    <description>Return the total mass of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="acceleration" >
    <description>Return the gravitational acceleration due to the distribution at the given coordinates.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
   </method>
   <method name="tidalTensor" >
    <description>Return the gravitational tidal tensor due to the distribution at the given coordinates.</description>
    <type>type(tensorRank2Dimension3Symmetric)</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
   </method>
   <method name="density" >
    <description>Return the density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
   </method>
   <method name="densityGradientRadial" >
    <description>Return the radial gradient of density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   )           :: coordinates</argument>
    <argument>logical          , intent(in   ), optional :: logarithmic</argument>
   </method>
   <method name="potential" >
    <description>Return the gravitational potential of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
   </method>
   <method name="massEnclosedBySphere" >
    <description>Return the mass enclosed in the distribution by a sphere of given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ) :: radius</argument>
   </method>
   <method name="densityRadialMoment" >
    <description>Return the radial moment of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: moment</argument>
    <argument>double precision, intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
    <argument>logical         , intent(  out), optional :: isInfinite</argument>
   </method>
   <method name="positionSample" >
    <description>Return a position sampled from the distribution.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(randomNumberGeneratorClass), intent(inout) :: randomNumberGenerator_</argument>
   </method>
   <data>logical :: dimensionless</data>
  </functionClass>
  !!]

  ! Enumeration of mass distribution symmetries.
  !![
  <enumeration>
   <name>massDistributionSymmetry</name>
   <description>Specifies the symmetry of {\normalfont \ttfamily massDistribution} objects.</description>
   <visibility>public</visibility>
   <entry label="none"        />
   <entry label="cylindrical" />
   <entry label="spherical"   />
  </enumeration>
  !!]

end module Mass_Distributions
