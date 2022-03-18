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
Contains a module which provides a class that implements galactic structure functions.
!!}

module Galactic_Structure
  !!{
  Provides an object that implements galactic structure functions.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>galacticStructure</name>
   <descriptiveName>Galactic Structure</descriptiveName>
   <description>Object providing galactic structure functions.</description>
   <default>standard</default>
   <method name="density">
    <description>Return the mass enclosed at the given {\normalfont \ttfamily position} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)               :: node</argument>
    <argument>double precision          , intent(in   ), dimension(3) :: position</argument>
    <argument>integer                   , intent(in   ), optional     :: coordinateSystem, componentType, massType, weightBy, weightIndex</argument>
   </method>
   <method name="massEnclosed">
    <description>Return the mass enclosed within the given {\normalfont \ttfamily radius} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   ), optional :: radius</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType, weightBy, weightIndex</argument>
   </method>
   <method name="radiusEnclosingMass">
    <description>Return the radius enclosing a given mass (or fractional mass) in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target   :: node</argument>
    <argument>double precision          , intent(in   ), optional :: mass         , massFractional</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType      , weightBy, weightIndex</argument>
   </method>
   <method name="velocityRotation">
    <description>Return the rotation velocity for a circular orbit at the given {\normalfont \ttfamily radius} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   )           :: radius</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType</argument>
   </method>
   <method name="velocityRotationGradient">
    <description>Return the gradient of the rotation velocity for a circular orbit at the given {\normalfont \ttfamily radius} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   )           :: radius</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType</argument>
   </method>
   <method name="potential">
    <description>Return the gravitational potential at the given {\normalfont \ttfamily radius} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   )           :: radius</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType</argument>
    <argument>integer                   , intent(  out), optional :: status</argument>
   </method>
   <method name="surfaceDensity">
    <description>Return the surface density of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)               :: node</argument>
    <argument>double precision          , intent(in   ), dimension(3) :: position</argument>
    <argument>integer                   , intent(in   ), optional     :: coordinateSystem, componentType, massType, weightBy, weightIndex</argument>
   </method>
   <method name="radiusEnclosingSurfaceDensity">
    <description>Return the radius enclosing a given surface density in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target   :: node</argument>
    <argument>double precision          , intent(in   )           :: surfaceDensity</argument>
    <argument>integer                   , intent(in   ), optional :: componentType, massType, weightBy, weightIndex</argument>
   </method>
   <method name="acceleration">
    <description>Compute the gravitational acceleration at a given position in {\normalfont \ttfamily node}.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)               :: node</argument>
    <argument>double precision          , intent(in   ), dimension(3) :: positionCartesian</argument>
    <argument>integer                   , intent(in   ), optional     :: componentType, massType</argument>
   </method>
   <method name="tidalTensor">
    <description>Compute the gravitational tidal tensor at a given position in {\normalfont \ttfamily node}.</description>
    <type>type(tensorRank2Dimension3Symmetric)</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)               :: node</argument>
    <argument>double precision          , intent(in   ), dimension(3) :: positionCartesian</argument>
    <argument>integer                   , intent(in   ), optional     :: componentType, massType</argument>
   </method>
   <method name="chandrasekharIntegral">
    <description>Compute the integral appearing in the \cite{chandrasekhar_dynamical_1943} dynamical friction model in {\normalfont \ttfamily node}.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)               :: node</argument>
    <argument>double precision          , intent(in   ), dimension(3) :: positionCartesian, velocityCartesian</argument>
    <argument>integer                   , intent(in   ), optional     :: componentType    , massType</argument>
   </method>
   <method name="velocityDispersion">
    <description>Return the velocity dispersion at the given {\normalfont \ttfamily radius} in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: radius       , radiusOuter</argument>
    <argument>integer                   , intent(in   )         :: componentType, massType</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Structure
