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
Contains a module which implements a merger tree branching probability class.
!!}

module Merger_Tree_Branching
  !!{
  Implements a merger tree branching probability class.
  !!}
  use :: Galacticus_Nodes        , only : treeNode
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  implicit none
  private

  !![
  <enumeration>
   <name>mergerTreeBranchingBound</name>
   <description>Upper/lower bound labels used in merger tree branching calculations.</description>
   <entry label="lower"/>
   <entry label="upper"/>
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>mergerTreeBranchingProbability</name>
   <descriptiveName>Merger Tree Branching Probabilities</descriptiveName>
   <description>
    Class providing merger tree branching probabilities.
   </description>
   <default>parkinsonColeHelly</default>
   <method name="rate" >
    <description>Computes the probability per unit ``time'' of branching at the given mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )         :: mass, deltaCritical, time, massBranch</argument>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="probability" >
    <description>Computes the probability per unit ``time'' that a branching event occurs.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )         :: haloMass, deltaCritical, time, massResolution</argument>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="probabilityBound" >
    <description>Computes a bound (upper or lower) to the probability per unit ``time'' that a branching event occurs.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                                         , intent(in   )         :: haloMass, deltaCritical, time, massResolution</argument>
    <argument>type            (enumerationMergerTreeBranchingBoundType), intent(in   )         :: bound</argument>
    <argument>type            (treeNode                               ), intent(inout), target :: node</argument>
   </method>
   <method name="fractionSubresolution" >
    <description>Computes the fraction of subresolution mass accreted per unit ``time''</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )         :: haloMass, deltaCritical, time, massResolution</argument>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="massBranch" >
    <description>Returns the mass of a new halo created by a branching event.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                            , intent(in   )         :: haloMass             , deltaCritical, time, massResolution, probabilityFraction</argument>
    <argument>class           (randomNumberGeneratorClass), intent(inout)         :: randomNumberGenerator_</argument>
    <argument>type            (treeNode                  ), intent(inout), target :: node</argument>
   </method>
   <method name="stepMaximum" >
    <description>Returns the maximum step in ``time'' allowed by this algorithm.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: haloMass, deltaCritical, time, massResolution</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Branching
