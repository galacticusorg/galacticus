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
Contains a module which provides a class that operators on distributions used in on-the-fly output analyses.
!!}

module Output_Analysis_Distribution_Operators
  !!{
  Provides a class that operators on distributions used in on-the-fly output analyses.
  !!}
  use            :: Galacticus_Nodes       , only : treeNode
  use            :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyTypeType
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  private

  !![
  <functionClass>
   <name>outputAnalysisDistributionOperator</name>
   <descriptiveName>Output Analysis Distribution Operator</descriptiveName>
   <description>Class providing operators on distributions for on-the-fly output analysis---transformations
    applied to the binned histogram of galaxy properties to convert it from an ideal theoretical distribution
    to the distribution that would be observed. Typical operations include convolution with a Gaussian
    measurement error kernel (to model photometric scatter), projection effects, or survey selection
    functions. Implementations can operate on a scalar property value (mapping it to a smeared distribution)
    or on an existing distribution (convolving it with a response function).</description>
   <default>identity</default>
   <method name="operateScalar" >
    <description>Map a single scalar galaxy property value to a smeared distribution across bins, applying the operator's response function (e.g.\ a Gaussian measurement error kernel) to produce a probability-weighted contribution to each bin.</description>
    <type>double precision, dimension(size(propertyValueMinimum))</type>
    <pass>yes</pass>
    <argument>double precision                                           , intent(in   )               :: propertyValue</argument>
    <argument>type            (enumerationOutputAnalysisPropertyTypeType), intent(in   )               :: propertyType</argument>
    <argument>double precision                                           , intent(in   ), dimension(:) :: propertyValueMinimum, propertyValueMaximum</argument>
    <argument>integer         (c_size_t                                 ), intent(in   )               :: outputIndex</argument>
    <argument>type            (treeNode                                 ), intent(inout)               :: node</argument>
   </method>
   <method name="operateDistribution" >
    <description>Transform an existing binned distribution by applying the operator's response function (e.g.\ convolution with a measurement error kernel or survey selection function), returning the modified distribution.</description>
    <type>double precision, dimension(size(propertyValueMinimum))</type>
    <pass>yes</pass>
    <argument>double precision                                           , intent(in   ), dimension(:) :: distribution</argument>
    <argument>type            (enumerationOutputAnalysisPropertyTypeType), intent(in   )               :: propertyType</argument>
    <argument>double precision                                           , intent(in   ), dimension(:) :: propertyValueMinimum, propertyValueMaximum</argument>
    <argument>integer         (c_size_t                                 ), intent(in   )               :: outputIndex</argument>
    <argument>type            (treeNode                                 ), intent(inout)               :: node</argument>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Distribution_Operators
