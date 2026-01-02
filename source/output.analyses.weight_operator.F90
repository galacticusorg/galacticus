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
Contains a module which provides a class that implements weight on operators for on-the-fly analyses.
!!}

module Output_Analysis_Weight_Operators
  !!{
  Provides a class that implements operators on weights for on-the-fly analyses.
  !!}
  use            :: Galacticus_Nodes       , only : treeNode
  use            :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyTypeType, enumerationOutputAnalysisPropertyQuantityType
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  private

  !![
  <functionClass>
   <name>outputAnalysisWeightOperator</name>
   <descriptiveName>Output Analysis Weight Operator</descriptiveName>
   <description>Class providing operators on weights for on-the-fly analysis of outputs.</description>
   <default>identity</default>
   <method name="operate" >
    <description>Operate on the given weight.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                                               , intent(in   ) :: weightValue                             </argument>
    <argument>type            (treeNode                                     ), intent(inout) :: node                                    </argument>
    <argument>double precision                                               , intent(in   ) :: propertyValue   , propertyValueIntrinsic</argument>
    <argument>type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType                            </argument>
    <argument>type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity                        </argument>
    <argument>integer         (c_size_t                                     ), intent(in   ) :: outputIndex                             </argument>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Weight_Operators
