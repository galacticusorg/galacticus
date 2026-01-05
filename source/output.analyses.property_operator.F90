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
Contains a module which provides a class that implements on-the-fly analyses.
!!}

module Output_Analysis_Property_Operators
  !!{
  Provides a class that implements operators on properties for on-the-fly analyses.
  !!}
  use            :: Galacticus_Nodes       , only : treeNode
  use            :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyTypeType
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  private

  !![
  <functionClass>
   <name>outputAnalysisPropertyOperator</name>
   <descriptiveName>Output Analysis Property Operator</descriptiveName>
   <description>Class providing operators on properties for on-the-fly analysis of outputs.</description>
   <default>identity</default>
   <method name="operate" >
    <description>Operate on the given property.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                                           , intent(in   )           :: propertyValue</argument>
    <argument>type            (treeNode                                 ), intent(inout), optional :: node         </argument>
    <argument>type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType </argument>
    <argument>integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex  </argument>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Property_Operators
