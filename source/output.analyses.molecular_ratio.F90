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

module Output_Analysis_Molecular_Ratios
  !!{
  Provides a class that implements operators on properties for on-the-fly analyses.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>outputAnalysisMolecularRatio</name>
   <descriptiveName>Output Analysis Molecular Ratio</descriptiveName>
   <description>Class providing H$_2$ molecular ratios for on-the-fly analysis of outputs.</description>
   <default>obreschkow2009</default>
   <method name="ratio" >
    <description>Return the molecular ratio, $R_\mathrm{mol}=M_\mathrm{H_2}/M_\mathrm{HI}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: massISM</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="ratioScatter" >
    <description>Return the scatter in logarithmic molecular ratio, $\log_{10}R_\mathrm{mol}=\log_{10}(M_\mathrm{H_2}/M_\mathrm{HI})$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: massISM</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Molecular_Ratios
