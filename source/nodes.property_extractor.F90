!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which provides a class that implements extraction of properties from nodes.
!!}

module Node_Property_Extractors
  !!{
  Provides a class that implements extraction of properties from nodes.
  !!}
  use :: Galacticus_Nodes       , only : treeNode
  use :: Multi_Counters         , only : multiCounter
  use :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyQuantityType, enumerationOutputAnalysisPropertyTypeType, outputAnalysisPropertyQuantityUnknown, outputAnalysisPropertyTypeLinear
  private

  !![
  <functionClass>
   <name>nodePropertyExtractor</name>
   <descriptiveName>Node Property Extractor</descriptiveName>
   <description>Class providing extraction of properties from nodes.</description>
   <default>nodeIndices</default>
   <method name="type" >
    <description>Return the type of the extracted property.</description>
    <type>type(enumerationOutputAnalysisPropertyTypeType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorType=outputAnalysisPropertyTypeLinear
    </code>
   </method>
   <method name="quantity" >
    <description>Return the class of the extracted property.</description>
    <type>type(enumerationOutputAnalysisPropertyQuantityType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorQuantity=outputAnalysisPropertyQuantityUnknown
    </code>
   </method>
   <method name="addInstances" >
    <description>Add multiple instances of this property to a {\normalfont \ttfamily multiCounter} object.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode    ), intent(inout) :: node</argument>
    <argument>type(multiCounter), intent(inout) :: instance</argument>
    <code>
     !$GLC attributes unused :: self, node, instance
     ! Nothing to do.
    </code>
   </method>
  </functionClass>
  !!]

  ! Enumerations for galactic components.
  !![
  <enumeration>
   <name>galacticComponent</name>
   <description>Specifies the galactic component for various node property extractors.</description>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="disk"              />
   <entry label="spheroid"          />
   <entry label="nuclearStarCluster"/>
   <entry label="total"             />
  </enumeration>
  !!]

end module Node_Property_Extractors
