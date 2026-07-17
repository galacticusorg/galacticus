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

!!{RST
Contains a module which provides a class that implements extraction of properties from nodes.
!!}

module Node_Property_Extractors
  !!{RST
  Provides a class that implements extraction of properties from nodes.
  !!}
  use :: Galacticus_Nodes       , only : treeNode
  use :: Multi_Counters         , only : multiCounter
  use :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyQuantityType, enumerationOutputAnalysisPropertyTypeType, outputAnalysisPropertyQuantityUnknown, outputAnalysisPropertyTypeLinear
  private

  !![
  <functionClass docformat="rst">
   <name>nodePropertyExtractor</name>
   <descriptiveName>Node Property Extractor</descriptiveName>
   <description>
   Class providing extraction of scalar, 1D, or multi-D properties from merger tree nodes for output. Property extractors are used by output analysis classes to retrieve galaxy and halo properties (e.g. stellar mass, dark matter halo mass, star formation rate, positions, velocities) and convert them to formats suitable for comparison with observational data or for writing to the Galacticus output file.
   </description>
   <default>nodeIndices</default>
   <method name="type" >
    <description>
    Return the type of the extracted property.
    </description>
    <type>type(enumerationOutputAnalysisPropertyTypeType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorType=outputAnalysisPropertyTypeLinear
    </code>
   </method>
   <method name="quantity" >
    <description>
    Return the class of the extracted property.
    </description>
    <type>type(enumerationOutputAnalysisPropertyQuantityType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorQuantity=outputAnalysisPropertyQuantityUnknown
    </code>
   </method>
   <method name="addInstances" >
    <description>
    Add multiple instances of this property to a ``multiCounter`` object.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode    ), intent(inout) :: node</argument>
    <argument>type(multiCounter), intent(inout) :: instance</argument>
    <code>
     !$GLC attributes unused :: self, node, instance
     ! Nothing to do.
    </code>
   </method>
   <method name="extractScalar" >
    <description>
    Extract a scalar property from the given ``node``. This is a convenience method for callers (such as the
    ``libgalacticus`` library interface) that hold an extractor through the base class: it dispatches to the
    ``extract`` method of the scalar extractor subclass, and reports an error for extractors of any other class
    (whose ``extract`` methods have different signatures).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <modules>Error</modules>
    <argument>type(treeNode), intent(inout), target :: node</argument>
    <code>
     select type (self)
     class is (nodePropertyExtractorScalar)
        nodePropertyExtractorExtractScalar=self%extract(node)
     class default
        nodePropertyExtractorExtractScalar=0.0d0
        call Error_Report('extractScalar requires an extractor of the scalar class'//{introspection:location})
     end select
    </code>
   </method>
  </functionClass>
  !!]

  ! Enumerations for galactic components.
  !![
  <enumeration docformat="rst">
   <name>galacticComponent</name>
   <description>
   Specifies the galactic component for various node property extractors.
   </description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <entry label="disk"              />
   <entry label="spheroid"          />
   <entry label="nuclearStarCluster"/>
   <entry label="total"             />
  </enumeration>
  !!]

end module Node_Property_Extractors
