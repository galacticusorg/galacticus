!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides a class that implements on-the-fly analyses.

module Output_Analysis_Property_Extractions
  !% Provides a class that implements extraction of properties for on-the-fly analyses.
  use Galacticus_Nodes
  private
  
  !# <functionClass>
  !#  <name>outputAnalysisPropertyExtractor</name>
  !#  <descriptiveName>Output Analysis Property Extractor</descriptiveName>
  !#  <description>Class providing extraction of properties for on-the-fly analysis of outputs.</description>
  !#  <default>null</default>
  !#  <method name="extract" >
  !#   <description>Extract the property from the given {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !# </functionClass>

end module Output_Analysis_Property_Extractions
