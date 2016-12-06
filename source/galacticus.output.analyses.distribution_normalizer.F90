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

!% Contains a module which provides a class that normalizors on distributions used in on-the-fly output analyses.

module Output_Analysis_Distribution_Normalizers
  !% Provides a class that normalizors on distributions used in on-the-fly output analyses.
  use Galacticus_Nodes
  private
  
  !# <functionClass>
  !#  <name>outputAnalysisDistributionNormalizer</name>
  !#  <descriptiveName>Output Analysis Distribution Normalizer</descriptiveName>
  !#  <description>Class providing normalizers on distributions for on-the-fly analysis of outputs.</description>
  !#  <default>identity</default>
  !#  <method name="normalize" >
  !#   <description>Normalize a distribution.</description>
  !#   <type>double precision, dimension(size(propertyValueMinimum))</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), dimension(:) :: distribution</argument>
  !#   <argument>double precision, intent(in   ), dimension(:) :: propertyValueMinimum, propertyValueMaximum</argument>
  !#  </method>
  !# </functionClass>

end module Output_Analysis_Distribution_Normalizers
