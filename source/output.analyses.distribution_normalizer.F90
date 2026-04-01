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
Contains a module which provides a class that normalizers on distributions used in on-the-fly output analyses.
!!}

module Output_Analysis_Distribution_Normalizers
  !!{
  Provides a class that normalizers on distributions used in on-the-fly output analyses.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>outputAnalysisDistributionNormalizer</name>
   <descriptiveName>Output Analysis Distribution Normalizer</descriptiveName>
   <description>Class providing normalizers for binned distributions in on-the-fly output analysis---operations
    that convert a raw histogram (e.g.\ galaxy counts per bin) into the desired normalized quantity
    (e.g.\ a number density per unit dex or a probability distribution). The \mono{normalize} method
    accepts the distribution array and its covariance matrix and modifies them in place, applying
    bin widths, volume normalizations, or other scale factors. Implementations include identity
    (no normalization), bin-width division, and survey-volume normalization for luminosity and mass
    functions.</description>
   <default>identity</default>
   <method name="normalize" >
    <description>Normalize the supplied binned distribution array and its covariance matrix in place, applying the bin-width or volume factors appropriate for this normalizer class.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(inout), dimension(:  ), optional :: distribution</argument>
    <argument>double precision, intent(inout), dimension(:,:), optional :: covariance</argument>
    <argument>double precision, intent(in   ), dimension(:  )           :: propertyValueMinimum, propertyValueMaximum</argument>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Distribution_Normalizers
