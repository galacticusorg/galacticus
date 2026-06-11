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
Contains a module which provides a class that packages the shared "target data" arguments used by 1D function output analyses (axis labels, log-scale flags, and target value/covariance arrays).
!!}

module Output_Analysis_Target_Data
  !!{RST
  Provides a class that packages the shared "target data" arguments used by the various 1D function output analyses (``outputAnalysisMeanFunction1D``, ``outputAnalysisScatterFunction1D``, ``outputAnalysisVolumeFunction1D``).  These fields are conceptually coupled --- a comparison dataset for a 1D function output --- so wrapping them into a single object both clarifies intent at call sites and collapses what would otherwise be a :math:`2^N` presence-combination explosion in the Python wrapper's optional-argument branching down to a single optional object argument on each outer constructor.
  !!}
  use :: ISO_Varying_String, only : varying_string
  private

  !![
  <functionClass docformat="rst">
   <name>outputAnalysisTargetData</name>
   <descriptiveName>Output Analysis Target Data</descriptiveName>
   <description>
   Class packaging the shared "target data" arguments used by 1D function output analyses (axis labels, log-scale flags, and target value/covariance arrays).  Bundling these into a single object lets downstream output-analysis constructors expose one optional argument instead of seven, which avoids a :math:`2^N` combinatorial explosion in the Python wrapper's optional-argument branching. All fields are themselves optional at construction; omitted labels default to empty (or ``'x'``/``'y'`` for the axis labels), omitted log-scale flags default to false, and the target arrays default to unallocated.
   </description>
   <default>standard</default>
   <method name="hasTarget" >
    <description>
    Return whether the value and covariance target arrays are both allocated (i.e.\ a comparison dataset is present).
    </description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Output_Analysis_Target_Data
