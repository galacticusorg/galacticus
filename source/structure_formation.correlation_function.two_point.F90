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
Contains a module which implements two-point correlation functions.
!!}

module Correlation_Functions_Two_Point
  !!{
  Implements two-point correlation functions.
  !!}
  private

  !![
  <functionClass>
   <name>correlationFunctionTwoPoint</name>
   <descriptiveName>Two-point Correlation Functions</descriptiveName>
   <description>Class providing two-point correlation functions.</description>
   <default>powerSpectrumTransform</default>
   <method name="correlation" >
    <description>Return the two-point correlation function for $r=${\normalfont \ttfamily separation} [Mpc].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: separation, time</argument>
   </method>
   <method name="correlationVolumeAveraged" >
    <description>Return the volume-averaged two-point correlation function for $r=${\normalfont \ttfamily separation} [Mpc].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: separation, time</argument>
   </method>
  </functionClass>
  !!]

end module Correlation_Functions_Two_Point
