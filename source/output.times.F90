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
Contains a module which provides a class that implements output times for \glc.
!!}

module Output_Times
  !!{
  Provides a class that implements output times for \glc.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  private

  !![
  <functionClass>
   <name>outputTimes</name>
   <descriptiveName>Output Times</descriptiveName>
   <description>Class providing output times for \glc.</description>
   <default>list</default>
   <method name="count" >
    <description>Return the number of output times.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
   </method>
   <method name="time" >
    <description>Return the output time index by {\normalfont \ttfamily indexOutput}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: indexOutput</argument>
   </method>
   <method name="redshift" >
    <description>Return the output redshift index by {\normalfont \ttfamily indexOutput}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: indexOutput</argument>
   </method>
   <method name="index" >
    <description>Return the index of the output at the given {\normalfont \ttfamily time}. If {\normalfont \ttfamily findClosest} is given and is true then the closest matching output is returned.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time       </argument>
    <argument>logical         , intent(in   ), optional :: findClosest</argument>
   </method>
   <method name="timeNext" >
    <description>Given a {\normalfont \ttfamily time}, return the time of the next output, and (optionally) the index of that output.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   )           :: timeCurrent</argument>
    <argument>integer         (c_size_t), intent(  out), optional :: indexOutput</argument>
   </method>
   <method name="timePrevious" >
    <description>Given a {\normalfont \ttfamily time}, return the time of the previous output, and (optionally) the index of that output.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   )           :: timeCurrent</argument>
    <argument>integer         (c_size_t), intent(  out), optional :: indexOutput</argument>
   </method>
  </functionClass>
  !!]

end module Output_Times
