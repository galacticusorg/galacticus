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
Contains a module which implements a class for calculations of black hole binary mergers.
!!}

module Black_Hole_Binary_Mergers
  !!{
  Implements a class for calculations of black hole binary mergers.
  !!}
  implicit none
  private

  !![
  <functionClass>
   <name>blackHoleBinaryMerger</name>
   <descriptiveName>Black Hole Binaries Merger</descriptiveName>
   <description>
    Class providing models of black hole binary mergers.
   </description>
   <default>rezzolla2008</default>
   <method name="merge" >
    <description>The properties of the black hole resulting from a binary merger.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massBlackHoleA    , massBlackHoleB, spinBlackHoleA    , spinBlackHoleB</argument>
    <argument>double precision, intent(  out) :: massBlackHoleFinal                , spinBlackHoleFinal</argument>
   </method>
  </functionClass>
  !!]

end module Black_Hole_Binary_Mergers
