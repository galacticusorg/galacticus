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
Contains a module which provides a class that implements surface density rates of star formation in disks.
!!}

module Star_Formation_Rate_Surface_Density_Disks
  !!{
  Provides a class that implements surface density rates of star formation in disks.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>starFormationRateSurfaceDensityDisks</name>
   <descriptiveName>Surface density rates of star formation in disks.</descriptiveName>
   <description>
    Class providing models of the surface density rate of star formation in disks.
   </description>
   <default>krumholz2009</default>
   <method name="intervals" >
    <description>Return a set of integration intervals to use when integrating over the surface density of star formation rate.</description>
    <type>double precision, allocatable, dimension(:,:)</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target                    :: node                           </argument>
    <argument>double precision          , intent(in   )                            :: radiusInner       , radiusOuter</argument>
    <argument>logical                   , intent(inout), allocatable, dimension(:) :: intervalIsAnalytic             </argument>
    <argument>double precision          , intent(inout), allocatable, dimension(:) :: integralsAnalytic              </argument>
    <code>
     !$GLC attributes unused :: self, node
     allocate(starFormationRateSurfaceDensityDisksIntervals(2,1))
     allocate(intervalIsAnalytic                           (  1))
     intervalIsAnalytic                           =.false.
     starFormationRateSurfaceDensityDisksIntervals=reshape([radiusInner,radiusOuter],[2,1])
    </code>
   </method>
   <method name="unchanged" >
    <description>Return true if the surface density rate of star formation is unchanged since the previous evaluation.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
    <code>
     !$GLC attributes unused :: self, node
     starFormationRateSurfaceDensityDisksUnchanged=.false.
    </code>
   </method>
   <method name="rate" >
    <description>Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) in the disk component of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
  </functionClass>
  !!]

end module Star_Formation_Rate_Surface_Density_Disks
