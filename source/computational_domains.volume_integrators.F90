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
Contains a module which provides computational domains.
!!}

module Computational_Domain_Volume_Integrators
  !!{
  Provides a class that implements computational domains.
  !!}
  use :: Coordinates, only : coordinate
  private
  
  !![
  <functionClass>
   <name>computationalDomainVolumeIntegrator</name>
   <descriptiveName>Computational Domain Volume Integrators</descriptiveName>
   <description>Class providing volume integrators over computational domains.</description>
   <default>cartesian3D</default>
   <method name="volume" >
    <description>Gives the volume of the integrated region.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="integrate" >
    <description>Integrate a function over the computational domain.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>procedure(computationalDomainVolumeIntegrand) :: integrand</argument>
   </method>
  </functionClass>
  !!]

  abstract interface
     double precision function computationalDomainVolumeIntegrand(coordinates)
       !!{
       Interface for integrands used by computational domain volume integrators.
       !!}
       import coordinate
       class(coordinate), intent(in   ) :: coordinates
     end function computationalDomainVolumeIntegrand
  end interface
  
end module Computational_Domain_Volume_Integrators
