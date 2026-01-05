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
Contains a module which provieds a class that implements importing of data from N-body simulations.
!!}

module NBody_Importers
  !!{
  Provides a class that implements importing of data from N-body simulations.
  !!}
  use :: NBody_Simulation_Data, only : nBodyData
  private

  !![
  <functionClass>
   <name>nbodyImporter</name>
   <descriptiveName>N-Body Simulation Data Importer</descriptiveName>
   <description>Class providing importing of data from N-body simulations.</description>
   <default>gadgetHDF5</default>
   <method name="import" >
    <description>Import position and velocity data from the named N-body data file.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(nBodyData), intent(  out), allocatable, dimension(:) :: simulations</argument>
   </method>
   <method name="isHDF5" >
    <description>Return true if the imported data is from an HDF5 file (to which new data can be written).</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

  type, public :: nbodyImporterList
     !!{
     Class used to build linked list of N-body data importers.
     !!}
     class(nbodyImporterClass), pointer                   :: importer_   => null()
     type (nbodyImporterList ), pointer                   :: next        => null()
     type (nBodyData         ), allocatable, dimension(:) :: simulations
  end type nbodyImporterList

  type, public :: nbodyPropertiesRealList
     !!{
     Class used to construct lists of N-body data scalar real properties.
     !!}
     double precision, pointer, dimension(:) :: property
  end type nbodyPropertiesRealList

  type, public :: nbodyPropertiesRealRank1List
     !!{
     Class used to construct lists of N-body data rank-1 real properties.
     !!}
     double precision, pointer, dimension(:,:) :: property
  end type nbodyPropertiesRealRank1List

  type, public :: nbodyPropertiesIntegerList
     !!{
     Class used to construct lists of N-body data scalar integer properties.
     !!}
     integer(c_size_t), pointer, dimension(:) :: property
  end type nbodyPropertiesIntegerList

end module NBody_Importers
