!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which provides a class that implements importing of data from N-body simulations.

module NBody_Importers
  !% Provides a class that implements importing of data from N-body simulations.
  use NBody_Simulation_Data
  private
  
  !# <functionClass>
  !#  <name>nbodyImporter</name>
  !#  <descriptiveName>N-Body Simulation Data Importer</descriptiveName>
  !#  <description>Class providing importing of data from N-body simulations.</description>
  !#  <default>gadgetHDF5</default>
  !#  <method name="import" >
  !#   <description>Import position and velocity data from the named N-body data file.</description>
  !#   <type>type(nBodyData)</type>
  !#   <pass>yes</pass>
  !#   <argument>character(len=*), intent(in   )           :: fileName</argument>
  !#   <argument>character(len=*), intent(in   ), optional :: fileNamePrevious</argument>
  !#  </method>
  !# </functionClass>
  
end module NBody_Importers
