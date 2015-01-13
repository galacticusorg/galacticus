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

!% Contains a module that implements calculations of accretion disk spectra.

module Accretion_Disk_Spectra
  !% Implements calculations of accretion disk spectra.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="accretionDiskSpectra" type="functionModules" >
  include 'accretionDiskSpectra.functionModules.inc'
  !# </include>
  private

  !# <include directive="accretionDiskSpectra" type="function" >
  !#  <descriptiveName>Accretion Disk Spectra</descriptiveName>
  !#  <description>Class providing spectra of accretion disks.</description>
  !#  <default>hopkins2007</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="spectrum" >
  !#   <description>Returns the spectrum (in units of $L_\odot$~Hz$^{-1}$) of the accretion disk at the given wavelength (in units of \AA) for {\normalfont \ttfamily node} .</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: wavelength</argument>
  !#  </method>
  include 'accretionDiskSpectra.type.inc'
  !# </include>

end module Accretion_Disk_Spectra
