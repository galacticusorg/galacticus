!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module that implements calculations of accretion disk spectra.

module Accretion_Disk_Spectra
  !% Implements calculations of accretion disk spectra.
  use :: Galacticus_Nodes, only : treeNode
  private

  !# <functionClass>
  !#  <name>accretionDiskSpectra</name>
  !#  <descriptiveName>Accretion Disk Spectra</descriptiveName>
  !#  <description>Class providing spectra of accretion disks.</description>
  !#  <default>hopkins2007</default>
  !#  <method name="spectrumNode" >
  !#   <description>Returns the spectrum (in units of $L_\odot$~Hz$^{-1}$) of the accretion disk at the given wavelength (in units of \AA) for {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: wavelength</argument>
  !#  </method>
  !#  <method name="spectrumMassRate" >
  !#   <description>Returns the spectrum (in units of $L_\odot$~Hz$^{-1}$) of the accretion disk at the given wavelength (in units of \AA) for a specified accretion rate, and radiative efficiency.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: accretionRate, efficiencyRadiative, wavelength</argument>
  !#  </method>
  !#  <generic name="spectrum">
  !#   <method>spectrumNode</method>
  !#   <method>spectrumMassRate</method>
  !#  </generic>
  !# </functionClass>

end module Accretion_Disk_Spectra
