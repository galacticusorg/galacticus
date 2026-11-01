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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides a class that implements heating of galactic disks.
!!}

module Galactic_Disk_Heatings
  !!{RST
  Provides a class that implements heating of galactic disks.
  !!}
  private

  !![
  <functionClass docformat="rst">
   <name>galacticDiskHeating</name>
   <descriptiveName>Galactic Disk Heatings</descriptiveName>
   <description>
   Class providing models of the heating of galactic disks by tidal shocks. The degree of heating is described by a dimensionless
   parameter, :math:`s`, such that a ring of the disk at radius :math:`R` is spread in radius with a dispersion :math:`s
   R^2/R_\mathrm{d}`, where :math:`R_\mathrm{d}` is the scale length of the disk (see
   :galacticus-class:`massDistributionExponentialDiskTidal`). Heating is applied after any truncation of the disk (see
   :galacticus-class:`galacticComponentTruncationClass`), and only to satellites.
   </description>
   <default>none</default>
   <method name="isActive" >
    <description>
    Returns true if this class heats disks.
    </description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="heating" >
    <description>
    Returns the dimensionless heating parameter, :math:`s`, of a disk of scale length ``radius``, and circular velocity
    ``velocity`` at that scale length, in a satellite which has experienced tidal heating ``tidalHeatingNormalized`` (see
    :galacticus-class:`nodeOperatorSatelliteTidalHeating`).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: tidalHeatingNormalized, radius, velocity</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Disk_Heatings
