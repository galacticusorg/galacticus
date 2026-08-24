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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides a class implementing the inclination of a galaxy to the line of sight.
!!}

module Galactic_Inclinations
  !!{RST
  Provides a class implementing the inclination of a galaxy to the line of sight.

  Inclination is the angle between a galaxy's symmetry axis and the direction to the observer, so that zero is
  face-on and :math:`\pi/2` edge-on. It is needed by any model of an anisotropic process---dust attenuation through a
  flattened dust distribution above all, but also H I line widths, disk size distributions, and surface brightness
  selection---and putting it behind one class guarantees that every such consumer sees the same angle for a given
  galaxy. That matters: if two dust attenuators drew their own inclinations, a galaxy's colors would be formed from
  luminosities computed at different orientations, and would be meaningless.

  Implementations differ in where the angle comes from, and in whether they need to store anything:

  * :galacticus-class:`galacticInclinationNull` supplies no inclination at all, and is the default. Consumers must
    test ``isAvailable`` and either refuse to run or average over orientation instead---see
    :galacticus-class:`dustAttenuationInclinationAveraged`. Costs nothing.
  * :galacticus-class:`galacticInclinationSpinVector` derives the angle from the halo angular momentum vector, and so
    needs no storage of its own. It is the physically motivated choice, and gives orientations correlated with
    large-scale structure, which matters for lightcones and mock catalogs.
  * :galacticus-class:`galacticInclinationRandom` assigns each galaxy an isotropically distributed orientation,
    stored in a meta-property of the ``basic`` component by a ``nodeOperator``, and so costs 8 bytes per node---but
    only if it is used.
  * :galacticus-class:`galacticInclinationFixed` returns a single angle for every galaxy, for testing and for
    face-on/edge-on comparisons.

  The angle is returned in radians. The tabulated atlases are in degrees, so a conversion has to happen somewhere;
  radians are the natural unit for the class itself, and ``degreesToRadians`` is available for the tables.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass docformat="rst">
   <name>galacticInclination</name>
   <descriptiveName>Galactic Inclinations</descriptiveName>
   <description>
   Class providing the inclination of a galaxy to the line of sight.
   </description>
   <default>null</default>
   <method name="inclination" >
    <description>
    Return the inclination of the galaxy in ``node`` to the line of sight, in radians, in the range
    :math:`[0,\pi/2]`. Zero is face-on---the symmetry axis pointing at the observer---and :math:`\pi/2` is edge-on.

    Only the angle between the axis and the line of sight is meaningful, not its sign, so implementations which
    derive the angle from a vector fold it into the first quadrant.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="isAvailable" >
    <description>
    Return true if this class supplies a per-galaxy inclination.

    Consumers which depend on orientation must test this at construction, and report a clear error if it is false,
    rather than silently using a meaningless angle. The default is true; only
    :galacticus-class:`galacticInclinationNull` returns false.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     galacticInclinationIsAvailable=.true.
    </code>
   </method>
  </functionClass>
  !!]

end module Galactic_Inclinations
