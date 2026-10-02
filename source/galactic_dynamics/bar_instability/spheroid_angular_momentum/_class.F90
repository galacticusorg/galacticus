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

!+    Contributions to this file made by: Andrew Robertson, Codex.

!!{RST
Contains a module providing models for the pseudo-angular momentum assigned to spheroids during bar instabilities.
!!}

module Bar_Instability_Spheroid_Angular_Momenta
  !!{RST
  Provides models for the pseudo-angular momentum assigned to spheroids during bar instabilities.
  !!}
  use :: Galactic_Structure_Solvers, only : galacticStructureSolverClass
  use :: Galacticus_Nodes          , only : treeNode
  implicit none
  private

  !![
  <functionClass docformat="rst">
   <name>barInstabilitySpheroidAngularMomentum</name>
   <descriptiveName>Bar-instability spheroid pseudo-angular momentum</descriptiveName>
   <description>
   Class providing models for the rate at which pseudo-angular momentum is assigned to a spheroid as disk material is transferred to it by a bar instability.
   </description>
   <default>retained</default>
   <method name="rate">
    <description>
    Return the bar-instability contribution to the spheroid pseudo-angular momentum rate.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: rateMassGasTransfer, rateMassStellarTransfer, rateAngularMomentumDisk, rateAngularMomentumTransfer, fractionAngularMomentumRetainedSpheroid</argument>
   </method>
   <method name="galacticStructureSolverSet">
    <description>
    Provide the galactic structure solver used during node evolution.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(galacticStructureSolverClass), intent(inout), target :: galacticStructureSolver</argument>
    <code>
      !$GLC attributes unused :: self, galacticStructureSolver
    </code>
   </method>
  </functionClass>
  !!]

end module Bar_Instability_Spheroid_Angular_Momenta
