!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the standard black hole component.

double precision function Node_Component_Black_Hole_Standard_Seed_Mass(self)
  !% Return the seed mass for standard black holes.
  use Input_Parameters
  implicit none
  class           (nodeComponentBlackHoleStandard), intent(inout) :: self   
  logical                                         , save          :: functionInitialized=.false.
  double precision                                , save          :: blackHoleSeedMass

  if (.not.functionInitialized) then
     !$omp critical (Node_Component_Black_Hole_Standard_Seed_Mass)
     if (.not.functionInitialized) then
        ! Get the black hole seed mass.
        !@ <inputParameter>
        !@   <name>blackHoleSeedMass</name>
        !@   <defaultValue>100</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     The mass of the seed black hole placed at the center of each newly formed galaxy.
        !@   </description>
        !@   <type>real</type>
        !@   <cardinality>1</cardinality>
        !@   <group>blackHoles</group>
        !@ </inputParameter>
        call Get_Input_Parameter("blackHoleSeedMass",blackHoleSeedMass,defaultValue=100.0d0)
        ! Record that the function is now initialized.
        functionInitialized=.true.
     end if
     !$omp end critical (Node_Component_Black_Hole_Standard_Seed_Mass)
  end if
  ! Return the seed mass.
  Node_Component_Black_Hole_Standard_Seed_Mass=blackHoleSeedMass
  return
end function Node_Component_Black_Hole_Standard_Seed_Mass

double precision function Node_Component_Black_Hole_Standard_Seed_Spin(self)
  !% Return the seed spin for standard black holes.
  use Input_Parameters
  implicit none
  class(nodeComponentBlackHoleStandard), intent(inout) :: self   

  Node_Component_Black_Hole_Standard_Seed_Spin=1.0d-3
  return
end function Node_Component_Black_Hole_Standard_Seed_Spin

double precision function Node_Component_Black_Hole_Standard_Spin(self)
  !% Return the spin of a standard black hole.
  implicit none
  class           (nodeComponentBlackHoleStandard), intent(inout) :: self   
  double precision                                , parameter     :: blackHoleSpinMaximum=0.9999d0 ! Maximum allowed spin (useful to avoid infinities).

  Node_Component_Black_Hole_Standard_Spin=self%spinData%value
  ! Ensure that the spin is within range.
  if (Node_Component_Black_Hole_Standard_Spin < 0.0d0) then
     Node_Component_Black_Hole_Standard_Spin=0.0d0
     call self%spinSet(0.0d0)
  end if
  if (Node_Component_Black_Hole_Standard_Spin > blackHoleSpinMaximum) then
     Node_Component_Black_Hole_Standard_Spin=blackHoleSpinMaximum
     call self%spinSet(blackHoleSpinMaximum)
  end if
  return
end function Node_Component_Black_Hole_Standard_Spin
