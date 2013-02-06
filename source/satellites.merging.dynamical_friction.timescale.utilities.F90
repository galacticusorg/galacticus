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

!% Contains a module which implements utiilities for dynamical friction timescale calculations.

module Dynamical_Friction_Timescale_Utilities
  !% Implements utiilities for dynamical friction timescale calculations.
  implicit none
  private
  public :: Dynamical_Friction_Timescale_Multiplier

  ! Multiplier for the merging time.
  double precision :: mergingTimescaleMultiplier
  
  ! Flag indicating if the module is initialized.
  logical          :: dynamicalFrictionMultiplierInitialized=.false.

contains

  double precision function Dynamical_Friction_Timescale_Multiplier()
    !% Returns a multiplicative factor for scaling of dynamical friction timescales for satellite merging time calculations.
    use Input_Parameters
    implicit none

    if (.not.dynamicalFrictionMultiplierInitialized) then
       !$omp critical (Dynamical_Friction_Timescale_Multiplier_Initialize)
       if (.not.dynamicalFrictionMultiplierInitialized) then
          !@ <inputParameter>
          !@   <name>mergingTimescaleMultiplier</name>
          !@   <defaultValue>1.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     A multiplier for the merging timescale in dynamical friction timescale calculations.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergingTimescaleMultiplier',mergingTimescaleMultiplier,defaultValue=1.0d0)
          ! Flag that the module is now initialized. 
          dynamicalFrictionMultiplierInitialized=.true.
       end if
       !$omp end critical (Dynamical_Friction_Timescale_Multiplier_Initialize)
    end if

    ! Return the stored multiplier.
    Dynamical_Friction_Timescale_Multiplier=mergingTimescaleMultiplier
    return
  end function Dynamical_Friction_Timescale_Multiplier
  
end module Dynamical_Friction_Timescale_Utilities
