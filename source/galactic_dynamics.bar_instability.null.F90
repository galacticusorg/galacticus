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

!% Contains a module which implements a null calculation of bar instability.

module Galactic_Dynamics_Bar_Instabilities_Null
  !% Implements a null calculation of bar instability.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galactic_Dynamics_Bar_Instabilities_Null_Initialize

contains

  !# <barInstabilityMethod>
  !#  <unitName>Galactic_Dynamics_Bar_Instabilities_Null_Initialize</unitName>
  !# </barInstabilityMethod>
  subroutine Galactic_Dynamics_Bar_Instabilities_Null_Initialize(barInstabilityMethod,Bar_Instability_Timescale_Get)
    !% Initializes the ``Null'' bar instability module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                ), intent(in   )          :: barInstabilityMethod
    procedure(Bar_Instability_Timescale_Null), intent(inout), pointer :: Bar_Instability_Timescale_Get

    if (barInstabilityMethod == 'null') Bar_Instability_Timescale_Get => Bar_Instability_Timescale_Null

    return
  end subroutine Galactic_Dynamics_Bar_Instabilities_Null_Initialize

  double precision function Bar_Instability_Timescale_Null(thisNode)
    !% Assumes that disks are never bar unstable and so returns an infinite timescale for bar instability.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Assume infinite timescale (i.e. no instability).
    Bar_Instability_Timescale_Null=-1.0d0

    return
  end function Bar_Instability_Timescale_Null

end module Galactic_Dynamics_Bar_Instabilities_Null
