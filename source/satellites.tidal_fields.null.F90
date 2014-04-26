!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null tidal field for satellites.

module Satellites_Tidal_Fields_Null
  !% Implements a null tidal field for satellites.
  use Galacticus_Nodes
  implicit none
  private
  public :: Satellites_Tidal_Fields_Null_Initialize

contains

  !# <satellitesTidalFieldMethod>
  !#  <unitName>Satellites_Tidal_Fields_Null_Initialize</unitName>
  !# </satellitesTidalFieldMethod>
  subroutine Satellites_Tidal_Fields_Null_Initialize(satellitesTidalFieldMethod,Satellites_Tidal_Field_Get)
    !% Initializes the ``Null'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in   ) :: satellitesTidalFieldMethod
    procedure(Satellites_Tidal_Field_Null_Get), pointer, intent(inout) :: Satellites_Tidal_Field_Get
    
    if (satellitesTidalFieldMethod == 'null') Satellites_Tidal_Field_Get => Satellites_Tidal_Field_Null_Get
    return
  end subroutine Satellites_Tidal_Fields_Null_Initialize

  double precision function Satellites_Tidal_Field_Null_Get(thisNode)
    !% Computes the tidal field acting on a satellite in the {\tt null} implementation. Always returns zero.
    use Galacticus_Nodes
    implicit none
    type (treeNode), intent(inout), pointer :: thisNode

    Satellites_Tidal_Field_Null_Get=0.0d0
    return
  end function Satellites_Tidal_Field_Null_Get
  
end module Satellites_Tidal_Fields_Null
