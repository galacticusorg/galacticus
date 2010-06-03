!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module that stores data that is shared between modules for the standard hot halo method.

module Tree_Node_Methods_Hot_Halo_Data
  !% Stores data that is shared between modules for the standard hot halo method.
  public

  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

end module Tree_Node_Methods_Hot_Halo_Data
