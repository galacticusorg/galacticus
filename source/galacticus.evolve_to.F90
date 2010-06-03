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










!% Contains a module which handles the current evolve to point in {\sc Galacticus}.

module Galacticus_Evolve_To_Module
  !% Handles the current evolve to point in {\sc Galacticus}.
  private
  public :: Galacticus_Evolve_to_Point_Next

contains

  logical function Galacticus_Evolve_to_Point_Next()
    !% Move to the next evolve to point in {\sc Galacticus}. Return false if no more evolve to points exist.
    implicit none
    
    Galacticus_Evolve_to_Point_Next=.false.
    return
  end function Galacticus_Evolve_to_Point_Next

end module Galacticus_Evolve_To_Module
