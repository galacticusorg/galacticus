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

!% Contains a module which provides options for halo accretion.

module Accretion_Halos_Options
  !% Provides options for halo accretion.
  public

  ! Options controlling whether hot, cold, or total accretion is required.
  integer                                         , parameter :: accretionModeTotal=0
  integer                                         , parameter :: accretionModeHot  =1
  integer                                         , parameter :: accretionModeCold =2

end module Accretion_Halos_Options
