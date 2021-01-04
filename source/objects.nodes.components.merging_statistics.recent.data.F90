!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which stores data for the recent merging statistics component.

module Node_Component_Merging_Statistics_Recent_Data
  !% Stores data for the recent merging statistics component.
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Recent_Count

  integer(c_size_t), public :: mergingStatisticsRecentCount

contains

  function Node_Component_Merging_Statistics_Recent_Count()
    !% Return the size of the merging statistics property
    implicit none
    integer(c_size_t) :: Node_Component_Merging_Statistics_Recent_Count

    Node_Component_Merging_Statistics_Recent_Count=mergingStatisticsRecentCount
    return
  end function Node_Component_Merging_Statistics_Recent_Count

end module Node_Component_Merging_Statistics_Recent_Data
