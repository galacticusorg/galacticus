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






module IO_HDF5
  use H5lt
  use HDF5
  use Galacticus_Error
  private
  public :: IO_HDF5_Initialize, IO_HDF5_Uninitialize

  logical :: hdf5IsInitalized=.false.
  integer :: initializationsCount=0

contains

  subroutine IO_HDF5_Initialize
    !% Initialize the HDF5 subsystem.
    implicit none
    integer :: errorCode

    if (.not.hdf5IsInitalized) then
       call h5open_f(errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('IO_HDF5_Initialize','failed to initialize HDF5 subsystem')
       hdf5IsInitalized=.true.
    end if
    initializationsCount=initializationsCount+1
    return
  end subroutine IO_HDF5_Initialize

  subroutine IO_HDF5_Uninitialize
    !% Uninitialize the HDF5 subsystem.
    implicit none
    integer :: errorCode

    if (hdf5IsInitalized) then
       initializationsCount=initializationsCount-1
       if (initializationsCount == 0) then
          call h5close_f(errorCode)
          if (errorCode < 0) call Galacticus_Error_Report('IO_HDF5_Uninitialize','failed to uninitialize HDF5 subsystem')
          hdf5IsInitalized=.false.
       end if
    end if
    return
  end subroutine IO_HDF5_Uninitialize

end module IO_HDF5
