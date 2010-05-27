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
