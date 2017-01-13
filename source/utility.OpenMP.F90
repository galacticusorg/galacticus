!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module that implements useful OpenMP utilities.

module OpenMP_Utilities
  !% Implements useful OpenMP utilities.
  private
  public :: OpenMP_Critical_Wait_Times

contains
  
  !# <hdfPreCloseTask>
  !#  <unitName>OpenMP_Critical_Wait_Times</unitName>
  !# </hdfPreCloseTask>
  subroutine OpenMP_Critical_Wait_Times()
    !% Outputs collected data on OpenMP critical section wait times.
    use OpenMP_Utilities_Data
    use Galacticus_HDF5
    use ISO_Varying_String
    use IO_HDF5
    implicit none
#ifdef OMPPROFILE
    type(hdf5Object) :: waitTimeGroup, waitTimeDataset, metaDataGroup
    include 'openMPCriticalSections.enumerate.inc'
    
    ! If no data was collected, simply return.
    if (all(criticalSectionWaitTime == 0.0d0)) return
    ! Open output group.
    !$omp critical(HDF5_Access)
    metaDataGroup=galacticusOutputFile%openGroup('metaData','Galacticus meta data.'           )
    waitTimeGroup=metaDataGroup       %openGroup('openMP'  ,'Meta-data on OpenMP performance.')
    ! Write wait time data.
    call waitTimeGroup%writeDataset(criticalSectionNames   ,"criticalSectionNames"    ,"Names of OpenMP critical sections"                                                   )
    call waitTimeGroup%writeDataset(criticalSectionWaitTime,"criticalSectionWaitTimes","Total time spent waiting at OpenMP critical sections",datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(1.0d0,"unitsInSI")
    call waitTimeDataset%close()
    ! Close output groups.
    call waitTimeGroup%close()
    call metaDataGroup%close()
    !$omp end critical(HDF5_Access)
#endif
    return
  end subroutine OpenMP_Critical_Wait_Times

end module OpenMP_Utilities
