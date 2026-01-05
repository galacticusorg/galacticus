!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module that implements useful OpenMP utilities.
!!}

module OpenMP_Utilities
  !!{
  Implements useful OpenMP utilities.
  !!}
  private

contains

  !![
  <outputFileClose function="OpenMP_Critical_Wait_Times"/>
  !!]
  subroutine OpenMP_Critical_Wait_Times()
    !!{
    Outputs collected data on OpenMP critical section wait times.
    !!}
#ifdef OMPPROFILE
    use :: Output_HDF5          , only : outputFile
    use :: HDF5_Access          , only : hdf5Access
    use :: IO_HDF5              , only : hdf5Object
    use :: ISO_Varying_String   , only : varying_string         , var_str
    use :: OpenMP_Utilities_Data, only : criticalSectionWaitTime, criticalSectionCount
#endif
    implicit none
#ifdef OMPPROFILE
    type(hdf5Object) :: waitTimeGroup, waitTimeDataset, metaDataGroup
    include 'openMPCriticalSections.enumerate.inc'

    ! If no data was collected, simply return.
    if (all(criticalSectionWaitTime == 0.0d0)) return
    ! Open output group.
    !$ call hdf5Access%set()
    metaDataGroup=outputFile   %openGroup('metaData','Galacticus meta data.'           )
    waitTimeGroup=metaDataGroup%openGroup('openMP'  ,'Meta-data on OpenMP performance.')
    ! Write wait time data.
    call waitTimeGroup%writeDataset(criticalSectionNames   ,"criticalSectionNames"    ,"Names of OpenMP critical sections"                                                   )
    call waitTimeGroup%writeDataset(criticalSectionWaitTime,"criticalSectionWaitTimes","Total time spent waiting at OpenMP critical sections",datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(1.0d0,"unitsInSI")
    call waitTimeDataset%close()
    ! Close output groups.
    call waitTimeGroup%close()
    call metaDataGroup%close()
    !$ call hdf5Access%unset()
#endif
    return
  end subroutine OpenMP_Critical_Wait_Times

end module OpenMP_Utilities
