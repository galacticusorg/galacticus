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

module OpenMP_Utilities_Data
  !!{
  Implements data for useful OpenMP utilities.
  !!}
#ifdef OMPPROFILE
  use :: ISO_Varying_String, only : varying_string, var_str
#endif
  private

  ! Include auto-generated content describing number of OpenMP critical sections in the code.
#ifdef OMPPROFILE
  include 'openMPCriticalSections.count.inc'

  ! Define arrays which will be used to accumulate the time spent waiting for OpenMP critical sections.
  double precision, public, dimension(criticalSectionCount) :: criticalSectionWaitTime=0.0d0

  ! Variables used for storing timing intervals.
  double precision, public                                  :: ompProfileTimeWaitStart      , ompProfileTimeWaitEnd
  !$omp threadprivate(ompProfileTimeWaitStart,ompProfileTimeWaitEnd)
#endif

end module OpenMP_Utilities_Data
