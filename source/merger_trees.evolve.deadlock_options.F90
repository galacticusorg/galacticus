!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which provides an enumeration for tree deadlock statuses.
!!}

module Merger_Trees_Evolve_Deadlock_Status
  !!{
  Provides an enumeration for tree deadlock statuses.
  !!}
  public

  !![
  <enumeration>
   <name>deadlockStatus</name>
   <description>Specifies deadlock status during merger tree evolution.</description>
   <entry label="isNotDeadlocked" />
   <entry label="isReporting"     />
   <entry label="isSuspendable"   />
   <entry label="isDeadlocked"    />
  </enumeration>
  !!]
end module Merger_Trees_Evolve_Deadlock_Status
