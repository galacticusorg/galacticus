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
Contains a module which provides a class that implements general tasks to be performed by \glc.
!!}

module Tasks
  !!{
  Provides a class that implements general tasks to be performed by \glc.
  !!}
  private

  !![
  <functionClass>
   <name>task</name>
   <descriptiveName>Tasks</descriptiveName>
   <description>Class providing general top-level tasks to be performed by \glc\---the primary unit of
    computation that the code executes when run. Each task implementation defines a self-contained
    operation, such as evolving a forest of merger trees to produce a galaxy catalogue, running a
    Bayesian parameter estimation, performing N-body analysis, or executing a radiative transfer
    calculation. The \mono{perform} method carries out the task and optionally returns an exit status,
    while \mono{requiresOutputFile} indicates whether HDF5 output should be opened beforehand.</description>
   <default>evolveForests</default>
   <method name="perform" >
    <description>Perform the task.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>integer, intent(  out), optional :: status</argument>
   </method>
   <method name="requiresOutputFile" >
    <description>Should return true if the task requires the main output file to be open.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     taskRequiresOutputFile=.true.
    </code>
   </method>
  </functionClass>
  !!]

end module Tasks
