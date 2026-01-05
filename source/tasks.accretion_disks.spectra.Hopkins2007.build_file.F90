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

  !![
  <task name="taskAGNSpectraHopkins2008BuildFile">
   <description>A task which evolves galaxies within a set of merger tree forests.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskAGNSpectraHopkins2008BuildFile
     !!{
     Implementation of a task which builds a file containing a tabulation of AGN spectra from the model of \cite{hopkins_observational_2007}.
     !!}
     private
   contains
     procedure :: perform            => agnSpectraHopkins2008BuildFilePerform
     procedure :: requiresOutputFile => agnSpectraHopkins2008BuildFileRequiresOutputFile
  end type taskAGNSpectraHopkins2008BuildFile

  interface taskAGNSpectraHopkins2008BuildFile
     !!{
     Constructors for the \refClass{taskAGNSpectraHopkins2008BuildFile} task.
     !!}
     module procedure agnSpectraHopkins2008BuildFileParameters
  end interface taskAGNSpectraHopkins2008BuildFile

contains

  function agnSpectraHopkins2008BuildFileParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskAGNSpectraHopkins2008BuildFile} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskAGNSpectraHopkins2008BuildFile)                :: self
    type(inputParameters                   ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskAGNSpectraHopkins2008BuildFile()
    return
  end function agnSpectraHopkins2008BuildFileParameters

  subroutine agnSpectraHopkins2008BuildFilePerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraHopkins2007
    use :: Display               , only : displayIndent                  , displayUnindent
    use :: Error      , only : errorStatusSuccess
    implicit none
    class  (taskAGNSpectraHopkins2008BuildFile), intent(inout), target   :: self
    integer                                    , intent(  out), optional :: status
    type   (accretionDiskSpectraHopkins2007   )                          :: accretionDiskSpectra_
    !$GLC attributes unused :: self

    call displayIndent  ('Begin task: hopkins2007 AGN spectra file build')
    accretionDiskSpectra_=accretionDiskSpectraHopkins2007()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: hopkins2007 AGN spectra file build')
    return
  end subroutine agnSpectraHopkins2008BuildFilePerform

  logical function agnSpectraHopkins2008BuildFileRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskAGNSpectraHopkins2008BuildFile), intent(inout) :: self
    !$GLC attributes unused :: self

    agnSpectraHopkins2008BuildFileRequiresOutputFile=.false.
    return
  end function agnSpectraHopkins2008BuildFileRequiresOutputFile
