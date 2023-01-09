!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which handles setting of error wait times.
!!}

module Error_Wait
  !!{
  Handle setting of error wait times.
  !!}
  implicit none
  private
  public :: Error_Wait_Set_From_Parameters

contains

  subroutine Error_Wait_Set_From_Parameters(parameters)
    !!{
    Read the parameter that controls the verbosity level, and set that level.
    !!}
    use :: Error           , only : Error_Wait_Set
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer                                 :: errorWaitTime

    ! Get the verbosity level parameter.
    !![
    <inputParameter>
      <name>errorWaitTime</name>
      <defaultValue>86400</defaultValue>
      <description>The time, in seconds, for which \glc\ should sleep after a fatal error when running under MPI.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call Error_Wait_Set(errorWaitTime)
    return
  end subroutine Error_Wait_Set_From_Parameters

end module Error_Wait
