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
Contains a module which handles setting of verbosity.
!!}

module Display_Verbosity
  !!{
  Handle setting of verbosity.
  !!}
  implicit none
  private
  public :: displayVerbositySetFromParameters

contains

  subroutine displayVerbositySetFromParameters(parameters)
    !!{
    Read the parameter that controls the verbosity level, and set that level.
    !!}
    use :: Display           , only : displayVerbositySet, enumerationVerbosityLevelEncode
    use :: ISO_Varying_String, only : char               , var_str                        , varying_string
    use :: Input_Parameters  , only : inputParameter     , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(varying_string )               :: verbosityLevel

    ! Get the verbosity level parameter.
    !![
    <inputParameter>
      <name>verbosityLevel</name>
      <defaultValue>var_str('standard')</defaultValue>
      <description>The level of verbosity for \glc\ (higher values give more verbosity).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call displayVerbositySet(enumerationVerbosityLevelEncode(char(verbosityLevel),includesPrefix=.false.))
    return
  end subroutine displayVerbositySetFromParameters

end module Display_Verbosity
