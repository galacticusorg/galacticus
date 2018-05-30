!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which handles setting of verbosity.

module Galacticus_Display_Verbosity
  !% Handle setting of verbosity.
  implicit none
  private
  public :: Galacticus_Verbosity_Set_From_Parameters

contains

  subroutine Galacticus_Verbosity_Set_From_Parameters()
    !% Read the parameter that controls the verbosity level, and set that level.
    use Input_Parameters
    use Galacticus_Display
    implicit none
    integer :: verbosityLevel

    ! Get the verbosity level parameter.
    !# <inputParameter>
    !#   <name>verbosityLevel</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1</defaultValue>
    !#   <description>The level of verbosity for \glc\ (higher values give more verbosity).</description>
    !#   <source>globalParameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    call Galacticus_Verbosity_Level_Set(verbosityLevel)
    return
  end subroutine Galacticus_Verbosity_Set_From_Parameters
  
end module Galacticus_Display_Verbosity
