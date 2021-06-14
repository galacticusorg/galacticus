!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which handles destruction of default {\normalfont \ttfamily functionClass} objects.

module Galacticus_Function_Classes_Destroys
  !% Handles resetting of calculations before a new or updated node is processed.
  implicit none
  private
  public :: Galacticus_Function_Classes_Destroy

contains

  subroutine Galacticus_Function_Classes_Destroy()
    !% Calls any functions required to destroy all default {\normalfont \ttfamily functionClass} objects.
    !# <include directive="functionClassDestroyTask" type="moduleUse">
    include 'galacticus.function_class_destroy.tasks.modules.inc'
    !# </include>
    implicit none

    !# <include directive="functionClassDestroyTask" type="functionCall" functionType="void">
    include 'galacticus.function_class_destroy.tasks.inc'
    !# </include>
    return
  end subroutine Galacticus_Function_Classes_Destroy

end module Galacticus_Function_Classes_Destroys
