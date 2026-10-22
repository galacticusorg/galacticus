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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides utility functions for halo environments.
!!}

module Halo_Environment_Utilities
  implicit none
  private
  public :: Halo_Environment_Time

contains

  double precision function Halo_Environment_Time(cosmologyFunctions_,node,time,expansionFactor)
    !!{RST
    Return the cosmic time at which to evaluate a halo environment: the given ``time`` or ``expansionFactor`` if present,
    otherwise the cosmic time of the ``node``.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Error              , only : Error_Report
    use :: Galacticus_Nodes   , only : nodeComponentBasic     , treeNode
    implicit none
    class           (cosmologyFunctionsClass), intent(inout)           :: cosmologyFunctions_
    type            (treeNode               ), intent(inout), optional :: node
    double precision                         , intent(in   ), optional :: time               , expansionFactor
    class           (nodeComponentBasic     ), pointer                 :: basic

    if      (present(time           )) then
       Halo_Environment_Time =  time
    else if (present(expansionFactor)) then
       Halo_Environment_Time =  cosmologyFunctions_%cosmicTime(expansionFactor)
    else if (present(node           )) then
       basic                 => node %basic()
       Halo_Environment_Time =  basic%time ()
    else
       Halo_Environment_Time =  0.0d0
       call Error_Report('one of `time`, `expansionFactor`, or `node` must be provided'//{introspection:location})
    end if
    return
  end function Halo_Environment_Time

end module Halo_Environment_Utilities
