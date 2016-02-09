!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides enumerations for merger tree output.

module Galacticus_Output_Merger_Tree_Data
  !% Provides enumerations for merger tree output.
  implicit none
  public

  !@ <enumeration>
  !@  <name>nodeStatus</name>
  !@  <description>Used to specify the status (first, last, neither) of a node in a merger tree during output.</description>
  !@  <entry label="nodeStatusFirst"/>
  !@  <entry label="nodeStatusLast" />
  !@  <entry label="nodeStatusNull" />
  !@  <entry label="nodeStatusFinal"/>
  !@ </enumeration>
  integer, parameter, public :: nodeStatusNull = 0
  integer, parameter, public :: nodeStatusFirst=+1
  integer, parameter, public :: nodeStatusLast =-1
  integer, parameter, public :: nodeStatusFinal=-2
 
end module Galacticus_Output_Merger_Tree_Data
