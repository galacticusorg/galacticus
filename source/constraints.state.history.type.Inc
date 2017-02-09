!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  type, public, extends(stateSimple) :: stateHistory
     !% Defines a ``history'' state class.
     private
     double precision, allocatable, dimension(:) :: stateSum, stateSquaredSum
   contains
     final     ::             stateHistoryDestroy
     procedure :: update   => stateHistoryUpdate
     procedure :: mean     => stateHistoryMean
     procedure :: variance => stateHistoryVariance
     procedure :: reset    => stateHistoryReset
     procedure :: restore  => stateHistoryRestore
  end type stateHistory
  
  interface stateHistory
     !% Constructors for the ``history'' state class.
     module procedure stateHistoryConstructor
  end interface stateHistory
