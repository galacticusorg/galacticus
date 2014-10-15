!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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


  type, public, extends(stateHistory) :: stateCorrelation
     !% Defines a ``correlation'' state class.
     private
     integer                                       :: storedStateCount               , stepComputePrevious, &
          &                                           convergedCorrelationLengthCount
     double precision, allocatable, dimension(:,:) :: states
     integer         , allocatable, dimension(:  ) :: correlationLengths
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     ::             stateCorrelationDestroy
     !# </workaround>
     procedure :: update                          => stateCorrelationUpdate
     procedure :: reset                           => stateCorrelationReset
     procedure :: postConvergenceCorrelationCount => stateCorrelationPostConvergenceCorrelationCount
     procedure :: correlationLength              => stateCorrelationCorrelationLength
  end type stateCorrelation
  
  interface stateCorrelation
     !% Constructors for the ``correlation'' state class.
     module procedure stateCorrelationConstructor
  end interface stateCorrelation
