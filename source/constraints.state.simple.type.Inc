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


  type, public, extends(state) :: stateSimple
     !% Defines a simple state class.
     private
     double precision, allocatable, dimension(:) :: current
     integer         , allocatable, dimension(:) :: accepted
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     ::             stateSimpleDestroy
     !# </workaround>      
     procedure :: get            => stateSimpleGet
     procedure :: update         => stateSimpleUpdate
     procedure :: mean           => stateSimpleMean
     procedure :: variance       => stateSimpleVariance
     procedure :: acceptanceRate => stateSimpleAcceptanceRate
     procedure :: reset          => stateSimpleReset
  end type stateSimple
  
  interface stateSimple
     !% Constructors for the ``simple'' state class.
     module procedure stateSimpleConstructor
  end interface stateSimple
