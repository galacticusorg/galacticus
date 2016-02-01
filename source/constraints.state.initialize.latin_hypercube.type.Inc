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

  !% Defines a state initializor that draws from the parameter priors using a Latin hypercube sampling.

  type, public, extends(stateInitializor) :: stateInitializorLatinHypercube
     private
     integer :: maximinTrialCount
   contains
     procedure :: initialize => stateInitializorLatinHypercubeInitialize
  end type stateInitializorLatinHypercube
  
  interface stateInitializorLatinHypercube
     !% Constructors for the prior random state initializor class.
     module procedure stateInitializorLatinHypercubeConstructor
  end interface stateInitializorLatinHypercube
