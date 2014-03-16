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

  !% Defines a \glc\ likelihood function.

  type, public, extends(likelihood) :: likelihoodGalacticus
     private
     type(varying_string) :: configFileName
   contains
     procedure :: evaluate => likelihoodGalacticusEvaluate
  end type likelihoodGalacticus
  
  interface likelihoodGalacticus
     !% Constructors for the \glc\ likelihood object.
     module procedure likelihoodGalacticusConstructor
  end interface likelihoodGalacticus
