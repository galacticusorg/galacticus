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


  type, public, extends(stoppingCriterion) :: stoppingCriterionCorrelationLengthCount
     !% Defines a stopping criterion class which stops after a given number of post-convergence correlation lengths.
     private
     integer :: stopAfterCount
   contains
     procedure :: stop => stoppingCriterionCorrelationLengthCountStop
  end type stoppingCriterionCorrelationLengthCount
  
  interface stoppingCriterionCorrelationLengthCount
     !% Constructors for the ``correlationLengthCount'' stoppingCriterion class.
     module procedure stoppingCriterionCorrelationLengthCountConstructor
  end interface stoppingCriterionCorrelationLengthCount
