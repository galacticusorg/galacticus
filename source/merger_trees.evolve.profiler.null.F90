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

!% Implements a merger tree evolve profiler that does nothing.

  !# <mergerTreeEvolveProfiler name="mergerTreeEvolveProfilerNull">
  !#  <description>A merger tree evolve profiler that does nothing.</description>
  !# </mergerTreeEvolveProfiler>
  type, extends(mergerTreeEvolveProfilerClass) :: mergerTreeEvolveProfilerNull
     !% A merger tree evolve profiler that does nothing.
     private
   contains
     procedure :: profile => nullProfile
  end type mergerTreeEvolveProfilerNull

  interface mergerTreeEvolveProfilerNull
     !% Constructors for the {\normalfont \ttfamily null} merger tree evolve profiler class.
     module procedure nullConstructorParameters
  end interface mergerTreeEvolveProfilerNull

contains

  function nullConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily null} merger tree evolve profiler class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeEvolveProfilerNull)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters
    
    self=mergerTreeEvolveProfilerNull()
    return
  end function nullConstructorParameters

  subroutine nullProfile(self,timestep,propertyName)
    !% Profile the differential evolution step.
    implicit none
    class           (mergerTreeEvolveProfilerNull), intent(inout) :: self
    double precision                              , intent(in   ) :: timeStep
    type            (varying_string              ), intent(in   ) :: propertyName
    !$GLC attributes unused :: self, timeStep, propertyName

    return
  end subroutine nullProfile
