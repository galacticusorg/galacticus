!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a merger tree evolve profiler that does nothing.
!!}

  !![
  <mergerTreeEvolveProfiler name="mergerTreeEvolveProfilerNull">
   <description>A merger tree evolve profiler that does nothing.</description>
  </mergerTreeEvolveProfiler>
  !!]
  type, extends(mergerTreeEvolveProfilerClass) :: mergerTreeEvolveProfilerNull
     !!{
     A merger tree evolve profiler that does nothing.
     !!}
     private
   contains
     procedure :: stepDescriptor => nullStepDescriptor
     procedure :: profile        => nullProfile
  end type mergerTreeEvolveProfilerNull

  interface mergerTreeEvolveProfilerNull
     !!{
     Constructors for the \refClass{mergerTreeEvolveProfilerNull} merger tree evolve profiler class.
     !!}
     module procedure nullConstructorParameters
  end interface mergerTreeEvolveProfilerNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveProfilerNull} merger tree evolve profiler class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeEvolveProfilerNull)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    
    self=mergerTreeEvolveProfilerNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  subroutine nullStepDescriptor(self,descriptor)
    !!{
    Set the descriptor for the current step.
    !!}
    implicit none
    class(mergerTreeEvolveProfilerNull), intent(inout) :: self
    type (varying_string              ), intent(in   ) :: descriptor
    !$GLC attributes unused :: self, descriptor
    
    return
  end subroutine nullStepDescriptor
  
  subroutine nullProfile(self,node,time,timeStart,timeEnd,timestep,countEvaluations,interrupted,propertyIndex,propertyName,propertyValue,propertyRate,propertyScale,propertyError,timeCPU)
    !!{
    Profile the differential evolution step.
    !!}
    implicit none
    class           (mergerTreeEvolveProfilerNull), intent(inout)               :: self
    type            (treeNode                    ), intent(in   )               :: node
    double precision                              , intent(in   )               :: time            , timeStep     , &
         &                                                                         timeStart       , timeEnd      , &
         &                                                                         timeCPU
    integer         (c_size_t                    ), intent(in   )               :: countEvaluations
    logical                                       , intent(in   )               :: interrupted
    integer         (c_size_t                    ), intent(in   )               :: propertyIndex
    type            (varying_string              ), intent(in   )               :: propertyName
    double precision                              , intent(in   ), dimension(:) :: propertyValue   , propertyScale, &
         &                                                                         propertyError   , propertyRate
    !$GLC attributes unused :: self, node, time, timeStart, timeEnd, timeStep, countEvaluations, interrupted, propertyIndex, propertyName, propertyValue, propertyRate, propertyScale, propertyError, timeCPU

    return
  end subroutine nullProfile
