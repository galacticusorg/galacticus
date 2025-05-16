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
  Implements a merger tree evolve profiler that warns about tiny steps.
  !!}
  
  !![
  <mergerTreeEvolveProfiler name="mergerTreeEvolveProfilerTinySteps">
   <description>
    A merger tree evolve profiler that warns about tiny steps.
   </description>
  </mergerTreeEvolveProfiler>
  !!]
  type, extends(mergerTreeEvolveProfilerClass) :: mergerTreeEvolveProfilerTinySteps
     !!{
     A merger tree evolve profiler that warns about tiny steps.
     !!}
     private
     double precision                 :: timeStepTiny   , timeStepMinimum
     type            (varying_string) :: stepDescriptor_
   contains
     procedure :: stepDescriptor => tinyStepsStepDescriptor
     procedure :: profile        => tinyStepsProfile
  end type mergerTreeEvolveProfilerTinySteps

  interface mergerTreeEvolveProfilerTinySteps
     !!{
     Constructors for the \refClass{mergerTreeEvolveProfilerTinySteps} merger tree evolve profiler class.
     !!}
     module procedure tinyStepsConstructorParameters
     module procedure tinyStepsConstructorInternal
  end interface mergerTreeEvolveProfilerTinySteps

contains
  
  function tinyStepsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveProfilerTinySteps} merger tree evolve profiler class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeEvolveProfilerTinySteps)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: timeStepTiny

    !![
    <inputParameter>
      <name>timeStepTiny</name>
      <defaultValue>1.0d-9</defaultValue>
      <description>The time step below which warnings will be issued.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeEvolveProfilerTinySteps(timeStepTiny)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tinyStepsConstructorParameters

  function tinyStepsConstructorInternal(timeStepTiny) result(self)
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type            (mergerTreeEvolveProfilerTinySteps)                :: self
    double precision                                   , intent(in   ) :: timeStepTiny
    !![
    <constructorAssign variables="timeStepTiny"/>
    !!]

    self%timeStepMinimum=huge(0.0d0)
    self%stepDescriptor_="unknown"
    return
  end function tinyStepsConstructorInternal

  subroutine tinyStepsStepDescriptor(self,descriptor)
    !!{
    Set the descriptor for the current step.
    !!}
    implicit none
    class(mergerTreeEvolveProfilerTinySteps), intent(inout) :: self
    type (varying_string                   ), intent(in   ) :: descriptor

    self%stepDescriptor_=descriptor
    return
  end subroutine tinyStepsStepDescriptor
  
  subroutine tinyStepsProfile(self,node,time,timeStart,timeEnd,timestep,countEvaluations,interrupted,propertyIndex,propertyName,propertyValue,propertyRate,propertyScale,propertyError,timeCPU)
    !!{
    Profile the differential evolution step.
    !!}
    use :: Display           , only : displayIndent, displayUnindent, displayMessage, displayMagenta, &
         &                            displayReset , displayRed
    use :: ISO_Varying_String, only : var_str      , assignment(=)  , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (mergerTreeEvolveProfilerTinySteps), intent(inout)               :: self
    type            (treeNode                         ), intent(in   )               :: node
    double precision                                   , intent(in   )               :: time                   , timeStep     , &
         &                                                                              timeStart              , timeEnd      , &
         &                                                                              timeCPU
    integer         (c_size_t                         ), intent(in   )               :: countEvaluations
    logical                                            , intent(in   )               :: interrupted
    integer         (c_size_t                         ), intent(in   )               :: propertyIndex
    type            (varying_string                   ), intent(in   )               :: propertyName
    double precision                                   , intent(in   ), dimension(:) :: propertyValue          , propertyScale, &
         &                                                                              propertyError          , propertyRate
    double precision                                   , parameter                   :: tolerance       =1.0d-2
    character       (len=24                           )                              :: label
    type            (varying_string                   )                              :: report

    ! Return immediately if the timestep is not tiny.
    if     (                                             &
         &   timeStep               >= self%timeStepTiny &
         & ) return
    ! Check for cases where the timestep is small, but was likely set to hit the end time precisely.
    if     (                                              &
         &   timeEnd     -timeStart >= self%timeStepTiny  &
         &  .and.                                         &
         &   timeEnd-time-timeStep  <  tolerance*timeStep &
         & ) return
    ! Timestep is tiny - issue a report.
    write (label,'(e12.6)') timeStep
    report=displayMagenta()//'Warning:'//displayReset()//' tiny timestep of '//trim(adjustl(label))//' Gyr taken'
    call displayIndent(report)
    report=var_str('     tree index = ')//node%hostTree%index
    call displayMessage(report)
    report=var_str('     node index = ')//node         %index()
    call displayMessage(report)
    write (label,'(e12.6)') time
    report=        '           time = '//trim(adjustl(label))//' Gyr'
    call displayMessage(report)
    call displayMessage('step descriptor = "'//self%stepDescriptor_//'"')
    if (timeStep < self%timeStepMinimum) then
       call displayMessage(displayRed()//'this is the smallest timestep seen so far'//displayReset())
       self%timeStepMinimum=timeStep
    end if
    ! Determine cause of the tiny timestep.
    if (timeEnd-timeStart < self%timeStepTiny) then
       ! Timestep was limited by the end time.
       report='Cause: imposed evolution interval'
       call displayIndent(report)
       write (label,'(e24.16)')         timeStart
       report='t₀ = '//label//' Gyr'
       call displayMessage(report)
       write (label,'(e24.16)') timeEnd
       report='t₁ = '//label//' Gyr'
       call displayMessage(report)
       write (label,'(e24.16)') timeEnd-timeStart
       report='Δt = '//label//' Gyr'
       call displayMessage(report)
       call displayUnindent('done')
    else
       ! Timestep was limited by the ODE solver.
       report='Cause: ODE solver'
       call displayIndent(report)
       report=var_str('name of limiting property    = ')//propertyName
       call displayMessage(report)
       report=var_str('number of failed evaluations = ')//countEvaluations
       call displayMessage(report)
       if (interrupted) then
          report='evolution was interrupted'
       else
          report='evolution was not interrupted'
       end if
       call displayMessage(report)
       report=var_str('value / rate / scale / error = ')
       write (label,'(e24.16)') propertyValue(propertyIndex)
       report=report//trim(label)//' / '
       write (label,'(e24.16)') propertyRate (propertyIndex)
       report=report//trim(label)//' / '
       write (label,'(e24.16)') propertyScale(propertyIndex)
       report=report//trim(label)//' / '
       write (label,'(e24.16)') propertyError(propertyIndex)
       report=report//trim(label)
       call displayMessage(report)
       call displayUnindent('done')
    end if
    call displayUnindent('done')
    return
  end subroutine tinyStepsProfile
