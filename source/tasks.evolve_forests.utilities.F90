!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module of globally-accessible functions supporting the {\normalfont \ttfamily evolveForests} task class.
!!}

module Tasks_Evolve_Forests_Utilities
  !!{
  Provides globally-accessible functions supporting the {\normalfont \ttfamily evolveForests} task class.
  !!}
  private
  public :: Tasks_Evolve_Forest_Construct, Tasks_Evolve_Forest_Perform, &
       &    Tasks_Evolve_Forest_Destruct

  ! Module-scope pointer to our task object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: task__
  !$omp threadprivate(task__)

contains

  !![
  <functionGlobal>
   <unitName>Tasks_Evolve_Forest_Construct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout)          :: parameters</arguments>
   <arguments>class(*              ), intent(  out), pointer :: task_</arguments>
  </functionGlobal>
  !!]
  subroutine Tasks_Evolve_Forest_Construct(parameters,task_)
    !!{
    Build a {\normalfont \ttfamily taskEvolveForests} object from a given parameter set. This is a globally-callable function
    to allow us to subvert the class/module hierarchy.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Tasks           , only : task          , taskEvolveForests
    implicit none
    type (inputParameters), intent(inout)          :: parameters
    class(*              ), intent(  out), pointer :: task_

    task__ => task(parameters)
    select type (task__)
    class is (taskEvolveForests)
       !![
       <referenceCountIncrement object="task__"/>
       !!]
       call task__%autoHook()
    class default
       call Error_Report('task must be of the "taskEvolveForests" class'//{introspection:location})
    end select
    task_ => task__
    return
  end subroutine Tasks_Evolve_Forest_Construct

  !![
  <functionGlobal>
   <unitName>Tasks_Evolve_Forest_Perform</unitName>
   <type>void</type>
   <arguments>class  (*), intent(inout)           :: task_</arguments>
   <arguments>integer   , intent(  out), optional :: status</arguments>
  </functionGlobal>
  !!]
  subroutine Tasks_Evolve_Forest_Perform(task_,status)
    !!{
    Perform the task for a {\normalfont \ttfamily taskEvolveForests} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error, only : Error_Report
    use :: Tasks, only : task        , taskEvolveForests
    implicit none
    class  (*), intent(inout)           :: task_
    integer   , intent(  out), optional :: status

    select type (task_)
    class is (taskEvolveForests)
       call task_%perform(status)
    class default
       call Error_Report('task must be of the "taskEvolveForests" class'//{introspection:location})
    end select
    return
  end subroutine Tasks_Evolve_Forest_Perform

  !![
  <functionGlobal>
   <unitName>Tasks_Evolve_Forest_Destruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: task_</arguments>
  </functionGlobal>
  !!]
  subroutine Tasks_Evolve_Forest_Destruct(task_)
    !!{
    Destruct a {\normalfont \ttfamily taskEvolveForests} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error, only : Error_Report
    use :: Tasks, only : task        , taskEvolveForests
    implicit none
    class(*), intent(inout), pointer :: task_

    task__ => task_
    select type (task__)
    class is (taskEvolveForests)
       !![
       <objectDestructor name="task__"/>
       !!]
    class default
       call Error_Report('task must be of the "taskEvolveForests" class'//{introspection:location})
    end select
    return
  end subroutine Tasks_Evolve_Forest_Destruct

end module Tasks_Evolve_Forests_Utilities
