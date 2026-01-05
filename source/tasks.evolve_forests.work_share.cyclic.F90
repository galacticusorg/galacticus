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

  !$ use :: OMP_Lib, only : omp_lock_kind

  !![
  <evolveForestsWorkShare name="evolveForestsWorkShareCyclic">
   <description>A forest evolution work sharing class in which forests are assigned by cycling through processes.</description>
  </evolveForestsWorkShare>
  !!]
  type, extends(evolveForestsWorkShareClass) :: evolveForestsWorkShareCyclic
     !!{
     Implementation of a forest evolution work sharing class in which forests are assigned by cycling through processes.
     !!}
     private
     !$ integer(omp_lock_kind)                            :: lock
     integer   (c_size_t     ), allocatable, dimension(:) :: treeNumber_
     logical                                              :: utilizeOpenMPThreads, first
   contains
     final     ::                 cyclicDestructor
     procedure :: forestNumber => cyclicForestNumber
  end type evolveForestsWorkShareCyclic

  interface evolveForestsWorkShareCyclic
     !!{
     Constructors for the \refClass{evolveForestsWorkShareCyclic} forest evolution work sharing class.
     !!}
     module procedure cyclicConstructorParameters
     module procedure cyclicConstructorInternal
  end interface evolveForestsWorkShareCyclic

contains

  function cyclicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{evolveForestsWorkShareCyclic} forest evolution work sharing class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(evolveForestsWorkShareCyclic)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=evolveForestsWorkShareCyclic()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cyclicConstructorParameters

  function cyclicConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{evolveForestsWorkShareCyclic} forest evolution work sharing class.
    !!}
    implicit none
    type(evolveForestsWorkShareCyclic) :: self

    self%first               =.true.
    self%utilizeOpenMPThreads=.true.
    !$ call OMP_Init_Lock(self%lock)
    return
  end function cyclicConstructorInternal

  subroutine cyclicDestructor(self)
    !!{
    Destructor for the \refClass{evolveForestsWorkShareCyclic} forest evolution work sharing class.
    !!}
    implicit none
    type(evolveForestsWorkShareCyclic), intent(inout) :: self

    !$ call OMP_Destroy_Lock(self%lock)
    return
  end subroutine cyclicDestructor

  function cyclicForestNumber(self,utilizeOpenMPThreads)
    !!{
    Return the number of the next forest to process.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer(c_size_t                    )                :: cyclicForestNumber
    class  (evolveForestsWorkShareCyclic), intent(inout) :: self
    logical                              , intent(in   ) :: utilizeOpenMPThreads
    integer(c_size_t                    )                :: i

    !$ call OMP_Set_Lock(self%lock)
    if (self%first) then
       self%utilizeOpenMPThreads=utilizeOpenMPThreads
       self%first               =.false.
       allocate(self%treeNumber_(0:self%workerCount(utilizeOpenMPThreads)-1))
       do i=0,self%workerCount(utilizeOpenMPThreads)-1
          self%treeNumber_(i)=+i                                      &
               &              -self%workerCount(utilizeOpenMPThreads) &
               &              +1_c_size_t
       end do
    else
       if (self%utilizeOpenMPThreads .neqv. utilizeOpenMPThreads) call Error_Report('"cyclic" work share can not support transitions between utilizing/not utilizing OpenMP threads'//{introspection:location})
    end if
    !$ call OMP_Unset_Lock(self%lock)
    self%treeNumber_(self%workerID(utilizeOpenMPThreads))=+self%treeNumber_(self%workerID   (utilizeOpenMPThreads)) &
         &                                                +                 self%workerCount(utilizeOpenMPThreads)
    cyclicForestNumber                                   =+self%treeNumber_(self%workerID   (utilizeOpenMPThreads))
    return
  end function cyclicForestNumber
