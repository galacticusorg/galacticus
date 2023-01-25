!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

  use :: MPI_Utilities, only : mpiCounter

  !![
  <evolveForestsWorkShare name="evolveForestsWorkShareFCFS">
   <description>A forest evolution work sharing class in which forests are assigned on a first-come-first-served basis.</description>
  </evolveForestsWorkShare>
  !!]
  type, extends(evolveForestsWorkShareClass) :: evolveForestsWorkShareFCFS
     !!{
     Implementation of a forest evolution work sharing class in which forests are assigned on a first-come-first-served basis.
     !!}
     private
     integer(c_size_t), allocatable, dimension(:) :: activeProcessRanks
   contains
     procedure :: forestNumber => fcfsForestNumber
  end type evolveForestsWorkShareFCFS

  interface evolveForestsWorkShareFCFS
     !!{
     Constructors for the {\normalfont \ttfamily fcfs} forest evolution work sharing class.
     !!}
     module procedure fcfsConstructorParameters
     module procedure fcfsConstructorInternal
  end interface evolveForestsWorkShareFCFS

  ! Global counter of forests assigned.
  type   (mpiCounter) :: forestCounter
  logical             :: forestCounterInitialized=.false.

contains

  function fcfsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily fcfs} forest evolution work sharing class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (evolveForestsWorkShareFCFS)                              :: self
    type   (inputParameters           ), intent(inout)               :: parameters
    integer(c_size_t                  ), allocatable  , dimension(:) :: activeProcessRanks

    if (parameters%isPresent('activeProcessRanks')) then
       allocate(activeProcessRanks(parameters%count('activeProcessRanks')))
       !![
       <inputParameter>
         <name>activeProcessRanks</name>
         <description>A list of MPI process ranks which will be allowed to process trees---all other ranks are given no work.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       self=evolveForestsWorkShareFCFS(activeProcessRanks)
    else
       self=evolveForestsWorkShareFCFS(                  )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fcfsConstructorParameters

  function fcfsConstructorInternal(activeProcessRanks) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily fcfs} forest evolution work sharing class.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    type   (evolveForestsWorkShareFCFS)                                        :: self
    integer(c_size_t                  ), intent(in   ), dimension(:), optional :: activeProcessRanks
#ifdef USEMPI
    integer                                                                    :: i
#endif
    
    if (.not.forestCounterInitialized) then
       forestCounter=mpiCounter()
       forestCounterInitialized=.true.
    end if
    call forestCounter%reset()
#ifdef USEMPI
    if (present(activeProcessRanks)) then
       allocate(self%activeProcessRanks(size(activeProcessRanks)))
       self%activeProcessRanks=activeProcessRanks
    else
       allocate(self%activeProcessRanks(0:mpiSelf%count()-1))
       do i=0,mpiSelf%count()-1
          self%activeProcessRanks(i)=i
       end do
    end if
#endif
    return
  end function fcfsConstructorInternal

  function fcfsForestNumber(self,utilizeOpenMPThreads)
    !!{
    Return the number of the next forest to process.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    integer(c_size_t                  )                :: fcfsForestNumber
    class  (evolveForestsWorkShareFCFS), intent(inout) :: self
    logical                            , intent(in   ) :: utilizeOpenMPThreads
    !$GLC attributes unused :: self, utilizeOpenMPThreads

#ifdef USEMPI
    if (any(self%activeProcessRanks == mpiSelf%rank())) then
#endif
       fcfsForestNumber=forestCounter%increment()+1_c_size_t
#ifdef USEMPI
    else
       fcfsForestNumber=huge(0_c_size_t)
    end if
#endif
    return
  end function fcfsForestNumber
