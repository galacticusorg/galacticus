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
     logical                                      :: doPing            =.false., reportWaitTime
   contains
     final     ::                 fcfsDestructor
     procedure :: forestNumber => fcfsForestNumber
  end type evolveForestsWorkShareFCFS

  interface evolveForestsWorkShareFCFS
     !!{
     Constructors for the \refClass{evolveForestsWorkShareFCFS} forest evolution work sharing class.
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
    Constructor for the \refClass{evolveForestsWorkShareFCFS} forest evolution work sharing class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (evolveForestsWorkShareFCFS)                              :: self
    type   (inputParameters           ), intent(inout)               :: parameters
    integer(c_size_t                  ), allocatable  , dimension(:) :: activeProcessRanks
    logical                                                          :: doPing            , reportWaitTime

    !![
    <inputParameter>
      <name>doPing</name>
      <defaultValue>.false.</defaultValue>
      <description>
        If true, the master MPI process will attach to the {\normalfont \ttfamily calculationReset} event and ping the MPI
        counter. This can help to ensure that the counter updates regularly.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportWaitTime</name>
      <defaultValue>.false.</defaultValue>
      <description>
        If true, the time spent waiting to increment the counter will be reported.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=evolveForestsWorkShareFCFS(doPing,reportWaitTime,activeProcessRanks)
    if (parameters%isPresent('activeProcessRanks')) then
       allocate(activeProcessRanks(parameters%count('activeProcessRanks')))
       !![
       <inputParameter>
         <name>activeProcessRanks</name>
         <description>A list of MPI process ranks which will be allowed to process trees---all other ranks are given no work.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       self=evolveForestsWorkShareFCFS(doPing,reportWaitTime,activeProcessRanks)
    else
       self=evolveForestsWorkShareFCFS(doPing,reportWaitTime                   )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fcfsConstructorParameters

  function fcfsConstructorInternal(doPing,reportWaitTime,activeProcessRanks) result(self)
    !!{
    Internal constructor for the \refClass{evolveForestsWorkShareFCFS} forest evolution work sharing class.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
    use :: Events_Hooks , only : calculationResetEvent, openMPThreadBindingAllLevels
#endif    
    implicit none
    type   (evolveForestsWorkShareFCFS)                                        :: self
    integer(c_size_t                  ), intent(in   ), dimension(:), optional :: activeProcessRanks
    logical                            , intent(in   )                         :: doPing            , reportWaitTime
#ifdef USEMPI
    integer                                                                    :: i
#endif
    !![
    <constructorAssign variables="doPing, reportWaitTime"/>
    !!]
    
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
    if (self%doPing .and. mpiSelf%isMaster()) call calculationResetEvent%attach(self,fcfsPing,openMPThreadBindingAllLevels,label='fcfsPing')
#endif
    return
  end function fcfsConstructorInternal

  subroutine fcfsDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily evolveForestsWorkShareFCFS} object.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
    use :: Events_Hooks , only : calculationResetEvent
#endif
    implicit none
    type(evolveForestsWorkShareFCFS), intent(inout) :: self

#ifdef USEMPI
    if (self%doPing.and.mpiSelf%isMaster().and.calculationResetEvent%isAttached(self,fcfsPing)) call calculationResetEvent%detach(self,fcfsPing)
#endif
    return
  end subroutine fcfsDestructor
  
  function fcfsForestNumber(self,utilizeOpenMPThreads)
    !!{
    Return the number of the next forest to process.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    use :: Display      , only : displayMessage
    use :: Timers       , only : timer
    implicit none
    integer(c_size_t                  )                :: fcfsForestNumber
    class  (evolveForestsWorkShareFCFS), intent(inout) :: self
    logical                            , intent(in   ) :: utilizeOpenMPThreads
    type   (timer                     )                :: timer_
    !$GLC attributes unused :: self, utilizeOpenMPThreads

#ifdef USEMPI
    if (any(self%activeProcessRanks == mpiSelf%rank())) then
#endif
       if (self%reportWaitTime) then
          timer_=timer()
          call timer_%start()
       end if
       fcfsForestNumber=forestCounter%increment()+1_c_size_t
       if (self%reportWaitTime) then
          call timer_%stop()
          call displayMessage('evolveForestsWorkShareFCFS wait time = '//trim(timer_%reportText()))
          end if
#ifdef USEMPI
    else
       fcfsForestNumber=huge(0_c_size_t)
    end if
#endif
    return
  end function fcfsForestNumber

  subroutine fcfsPing(self,node,uniqueID)
    !!{
    Return the number of the next forest to process.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities   , only : mpiSelf
#endif
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class  (*        ), intent(inout) :: self
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
#ifdef USEMPI
    integer(c_size_t )                :: forestNumber
#endif
    !$GLC attributes unused :: self, node, uniqueID

#ifdef USEMPI
    !$omp master
    if (mpiSelf%isMaster()) forestNumber=forestCounter%get()
    !$omp end master
#endif
    return
  end subroutine fcfsPing
