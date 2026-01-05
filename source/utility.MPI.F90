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
Contains a module that implements useful MPI utilities.
!!}

module MPI_Utilities
  !!{
  Implements useful MPI utilities.
  !!}
#ifdef USEMPI
  use               :: MPI_F08           , only : MPI_Win        , MPI_Datatype, MPI_Comm
#endif
  !$ use            :: Locks             , only : ompLock
  use   , intrinsic :: ISO_C_Binding     , only : c_size_t       , c_ptr
  use               :: ISO_Varying_String, only : varying_string
  use               :: Resource_Manager  , only : resourceManager
  private
  public :: mpiInitialize, mpiFinalize, mpiBarrier, mpiSelf, mpiCounter
 
  ! Define a type for interacting with MPI.
  type :: mpiObject
     private
     integer                                            :: rankValue        , countValue      , &
          &                                                nodeCountValue   , rankOnNodeValue , &
          &                                                iCommunicator    , countOnNodeValue
#ifdef USEMPI
     type   (MPI_Comm      )                            :: communicator
     type   (MPI_Comm      ), allocatable, dimension(:) :: communicatorStack
#endif
     type   (varying_string)                            :: hostName
     integer                , allocatable, dimension(:) :: allRanks         , nodeAffinities
   contains
     !![
     <methods>
       <method description="Return true if this is the master process (i.e. rank-0 process)." method="isMaster" />
       <method description="Return true if MPI is active." method="isActive" />
       <method description="Return the rank of this process." method="rank" />
       <method description="Return the total number of processes." method="count" />
       <method description="Return the rank of this process on its node." method="rankOnNode" />
       <method description="Return the total number of processes on this processes' node." method="countOnNode" />
       <method description="Return a label containing the rank of the process." method="rankLabel" />
       <method description="Return the number of nodes on which this MPI job is running." method="nodeCount" />
       <method description="Return the index of the node on which the MPI process of the given rank (or this process if no rank is given) is running." method="nodeAffinity" />
       <method description="Return the name of the host on which this MPI process is running." method="hostAffinity" />
       <method description="Request the content of {\normalfont \ttfamily array} from each processes listed in {\normalfont \ttfamily requestFrom}." method="requestData" />
       <method description="Broadcast the content of {\normalfont \ttfamily array} from the {\normalfont \ttfamily sendFrom} processes to all other processes." method="broadcastData" />
       <method description="Return true if a message is waiting, optionally from the specified process and with the specified tag." method="messageWaiting" />
       <method description="Return the average of {\normalfont \ttfamily array} over all processes." method="average" />
       <method description="Return the median of {\normalfont \ttfamily array} over all processes." method="median" />
       <method description="Return the sum of {\normalfont \ttfamily array} over all processes." method="sum" />
       <method description="Return the maximum value of {\normalfont \ttfamily array} over all processes." method="maxval" />
       <method description="Return the rank of the process with the maximum value of {\normalfont \ttfamily array} over all processes." method="maxloc" />
       <method description="Return the minimum value of {\normalfont \ttfamily array} over all processes." method="minval" />
       <method description="Return true if any of {\normalfont \ttfamily scalar} is true over all processes." method="any" />
       <method description="Return true if every {\normalfont \ttfamily scalar} is true over all processes." method="all" />
       <method description="Return the rank of the process with the minimum value of {\normalfont \ttfamily array} over all processes." method="minloc" />
       <method description="Gather arrays from all processes into an array of rank one higher." method="gather" />
       <method description="Create a new communicator and push onto the stack." method="communicatorPush" />
       <method description="Pop a communicator off of the stack, restoring the previous communicator." method="communicatorPop" />
       <method description="Initialize the state of the MPI group." method="stateInitialize" />
     </methods>
     !!]
     procedure :: isMaster         => mpiIsMaster
     procedure :: isActive         => mpiIsActive
     procedure :: rank             => mpiGetRank
     procedure :: rankOnNode       => mpiGetRankOnNode
     procedure :: rankLabel        => mpiGetRankLabel
     procedure :: count            => mpiGetCount
     procedure :: countOnNode      => mpiGetCountOnNode
     procedure :: nodeCount        => mpiGetNodeCount
     procedure :: nodeAffinity     => mpiGetNodeAffinity
     procedure :: hostAffinity     => mpiGetHostAffinity
     procedure ::                     mpiRequestData1D           , mpiRequestData2D       , &
          &                           mpiRequestDataInt1D        , mpiRequestDataLogical1D
     generic   :: requestData      => mpiRequestData1D           , mpiRequestData2D       , &
          &                           mpiRequestDataInt1D        , mpiRequestDataLogical1D
     procedure ::                     mpiBroadcastData1D         , mpiBroadcastData2D     , &
          &                           mpiBroadcastData3D         , mpiBroadcastDataScalar , &
          &                           mpiBroadcastDataSizeTScalar
     generic   :: broadcastData    => mpiBroadcastData1D         , mpiBroadcastData2D     , &
          &                           mpiBroadcastData3D         , mpiBroadcastDataScalar , &
          &                           mpiBroadcastDataSizeTScalar
     procedure :: messageWaiting   => mpiMessageWaiting
     procedure ::                     mpiAverageScalar           , mpiAverageArray
     generic   :: average          => mpiAverageScalar           , mpiAverageArray
     procedure ::                     mpiMedianArray
     generic   :: median           => mpiMedianArray
     procedure ::                     mpiSumScalarInt            , mpiSumArrayInt
     procedure ::                     mpiSumScalarSizeT          , mpiSumArraySizeT       , &
          &                           mpiSumArrayTwoSizeT        , mpiSumArrayThreeSizeT
     procedure ::                     mpiSumScalarDouble         , mpiSumArrayDouble      , &
          &                           mpiSumArrayTwoDouble       , mpiSumArrayThreeDouble
     procedure ::                     mpiSumArrayReal
     generic   :: sum              => mpiSumScalarInt            , mpiSumArrayInt         , &
          &                           mpiSumScalarDouble         , mpiSumArrayDouble      , &
          &                           mpiSumArrayTwoDouble       , mpiSumArrayThreeDouble , &
          &                           mpiSumScalarSizeT          , mpiSumArraySizeT       , &
          &                           mpiSumArrayTwoSizeT        , mpiSumArrayThreeSizeT  , &
          &                           mpiSumArrayReal
     procedure ::                     mpiAnyLogicalScalar
     generic   :: any              => mpiAnyLogicalScalar
     procedure ::                     mpiAllLogicalScalar
     generic   :: all              => mpiAllLogicalScalar
     procedure :: maxloc           => mpiMaxloc
     procedure ::                     mpiMaxvalScalar            , mpiMaxvalArray         , &
          &                           mpiMaxvalScalarSizeT       , mpiMaxvalArraySizeT
     generic   :: maxval           => mpiMaxvalScalar            , mpiMaxvalArray         , &
          &                           mpiMaxvalScalarSizeT       , mpiMaxvalArraySizeT
     procedure :: minloc           => mpiMinloc
     procedure ::                     mpiMinvalScalar            , mpiMinvalArray         , &
          &                           mpiMinValIntScalar         , mpiMinvalIntArray
     generic   :: minval           => mpiMinvalScalar            , mpiMinvalArray         , &
          &                           mpiMinValIntScalar         , mpiMinvalIntArray
     procedure ::                     mpiGather1D                , mpiGather2D            , &
          &                           mpiGatherScalar            , mpiGatherInt1D         , &
          &                           mpiGatherIntScalar         , mpiGatherLogicalScalar
     generic   :: gather           => mpiGather1D                , mpiGather2D            , &
          &                           mpiGatherScalar            , mpiGatherInt1D         , &
          &                           mpiGatherIntScalar         , mpiGatherLogicalScalar
     procedure :: communicatorPush => mpiCommunicatorPush
     procedure :: communicatorPop  => mpiCommunicatorPop
     procedure :: stateInitialize  => mpiStateInitialize
  end type mpiObject

  ! Declare an object for interaction with MPI.
  type(mpiObject) :: mpiSelf

  ! Wrapper types for MPI shared resources.
  type :: mpiMemory
     !!{
     A wrapper type for MPI shared memory.
     !!}
     type(c_ptr) :: memory
   contains
     final :: mpiMemoryDestructor
  end type mpiMemory

  type :: mpiWindow
     !!{
     A wrapper type for MPI windows.
     !!}
#ifdef USEMPI
     type(MPI_Win) :: window
#endif
   contains
     final :: mpiWindowDestructor
  end type mpiWindow
  
  ! Define an MPI counter type.
  type :: mpiCounter
     !!{
     An MPI-global counter class. The counter can be incremented and will return a globally unique integer, beginning at 0.
     !!}
#ifdef USEMPI
     type   (mpiWindow      ), pointer :: window              => null()
     type   (MPI_Datatype   )          :: typeClass
     type   (mpiMemory      ), pointer :: counter             => null()
     type   (resourceManager)          :: windowManager                , counterManager
     integer(c_size_t       )          :: countThreadsWaiting          , countIncrementsHeld, &
          &                               counterHeld
#else
     integer(c_size_t       )          :: counter
#endif
     !$ type(ompLock        )          :: ompLock_
   contains
     !![
     <methods>
       <method description="Increment the counter and return the new value." method="increment"/>
       <method description="Decrement the counter and return the new value." method="decrement"/>
       <method description="Get the current value of the counter."           method="get"      />
       <method description="Reset the counter."                              method="reset"    />
     </methods>
     !!]
     procedure :: increment => counterIncrement
     procedure :: decrement => counterDecrement
     procedure :: get       => counterGet
     procedure :: reset     => counterReset
  end type mpiCounter

  interface mpiCounter
     module procedure counterConstructor
  end interface mpiCounter

  ! Record of whether we're running under MPI or not.
  logical            :: mpiIsActiveValue=.false.

  ! Tags.
  integer, parameter :: tagRequestForData= 1, tagState=2
  integer, parameter :: nullRequester    =-1

contains

  subroutine mpiInitialize(mpiThreadingRequired)
    !!{
    Initialize MPI.
    !!}
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Thread_Funneled, MPI_Thread_Single, MPI_Comm_World
    use :: Error  , only : Error_Report
#endif
    implicit none
    integer, optional, intent(in   ) :: mpiThreadingRequired
#ifdef USEMPI
    integer                          :: iError             , mpiThreadingProvided
    !![
    <optionalArgument name="mpiThreadingRequired" defaultsTo="MPI_Thread_Funneled" if="USEMPI"/>
    !!]

    if (mpiThreadingRequired_ == MPI_Thread_Single) then
       call MPI_Init              (                                           iError)
       if (iError               /= 0                    ) call Error_Report('failed to initialize MPI'                                        //{introspection:location})
    else
       call MPI_Init_Thread       (mpiThreadingRequired_,mpiThreadingProvided,iError)
       if (iError               /= 0                    ) call Error_Report('failed to initialize MPI'                                        //{introspection:location})
       if (mpiThreadingProvided <  mpiThreadingRequired_) call Error_Report('MPI library does not provide required level of threading support'//{introspection:location})
    end if
    ! Initialize our communicator to "world".
    allocate(mpiSelf%communicatorStack(1))
    mpiSelf%iCommunicator                           =1
    mpiSelf%communicatorStack(mpiSelf%iCommunicator)=MPI_Comm_World
    mpiSelf%communicator                            =mpiSelf%communicatorStack(mpiSelf%iCommunicator)
    call mpiSelf%stateInitialize()
    ! Record that MPI is active.
    mpiIsActiveValue=.true.
#else
    !$GLC attributes unused :: mpiThreadingRequired
#endif
    return
  end subroutine mpiInitialize

  subroutine mpiFinalize()
    !!{
    Finalize MPI.
    !!}
#ifdef USEMPI
    use Error, only : Error_Report
    implicit none
    integer :: iError

    call MPI_Finalize(iError)
    if (iError /= 0) call Error_Report('failed to finalize MPI'//{introspection:location})
    ! Record that MPI is inactive.
    mpiIsActiveValue=.false.
#endif
    return
  end subroutine mpiFinalize

  subroutine mpiBarrier()
    !!{
    Block until all MPI processes are synchronized.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Barrier
    implicit none
    integer :: iError

    call MPI_Barrier(mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI barrier failed'//{introspection:location})
#endif
    return
  end subroutine mpiBarrier

  logical function mpiIsActive(self)
    !!{
    Return true if MPI is active.
    !!}
    implicit none
    class(mpiObject), intent(in   ) :: self
    !$GLC attributes unused :: self

    mpiIsActive=mpiIsActiveValue
    return
  end function mpiIsActive

  logical function mpiIsMaster(self)
    !!{
    Return true if this is the master process.
    !!}
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiIsMaster=(.not.self%isActive() .or. self%rank() == 0)
#else
    !$GLC attributes unused :: self
    mpiIsMaster=.true.
#endif
    return
  end function mpiIsMaster

  integer function mpiGetRank(self)
    !!{
    Return MPI rank.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiGetRank=self%rankValue
#else
    !$GLC attributes unused :: self
    mpiGetRank=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetRank

  integer function mpiGetRankOnNode(self)
    !!{
    Return MPI rank on the node.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiGetRankOnNode=self%rankOnNodeValue
#else
    !$GLC attributes unused :: self
    mpiGetRankOnNode=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetRankOnNode

  function mpiGetRankLabel(self)
    !!{
    Return MPI rank label.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
#ifndef USEMPI
    use :: Error             , only : Error_Report
#endif
    implicit none
    type     (varying_string)                :: mpiGetRankLabel
    class    (mpiObject     ), intent(in   ) :: self
#ifdef USEMPI
    character(len=4         )                :: label
#endif

#ifdef USEMPI
    if (self%isActive()) then
       write (label,'(i4.4)') self%rankValue
       mpiGetRankLabel=label
    else
       mpiGetRankLabel=''
    end if
#else
    !$GLC attributes unused :: self
    mpiGetRankLabel=''
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetRankLabel

  integer function mpiGetCount(self)
    !!{
    Return MPI count.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiGetCount=self%countValue
#else
    !$GLC attributes unused :: self
    mpiGetCount=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetCount

  integer function mpiGetNodeCount(self)
    !!{
    Return count of nodes used by MPI.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiGetNodeCount=self%nodeCountValue
#else
    !$GLC attributes unused :: self
    mpiGetNodeCount=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetNodeCount

  integer function mpiGetCountOnNode(self)
    !!{
    Return count of nodes used by MPI.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class(mpiObject), intent(in   ) :: self

#ifdef USEMPI
    mpiGetCountOnNode=self%countOnNodeValue
#else
    !$GLC attributes unused :: self
    mpiGetCountOnNode=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetCountOnNode

  integer function mpiGetNodeAffinity(self,rank)
    !!{
    Return node affinity of given MPI process.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class  (mpiObject), intent(in   )           :: self
    integer           , intent(in   ), optional :: rank
#ifdef USEMPI
    integer                                     :: rankActual
#endif

#ifdef USEMPI
    rankActual=self%rank()
    if (present(rank)) rankActual=rank
    mpiGetNodeAffinity=self%nodeAffinities(rankActual)
#else
    !$GLC attributes unused :: self, rank
    mpiGetNodeAffinity=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetNodeAffinity

  function mpiGetHostAffinity(self)
    !!{
    Return host affinity of this MPI process.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
#ifndef USEMPI
    use :: Error             , only : Error_Report
#endif
    implicit none
    type (varying_string)                :: mpiGetHostAffinity
    class(mpiObject     ), intent(in   ) :: self

#ifdef USEMPI
    mpiGetHostAffinity=self%hostName
#else
    !$GLC attributes unused :: self
    mpiGetHostAffinity=""
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGetHostAffinity

  logical function mpiMessageWaiting(self,from,tag)
    !!{
    Return true if an MPI message (matching the optional {\normalfont \ttfamily from} and {\normalfont \ttfamily tag} if given) is waiting for receipt.
    !!}
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Status  , MPI_Any_Source, MPI_Any_Tag, MPI_IProbe
#endif
    use :: Error  , only : Error_Report
    implicit none
    class  (mpiObject ), intent(in   )           :: self
    integer            , intent(in   ), optional :: from         , tag
#ifdef USEMPI
    type   (MPI_Status)                          :: messageStatus
    integer                                      :: fromActual   , tagActual, iError
#endif

#ifdef USEMPI
    fromActual=MPI_Any_Source
    tagActual =MPI_Any_Tag
    if (present(from)) fromActual=from
    if (present(tag ) ) tagActual=tag
    call MPI_IProbe(fromActual,tagActual,mpiSelf%communicator,mpiMessageWaiting,messageStatus,iError)
    if (iError /= 0) call Error_Report('failed to probe for waiting messages'//{introspection:location})
#else
    !$GLC attributes unused :: self, from, tag
    mpiMessageWaiting=.false.
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMessageWaiting

  function mpiRequestData1D(self,requestFrom,array)
    !!{
    Request and receive data from other MPI processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_Request , MPI_Status , MPI_Wait      , MPI_ISend           , &
         &                 MPI_Recv    , MPI_Integer, MPI_Any_Source, MPI_Double_Precision
#endif
    implicit none
    class           (mpiObject), intent(in   )                                           :: self
    integer                    , intent(in   ), dimension(                            :) :: requestFrom
    double precision           , intent(in   ), dimension(          :                  ) :: array
    double precision                          , dimension(size(array),size(requestFrom)) :: mpiRequestData1D
#ifdef USEMPI
    double precision                          , dimension(size(array)                  ) :: receivedData
    integer                                   , dimension(                            1) :: requester       , requestedBy
    type            (MPI_Request)             , dimension(         0: self%countValue-1) :: requestFromID
    type            (MPI_Request), allocatable, dimension(                            :) :: requestID       , requestIDtemp
    type            (MPI_Status )                                                        :: messageStatus
    integer                                                                              :: i               , iError       , &
         &                                                                                  iRequest        , receivedFrom , &
         &                                                                                  j
#endif

#ifdef USEMPI
    ! Record our own rank as the requester.
    requester=self%rank()
    ! Send requests.
    call mpiBarrier()
    do i=0,self%count()-1
       if (any(requestFrom == i)) then
          call MPI_ISend(requester    ,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       else
          call MPI_ISend(nullRequester,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       end if
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    allocate(requestID(size(requestFrom)))
    iRequest=0
    do i=0,self%count()-1
       ! Receive the request.
       call MPI_Recv(requestedBy,1,MPI_Integer,MPI_Any_Source,tagRequestForData,mpiSelf%communicator,messageStatus,iError)
       ! Check for a non-null request.
       if (requestedBy(1) /= nullRequester) then
          ! Expand the requestID buffer as required.
          iRequest=iRequest+1
          if (iRequest > size(requestID)) then
             call Move_Alloc(requestID,requestIDtemp)
             allocate(requestID(2*iRequest))
             requestID(1:size(requestIDtemp))=requestIDtemp
             deallocate(requestIDtemp)
          end if
          ! Send our data in reply.
          call MPI_ISend(array,size(array),MPI_Double_Precision,requestedBy(1),tagState,mpiSelf%communicator,requestID(iRequest),iError)
       end if
    end do
    call mpiBarrier()
    ! Wait until all of our sends have been received.
    do i=0,self%count()-1
       call MPI_Wait(requestFromID(i),messageStatus,iError)
    end do
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(receivedData,size(array),MPI_Double_Precision,MPI_Any_Source,tagState,mpiSelf%communicator,messageStatus,iError)
       ! Find who sent this data and apply to the relevant part of the results array.
       receivedFrom=messageStatus%MPI_Source
       do j=1,size(requestFrom)
          if (requestFrom(j) == receivedFrom) mpiRequestData1D(:,j)=receivedData
       end do
    end do
    ! Wait until all of our sends have been received.
    do i=1,iRequest
       call MPI_Wait(requestID(i),messageStatus,iError)
    end do
    call mpiBarrier()
    ! Deallocate request ID workspace.
    deallocate(requestID)
#else
    !$GLC attributes unused :: self, requestFrom, array
    mpiRequestData1D=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiRequestData1D

  function mpiRequestData2D(self,requestFrom,array)
    !!{
    Request and receive data from other MPI processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_Request , MPI_Status , MPI_Wait      , MPI_ISend           , &
         &                 MPI_Recv    , MPI_Integer, MPI_Any_Source, MPI_Double_Precision
#endif
    implicit none
    class           (mpiObject  ), intent(in   )                                                                   :: self
    integer                      , intent(in   ), dimension(                                                    :) :: requestFrom
    double precision             , intent(in   ), dimension(                :,                :                  ) :: array
    double precision                            , dimension(size(array,dim=1),size(array,dim=2),size(requestFrom)) :: mpiRequestData2D
#ifdef USEMPI
    double precision                            , dimension(size(array,dim=1),size(array,dim=2)                  ) :: receivedData
    integer                                     , dimension(                                                    1) :: requester       , requestedBy
    type            (MPI_Request)               , dimension(                                 0: self%countValue-1) :: requestFromID
    type            (MPI_Request), allocatable  , dimension(                                                    :) :: requestID       , requestIDtemp
    type            (MPI_Status )                                                                                  :: messageStatus
    integer                                                                                                        :: i               , iError       , &
         &                                                                                                            iRequest        , j            , &
         &                                                                                                            receivedFrom
#endif

#ifdef USEMPI
    ! Record our own rank as the requester.
    requester=self%rank()
    ! Send requests.
    call mpiBarrier()
    do i=0,self%count()-1
       if (any(requestFrom == i)) then
          call MPI_ISend(requester    ,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       else
          call MPI_ISend(nullRequester,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       end if
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    allocate(requestID(size(requestFrom)))
    iRequest=0
    do i=0,self%count()-1
       ! Receive the request.
       call MPI_Recv(requestedBy,1,MPI_Integer,MPI_Any_Source,tagRequestForData,mpiSelf%communicator,messageStatus,iError)
       ! Check for a non-null request.
       if (requestedBy(1) /= nullRequester) then
          ! Expand the requestID buffer as required.
          iRequest=iRequest+1
          if (iRequest > size(requestID)) then
             call Move_Alloc(requestID,requestIDtemp)
             allocate(requestID(2*iRequest))
             requestID(1:size(requestIDtemp))=requestIDtemp
             deallocate(requestIDtemp)
          end if
          ! Send our data in reply.
          call MPI_ISend(array,product(shape(array)),MPI_Double_Precision,requestedBy(1),tagState,mpiSelf%communicator,requestID(iRequest),iError)
       end if
    end do
    call mpiBarrier()
    ! Wait until all of our sends have been received.
    do i=0,self%count()-1
       call MPI_Wait(requestFromID(i),messageStatus,iError)
    end do
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(receivedData,product(shape(array)),MPI_Double_Precision,requestFrom(i),tagState,mpiSelf%communicator,messageStatus,iError)
       ! Find who sent this data and apply to the relevant part of the results array.
       receivedFrom=messageStatus%MPI_Source
       do j=1,size(requestFrom)
          if (requestFrom(j) == receivedFrom) mpiRequestData2D(:,:,j)=receivedData
       end do
    end do
    ! Wait until all of our sends have been received.
    do i=1,iRequest
       call MPI_Wait(requestID(i),messageStatus,iError)
    end do
    call mpiBarrier()
    ! Deallocate request ID workspace.
    deallocate(requestID)
#else
    !$GLC attributes unused :: self, requestFrom, array
    mpiRequestData2D=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiRequestData2D

  function mpiRequestDataInt1D(self,requestFrom,array)
    !!{
    Request and receive data from other MPI processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_Request              , MPI_Status    , MPI_Wait, MPI_ISend, &
         &                 MPI_Recv    , MPI_Integer, MPI_Any_Source
#endif
    implicit none
    class  (mpiObject  ), intent(in   )                                           :: self
    integer             , intent(in   ), dimension(                            :) :: requestFrom
    integer             , intent(in   ), dimension(          :                  ) :: array
    integer                            , dimension(size(array),size(requestFrom)) :: mpiRequestDataInt1D
#ifdef USEMPI
    integer                            , dimension(size(array)                  ) :: receivedData
    integer                            , dimension(                            1) :: requester       , requestedBy
    type   (MPI_Request)               , dimension(         0: self%countValue-1) :: requestFromID
    type   (MPI_Request), allocatable  , dimension(                            :) :: requestID       , requestIDtemp
    type   (MPI_Status )                                                          :: messageStatus
    integer                                                                       :: i               , iError       , &
         &                                                                           iRequest        , j            , &
         &                                                                           receivedFrom
#endif

#ifdef USEMPI
    ! Record our own rank as the requester.
    requester=self%rank()
    ! Send requests.
    call mpiBarrier()
    do i=0,self%count()-1
       if (any(requestFrom == i)) then
          call MPI_ISend(requester    ,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       else
          call MPI_ISend(nullRequester,1,MPI_Integer,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       end if
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    allocate(requestID(size(requestFrom)))
    iRequest=0
    do i=0,self%count()-1
       ! Receive the request.
       call MPI_Recv(requestedBy,1,MPI_Integer,MPI_Any_Source,tagRequestForData,mpiSelf%communicator,messageStatus,iError)
       ! Check for a non-null request.
       if (requestedBy(1) /= nullRequester) then
          ! Expand the requestID buffer as required.
          iRequest=iRequest+1
          if (iRequest > size(requestID)) then
             call Move_Alloc(requestID,requestIDtemp)
             allocate(requestID(2*iRequest))
             requestID(1:size(requestIDtemp))=requestIDtemp
             deallocate(requestIDtemp)
          end if
          ! Send our data in reply.
          call MPI_ISend(array,size(array),MPI_Integer,requestedBy(1),tagState,mpiSelf%communicator,requestID(iRequest),iError)
       end if
    end do
    call mpiBarrier()
    ! Wait until all of our sends have been received.
    do i=0,self%count()-1
       call MPI_Wait(requestFromID(i),messageStatus,iError)
    end do
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(receivedData,size(array),MPI_Integer,requestFrom(i),tagState,mpiSelf%communicator,messageStatus,iError)
       ! Find who sent this data and apply to the relevant part of the results array.
       receivedFrom=messageStatus%MPI_Source
       do j=1,size(requestFrom)
          if (requestFrom(j) == receivedFrom) mpiRequestDataInt1D(:,j)=receivedData
       end do
    end do
    ! Wait until all of our sends have been received.
    do i=1,iRequest
       call MPI_Wait(requestID(i),messageStatus,iError)
    end do
    call mpiBarrier()
    ! Deallocate request ID workspace.
    deallocate(requestID)
#else
    !$GLC attributes unused :: self, requestFrom, array
    mpiRequestDataInt1D=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiRequestDataInt1D

  function mpiRequestDataLogical1D(self,requestFrom,array)
    !!{
    Request and receive data from other MPI processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_Request , MPI_Status , MPI_Wait      , MPI_ISend, &
         &                 MPI_Recv    , MPI_Logical, MPI_Any_Source
#endif
    implicit none
    class  (mpiObject  ), intent(in   )                                           :: self
    integer             , intent(in   ), dimension(                            :) :: requestFrom
    logical             , intent(in   ), dimension(          :                  ) :: array
    logical                            , dimension(size(array),size(requestFrom)) :: mpiRequestDataLogical1D
#ifdef USEMPI
    logical                            , dimension(size(array)                  ) :: receivedData
    integer                            , dimension(                            1) :: requester              , requestedBy
    type   (MPI_Request)               , dimension(         0: self%countValue-1) :: requestFromID
    type   (MPI_Request), allocatable  , dimension(                            :) :: requestID              , requestIDtemp
    type   (MPI_Status )                                                          :: messageStatus
    integer                                                                       :: i                      , iError       , &
         &                                                                           iRequest               , j            , &
         &                                                                           receivedFrom
#endif

#ifdef USEMPI
    ! Record our own rank as the requester.
    requester=self%rank()
    ! Send requests.
    call mpiBarrier()
    do i=0,self%count()-1
       if (any(requestFrom == i)) then
          call MPI_ISend(requester    ,1,MPI_Logical,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       else
          call MPI_ISend(nullRequester,1,MPI_Logical,i,tagRequestForData,mpiSelf%communicator,requestFromID(i),iError)
       end if
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    allocate(requestID(size(requestFrom)))
    iRequest=0
    do i=0,self%count()-1
       ! Receive the request.
       call MPI_Recv(requestedBy,1,MPI_Logical,MPI_Any_Source,tagRequestForData,mpiSelf%communicator,messageStatus,iError)
       ! Check for a non-null request.
       if (requestedBy(1) /= nullRequester) then
          ! Expand the requestID buffer as required.
          iRequest=iRequest+1
          if (iRequest > size(requestID)) then
             call Move_Alloc(requestID,requestIDtemp)
             allocate(requestID(2*iRequest))
             requestID(1:size(requestIDtemp))=requestIDtemp
             deallocate(requestIDtemp)
          end if
          ! Send our data in reply.
          call MPI_ISend(array,size(array),MPI_Logical,requestedBy(1),tagState,mpiSelf%communicator,requestID(iRequest),iError)
       end if
    end do
    call mpiBarrier()
    ! Wait until all of our sends have been received.
    do i=0,self%count()-1
       call MPI_Wait(requestFromID(i),messageStatus,iError)
    end do
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(receivedData,size(array),MPI_Logical,requestFrom(i),tagState,mpiSelf%communicator,messageStatus,iError)
       ! Find who sent this data and apply to the relevant part of the results array.
       receivedFrom=messageStatus%MPI_Source
       do j=1,size(requestFrom)
          if (requestFrom(j) == receivedFrom) mpiRequestDataLogical1D(:,j)=receivedData
       end do
    end do
    ! Wait until all of our sends have been received.
    do i=1,iRequest
       call MPI_Wait(requestID(i),messageStatus,iError)
    end do
    call mpiBarrier()
    ! Deallocate request ID workspace.
    deallocate(requestID)
#else
    !$GLC attributes unused :: self, requestFrom, array
    mpiRequestDataLogical1D=.false.
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiRequestDataLogical1D

  subroutine mpiBroadcastDataScalar(self,sendFrom,scalar)
    !!{
    Broadcast data to all other MPI processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Double_Precision, MPI_Bcast
#endif
    implicit none
    class           (mpiObject), intent(in   ) :: self
    integer                    , intent(in   ) :: sendFrom
    double precision           , intent(inout) :: scalar
#ifdef USEMPI
    integer                                    :: status
    !$GLC attributes unused :: self
    
    call MPI_Bcast(scalar,1,MPI_Double_Precision,sendFrom,mpiSelf%communicator,status)
    if (status /= 0) call Error_Report('failed to broadcast data'//{introspection:location})
#else
    !$GLC attributes unused :: self, sendFrom, scalar
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiBroadcastDataScalar

  subroutine mpiBroadcastDataSizeTScalar(self,sendFrom,scalar)
    !!{
    Broadcast data to all other MPI processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Integer8, MPI_Bcast
#endif
    implicit none
    class           (mpiObject), intent(in   ) :: self
    integer                    , intent(in   ) :: sendFrom
    integer         (c_size_t ), intent(inout) :: scalar
#ifdef USEMPI
    integer                                    :: status
    !$GLC attributes unused :: self
    
    call MPI_Bcast(scalar,1,MPI_Integer8,sendFrom,mpiSelf%communicator,status)
    if (status /= 0) call Error_Report('failed to broadcast data'//{introspection:location})
#else
    !$GLC attributes unused :: self, sendFrom, scalar
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiBroadcastDataSizeTScalar
  
  subroutine mpiBroadcastData1D(self,sendFrom,array)
    !!{
    Broadcast data to all other MPI processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Double_Precision, MPI_Bcast
#endif
    implicit none
    class           (mpiObject), intent(in   )               :: self
    integer                    , intent(in   )               :: sendFrom
    double precision           , intent(inout), dimension(:) :: array
#ifdef USEMPI
    integer                                                  :: status
    !$GLC attributes unused :: self
    
    call MPI_Bcast(array,size(array),MPI_Double_Precision,sendFrom,mpiSelf%communicator,status)
    if (status /= 0) call Error_Report('failed to broadcast data'//{introspection:location})
#else
    !$GLC attributes unused :: self, sendFrom, array
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiBroadcastData1D
  
  subroutine mpiBroadcastData2D(self,sendFrom,array)
    !!{
    Broadcast data to all other MPI processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Double_Precision, MPI_Bcast
#endif
    implicit none
    class           (mpiObject), intent(in   )                 :: self
    integer                    , intent(in   )                 :: sendFrom
    double precision           , intent(inout), dimension(:,:) :: array
#ifdef USEMPI
    integer                                                    :: status
    !$GLC attributes unused :: self
    
    call MPI_Bcast(array,size(array),MPI_Double_Precision,sendFrom,mpiSelf%communicator,status)
    if (status /= 0) call Error_Report('failed to broadcast data'//{introspection:location})
#else
    !$GLC attributes unused :: self, sendFrom, array
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiBroadcastData2D
  
  subroutine mpiBroadcastData3D(self,sendFrom,array)
    !!{
    Broadcast data to all other MPI processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Double_Precision, MPI_Bcast
#endif
    implicit none
    class           (mpiObject), intent(in   )                   :: self
    integer                    , intent(in   )                   :: sendFrom
    double precision           , intent(inout), dimension(:,:,:) :: array
#ifdef USEMPI
    integer                                                      :: status
    !$GLC attributes unused :: self

    call MPI_Bcast(array,size(array),MPI_Double_Precision,sendFrom,mpiSelf%communicator,status)
    if (status /= 0) call Error_Report('failed to broadcast data'//{introspection:location})
#else
    !$GLC attributes unused :: self, sendFrom, array
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiBroadcastData3D
  
  function mpiSumArrayInt(self,array,mask)
    !!{
    Sum an integer array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer, MPI_Sum
#endif
    implicit none
    class  (mpiObject), intent(in   )                                    :: self
    integer           , intent(in   ), dimension( :          )           :: array
    logical           , intent(in   ), dimension(0:          ), optional :: mask
    integer                          , dimension(size(array))            :: mpiSumArrayInt
#ifdef USEMPI
    integer                          , dimension(size(array))            :: maskedArray
    integer                                                              :: iError        , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayInt,size(array),MPI_Integer,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayInt=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayInt

  function mpiSumArraySizeT(self,array,mask)
    !!{
    Sum an integer array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer8, MPI_Sum
#endif
    implicit none
    class  (mpiObject), intent(in   )                                    :: self
    integer(c_size_t ), intent(in   ), dimension( :          )           :: array
    logical           , intent(in   ), dimension(0:          ), optional :: mask
    integer(c_size_t )               , dimension(size(array))            :: mpiSumArraySizeT
#ifdef USEMPI
    integer(c_size_t )               , dimension(size(array))            :: maskedArray
    integer                                                              :: iError        , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0_c_size_t
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArraySizeT,size(array),MPI_Integer8,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArraySizeT=0_c_size_t
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArraySizeT

  function mpiSumArrayTwoSizeT(self,array,mask)
    !!{
    Sum a rank-2 integer array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer8, MPI_Sum
#endif
    implicit none
    class  (mpiObject), intent(in   )                                                           :: self
    integer(c_size_t ), intent(in   ), dimension( :               , :               )           :: array
    logical           , intent(in   ), dimension(0:                                 ), optional :: mask
    integer(c_size_t )               , dimension(size(array,dim=1),size(array,dim=2))           :: mpiSumArrayTwoSizeT
#ifdef USEMPI
    integer(c_size_t )               , dimension(size(array,dim=1),size(array,dim=2))           :: maskedArray
    integer                                                                                     :: iError             , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0_c_size_t
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayTwoSizeT,size(array),MPI_Integer8,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayTwoSizeT=0_c_size_t
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayTwoSizeT

  function mpiSumArrayThreeSizeT(self,array,mask)
    !!{
    Sum a rank-3 integer array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer8, MPI_Sum
#endif
    implicit none
    class  (mpiObject), intent(in   )                                                                             :: self
    integer(c_size_t ), intent(in   ), dimension( :               , :               , :               )           :: array
    logical           , intent(in   ), dimension(0:                                                   ), optional :: mask
    integer(c_size_t )               , dimension(size(array,dim=1),size(array,dim=2),size(array,dim=3))           :: mpiSumArrayThreeSizeT
#ifdef USEMPI
    integer(c_size_t )               , dimension(size(array,dim=1),size(array,dim=2),size(array,dim=3))           :: maskedArray
    integer                                                                                                       :: iError               , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0_c_size_t
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayThreeSizeT,size(array),MPI_Integer8,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayThreeSizeT=0_c_size_t
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayThreeSizeT
  
  integer function mpiSumScalarInt(self,scalar,mask)
    !!{
    Sum an integer scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class  (mpiObject), intent(in   )                         :: self
    integer           , intent(in   )                         :: scalar
    logical           , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%sum([scalar],mask)
    mpiSumScalarInt=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiSumScalarInt=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumScalarInt

  function mpiSumScalarSizeT(self,scalar,mask)
    !!{
    Sum a {\normalfont \ttfamily size\_t} scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    integer(c_size_t )                                        :: mpiSumScalarSizeT
    class  (mpiObject), intent(in   )                         :: self
    integer(c_size_t ), intent(in   )                         :: scalar
    logical           , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer(c_size_t )               , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%sum([scalar],mask)
    mpiSumScalarSizeT=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiSumScalarSizeT=0_c_size_t
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumScalarSizeT

  function mpiSumArrayReal(self,array,mask)
    !!{
    Sum a real array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Real, MPI_Sum
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    real                       , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    real                                      , dimension(size(array))            :: mpiSumArrayReal
#ifdef USEMPI
    real                                      , dimension(size(array))            :: maskedArray
    integer                                                                       :: iError          , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayReal,size(array),MPI_Real,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayReal=0.0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayReal

  function mpiSumArrayDouble(self,array,mask)
    !!{
    Sum an integer array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Sum
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array))            :: mpiSumArrayDouble
#ifdef USEMPI
    double precision                          , dimension(size(array))            :: maskedArray
    integer                                                                       :: iError           , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0d0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayDouble,size(array),MPI_Double_Precision,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayDouble=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayDouble

  function mpiSumArrayTwoDouble(self,array,mask)
    !!{
    Sum an rank-2 double array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Sum
#endif
    implicit none
    class           (mpiObject), intent(in   )                                                           :: self
    double precision           , intent(in   ), dimension( :               , :               )           :: array
    logical                    , intent(in   ), dimension(0:                                 ), optional :: mask
    double precision                          , dimension(size(array,dim=1),size(array,dim=2))           :: mpiSumArrayTwoDouble
#ifdef USEMPI
    double precision                          , dimension(size(array,dim=1),size(array,dim=2))            :: maskedArray
    integer                                                                                               :: iError           , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0d0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayTwoDouble,size(array),MPI_Double_Precision,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayTwoDouble=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayTwoDouble

  function mpiSumArrayThreeDouble(self,array,mask)
    !!{
    Sum an rank-3 double array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Sum
#endif
    implicit none
    class           (mpiObject), intent(in   )                                                                             :: self
    double precision           , intent(in   ), dimension( :               , :               , :               )           :: array
    logical                    , intent(in   ), dimension(0:                                                   ), optional :: mask
    double precision                          , dimension(size(array,dim=1),size(array,dim=2),size(array,dim=3))           :: mpiSumArrayThreeDouble
#ifdef USEMPI
    double precision                          , dimension(size(array,dim=1),size(array,dim=2),size(array,dim=3))           :: maskedArray
    integer                                                                                                                :: iError                , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0d0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayThreeDouble,size(array),MPI_Double_Precision,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiSumArrayThreeDouble=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumArrayThreeDouble

  double precision function mpiSumScalarDouble(self,scalar,mask)
    !!{
    Sum an integer scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    double precision                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%sum([scalar],mask)
    mpiSumScalarDouble=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiSumScalarDouble=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiSumScalarDouble

  function mpiAverageArray(self,array,mask)
    !!{
    Average an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Sum
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiAverageArray
#ifdef USEMPI
    double precision                          , dimension(size(array) )           :: maskedArray
    integer                                                                       :: iError         , activeCount
#endif

#ifdef USEMPI
    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0d0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiAverageArray,size(array),MPI_Double_Precision,MPI_Sum,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
    ! Convert the sum into an average.
    mpiAverageArray=mpiAverageArray/dble(activeCount)
#else
    !$GLC attributes unused :: self, array, mask
    mpiAverageArray=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiAverageArray

  function mpiMedianArray(self,array,mask)
    !!{
    Find the median of an array over all processes, returning it to all processes.
    !!}
#ifdef USEMPI
    use :: Sorting, only : sort
#else
    use :: Error  , only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                                                   :: self
    integer                    , intent(in   ), dimension(:                          )           :: array
    logical                    , intent(in   ), dimension(:                          ), optional :: mask
    integer                                   , dimension(size(array)                )           :: mpiMedianArray
#ifdef USEMPI
    integer                                   , dimension(size(array),self%countValue)           :: allArray
    integer                                   , dimension(1:2                        )           :: indexMedian
    integer                                                                                      :: i             , activeCount
#endif

#ifdef USEMPI
    ! Get count of active process.
    if (present(mask)) then
       activeCount=self%countValue-count(mask)
    else
       activeCount=self%countValue
    end if
    ! Find the indices corresponding to the median.
    if (mod(activeCount,2) == 1) then
       indexMedian=               activeCount/2+1
    else
       indexMedian=[activeCount/2,activeCount/2+1]
    end if
    ! Gather the array from all processes.
    allArray=self%gather(array)
    ! Iterate over array index.
    do i=1,size(array)
       ! Set masked values to huge.
       if (present(mask)) then
          where (mask)
             allArray(i,:)=huge(1)
          end where
       end if
       ! Sort over processes.
       call sort(allArray(i,:))
       ! Compute the median.
       mpiMedianArray(i)=(allArray(i,indexMedian(1))+allArray(i,indexMedian(2)))/2
    end do
#else
    !$GLC attributes unused :: self, array, mask
    mpiMedianArray=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMedianArray

  double precision function mpiAverageScalar(self,scalar,mask)
    !!{
    Find the maximum values of a scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    double precision                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%average([scalar],mask)
    mpiAverageScalar=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiAverageScalar=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiAverageScalar

  function mpiMaxvalArray(self,array,mask)
    !!{
    Find the maximum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Max
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiMaxvalArray
#ifdef USEMPI
    double precision                          , dimension(size(array) )           :: maskedArray
    integer                                                                       :: iError
#endif

#ifdef USEMPI
    ! Find the maximum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-HUGE(1.0d0)
    end if
    call MPI_AllReduce(maskedArray,mpiMaxvalArray,size(array),MPI_Double_Precision,MPI_Max,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiMaxvalArray=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMaxvalArray

  double precision function mpiMaxvalScalar(self,scalar,mask)
    !!{
    Find the maximum values of a scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    double precision                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%maxval([scalar],mask)
    mpiMaxvalScalar=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiMaxvalScalar=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMaxvalScalar

  function mpiMaxvalArraySizeT(self,array,mask)
    !!{
    Find the maximum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer8, MPI_Max
#endif
    implicit none
    class  (mpiObject), intent(in   )                                    :: self
    integer(c_size_t ), intent(in   ), dimension( :          )           :: array
    logical           , intent(in   ), dimension(0:          ), optional :: mask
    integer(c_size_t )               , dimension(size(array) )           :: mpiMaxvalArraySizeT
#ifdef USEMPI
    integer(c_size_t )               , dimension(size(array) )           :: maskedArray
    integer                                                              :: iError
#endif

#ifdef USEMPI
    ! Find the maximum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-huge(1_c_size_t)
    end if
    call MPI_AllReduce(maskedArray,mpiMaxvalArraySizeT,size(array),MPI_Integer8,MPI_Max,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiMaxvalArraySizeT=0_c_size_t
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMaxvalArraySizeT

  function mpiMaxvalScalarSizeT(self,scalar,mask)
    !!{
    Find the maximum values of a {\normalfont \ttfamily size\_t} scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    integer(c_size_t )                                        :: mpiMaxvalScalarSizeT
    class  (mpiObject), intent(in   )                         :: self
    integer(c_size_t ), intent(in   )                         :: scalar
    logical           , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer(c_size_t )               , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%maxval([scalar],mask)
    mpiMaxvalScalarSizeT=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiMaxvalScalarSizeT=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMaxvalScalarSizeT

  function mpiMaxloc(self,array,mask)
    !!{
    Find the rank of the process having maximum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_2Double_Precision, MPI_MaxLoc
#endif
    implicit none
    class           (mpiObject), intent(in   )                                      :: self
    double precision           , intent(in   ), dimension( :            )           :: array
    logical                    , intent(in   ), dimension(0:            ), optional :: mask
    integer                                   , dimension(   size(array))           :: mpiMaxloc
#ifdef USEMPI
    double precision                          , dimension(2 ,size(array))           :: arrayIn  , arrayOut
    integer                                                                         :: iError
#endif

#ifdef USEMPI
    ! Find the maximum over all processes.
    arrayIn(1,:)=array
    if (present(mask)) then
       if (.not.mask(self%rank())) arrayIn(1,:)=-HUGE(1.0d0)
    end if
    arrayIn(2,:)=self%rank()
    call MPI_AllReduce(arrayIn,arrayOut,size(array),MPI_2Double_Precision,MPI_MaxLoc,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
    mpiMaxloc=int(arrayOut(2,:))
#else
    !$GLC attributes unused :: self, array, mask
    mpiMaxloc=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMaxloc

  function mpiMinvalArray(self,array,mask)
    !!{
    Find the minimum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Double_Precision, MPI_Min
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiMinvalArray
#ifdef USEMPI
    double precision                          , dimension(size(array) )           :: maskedArray
    integer                                                                       :: iError
#endif

#ifdef USEMPI
   ! Find the minimum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-HUGE(1.0d0)
    end if
    call MPI_AllReduce(maskedArray,mpiMinvalArray,size(array),MPI_Double_Precision,MPI_Min,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiMinvalArray=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMinvalArray

  function mpiMinvalIntArray(self,array,mask)
    !!{
    Find the minimum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_Integer, MPI_Min
#endif
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    integer                    , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    integer                                   , dimension(size(array) )           :: mpiMinvalIntArray
#ifdef USEMPI
    integer                                   , dimension(size(array) )           :: maskedArray
    integer                                                                       :: iError
#endif

#ifdef USEMPI
   ! Find the minimum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-huge(1)
    end if
    call MPI_AllReduce(maskedArray,mpiMinvalIntArray,size(array),MPI_Integer,MPI_Min,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
#else
    !$GLC attributes unused :: self, array, mask
    mpiMinvalIntArray=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMinvalIntArray

  double precision function mpiMinvalScalar(self,scalar,mask)
    !!{
    Find the minimum values of a scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    double precision                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%minval([scalar],mask)
    mpiMinvalScalar=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiMinvalScalar=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMinvalScalar

  integer function mpiMinvalIntScalar(self,scalar,mask)
    !!{
    Find the minimum values of a scalar over all processes, returning it to all processes.
    !!}
#ifndef USEMPI
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    integer                    , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer                                   , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=self%minval([scalar],mask)
    mpiMinvalIntScalar=array(1)
#else
    !$GLC attributes unused :: self, scalar, mask
    mpiMinvalIntScalar=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMinvalIntScalar

  function mpiMinloc(self,array,mask)
    !!{
    Find the rank of the process having minimum values of an array over all processes, returning it to all processes.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_AllReduce, MPI_2Double_Precision, MPI_MinLoc
#endif
    implicit none
    class           (mpiObject), intent(in   )                                      :: self
    double precision           , intent(in   ), dimension( :            )           :: array
    logical                    , intent(in   ), dimension(0:            ), optional :: mask
    integer                                   , dimension(   size(array))           :: mpiMinloc
#ifdef USEMPI
    double precision                          , dimension(2 ,size(array))           :: arrayIn  , arrayOut
    integer                                                                         :: iError
#endif

#ifdef USEMPI
    ! Find the minimum over all processes.
    arrayIn(1,:)=array
    if (present(mask)) then
       if (.not.mask(self%rank())) arrayIn(1,:)=-HUGE(1.0d0)
    end if
    arrayIn(2,:)=self%rank()
    call MPI_AllReduce(arrayIn,arrayOut,size(array),MPI_2Double_Precision,MPI_MinLoc,mpiSelf%communicator,iError)
    if (iError /= 0) call Error_Report('MPI all reduce failed'//{introspection:location})
    mpiMinloc=int(arrayOut(2,:))
#else
    !$GLC attributes unused :: self, array, mask
    mpiMinloc=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiMinloc

  logical function mpiAnyLogicalScalar(self,boolean,mask)
    !!{
    Return true if any of the given booleans is true over all processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_AllReduce, MPI_Logical, MPI_LOr
#endif
    implicit none
    class  (mpiObject), intent(in   )                         :: self
    logical           , intent(in   )                         :: boolean
    logical           , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer                                                   :: iError
    logical                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=boolean
    if (present(mask)) then
       if (.not.mask(self%rank())) array=.false.
    end if
    call MPI_AllReduce(array,mpiAnyLogicalScalar,size(array),MPI_Logical,MPI_LOr,mpiSelf%communicator,iError)
#else
    !$GLC attributes unused :: self, boolean, mask
    mpiAnyLogicalScalar=.false.
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiAnyLogicalScalar

  logical function mpiAllLogicalScalar(self,boolean,mask)
    !!{
    Return true if all of the given booleans are true over all processes.
    !!}
#ifndef USEMPI
    use :: Error  , only : Error_Report
#else
    use :: MPI_F08, only : MPI_AllReduce, MPI_Logical, MPI_LAnd
#endif
    implicit none
    class  (mpiObject), intent(in   )                         :: self
    logical           , intent(in   )                         :: boolean
    logical           , intent(in   ), dimension(:), optional :: mask
#ifdef USEMPI
    integer                                                   :: iError
    logical                          , dimension(1)           :: array
#endif

#ifdef USEMPI
    array=boolean
    if (present(mask)) then
       if (.not.mask(self%rank())) array=.false.
    end if
    call MPI_AllReduce(array,mpiAllLogicalScalar,size(array),MPI_Logical,MPI_LAnd,mpiSelf%communicator,iError)
#else
    !$GLC attributes unused :: self, boolean, mask
    mpiAllLogicalScalar=.false.
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiAllLogicalScalar

  function mpiGatherScalar(self,scalar)
    !!{
    Gather a scalar from all processes, returning it as a 1-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                :: self
    double precision           , intent(in   )                :: scalar
    double precision           , dimension(  self%countValue) :: mpiGatherScalar
#ifdef USEMPI
    double precision           , dimension(1,self%countValue) :: array
#endif

#ifdef USEMPI
    array=self%requestData(self%allRanks,[scalar])
    mpiGatherScalar=array(1,:)
#else
    !$GLC attributes unused :: self, scalar
    mpiGatherScalar=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGatherScalar

  function mpiGather1D(self,array)
    !!{
    Gather a 1-D array from all processes, returning it as a 2-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                                         :: self
    double precision           , intent(in   ), dimension(          :                ) :: array
    double precision           ,                dimension(size(array),self%countValue) :: mpiGather1D

#ifdef USEMPI
    mpiGather1D=self%requestData(self%allRanks,array)
#else
    !$GLC attributes unused :: self, array
    mpiGather1D=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGather1D

  function mpiGather2D(self,array)
    !!{
    Gather a 1-D array from all processes, returning it as a 2-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                                                                 :: self
    double precision           , intent(in   ), dimension(                :,                :                ) :: array
    double precision           ,                dimension(size(array,dim=1),size(array,dim=2),self%countValue) :: mpiGather2D

#ifdef USEMPI
    mpiGather2D=self%requestData(self%allRanks,array)
#else
    !$GLC attributes unused :: self, array
    mpiGather2D=0.0d0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGather2D

  function mpiGatherLogicalScalar(self,scalar)
    !!{
    Gather a logical scalar from all processes, returning it as a 1-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                :: self
    logical                    , intent(in   )                :: scalar
    logical                    , dimension(  self%countValue) :: mpiGatherLogicalScalar
#ifdef USEMPI
    logical                    , dimension(1,self%countValue) :: array
#endif

#ifdef USEMPI
    array=self%requestData(self%allRanks,[scalar])
    mpiGatherLogicalScalar=array(1,:)
#else
    !$GLC attributes unused :: self, scalar
    mpiGatherLogicalScalar=.false.
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGatherLogicalScalar

  function mpiGatherIntScalar(self,scalar)
    !!{
    Gather an integer scalar from all processes, returning it as a 1-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class           (mpiObject), intent(in   )                :: self
    integer                    , intent(in   )                :: scalar
    integer                    , dimension(  self%countValue) :: mpiGatherIntScalar
#ifdef USEMPI
    integer                    , dimension(1,self%countValue) :: array
#endif

#ifdef USEMPI
    array=self%requestData(self%allRanks,[scalar])
    mpiGatherIntScalar=array(1,:)
#else
    !$GLC attributes unused :: self, scalar
    mpiGatherIntScalar=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGatherIntScalar

  function mpiGatherInt1D(self,array)
    !!{
    Gather an integer 1-D array from all processes, returning it as a 2-D array.
    !!}
#ifndef USEMPI
    use Error, only : Error_Report
#endif
    implicit none
    class  (mpiObject), intent(in   )                                         :: self
    integer           , intent(in   ), dimension(          :                ) :: array
    integer           ,                dimension(size(array),self%countValue) :: mpiGatherInt1D

#ifdef USEMPI
    mpiGatherInt1D=self%requestData(self%allRanks,array)
#else
    !$GLC attributes unused :: self, array
    mpiGatherInt1D=0
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end function mpiGatherInt1D

  subroutine mpiCommunicatorPush(self,color)
    !!{
    Create a new communicator and push it onto the stack.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Comm_Split
#endif
    implicit none
    class  (mpiObject), intent(inout)               :: self
    integer           , intent(in   )               :: color
#ifdef USEMPI
    integer                                         :: iError
    type   (MPI_Comm ), allocatable  , dimension(:) :: communicatorStack
    
    self%iCommunicator=self%iCommunicator+1
    if (self%iCommunicator > size(self%communicatorStack)) then
       call move_alloc(self%communicatorStack,communicatorStack)
       allocate(self%communicatorStack(self%iCommunicator))
       self%communicatorStack(1:self%iCommunicator-1)=communicatorStack
       deallocate(communicatorStack)
    end if
    call MPI_Comm_Split(self%communicatorStack(self%iCommunicator-1),color,0,self%communicatorStack(self%iCommunicator),iError)
    if (iError /= 0) call Error_Report('failed to split communicator'//{introspection:location})
    self%communicator=self%communicatorStack(self%iCommunicator)
    call self%stateInitialize()
#else
    !$GLC attributes unused :: self, color
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiCommunicatorPush
  
  subroutine mpiCommunicatorPop(self)
    !!{
    Pop a communicator off of the stack and destroy it.
    !!}
    use :: Error  , only : Error_Report
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Comm_Free
#endif
    implicit none
    class  (mpiObject), intent(inout) :: self
#ifdef USEMPI
    integer                           :: iError
    
    if (self%iCommunicator == 1) call Error_Report('can not pop MPI_Comm_World off of the stack'//{introspection:location})
    call MPI_Comm_Free(self%communicator,iError)
    if (iError /= 0) call Error_Report('failed to free communicator'//{introspection:location})
    self%iCommunicator=self%iCommunicator-1
    self%communicator =self%communicatorStack(self%iCommunicator)
    call self%stateInitialize()
#else
    !$GLC attributes unused :: self
    call Error_Report('code was not compiled for MPI'//{introspection:location})
#endif
    return
  end subroutine mpiCommunicatorPop

  subroutine mpiStateInitialize(self)
    !!{
    Initialize the state (rank, group, node affinities, etc.) of the current MPI group.
    !!}
#ifdef USEMPI
    use :: MPI_F08           , only : MPI_Max_Processor_Name, MPI_Character
    use :: Error             , only : Error_Report
    use :: Hashes            , only : integerHash
    use :: ISO_Varying_String, only : assignment(=)         , operator(==), var_str, operator(//), char
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class    (mpiObject                 ), intent(inout)             :: self
#ifdef USEMPI
    character(len=MPI_Max_Processor_Name), dimension(1)              :: processorName
    character(len=MPI_Max_Processor_Name), dimension(:), allocatable :: processorNames
    type     (integerHash               )                            :: processCount
    type     (varying_string            )                            :: message
    integer                                                          :: i                  , iError  , &
         &                                                              processorNameLength, iProcess
    
    ! Determine ranks, counts, and processors.
    call    MPI_Comm_Size         (self%communicator,self%countValue  ,iError)
    if    (iError               /= 0                    ) call Error_Report('failed to determine MPI count'                                   //{introspection:location})
    call    MPI_Comm_Rank         (self%communicator,self% rankValue  ,iError)
    if    (iError               /= 0                    ) call Error_Report('failed to determine MPI rank'                                    //{introspection:location})
    call    MPI_Get_Processor_Name(processorName(1) ,processorNameLength ,iError)
    if    (iError               /= 0                    ) call Error_Report('failed to get MPI processor name'                                //{introspection:location})
    self%hostName=trim(processorName(1))
    call mpiBarrier()
    ! Construct an array containing all ranks.
    if (allocated(self%allRanks)) deallocate(self%allRanks)
    allocate(self%allRanks(0:self%countValue-1))
    forall(i=0:self%countValue-1)
       self%allRanks(i)=i
    end forall
    ! Get processor names from all processes.
    allocate(processorNames(0:self%countValue-1))
    call MPI_AllGather(processorName,MPI_Max_Processor_Name,MPI_Character,processorNames,MPI_Max_Processor_Name,MPI_Character,self%communicator,iError)
    ! Count processes per node.
    self%rankOnNodeValue =-1
    self%countOnNodeValue=+0
    call processCount%initialize()
    do iProcess=0,self%countValue-1
       if (processCount%exists(trim(processorNames(iProcess)))) then
          call processCount%set(trim(processorNames(iProcess)),processCount%value(trim(processorNames(iProcess)))+1)
       else
          call processCount%set(trim(processorNames(iProcess)),                                                   1)
       end if
       if (trim(processorNames(iProcess)) == trim(processorName(1))) then
          if (iProcess <= self%rankValue) &
               & self%rankOnNodeValue =self%rankOnNodeValue +1
          self       %countOnNodeValue=self%countOnNodeValue+1
       end if
    end do
    self%nodeCountValue=processCount%size()
    if (allocated(self%nodeAffinities)) deallocate(self%nodeAffinities)
    allocate(self%nodeAffinities(0:self%countValue-1))
    self%nodeAffinities=-1
    do iProcess=0,self%countValue-1
       do i=1,self%nodeCountValue
          if (trim(processorNames(iProcess)) == processCount%key(i)) self%nodeAffinities(iProcess)=i
       end do
       if (self%nodeAffinities(iProcess) < 0) then
          message=var_str('failed to determine node affinity for process ')//iProcess//' with processor name "'//processorNames(iProcess)//' - known processor names are:'
          do i=1,self%nodeCountValue
             message=message//char(10)//'   '//processCount%key(i)
          end do
          call Error_Report(message//{introspection:location})
       end if
    end do
    deallocate(processorNames)
#else
    !$GLC attributes unused :: self
#endif
    return
  end subroutine mpiStateInitialize
  
  function counterConstructor() result(self)
    !!{
    Constructor for MPI counter class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : C_Null_Ptr           , C_F_Pointer
#ifdef USEMPI
    use            :: Error        , only : Error_Report
    use            :: MPI_F08      , only : MPI_Win_Create       , MPI_Address_Kind, MPI_Info_Null      , MPI_Alloc_Mem, &
         &                                  MPI_TypeClass_Integer, MPI_SizeOf      , MPI_Type_Match_Size
#endif
    implicit none
    type   (mpiCounter)          :: self
#ifdef USEMPI
    integer                      :: mpiSize            , iError
    integer(c_size_t  ), pointer :: countInitialPointer
    class  (*         ), pointer :: dummyPointer_

    call MPI_SizeOf(0_c_size_t,mpiSize,iError)
    if (iError /= 0) call Error_Report('failed to get type size'//{introspection:location})
    call MPI_Type_Match_Size(MPI_TypeClass_Integer,mpiSize,self%typeClass,iError)
    if (iError /= 0) call Error_Report('failed to get type'     //{introspection:location})
    allocate(self%window)
    if (mpiSelf%rank() == 0) then
       ! The rank-0 process allocates space for the counter and creates its window.
       allocate(self%counter)
       call MPI_Alloc_Mem(int(mpiSize,kind=MPI_Address_Kind),MPI_Info_Null,self%counter%memory,iError)
       if (iError /= 0) call Error_Report('failed to allocate counter memory'//{introspection:location})
       countInitialPointer => null()
       call C_F_Pointer(self%counter%memory,countInitialPointer)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       !$ dummyPointer_       => self%counter
       !$ self%counterManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
       call MPI_Win_Create(countInitialPointer,int(mpiSize,kind=MPI_Address_Kind),mpiSize,MPI_Info_Null,mpiSelf%communicator,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to create RMA window'//{introspection:location})
       call mpiBarrier()
    else
       ! Other processes create a zero-size window.
       call MPI_Win_Create(C_Null_Ptr  ,               0_MPI_Address_Kind,mpiSize,MPI_Info_Null,mpiSelf%communicator,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to create RMA window'//{introspection:location})
       call mpiBarrier()
    end if
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    !$ dummyPointer_      => self%window
    !$ self%windowManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
#endif
    call self%reset()
    !$ self%ompLock_=ompLock()    
    return
  end function counterConstructor

  subroutine counterReset(self)
    !!{
    Reset an MPI counter.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Put     , MPI_Win_Unlock, MPI_Lock_Exclusive, MPI_Address_Kind, &
         &                 MPI_Win_Lock
#endif
    implicit none
    class  (mpiCounter), intent(inout) :: self
#ifdef USEMPI
    integer                            :: iError
    integer(c_size_t  ), dimension(1)  :: countInitial

    if (mpiSelf%rank() == 0) then
       ! The rank-0 process resets the counter.
       !$omp masked
       ! Reset the counter to zero.
       call MPI_Win_Lock(MPI_Lock_Exclusive,0,0,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to lock RMA window'  //{introspection:location})
       countInitial=0_c_size_t
       call MPI_Put(countInitial,1,self%typeClass,0,0_MPI_Address_Kind,1,self%typeClass,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to set MPI counter'  //{introspection:location})
       call MPI_Win_Unlock(0,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to unlock RMA window'//{introspection:location})
       !$omp end masked
    end if
    self%countThreadsWaiting=0_c_size_t
    self%countIncrementsHeld=0_c_size_t
    call mpiBarrier()
#else
    self%counter=0
#endif
    return
  end subroutine counterReset

  subroutine mpiWindowDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mpiWindow} class.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Win_Free
#endif
    implicit none
    type   (mpiWindow), intent(inout) :: self
#ifdef USEMPI
    integer                           :: iError

    call MPI_Win_Free(self%window,iError)
    if (iError /= 0) call Error_Report('failed to free RMA window'//{introspection:location})
#else
    !$GLC attributes unused :: self
#endif
    return
  end subroutine mpiWindowDestructor
  
  subroutine mpiMemoryDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mpiMemory} class.
    !!}
#ifdef USEMPI
    use :: MPI_F08, only : MPI_Free_Mem
#endif
    implicit none
    type   (mpiMemory), intent(inout) :: self
#ifdef USEMPI
    integer                           :: iError

    call MPI_Free_Mem(self%memory,iError)
#else
    !$GLC attributes unused :: self
#endif
    return
  end subroutine mpiMemoryDestructor

  function counterIncrement(self)
    !!{
    Increment an MPI counter.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Win_Lock    , MPI_Get_Accumulate, MPI_Win_Unlock, MPI_Lock_Shared, &
         &                 MPI_Address_Kind, MPI_Sum
#endif
    implicit none
    integer(c_size_t  )                :: counterIncrement
    class  (mpiCounter), intent(inout) :: self
#ifdef USEMPI
    integer(c_size_t  ), dimension(1)  :: counterIn          , counterOut
    integer(c_size_t  )                :: countThreadsWaiting
    integer                            :: iError
    
    ! Increment the count of the number of OpenMP threads waiting to obtain the lock.
    !$omp atomic
    self%countThreadsWaiting=self%countThreadsWaiting+1_c_size_t
    !$ call self%ompLock_%  set()
    ! If we are currently not holding any increments from the counter, we need to get some now.
    if (self%countIncrementsHeld == 0_c_size_t) then
       ! Begin a lock on the MPI shared window.
       call MPI_Win_Lock(MPI_Lock_Shared,0,0,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to lock RMA window'          //{introspection:location})
       ! Take a snapshot of the number of OpenMP threads that are currently waiting on the lock (plus the current thread). This is
       ! done without locking, so is subject to race conditions, but this does not matter - we just want a good estimate of the
       ! number of increments that we need to obtain.
       countThreadsWaiting=self%countThreadsWaiting
       counterIn          =countThreadsWaiting
       ! Increment the counter by the number of waiting OpenMP threads.
       call MPI_Get_Accumulate(counterIn,1,self%typeClass,counterOut,1,self%typeClass,0,0_MPI_Address_Kind,1,self%typeClass,MPI_Sum,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to accumulate to MPI counter'//{introspection:location})
       ! Unlock the MPI shared window to force synchronization to local variables.
       call MPI_Win_Unlock(0,self%window%window,iError)
       if (iError /= 0) call Error_Report('failed to unlock RMA window'        //{introspection:location})
       ! Record the current value of the counter for this MPI process, and the number of increments of it that we have held.
       self%countIncrementsHeld=countThreadsWaiting
       self%counterHeld        =counterOut(1)
    end if
    ! Return the current value of the counter on this process, increment the counter on this process, and decrement the number of
    ! held increments.
    counterIncrement        =self%counterHeld
    self%counterHeld        =self%counterHeld        +1_c_size_t
    self%countIncrementsHeld=self%countIncrementsHeld-1_c_size_t
    !$ call self%ompLock_%unset()
    ! Decrement the count of OpenMP threads waiting on the lock.
    !$omp atomic
    self%countThreadsWaiting=self%countThreadsWaiting-1_c_size_t
#else
    !$ call self%ompLock_%  set()
    counterIncrement=self%counter
    self%counter    =self%counter+1_c_size_t
    !$ call self%ompLock_%unset()
#endif
    return
  end function counterIncrement

  function counterDecrement(self)
    !!{
    Decrement an MPI counter.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Win_Lock    , MPI_Get_Accumulate, MPI_Win_Unlock, MPI_Lock_Exclusive, &
         &                 MPI_Address_Kind, MPI_Sum
#endif
    implicit none
    integer(c_size_t  )                :: counterDecrement
    class  (mpiCounter), intent(inout) :: self
#ifdef USEMPI
    integer(c_size_t  ), dimension(1)  :: counterIn       , counterOut
    integer                            :: iError

    counterIn=-1_c_size_t
    !$ call self%ompLock_%  set()
    call MPI_Win_Lock(MPI_Lock_Exclusive,0,0,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to lock RMA window'          //{introspection:location})
    call MPI_Get_Accumulate(counterIn,1,self%typeClass,counterOut,1,self%typeClass,0,0_MPI_Address_Kind,1,self%typeClass,MPI_Sum,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to accumulate to MPI counter'//{introspection:location})
    call MPI_Win_Unlock(0,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to unlock RMA window'        //{introspection:location})
    !$ call self%ompLock_%unset()
    counterDecrement=counterOut(1)
#else
    !$ call self%ompLock_%  set()
    counterDecrement=self%counter
    self%counter    =self%counter-1_c_size_t
    !$ call self%ompLock_%unset()
#endif
    return
  end function counterDecrement

  function counterGet(self)
    !!{
    Return the current value of an MPI counter.
    !!}
#ifdef USEMPI
    use :: Error  , only : Error_Report
    use :: MPI_F08, only : MPI_Win_Lock    , MPI_Get, MPI_Win_Unlock, MPI_Lock_Exclusive, &
    &                      MPI_Address_Kind
#endif
    implicit none
    integer(c_size_t  )                :: counterGet
    class  (mpiCounter), intent(inout) :: self
#ifdef USEMPI
    integer(c_size_t  ), dimension(1)  :: counterOut
    integer                            :: iError

    !$ call self%ompLock_%  set()
    call MPI_Win_Lock(MPI_Lock_Exclusive,0,0,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to lock RMA window'           //{introspection:location})
    call MPI_Get(counterOut,1,self%typeClass,0,0_MPI_Address_Kind,1,self%typeClass,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to get value from MPI counter'//{introspection:location})
    call MPI_Win_Unlock(0,self%window%window,iError)
    if (iError /= 0) call Error_Report('failed to unlock RMA window'         //{introspection:location})
    !$ call self%ompLock_%unset()
    counterGet=counterOut(1)-1_c_size_t
#else
    !$ call self%ompLock_%  set()
    counterGet=self%counter-1_c_size_t
    !$ call self%ompLock_%unset()
#endif
    return
  end function counterGet

end module MPI_Utilities
