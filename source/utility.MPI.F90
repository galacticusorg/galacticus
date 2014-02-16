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

!% Contains a module that implements useful MPI utilities.

module MPI_Utilities
  !% Implements useful MPI utilities.
  private
  public :: mpiInitialize, mpiFinalize, mpiBarrier, mpiSelf

  ! Define a type for interacting with MPI.
  type :: mpiObject
     private
     integer                            :: rankValue, countValue
     integer, allocatable, dimension(:) :: allRanks
   contains
     !@ <objectMethods>
     !@   <object>mpiObject</object>
     !@   <objectMethod>
     !@     <method>isMaster</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if this is the master process (i.e. rank-0 process).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rank</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the rank of this process.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>count</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the total number of processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rankLabel</method>
     !@     <type>\textcolor{red}{\textless type(varying\_string)\textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Return a label containing the rank of the process.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>requestData</method>
     !@     <type>\doubletwo</type>
     !@     <arguments>\intone requestFrom\argin, \doubleone array</arguments>
     !@     <description>Request the content of {\tt array} from each processes listed in {\tt requestFrom}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>messageWaiting</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\intzero\ [from]\argin, \intzero\ [tag]\argin</arguments>
     !@     <description>Return true if a message is waiting, optionally from the specified process and with the specified tag.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>average</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\doubleone array\argin</arguments>
     !@     <description>Return the average of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>sum</method>
     !@     <type>\intzero|\intone</type>
     !@     <arguments>(\intzero|\intone) array\argin</arguments>
     !@     <description>Return the sum of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>maxval</method>
     !@     <type>\doubleone</type>
     !@     <arguments>(\doublezero|\doubleone) array\argin</arguments>
     !@     <description>Return the maximum value of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>maxloc</method>
     !@     <type>\intone</type>
     !@     <arguments>\doubleone array\argin</arguments>
     !@     <description>Return the rank of the process with the maximum value of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>minval</method>
     !@     <type>\doubleone</type>
     !@     <arguments>(\doublezero|\doubleone) array\argin</arguments>
     !@     <description>Return the minimum value of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>minloc</method>
     !@     <type>\intone</type>
     !@     <arguments>\doubleone array\argin</arguments>
     !@     <description>Return the rank of the process with the minimum value of {\tt array} over all processes.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>gather</method>
     !@     <type>(\doubletwo|\doublethree)</type>
     !@     <arguments>(\doubleone|\doubletwo) array\argin</arguments>
     !@     <description>Gather arrays from all processes into an array of rank one higher.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: isMaster       => mpiIsMaster
     procedure :: rank           => mpiGetRank
     procedure :: rankLabel      => mpiGetRankLabel
     procedure :: count          => mpiGetCount
     procedure ::                   mpiRequestData1D , mpiRequestData2D
     generic   :: requestData    => mpiRequestData1D , mpiRequestData2D
     procedure :: messageWaiting => mpiMessageWaiting
     procedure ::                   mpiAverageScalar , mpiAverageArray
     generic   :: average        => mpiAverageScalar , mpiAverageArray
     procedure ::                   mpiSumScalarInt  , mpiSumArrayInt
     generic   :: sum            => mpiSumScalarInt  , mpiSumArrayInt
     procedure :: maxloc         => mpiMaxloc
     procedure ::                   mpiMaxvalScalar  , mpiMaxvalArray
     generic   :: maxval         => mpiMaxvalScalar  , mpiMaxvalArray
     procedure :: minloc         => mpiMinloc
     procedure ::                   mpiMinvalScalar  , mpiMinvalArray
     generic   :: minval         => mpiMinvalScalar  , mpiMinvalArray
     procedure ::                   mpiGather1D      , mpiGather2D
     generic   :: gather         => mpiGather1D      , mpiGather2D
  end type mpiObject

  ! Declare an object for interaction with MPI.
  type(mpiObject) :: mpiSelf

  ! Tags.
  integer, parameter :: tagRequestForData=1, tagState=2

contains

  subroutine mpiInitialize()
    !% Initialize MPI.
    use MPI
    use Memory_Management
    use Galacticus_Error
    implicit none
    integer :: i, iError
    
    call MPI_Init(iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiInitialize','failed to initialize MPI'     )
    call MPI_Comm_Size(MPI_Comm_World,mpiSelf%countValue,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiInitialize','failed to determine MPI count')
    call MPI_Comm_Rank(MPI_Comm_World,mpiSelf% rankValue,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiInitialize','failed to determine MPI rank' )
    ! Construct an array containing all ranks.
    call Alloc_Array(mpiSelf%allRanks,[mpiSelf%countValue],[0])
    forall(i=0:mpiSelf%countValue-1)
       mpiSelf%allRanks(i)=i
    end forall
    return
  end subroutine mpiInitialize
  
  subroutine mpiFinalize()
    !% Finalize MPI.
    use MPI
    use Galacticus_Error
    implicit none
    integer :: iError
    
    call MPI_Finalize(iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiFinalize','failed to finalize MPI')
    return
  end subroutine mpiFinalize
  
  subroutine mpiBarrier()
    !% Block until all MPI processes are synchronized.
    use MPI
    use Galacticus_Error
    implicit none
    integer :: iError
    
    call MPI_Barrier(MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiInitialize','MPI barrier failed')
    return
  end subroutine mpiBarrier
  
  logical function mpiIsMaster(self)
    !% Return true if this is the master process.
    implicit none
    class(mpiObject), intent(in   ) :: self
    
    mpiIsMaster=(self%rank() == 0)
    return
  end function mpiIsMaster
  
  integer function mpiGetRank(self)
    !% Return MPI rank.
    implicit none
    class(mpiObject), intent(in   ) :: self
    
    mpiGetRank=self%rankValue
    return
  end function mpiGetRank
  
  function mpiGetRankLabel(self)
    !% Return MPI rank label.
    use ISO_Varying_String
    implicit none
    type     (varying_string)                :: mpiGetRankLabel
    class    (mpiObject     ), intent(in   ) :: self
    character(len=4         )                :: label

    write (label,'(i4.4)') self%rankValue
    mpiGetRankLabel=label
    return
  end function mpiGetRankLabel
  
  integer function mpiGetCount(self)
    !% Return MPI count.
    implicit none
    class(mpiObject), intent(in   ) :: self
    
    mpiGetCount=self%countValue
    return
  end function mpiGetCount
  
  logical function mpiMessageWaiting(self,from,tag)
    !% Return true if an MPI message (matching the optional {\tt from} and {\tt tag} if given) is waiting for receipt.
    use MPI
    use Galacticus_Error
    implicit none
    class  (mpiObject), intent(in   )                        :: self
    integer           , intent(in   )             , optional :: from         , tag
    integer           , dimension(MPI_Status_Size)           :: messageStatus
    integer                                                  :: fromActual   , tagActual, iError
   
    fromActual=MPI_Any_Source
    tagActual =MPI_Any_Tag
    if (present(from)) fromActual=from
    if (present(tag ) ) tagActual=tag
    call MPI_IProbe(fromActual,tagActual,MPI_Comm_World,mpiMessageWaiting,messageStatus,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiMessageWaiting','failed to probe for waiting messages')
    return
  end function mpiMessageWaiting

  function mpiRequestData1D(self,requestFrom,array)
    !% Request and receive data from other MPI processes.
    use MPI
    implicit none
    class           (mpiObject), intent(in   )                                           :: self
    integer                    , intent(in   ), dimension(                            :) :: requestFrom
    double precision           , intent(in   ), dimension(          :                  ) :: array
    double precision                          , dimension(size(array),size(requestFrom)) :: mpiRequestData1D
    integer                                   , dimension(                            1) :: requester
    integer                                   , dimension(              MPI_Status_Size) :: messageStatus
    integer                                                                              :: i             , iError
    
    ! Record our own rank as the requester.
    requester=mpiSelf%rank()
    ! Send requests.
    do i=1,size(requestFrom)
       call MPI_Send(requester,1,MPI_Integer,requestFrom(i),tagRequestForData,MPI_Comm_World,iError)
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    do while (mpiSelf%messageWaiting(tag=tagRequestForData))
       ! Receive the request.
       call MPI_Recv(requester,1,MPI_Integer,MPI_Any_Source,tagRequestForData,MPI_Comm_World,messageStatus,iError)
       ! Send our data in reply.
       call MPI_Send(array,size(array),MPI_Double_Precision,requester(1),tagState,MPI_Comm_World,iError)
    end do
    call mpiBarrier()
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(mpiRequestData1D(:,i),size(array),MPI_Double_Precision,requestFrom(i),tagState,MPI_Comm_World,messageStatus,iError)
    end do
    return
  end function mpiRequestData1D

  function mpiRequestData2D(self,requestFrom,array)
    !% Request and receive data from other MPI processes.
    use MPI
    implicit none
    class           (mpiObject), intent(in   )                                                                   :: self
    integer                    , intent(in   ), dimension(                                                    :) :: requestFrom
    double precision           , intent(in   ), dimension(                :,                :                  ) :: array
    double precision                          , dimension(size(array,dim=1),size(array,dim=2),size(requestFrom)) :: mpiRequestData2D
    integer                                   , dimension(                                                    1) :: requester
    integer                                   , dimension(                                      MPI_Status_Size) :: messageStatus
    integer                                                                                                      :: i             , iError
    
    ! Record our own rank as the requester.
    requester=mpiSelf%rank()
    ! Send requests.
    do i=1,size(requestFrom)
       call MPI_Send(requester,1,MPI_Integer,requestFrom(i),tagRequestForData,MPI_Comm_World,iError)
    end do
    call mpiBarrier()
    ! Check for waiting requests.
    do while (mpiSelf%messageWaiting(tag=tagRequestForData))
       ! Receive the request.
       call MPI_Recv(requester,1,MPI_Integer,MPI_Any_Source,tagRequestForData,MPI_Comm_World,messageStatus,iError)
       ! Send our data in reply.
       call MPI_Send(array,product(shape(array)),MPI_Double_Precision,requester(1),tagState,MPI_Comm_World,iError)
    end do
    call mpiBarrier()
    ! Receive data.
    do i=1,size(requestFrom)
       call MPI_Recv(mpiRequestData2D(:,:,i),product(shape(array)),MPI_Double_Precision,requestFrom(i),tagState,MPI_Comm_World,messageStatus,iError)
    end do
    return
  end function mpiRequestData2D

  function mpiSumArrayInt(self,array,mask)
    !% Sum an integer array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class  (mpiObject), intent(in   )                                    :: self
    integer           , intent(in   ), dimension( :          )           :: array
    logical           , intent(in   ), dimension(0:          ), optional :: mask
    integer                          , dimension(size(array))            :: mpiSumArrayInt, maskedArray
    integer                                                              :: iError        , activeCount

    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiSumArrayInt,size(array),MPI_Integer,MPI_Sum,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiSumArrayInt','MPI all reduce failed')
    return
  end function mpiSumArrayInt

  integer function mpiSumScalarInt(self,scalar,mask)
    !% Sum an integer scalar over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class  (mpiObject), intent(in   )                         :: self
    integer           , intent(in   )                         :: scalar
    logical           , intent(in   ), dimension(:), optional :: mask
    integer                          , dimension(1)           :: array

    array=self%sum([scalar],mask)
    mpiSumScalarInt=array(1)
    return
  end function mpiSumScalarInt

  function mpiAverageArray(self,array,mask)
    !% Average an array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiAverageArray, maskedArray
    integer                                                                       :: iError         , activeCount

    ! Sum the array over all processes.
    maskedArray=array
    activeCount=self%count()
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=0.0d0
       activeCount=count(mask)
    end if
    call MPI_AllReduce(maskedArray,mpiAverageArray,size(array),MPI_Double_Precision,MPI_Sum,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiAverageArray','MPI all reduce failed')
    ! Convert the sum into an average.
    mpiAverageArray=mpiAverageArray/dble(activeCount)
    return
  end function mpiAverageArray

  double precision function mpiAverageScalar(self,scalar,mask)
    !% Find the maximum values of a scalar over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
    double precision                          , dimension(1)           :: array

    array=self%average([scalar],mask)
    mpiAverageScalar=array(1)
    return
  end function mpiAverageScalar

  function mpiMaxvalArray(self,array,mask)
    !% Find the maximum values of an array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiMaxvalArray, maskedArray
    integer                                                                       :: iError

    ! Find the maximum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-HUGE(1.0d0)
    end if
    call MPI_AllReduce(maskedArray,mpiMaxvalArray,size(array),MPI_Double_Precision,MPI_Max,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiMaxvalArray','MPI all reduce failed')
    return
  end function mpiMaxvalArray

  double precision function mpiMaxvalScalar(self,scalar,mask)
    !% Find the maximum values of a scalar over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
    double precision                          , dimension(1)           :: array

    array=self%maxval([scalar],mask)
    mpiMaxvalScalar=array(1)
    return
  end function mpiMaxvalScalar

  function mpiMaxloc(self,array,mask)
    !% Find the rank of the process having maximum values of an array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                      :: self
    double precision           , intent(in   ), dimension( :            )           :: array
    logical                    , intent(in   ), dimension(0:            ), optional :: mask
    integer                                   , dimension(   size(array))           :: mpiMaxloc
    double precision                          , dimension(2 ,size(array))           :: arrayIn  , arrayOut
    integer                                                                         :: iError

    ! Find the maximum over all processes.
    arrayIn(1,:)=array
    if (present(mask)) then
       if (.not.mask(self%rank())) arrayIn(1,:)=-HUGE(1.0d0)
    end if
    arrayIn(2,:)=self%rank()
    call MPI_AllReduce(arrayIn,arrayOut,size(array),MPI_2Double_Precision,MPI_MaxLoc,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiMax','MPI all reduce failed')
    mpiMaxloc=arrayOut(2,:)
    return
  end function mpiMaxloc

  function mpiMinvalArray(self,array,mask)
    !% Find the minimum values of an array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                    :: self
    double precision           , intent(in   ), dimension( :          )           :: array
    logical                    , intent(in   ), dimension(0:          ), optional :: mask
    double precision                          , dimension(size(array) )           :: mpiMinvalArray, maskedArray
    integer                                                                       :: iError

    ! Find the minimum over all processes.
    maskedArray=array
    if (present(mask)) then
       if (.not.mask(self%rank())) maskedArray=-HUGE(1.0d0)
    end if
    call MPI_AllReduce(maskedArray,mpiMinvalArray,size(array),MPI_Double_Precision,MPI_Min,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiMinvalArray','MPI all reduce failed')
    return
  end function mpiMinvalArray

  double precision function mpiMinvalScalar(self,scalar,mask)
    !% Find the minimum values of a scalar over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                         :: self
    double precision           , intent(in   )                         :: scalar
    logical                    , intent(in   ), dimension(:), optional :: mask
    double precision                          , dimension(1)           :: array

    array=self%minval([scalar],mask)
    mpiMinvalScalar=array(1)
    return
  end function mpiMinvalScalar

  function mpiMinloc(self,array,mask)
    !% Find the rank of the process having minimum values of an array over all processes, returning it to all processes.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                      :: self
    double precision           , intent(in   ), dimension( :            )           :: array
    logical                    , intent(in   ), dimension(0:            ), optional :: mask
    integer                                   , dimension(   size(array))           :: mpiMinloc
    double precision                          , dimension(2 ,size(array))           :: arrayIn  , arrayOut
    integer                                                                         :: iError

    ! Find the minimum over all processes.
    arrayIn(1,:)=array
    if (present(mask)) then
       if (.not.mask(self%rank())) arrayIn(1,:)=-HUGE(1.0d0)
    end if
    arrayIn(2,:)=self%rank()
    call MPI_AllReduce(arrayIn,arrayOut,size(array),MPI_2Double_Precision,MPI_MinLoc,MPI_Comm_World,iError)
    if (iError /= 0) call Galacticus_Error_Report('mpiMin','MPI all reduce failed')
    mpiMinloc=arrayOut(2,:)
    return
  end function mpiMinloc

  function mpiGather1D(self,array)
    !% Gather a 1-D array from all processes, returning it as a 2-D array.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                         :: self
    double precision           , intent(in   ), dimension(          :                ) :: array
    double precision           ,                dimension(size(array),self%countValue) :: mpiGather1D
    integer                                                                            :: iError

    mpiGather1D=self%requestData(self%allRanks,array)
    return
  end function mpiGather1D

  function mpiGather2D(self,array)
    !% Gather a 1-D array from all processes, returning it as a 2-D array.
    use MPI
    use Galacticus_Error
    implicit none
    class           (mpiObject), intent(in   )                                                                 :: self
    double precision           , intent(in   ), dimension(                :,                :                ) :: array
    double precision           ,                dimension(size(array,dim=1),size(array,dim=2),self%countValue) :: mpiGather2D
    integer                                                                                                    :: iError

    mpiGather2D=self%requestData(self%allRanks,array)
    return
  end function mpiGather2D

end module MPI_Utilities
