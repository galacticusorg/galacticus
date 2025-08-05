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

  use, intrinsic :: ISO_C_Binding                  , only : c_size_t
  use            :: Multi_Counters                 , only : multiCounter
  use            :: Radiative_Transfer_Convergences, only : radiativeTransferConvergenceClass
  use            :: Radiative_Transfer_Matters     , only : radiativeTransferMatterClass     , radiativeTransferPropertiesMatter

  type :: cartesian3DBoundaries
     !!{
     Type used to store boundaries of computational domain cells for 3D Cartesian domains.
     !!}
     private
     double precision, allocatable, dimension(:) :: boundary
  end type cartesian3DBoundaries
  
  !![
  <computationalDomain name="computationalDomainCartesian3D">
   <description>A computational domain using a 3D Cartesian grid.</description>
  </computationalDomain>
  !!]
  type, extends(computationalDomainClass) :: computationalDomainCartesian3D
     !!{
     Implementation of a computational domain using a 3D Cartesian grid.
     !!}
     private
     double precision                                                , dimension(  2  ) :: xBoundaries                            , yBoundaries               , &
          &                                                                                zBoundaries
     double precision                                                , dimension(3,2  ) :: boundaries
     integer         (c_size_t                         )             , dimension(3    ) :: countCells
     double precision                                                                   :: convergencePercentile                  , convergenceMeasurePrevious, &
          &                                                                                convergenceThreshold                   , convergenceRatioThreshold
     integer         (c_size_t                         )                                :: countCellsConvergence
     integer         (c_size_t                         ), allocatable, dimension(    :) :: sliceMinimum                           , sliceMaximum
     class           (radiativeTransferPropertiesMatter), allocatable, dimension(:,:,:) :: properties
     class           (radiativeTransferMatterClass     ), pointer                       :: radiativeTransferMatter_      => null()
     class           (radiativeTransferConvergenceClass), pointer                       :: radiativeTransferConvergence_ => null()
     type            (cartesian3DBoundaries            )             , dimension(3    ) :: boundariesCells
   contains
     final     ::                             cartesian3DDestructor
     procedure :: iterator                 => cartesian3DIterator
     procedure :: initialize               => cartesian3DInitialize
     procedure :: reset                    => cartesian3DReset
     procedure :: indicesFromPosition      => cartesian3DIndicesFromPosition
     procedure :: lengthToCellBoundary     => cartesian3DLengthToCellBoundary
     procedure :: absorptionCoefficient    => cartesian3DAbsorptionCoefficient
     procedure :: interactWithPhotonPacket => cartesian3DInteractWithPhotonPacket
     procedure :: accumulatePhotonPacket   => cartesian3DAccumulatePhotonPacket
     procedure :: stateSolve               => cartesian3DStateSolve
     procedure :: converged                => cartesian3DConverged
     procedure :: output                   => cartesian3DOutput
  end type computationalDomainCartesian3D

  interface computationalDomainCartesian3D
     !!{
     Constructors for the \refClass{computationalDomainCartesian3D} computational domain.
     !!}
     module procedure cartesian3DConstructorParameters
     module procedure cartesian3DConstructorInternal
  end interface computationalDomainCartesian3D

  type, extends(domainIterator) :: domainIteratorCartesian3D
     !!{
     An interactor for 3D Cartesian computational domains.
     !!}
     private
     type(multiCounter) :: counter_
   contains
     !![
     <methods>
       <method description="Move to the next cell in the domain." method="next" />
     </methods>
     !!]
     procedure :: next => domainIteratorCartesian3DNext
  end type domainIteratorCartesian3D
  
contains

  function cartesian3DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{computationalDomainCartesian3D} computational domain class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainCartesian3D   )                 :: self
    type            (inputParameters                  ), intent(inout)  :: parameters
    double precision                                   , dimension(  2) :: xBoundaries                  , yBoundaries         , &
         &                                                                 zBoundaries
    double precision                                   , dimension(3,2) :: boundaries
    integer         (c_size_t                         ), dimension(3  ) :: countCells
    class           (radiativeTransferMatterClass     ), pointer        :: radiativeTransferMatter_
    class           (radiativeTransferConvergenceClass), pointer        :: radiativeTransferConvergence_
    double precision                                                    :: convergencePercentile        , convergenceThreshold, &
         &                                                                 convergenceRatioThreshold
    
    !![
    <inputParameter>
      <name>xBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $x$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>yBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $y$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>zBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $z$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countCells</name>
      <defaultValue>[3_c_size_t,3_c_size_t,3_c_size_t]</defaultValue>
      <description>The number of cells in the domain in each dimension.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>convergencePercentile</name>
      <defaultValue>0.99d0</defaultValue>
      <description>The percentile used in the convergence criterion.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>convergenceThreshold</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The threshold for the convergence measure.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>convergenceRatioThreshold</name>
      <defaultValue>1.1d0</defaultValue>
      <description>The threshold for the change in convergence criterion.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="radiativeTransferMatter"      name="radiativeTransferMatter_"      source="parameters"/>
    <objectBuilder class="radiativeTransferConvergence" name="radiativeTransferConvergence_" source="parameters"/>
    !!]
    boundaries(1,:)=xBoundaries
    boundaries(2,:)=yBoundaries
    boundaries(3,:)=zBoundaries
    self=computationalDomainCartesian3D(boundaries,countCells,convergencePercentile,convergenceThreshold,convergenceRatioThreshold,radiativeTransferMatter_,radiativeTransferConvergence_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="radiativeTransferMatter_"     />
    <objectDestructor name="radiativeTransferConvergence_"/>
    !!]
    return
  end function cartesian3DConstructorParameters

  function cartesian3DConstructorInternal(boundaries,countCells,convergencePercentile,convergenceThreshold,convergenceRatioThreshold,radiativeTransferMatter_,radiativeTransferConvergence_) result(self)
    !!{
    Constructor for the \refClass{computationalDomainCartesian3D} computational domain class which takes a parameter set as input.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLinear
    implicit none
    type            (computationalDomainCartesian3D   )                                :: self
    double precision                                   , dimension(3,2), intent(in   ) :: boundaries
    integer         (c_size_t                         ), dimension(3  ), intent(in   ) :: countCells
    double precision                                                   , intent(in   ) :: convergencePercentile        , convergenceThreshold, &
         &                                                                                convergenceRatioThreshold
    class           (radiativeTransferMatterClass     ), target        , intent(in   ) :: radiativeTransferMatter_
    class           (radiativeTransferConvergenceClass), target        , intent(in   ) :: radiativeTransferConvergence_
    integer         (c_size_t                         )                                :: i
    !![
    <constructorAssign variables="boundaries, countCells, convergencePercentile, convergenceThreshold, convergenceRatioThreshold, *radiativeTransferMatter_, *radiativeTransferConvergence_"/>
    !!]
    
    ! Store boundaries for each axis.
    self%xBoundaries=boundaries(1,:)
    self%yBoundaries=boundaries(2,:)
    self%zBoundaries=boundaries(3,:)
    ! Construct cell boundaries.
    do i=1,3
       allocate(self%boundariesCells(i)%boundary(countCells(i)+1))
       self%boundariesCells(i)%boundary=Make_Range(boundaries(i,1),boundaries(i,2),int(countCells(i)+1),rangeTypeLinear)
    end do
     ! Compute the number of cells used in convergence criteria.
    self%countCellsConvergence     =min(int(dble(product(self%countCells))*(1.0d0-convergencePercentile),c_size_t)+1,product(self%countCells))
    self%convergenceMeasurePrevious=huge(0.0d0)
    return
  end function cartesian3DConstructorInternal

  subroutine cartesian3DDestructor(self)
    !!{
    Destructor for the \refClass{computationalDomainCartesian3D} computational domain class.
    !!}
    implicit none
    type(computationalDomainCartesian3D), intent(inout) :: self

    !![
    <objectDestructor name="self%radiativeTransferMatter_"     />
    <objectDestructor name="self%radiativeTransferConvergence_"/>
    !!]
    return
  end subroutine cartesian3DDestructor

  subroutine cartesian3DInitialize(self)
    !!{
    Initialize the computational domain.
    !!}
    use :: Computational_Domain_Volume_Integrators, only : computationalDomainVolumeIntegratorCartesian3D
    use :: Display                                , only : displayCounter                                , displayCounterClear  , displayIndent, displayUnindent, &
          &                                                verbosityLevelStandard                        , verbosityLevelWorking
    use :: MPI_Utilities                          , only : mpiBarrier                                    , mpiSelf
    use :: Timers                                 , only : timer
    implicit none
    class           (computationalDomainCartesian3D                ), intent(inout)  :: self
    class           (radiativeTransferPropertiesMatter             ), allocatable    :: properties
    integer         (c_size_t                                      )                 :: i             , j               , &
         &                                                                              k             , slicesPerProcess, &
         &                                                                              slicesExtra
    type            (computationalDomainVolumeIntegratorCartesian3D), allocatable    :: integrator
    double precision                                                , dimension(3,2) :: boundariesCell
#ifdef USEMPI
    integer                                                                          :: p
#endif
    type            (timer                                         )                 :: timer_

    ! Establish a timer.
    timer_=timer()
    ! Divide up the domain between MPI processes.
    allocate(self%sliceMinimum(0:mpiSelf%count()-1))
    allocate(self%sliceMaximum(0:mpiSelf%count()-1))
    slicesPerProcess=self%countCells(3)/mpiSelf%count()
    slicesExtra     =self%countCells(3)-mpiSelf%count()*slicesPerProcess
    do i=0,mpiSelf%count()-1
       if (i == 0) then
          self%sliceMinimum(i)=1
       else
          self%sliceMinimum(i)=self%slicemaximum(i-1)+1
       end if
       self%sliceMaximum(i)=self%sliceMinimum(i)+slicesPerProcess-1
       if (slicesExtra > 0) then
          self%sliceMaximum(i)=self%sliceMaximum(i)+1
          slicesExtra=slicesExtra-1
       end if
    end do
    ! Construct the domain.
    if (mpiSelf%isMaster()) call displayIndent('populating computational domain',verbosityLevelStandard)
    call timer_%start                                 (          )
    call self  %radiativeTransferMatter_%propertyClass(properties)
    allocate(self%properties(self%countCells(1),self%countCells(2),self%countCells(3)),mold=properties)
    deallocate(properties)
    do i    =1,self%countCells(1)
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(i-1_c_size_t)/dble(self%countCells(1))),isNew=i==1,verbosity=verbosityLevelWorking)
       boundariesCell      (1,:)=self%boundariesCells(1)%boundary(i:i+1)
       do j   =1,self%countCells(2)
          boundariesCell   (2,:)=self%boundariesCells(2)%boundary(j:j+1)
          do k=1,self%countCells(3)
             boundariesCell(3,:)=self%boundariesCells(3)%boundary(k:k+1)
             ! Build a volume integrator for this cell.
             allocate(integrator)
             integrator=computationalDomainVolumeIntegratorCartesian3D(boundariesCell)
             ! Populate this cell.
             call self%radiativeTransferMatter_%populateDomain(self%properties(i,j,k),integrator,onProcess=k >= self%sliceMinimum(mpiSelf%rank()) .and. k <= self%sliceMaximum(mpiSelf%rank()))
             ! Destroy the integrator.
             deallocate(integrator)
          end do
       end do
    end do
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayCounterClear(                                   verbosityLevelWorking )
    if (mpiSelf%isMaster()) call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
#ifdef USEMPI
    ! Broadcast populated domain cells.
    if (mpiSelf%isMaster()) call displayIndent('broadcasting computational domain',verbosityLevelStandard)
    call timer_%start                                 (          )
    do p=0,mpiSelf%count()-1
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(p)/dble(mpiSelf%count())),isNew=p==1,verbosity=verbosityLevelWorking)
       do k      =self%sliceMinimum(p),self%sliceMaximum(p)
          do j   =1                   ,self%countCells  (2)
             do i=1                   ,self%countCells  (1)
                call self%radiativeTransferMatter_%broadcastDomain(p,self%properties(i,j,k))
                call mpiBarrier()
             end do
          end do
       end do
    end do
    call mpiBarrier()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayCounterClear(                                   verbosityLevelWorking )
    if (mpiSelf%isMaster()) call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
#endif
    return
  end subroutine cartesian3DInitialize

  subroutine cartesian3DIterator(self,iterator)
    !!{
    Construct an iterator for this domain,
    !!}
    implicit none
    class(computationalDomainCartesian3D), intent(inout)              :: self
    class(domainIterator                ), intent(inout), allocatable :: iterator

    ! Build a multi-counter object for use in iterating over domain cells.
    allocate(domainIteratorCartesian3D :: iterator)
    select type (iterator)
    type is (domainIteratorCartesian3D)
       iterator%counter_=multiCounter(self%countCells)
    end select
    return
  end subroutine cartesian3DIterator

  logical function domainIteratorCartesian3DNext(self)
    implicit none
    class(domainIteratorCartesian3D), intent(inout) :: self

    domainIteratorCartesian3DNext=self%counter_%increment()
    return
  end function domainIteratorCartesian3DNext
  
  subroutine cartesian3DReset(self)
    !!{
    Reset the computational domain prior to a new iteration.
    !!}
    implicit none
    class  (computationalDomainCartesian3D), intent(inout) :: self
    integer(c_size_t                      )                :: i   , j, &
         &                                                    k

    do i=1,self%countCells(1)
       do j=1,self%countCells(2)
          do k=1,self%countCells(3)
             call self%radiativeTransferMatter_%reset(self%properties(i,j,k))
          end do
       end do
    end do
    return
  end subroutine cartesian3DReset

  subroutine cartesian3DIndicesFromPosition(self,position,indices)
    !!{
    Determine the indices of the cell containing the given point.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (computationalDomainCartesian3D), intent(inout)                            :: self
    double precision                                             , dimension(3), intent(in   ) :: position
    integer         (c_size_t                      ), allocatable, dimension(:), intent(inout) :: indices
    integer                                                                                    :: i
    
    ! Allocate indices to the correct size if necessary.
    if (allocated(indices)) then
       if (size(indices) /= 3) then
          deallocate(indices   )
          allocate  (indices(3))
       end if
    else
       allocate     (indices(3))
    end if
    ! Determine indices.
    do i=1,3
       if (position(i) < self%boundariesCells(i)%boundary(1) .or. position(i) >= self%boundariesCells(i)%boundary(self%countCells(i)+1)) then
          indices=-huge(0_c_size_t)
          exit
       else
          indices(i)=searchArray(self%boundariesCells(i)%boundary,position(i))
       end if
    end do
    return
  end subroutine cartesian3DIndicesFromPosition

  double precision function cartesian3DAbsorptionCoefficient(self,photonPacket,indices)
    !!{
    Return the absorption coefficient for the given photon packet in the given domain cell.
    !!}
    implicit none
    class  (computationalDomainCartesian3D    )              , intent(inout) :: self
    class  (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer(c_size_t                          ), dimension(:), intent(in   ) :: indices

    cartesian3DAbsorptionCoefficient=self%radiativeTransferMatter_%absorptionCoefficient(                              &
         &                                                                               self%properties  (            &
         &                                                                                                 indices(1), &
         &                                                                                                 indices(2), &
         &                                                                                                 indices(3)  &
         &                                                                                                )          , &
         &                                                                                    photonPacket             &
         &                                                                              )
    return
  end function cartesian3DAbsorptionCoefficient
  
  double precision function cartesian3DLengthToCellBoundary(self,photonPacket,indices,indicesNeighbor,positionBoundary)
    !!{
    Return the length to the first domain cell boundary intersected by the given photon packet.
    !!}
    implicit none
    class           (computationalDomainCartesian3D    )                           , intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass)                           , intent(inout) :: photonPacket
    integer         (c_size_t                          )             , dimension(:), intent(in   ) :: indices
    integer         (c_size_t                          ), allocatable, dimension(:), intent(inout) :: indicesNeighbor
    double precision                                                 , dimension(3), intent(  out) :: positionBoundary
    double precision                                                 , dimension(3)                :: position        , direction
    integer                                                                                        :: i               , j
    double precision                                                                               :: lengthToFace

    ! Allocate indices to the correct size if necessary.
    if (allocated(indicesNeighbor)) then
       if (size(indicesNeighbor) /= 3) then
          deallocate(indicesNeighbor   )
          allocate  (indicesNeighbor(3))
       end if
    else
       allocate     (indicesNeighbor(3))
    end if
    ! Consider all six faces of the cell, find the distance to the face, and look for the smallest non-negative distance.
    cartesian3DLengthToCellBoundary=huge(0.0d0)
    position                       =photonPacket%position ()
    direction                      =photonPacket%direction()
    do i=1,3    ! Axes.
       if (direction(i) == 0.0d0) cycle
       do j=0,1 ! Faces.
          lengthToFace=+(                                                &
               &         +self%boundariesCells(i)%boundary(indices(i)+j) &
               &         -                        position(        i   ) &
               &        )                                                &
               &       /                          direction(       i   )
          if     (                                                     &
               &   (                                                   &
               &       lengthToFace >  0.0d0                           & ! Distance to face is positive (i.e. in direction of travel).
               &    .or.                                               &
               &     (                                                 &
               &       lengthToFace == 0.0d0                           & ! Photon packet is precisely on the lower face of the cell, and is moving in the negative direction
               &      .and.                                            & ! - in this case we want to transition it to the cell below.
               &       direction(i) <  0.0d0                           &
               &      .and.                                            &
               &       j            == 0                               &
               &     )                                                 &
               &   )                                                   &
               &  .and.                                                &
               &       lengthToFace <  cartesian3DLengthToCellBoundary & ! This is the shortest distance to a face we've found.
               & ) then
             cartesian3DLengthToCellBoundary=lengthToFace
             indicesNeighbor    =+indices
             indicesNeighbor (i)=+indicesNeighbor(i)+(2*j-1)
             positionBoundary   =+position                   &
                  &              +direction                  &
                  &              *lengthToFace
             positionBoundary(i)=+self%boundariesCells(i)%boundary(indices(i)+j)
             if (indicesNeighbor(i) < 1 .or. indicesNeighbor(i) > self%countCells(i)) indicesNeighbor(i)=-huge(0_c_size_t)
          end if
       end do
    end do
    return
  end function cartesian3DLengthToCellBoundary

  subroutine cartesian3DAccumulatePhotonPacket(self,photonPacket,indices,absorptionCoefficient,lengthTraversed)
    !!{
    Accumulate ``absorptions'' from the photon packet as it traverses a cell of the computational domain.
    !!}
    implicit none
    class           (computationalDomainCartesian3D    )              , intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer         (c_size_t                          ), dimension(:), intent(in   ) :: indices
    double precision                                                  , intent(in   ) :: absorptionCoefficient, lengthTraversed

    call self%radiativeTransferMatter_%accumulatePhotonPacket(                                       &
         &                                                    self%properties           (            &
         &                                                                               indices(1), &
         &                                                                               indices(2), &
         &                                                                               indices(3)  &
         &                                                                              )          , &
         &                                                         photonPacket                    , &
         &                                                         absorptionCoefficient           , &
         &                                                         lengthTraversed                   &
         &                                                   )
    return
  end subroutine cartesian3DAccumulatePhotonPacket

  logical function cartesian3DInteractWithPhotonPacket(self,photonPacket,indices)
    !!{
    Allow matter in a domain cell to interact with the photon packet.
    !!}
    implicit none
    class  (computationalDomainCartesian3D    )              , intent(inout) :: self
    class  (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer(c_size_t                          ), dimension(:), intent(inout) :: indices

    cartesian3DInteractWithPhotonPacket=self%radiativeTransferMatter_%interactWithPhotonPacket(                              &
         &                                                                                     self%properties  (            &
         &                                                                                                       indices(1), &
         &                                                                                                       indices(2), &
         &                                                                                                       indices(3)  &
         &                                                                                                      )          , &
         &                                                                                          photonPacket             &
         &                                                                                    )
    return
  end function cartesian3DInteractWithPhotonPacket
  
  subroutine cartesian3DStateSolve(self)
    !!{
    Solve for the state of matter in the computational domain.
    !!}
    use :: Display         , only : displayCounter    , displayCounterClear   , displayIndent        , displayMessage, &
          &                         displayUnindent   , verbosityLevelStandard, verbosityLevelWorking
    use :: Error, only : errorStatusSuccess
    use :: MPI_Utilities   , only : mpiBarrier        , mpiSelf
    use :: Timers          , only : timer
    implicit none
    class    (computationalDomainCartesian3D), intent(inout) :: self
    integer  (c_size_t                      )                :: i      , j        , &
         &                                                      k      , countFail
    integer                                                  :: status
#ifdef USEMPI
    integer                                                  :: p
#endif
    type     (timer                         )                :: timer_
    type     (varying_string                )                :: message
    character(len=12                        )                :: label

    ! Establish a timer.
    timer_=timer()
#ifdef USEMPI
    if (mpiSelf%isMaster()) call displayIndent  ('accumulating absorptions across processes',verbosityLevelStandard)
    call timer_%start()
    ! Reduce accumulated properties across all MPI processes.
    do k=1,self%countCells(3)
       do j=1,self%countCells(2)
          do i=1,self%countCells(1)
             call self%radiativeTransferMatter_%accumulationReduction(self%properties(i,j,k))            
          end do
       end do
    end do
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayUnindent('done ['//timer_%reportText()//']'        ,verbosityLevelStandard)
#endif
    ! Solve for state in domain cells local to this process.
    if (mpiSelf%isMaster()) call displayIndent  ('solving for domain state',verbosityLevelStandard)
    call timer_%start()
    countFail=0_c_size_t
    do k=self%sliceMinimum(mpiSelf%rank()),self%sliceMaximum(mpiSelf%rank())
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(k-self%sliceMinimum(mpiSelf%rank()))/dble(self%sliceMaximum(mpiSelf%rank())-self%sliceMinimum(mpiSelf%rank())+1)),isNew=k==self%sliceMinimum(mpiSelf%rank()),verbosity=verbosityLevelWorking)
       do j=1,self%countCells(2)
          do i=1,self%countCells(1)
             call self%radiativeTransferMatter_%stateSolve(self%properties(i,j,k),status)
             if (status /= errorStatusSuccess) countFail=countFail+1_c_size_t
          end do
       end do
    end do
    call mpiBarrier     ()
    countFail=mpiSelf%sum(countFail)
    call timer_    %stop()
    if (mpiSelf%isMaster()) then
       call displayCounterClear(                                   verbosityLevelWorking )
       if (countFail == 0_c_size_t) then
          message='state solve succeeded for all cells'
       else
          write (label,'(i12)') countFail
          message='state solve failed for '//trim(adjustl(label))//' cells'
       end if
       call displayMessage(message,verbosityLevelStandard)
       call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
    end if
    ! Broadcast populated domain cell solutions.
#ifdef USEMPI
    if (mpiSelf%isMaster()) call displayIndent('broadcasting computational domain state',verbosityLevelStandard)
    call timer_%start()
    do p=0,mpiSelf%count()-1
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(p)/dble(mpiSelf%count())),isNew=p==1,verbosity=verbosityLevelWorking)
       do k      =self%sliceMinimum(p),self%sliceMaximum(p)
          do j   =1                   ,self%countCells  (2)
             do i=1                   ,self%countCells  (1)
                call self%radiativeTransferMatter_%broadcastState(p,self%properties(i,j,k))
                call mpiBarrier()
             end do
          end do
       end do
    end do
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayCounterClear(                                   verbosityLevelWorking )
    if (mpiSelf%isMaster()) call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
#endif
    return
  end subroutine cartesian3DStateSolve

  logical function cartesian3DConverged(self)
    !!{
    Return the convergence state of the computational domain.
    !!}
    use :: Arrays_Search                  , only : searchArray
    use :: Disparity_Ratios               , only : Disparity_Ratio
    use :: Display                        , only : displayIndent     , displayMessage        , displayUnindent, displayVerbosity         , &
          &                                        verbosityLevelInfo, verbosityLevelStandard
    use :: MPI_Utilities                  , only : mpiBarrier        , mpiSelf
    use :: Radiative_Transfer_Convergences, only : statusCellFirst   , statusCellLast        , statusCellOther, enumerationStatusCellType
    use :: Timers                         , only : timer
    implicit none
    class           (computationalDomainCartesian3D), intent(inout)                         :: self
    double precision                                , dimension(self%countCellsConvergence) :: convergenceMeasures
    integer         (c_size_t                      )                                        :: i                  , j                      , &
         &                                                                                     k                  , l                      , &
         &                                                                                     m
    double precision                                                                        :: convergenceMeasure , convergenceMeasureRatio
    character       (len=128                       )                                        :: message
    logical                                                                                 :: convergedCriteria
    type            (enumerationStatusCellType     )                                        :: statusCell
    type            (timer                         )                                        :: timer_

    ! Establish a timer.
    timer_=timer()
    if (mpiSelf%isMaster()) call displayIndent('computing convergence state',verbosityLevelStandard)
    call timer_%start()
    ! Find the top k(=self%countCellsConvergence) convergence measures across all cells.
    do i=1,self%countCellsConvergence
       convergenceMeasures(i)=-dble(self%countCellsConvergence+1_c_size_t-i)
    end do
    statusCell=statusCellFirst
    do i=1,self%countCells(1)
       do j=1,self%countCells(2)
          do k=1,self%countCells(3)
             convergenceMeasure=self%radiativeTransferMatter_%convergenceMeasure(self%properties(i,j,k))
             if (convergenceMeasure > convergenceMeasures(1)) then
                ! A new highest-k convergence measure has been found. Insert it into our list.
                if (convergenceMeasure > convergenceMeasures(self%countCellsConvergence)) then
                   l=self%countCellsConvergence
                else
                   l=searchArray(convergenceMeasures,convergenceMeasure)
                end if
                if (l > 1) then
                   do m=1,l-1
                      convergenceMeasures(m)=convergenceMeasures(m+1)
                   end do
                end if
                convergenceMeasures(l)=convergenceMeasure
             end if
             ! Test other convergence criteria.
             if (i == self%countCells(1) .and. j == self%countCells(2) .and. k == self%countCells(3)) statusCell=statusCellLast
             call self%radiativeTransferConvergence_%testConvergence(self%radiativeTransferMatter_,self%properties(i,j,k),statusCell,convergedCriteria)
             statusCell=statusCellOther
          end do
       end do
    end do
    ! Determine if the domain is converged.
    convergenceMeasureRatio=Disparity_Ratio(convergenceMeasures(1),self%convergenceMeasurePrevious)
    cartesian3DConverged   = convergenceMeasures(1)  < self%convergenceThreshold      &
         &                  .and.                                                     &
         &                   convergenceMeasureRatio < self%convergenceRatioThreshold &
         &                  .and.                                                     &
         &                   convergedCriteria
    ! Report on convergence.
    if (mpiSelf%isMaster()) then
       call displayIndent('computational domain convergence',verbosityLevelStandard)    
       write (message,'(a,f6.3,a,e12.6,a,f5.3,a)') '       disparity ratio at ',100.0d0*self%convergencePercentile,'% percentile = ',convergenceMeasures    (1)," [target = ",self%convergenceThreshold     ,"]"
       call displayMessage(trim(message),verbosityLevelStandard)
       write (message,'(a,f6.3,a,e12.6,a,f5.3,a)') 'change disparity ratio at ',100.0d0*self%convergencePercentile,'% percentile = ',convergenceMeasureRatio   ," [target = ",self%convergenceRatioThreshold,"]"
       call displayMessage(trim(message),verbosityLevelStandard)
       if (displayVerbosity() >= verbosityLevelInfo) then
          call displayIndent('convergence measures',verbosityLevelInfo)    
          call displayMessage('rank measure',verbosityLevelInfo)
          call displayMessage('―――――――――――――――――',verbosityLevelInfo)
          do i=1,self%countCellsConvergence
             write (message,'(i4,1x,e12.6)') i,convergenceMeasures(i)
             call displayMessage(trim(message),verbosityLevelInfo)
          end do
          call displayUnindent('done measures',verbosityLevelInfo)    
       end if
       call displayUnindent('done',verbosityLevelStandard)
    end if
    self%convergenceMeasurePrevious=convergenceMeasures(1)
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayUnindent('done ['//timer_%reportText()//']',verbosityLevelStandard)
    return
  end function cartesian3DConverged

  subroutine cartesian3DOutput(self,outputGroup)
    !!{
    Output the computational domain.
    !!}
    !$ use :: HDF5_Access                     , only : hdf5Access
    use    :: ISO_Varying_String              , only : char
    use    :: Numerical_Constants_Astronomical, only : megaparsec
    implicit none
    class           (computationalDomainCartesian3D), intent(inout)                   :: self
    type            (hdf5Object                    ), intent(inout)                   :: outputGroup
    integer         (c_size_t                      )                                  :: i             , j     , &
         &                                                                               k             , output, &
         &                                                                               countOutputs
    double precision                                , allocatable  , dimension(:,:,:) :: propertyScalar
    type            (hdf5Object                    )                                  :: dataset
    
    !$ call hdf5Access%set  ()
    call outputGroup%writeDataset  (self%boundariesCells(1)%boundary                                                         ,'domainBoundariesX',datasetReturned=dataset)
    call dataset    %writeAttribute(megaparsec                                                                               ,'unitsInSI'                                )
    call dataset    %writeAttribute('Mpc'                                                                                    ,'units'                                    )
    call dataset    %writeAttribute('boundaries of computational domain cells in the x direction in 3D Cartesian coordinates','description'                              )
    call outputGroup%writeDataset  (self%boundariesCells(2)%boundary                                                         ,'domainBoundariesY',datasetReturned=dataset)
    call dataset    %writeAttribute(megaparsec                                                                               ,'unitsInSI'                                )
    call dataset    %writeAttribute('Mpc'                                                                                    ,'units'                                    )
    call dataset    %writeAttribute('boundaries of computational domain cells in the y direction in 3D Cartesian coordinates','description'                              )
    call outputGroup%writeDataset  (self%boundariesCells(3)%boundary                                                         ,'domainBoundariesZ',datasetReturned=dataset)
    call dataset    %writeAttribute(megaparsec                                                                               ,'unitsInSI'                                )
    call dataset    %writeAttribute('Mpc'                                                                                    ,'units'                                    )
    call dataset    %writeAttribute('boundaries of computational domain cells in the z direction in 3D Cartesian coordinates','description'                              )
    !$ call hdf5Access%unset()
    countOutputs=self%radiativeTransferMatter_%countOutputs()
    allocate(propertyScalar(self%countCells(1),self%countCells(2),self%countCells(3)))
    do output=1,countOutputs
       do i=1,self%countCells(1)
          do j=1,self%countCells(2)
             do k=1,self%countCells(3)
                propertyScalar(i,j,k)=self%radiativeTransferMatter_%outputProperty(self%properties(i,j,k),output)
             end do
          end do
       end do
       !$ call hdf5Access%set  ()
       call outputGroup%writeDataset(propertyScalar,char(self%radiativeTransferMatter_%outputName(output)))
       !$ call hdf5Access%unset()
    end do
    return
  end subroutine cartesian3DOutput
