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

  use, intrinsic :: ISO_C_Binding                  , only : c_size_t
  use            :: Multi_Counters                 , only : multiCounter
  use            :: Radiative_Transfer_Convergences, only : radiativeTransferConvergenceClass
  use            :: Radiative_Transfer_Matters     , only : radiativeTransferMatterClass     , radiativeTransferPropertiesMatter

  type :: sphericalBoundaries
     !!{
     Type used to store boundaries of computational domain cells for spherical domains.
     !!}
     private
     double precision, allocatable, dimension(:) :: boundary
  end type sphericalBoundaries
  
  !![
  <computationalDomain name="computationalDomainSpherical">
   <description>A computational domain using a spherical grid.</description>
  </computationalDomain>
  !!]
  type, extends(computationalDomainClass) :: computationalDomainSpherical
     !!{
     Implementation of a computational domain using a spherical grid.
     !!}
     private
     double precision                                                , dimension(2) :: boundaries
     integer         (c_size_t                         )                            :: countCells
     double precision                                                               :: convergencePercentile                  , convergenceMeasurePrevious, &
          &                                                                            convergenceThreshold                   , convergenceRatioThreshold
     integer         (c_size_t                         )                            :: countCellsConvergence
     integer         (c_size_t                         ), allocatable, dimension(:) :: sliceMinimum                           , sliceMaximum
     class           (radiativeTransferPropertiesMatter), allocatable, dimension(:) :: properties
     class           (radiativeTransferMatterClass     ), pointer                   :: radiativeTransferMatter_      => null()
     class           (radiativeTransferConvergenceClass), pointer                   :: radiativeTransferConvergence_ => null()
     type            (sphericalBoundaries              )                            :: boundariesCells
   contains
     final     ::                             sphericalDestructor
     procedure :: iterator                 => sphericalIterator
     procedure :: initialize               => sphericalInitialize
     procedure :: reset                    => sphericalReset
     procedure :: indicesFromPosition      => sphericalIndicesFromPosition
     procedure :: lengthToCellBoundary     => sphericalLengthToCellBoundary
     procedure :: absorptionCoefficient    => sphericalAbsorptionCoefficient
     procedure :: interactWithPhotonPacket => sphericalInteractWithPhotonPacket
     procedure :: accumulatePhotonPacket   => sphericalAccumulatePhotonPacket
     procedure :: stateSolve               => sphericalStateSolve
     procedure :: converged                => sphericalConverged
     procedure :: output                   => sphericalOutput
  end type computationalDomainSpherical

  interface computationalDomainSpherical
     !!{
     Constructors for the \refClass{computationalDomainSpherical} computational domain.
     !!}
     module procedure sphericalConstructorParameters
     module procedure sphericalConstructorInternal
  end interface computationalDomainSpherical

  type, extends(domainIterator) :: domainIteratorSpherical
     !!{
     An interactor for spherical computational domains.
     !!}
     private
     type(multiCounter) :: counter_
   contains
     !![
     <methods>
       <method description="Move to the next cell in the domain." method="next" />
     </methods>
     !!]
     procedure :: next => domainIteratorSphericalNext
  end type domainIteratorSpherical
  
contains

  function sphericalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{computationalDomainSpherical} computational domain class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainSpherical     )                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                   , dimension(2)  :: boundaries
    integer         (c_size_t                         )                :: countCells
    class           (radiativeTransferMatterClass     ), pointer       :: radiativeTransferMatter_
    class           (radiativeTransferConvergenceClass), pointer       :: radiativeTransferConvergence_
    double precision                                                   :: convergencePercentile        , convergenceThreshold, &
         &                                                                convergenceRatioThreshold
    
    !![
    <inputParameter>
      <name>boundaries</name>
      <defaultValue>[+0.0d0,+1.0d0]</defaultValue>
      <description>The $r$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countCells</name>
      <defaultValue>3_c_size_t</defaultValue>
      <description>The number of cells in the domain in radius.</description>
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
    self=computationalDomainSpherical(boundaries,countCells,convergencePercentile,convergenceThreshold,convergenceRatioThreshold,radiativeTransferMatter_,radiativeTransferConvergence_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="radiativeTransferMatter_"     />
    <objectDestructor name="radiativeTransferConvergence_"/>
    !!]
    return
  end function sphericalConstructorParameters

  function sphericalConstructorInternal(boundaries,countCells,convergencePercentile,convergenceThreshold,convergenceRatioThreshold,radiativeTransferMatter_,radiativeTransferConvergence_) result(self)
    !!{
    Constructor for the \refClass{computationalDomainSpherical} computational domain class which takes a parameter set as input.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLinear
    implicit none
    type            (computationalDomainSpherical     )                              :: self
    double precision                                   , dimension(2), intent(in   ) :: boundaries
    integer         (c_size_t                         )              , intent(in   ) :: countCells
    double precision                                                 , intent(in   ) :: convergencePercentile        , convergenceThreshold, &
         &                                                                              convergenceRatioThreshold
    class           (radiativeTransferMatterClass     ), target      , intent(in   ) :: radiativeTransferMatter_
    class           (radiativeTransferConvergenceClass), target      , intent(in   ) :: radiativeTransferConvergence_
    !![
    <constructorAssign variables="boundaries, countCells, convergencePercentile, convergenceThreshold, convergenceRatioThreshold, *radiativeTransferMatter_, *radiativeTransferConvergence_"/>
    !!]
    
    ! Construct cell boundaries.
    allocate(self%boundariesCells%boundary(countCells+1))
    self%boundariesCells%boundary=Make_Range(boundaries(1),boundaries(2),int(countCells+1),rangeTypeLinear)
     ! Compute the number of cells used in convergence criteria.
    self%countCellsConvergence     =min(int(dble(self%countCells)*(1.0d0-convergencePercentile),c_size_t)+1,self%countCells)
    self%convergenceMeasurePrevious=huge(0.0d0)
    return
  end function sphericalConstructorInternal

  subroutine sphericalDestructor(self)
    !!{
    Destructor for the \refClass{computationalDomainSpherical} computational domain class.
    !!}
    implicit none
    type(computationalDomainSpherical), intent(inout) :: self

    !![
    <objectDestructor name="self%radiativeTransferMatter_"     />
    <objectDestructor name="self%radiativeTransferConvergence_"/>
    !!]
    return
  end subroutine sphericalDestructor

  subroutine sphericalInitialize(self)
    !!{
    Initialize the computational domain.
    !!}
    use :: Computational_Domain_Volume_Integrators, only : computationalDomainVolumeIntegratorSpherical
    use :: Display                                , only : displayCounter                              , displayCounterClear  , displayIndent, displayUnindent, &
          &                                                verbosityLevelStandard                      , verbosityLevelWorking
    use :: MPI_Utilities                          , only : mpiBarrier                                  , mpiSelf
    use :: Timers                                 , only : timer
    implicit none
    class           (computationalDomainSpherical                ), intent(inout) :: self
    class           (radiativeTransferPropertiesMatter           ), allocatable   :: properties
    integer         (c_size_t                                    )                :: i             , slicesPerProcess, &
         &                                                                           slicesExtra
    type            (computationalDomainVolumeIntegratorSpherical), allocatable   :: integrator
    double precision                                              , dimension(2)  :: boundariesCell
#ifdef USEMPI
    integer                                                                       :: p
#endif
    type            (timer                                       )                :: timer_

    ! Establish a timer.
    timer_=timer()
    ! Divide up the domain between MPI processes.
    allocate(self%sliceMinimum(0:mpiSelf%count()-1))
    allocate(self%sliceMaximum(0:mpiSelf%count()-1))
    slicesPerProcess=self%countCells/mpiSelf%count()
    slicesExtra     =self%countCells-mpiSelf%count()*slicesPerProcess
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
    allocate(self%properties(self%countCells),mold=properties)
    deallocate(properties)
    do i    =1,self%countCells
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(i-1_c_size_t)/dble(self%countCells)),isNew=i==1,verbosity=verbosityLevelWorking)
       boundariesCell   (:)=self%boundariesCells%boundary(i:i+1)
       ! Build a volume integrator for this cell.
       allocate(integrator)
       integrator=computationalDomainVolumeIntegratorSpherical(boundariesCell)
       ! Populate this cell.
       call self%radiativeTransferMatter_%populateDomain(self%properties(i),integrator,onProcess=i >= self%sliceMinimum(mpiSelf%rank()) .and. i <= self%sliceMaximum(mpiSelf%rank()))
       ! Destroy the integrator.
       deallocate(integrator)
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
       do i=self%sliceMinimum(p),self%sliceMaximum(p)
          call self%radiativeTransferMatter_%broadcastDomain(p,self%properties(i))
          call mpiBarrier()
       end do
    end do
    call mpiBarrier()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayCounterClear(                                   verbosityLevelWorking )
    if (mpiSelf%isMaster()) call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
#endif
    return
  end subroutine sphericalInitialize

  subroutine sphericalIterator(self,iterator)
    !!{
    Construct an iterator for this domain,
    !!}
    implicit none
    class(computationalDomainSpherical), intent(inout)              :: self
    class(domainIterator              ), intent(inout), allocatable :: iterator

    ! Build a multi-counter object for use in iterating over domain cells.
    allocate(domainIteratorSpherical :: iterator)
    select type (iterator)
    type is (domainIteratorSpherical)
       iterator%counter_=multiCounter([self%countCells])
    end select
    return
  end subroutine sphericalIterator

  logical function domainIteratorSphericalNext(self)
    implicit none
    class(domainIteratorSpherical), intent(inout) :: self

    domainIteratorSphericalNext=self%counter_%increment()
    return
  end function domainIteratorSphericalNext
  
  subroutine sphericalReset(self)
    !!{
    Reset the computational domain prior to a new iteration.
    !!}
    implicit none
    class  (computationalDomainSpherical), intent(inout) :: self
    integer(c_size_t                    )                :: i
    
    do i=1,self%countCells
       call self%radiativeTransferMatter_%reset(self%properties(i))
    end do
    return
  end subroutine sphericalReset

  subroutine sphericalIndicesFromPosition(self,position,indices)
    !!{
    Determine the indices of the cell containing the given point.
    !!}
    use :: Arrays_Search, only : searchArray

    implicit none
    class           (computationalDomainSpherical), intent(inout)                            :: self
    double precision                                           , dimension(3), intent(in   ) :: position
    integer         (c_size_t                    ), allocatable, dimension(:), intent(inout) :: indices
    double precision                                                                         :: radius
    
    ! Allocate indices to the correct size if necessary.
    if (allocated(indices)) then
       if (size(indices) /= 1) then
          deallocate(indices   )
          allocate  (indices(1))
       end if
    else
       allocate     (indices(1))
    end if
    ! Determine indices.
    radius=sqrt(sum(position**2))
    if (radius < self%boundariesCells%boundary(1) .or. radius >= self%boundariesCells%boundary(self%countCells+1)) then
       indices=-huge(0_c_size_t)
    else
       indices(1)=searchArray(self%boundariesCells%boundary,radius)
    end if
    return
  end subroutine sphericalIndicesFromPosition

  double precision function sphericalAbsorptionCoefficient(self,photonPacket,indices)
    !!{
    Return the absorption coefficient for the given photon packet in the given domain cell.
    !!}
    implicit none
    class  (computationalDomainSpherical      )              , intent(inout) :: self
    class  (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer(c_size_t                          ), dimension(:), intent(in   ) :: indices

    sphericalAbsorptionCoefficient=self%radiativeTransferMatter_%absorptionCoefficient(                              &
         &                                                                             self%properties  (            &
         &                                                                                               indices(1)  &
         &                                                                                              )          , &
         &                                                                                  photonPacket             &
         &                                                                            )
    return
  end function sphericalAbsorptionCoefficient
  
  double precision function sphericalLengthToCellBoundary(self,photonPacket,indices,indicesNeighbor,positionBoundary)
    !!{
    Return the length to the first domain cell boundary intersected by the given photon packet.
    !!}
    implicit none
    class           (computationalDomainSpherical      )                           , intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass)                           , intent(inout) :: photonPacket
    integer         (c_size_t                          )             , dimension(:), intent(in   ) :: indices
    integer         (c_size_t                          ), allocatable, dimension(:), intent(inout) :: indicesNeighbor
    double precision                                                 , dimension(3), intent(  out) :: positionBoundary
    double precision                                                 , dimension(3)                :: position            , direction
    integer                                                                                        :: j
    double precision                                                                               :: lengthToFace        , argument            , &
         &                                                                                            lengthToFacePositive, lengthToFaceNegative, &
         &                                                                                            radiusSquared       , directionRadius     , &
         &                                                                                            radiusBoundary

    ! Allocate indices to the correct size if necessary.
    if (allocated(indicesNeighbor)) then
       if (size(indicesNeighbor) /= 1) then
          deallocate(indicesNeighbor   )
          allocate  (indicesNeighbor(1))
       end if
    else
       allocate     (indicesNeighbor(1))
    end if
    ! Initialize to infinite distance.
    sphericalLengthToCellBoundary=huge(0.0d0)
    ! Get position and direction of the photon packet.
    position                     =photonPacket%position ()
    direction                    =photonPacket%direction()
    ! Consider boundaries along the r-axis.
    radiusSquared  =sum(position          **2)
    directionRadius=sum(position*direction   )
    do j=0,1 ! Faces.
       argument=+self%boundariesCells%boundary(indices(1)+j)**2-radiusSquared+directionRadius**2
       if (argument > 0.0d0) then ! Negative argument means that the photon never crosses the boundary.
          lengthToFacePositive=+sqrt(argument)-directionRadius
          lengthToFaceNegative=-sqrt(argument)-directionRadius
          if (lengthToFaceNegative > 0.0d0) then
             lengthToFace=min(lengthToFacePositive,lengthToFaceNegative)
          else
             lengthToFace=    lengthToFacePositive
          end if
          if     (                                                   &
               &   (                                                 &
               &       lengthToFace    >  0.0d0                      & ! Distance to face is positive (i.e. in direction of travel).
               &    .or.                                             &
               &     (                                               &
               &       lengthToFace    == 0.0d0                      & ! Photon packet is precisely on the lower face of the cell, and is moving in the negative direction
               &      .and.                                          & ! - in this case we want to transition it to the cell below.
               &       directionRadius <  0.0d0                      &
               &      .and.                                          &
               &       j            == 0                             &
               &     )                                               &
               &   )                                                 &
               &  .and.                                              &
               &       lengthToFace <  sphericalLengthToCellBoundary & ! This is the shortest distance to a face we've found.
               & ) then
             sphericalLengthToCellBoundary=lengthToFace
             indicesNeighbor    =+indices
             indicesNeighbor (1)=+indicesNeighbor(1)+(2*j-1)
             positionBoundary   =+position                   &
                  &              +direction                  &
                  &              *lengthToFace
             ! Ensure that the boundary position we return is actually into the neighboring cell - due to rounding errors we
             ! can't guarantee that the photon packet is precisely on the domain boundary, so we simply make small shifts in
             ! the position until it is.
             radiusBoundary=sqrt(sum(positionBoundary**2)) 
             do while ((j == 0 .and. radiusBoundary >= self%boundariesCells%boundary(indices(1)+j)) .or. (j == 1 .and. radiusBoundary <= self%boundariesCells%boundary(indices(1)+j)) )
                select case (j)
                case (0)
                   positionBoundary=positionBoundary/(1.0d0+epsilon(0.0d0))
                case (1)
                   positionBoundary=positionBoundary*(1.0d0+epsilon(0.0d0))
                end select
                radiusBoundary=sqrt(sum(positionBoundary**2)) 
             end do
             ! If the neighbor cell is outside of the domain, mark it as such.
             if (indicesNeighbor(1) < 1 .or. indicesNeighbor(1) > self%countCells) &
                  & indicesNeighbor(1)=-huge(0_c_size_t)
          end if
       end if
    end do
    return
  end function sphericalLengthToCellBoundary

  subroutine sphericalAccumulatePhotonPacket(self,photonPacket,indices,absorptionCoefficient,lengthTraversed)
    !!{
    Accumulate ``absorptions'' from the photon packet as it traverses a cell of the computational domain.
    !!}
    implicit none
    class           (computationalDomainSpherical      )              , intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer         (c_size_t                          ), dimension(:), intent(in   ) :: indices
    double precision                                                  , intent(in   ) :: absorptionCoefficient, lengthTraversed

    call self%radiativeTransferMatter_%accumulatePhotonPacket(                                       &
         &                                                    self%properties           (            &
         &                                                                               indices(1)  &
         &                                                                              )          , &
         &                                                         photonPacket                    , &
         &                                                         absorptionCoefficient           , &
         &                                                         lengthTraversed                   &
         &                                                   )
    return
  end subroutine sphericalAccumulatePhotonPacket

  logical function sphericalInteractWithPhotonPacket(self,photonPacket,indices)
    !!{
    Allow matter in a domain cell to interact with the photon packet.
    !!}
    implicit none
    class  (computationalDomainSpherical      )              , intent(inout) :: self
    class  (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket
    integer(c_size_t                          ), dimension(:), intent(inout) :: indices

    sphericalInteractWithPhotonPacket=self%radiativeTransferMatter_%interactWithPhotonPacket(                              &
         &                                                                                   self%properties  (            &
         &                                                                                                     indices(1)  &
         &                                                                                                    )          , &
         &                                                                                        photonPacket             &
         &                                                                                  )
    ! Reset the second index to indicate that no boundary was just crossed and so we should consider both boundaries on the next
    ! "length-to-boundary" calculation.
    indices(2)=-1_c_size_t
    return
  end function sphericalInteractWithPhotonPacket
  
  subroutine sphericalStateSolve(self)
    !!{
    Solve for the state of matter in the computational domain.
    !!}
    use :: Display         , only : displayCounter    , displayCounterClear   , displayIndent        , displayMessage, &
          &                         displayUnindent   , verbosityLevelStandard, verbosityLevelWorking
    use :: Error, only : errorStatusSuccess
    use :: MPI_Utilities   , only : mpiBarrier        , mpiSelf
    use :: Timers          , only : timer
    implicit none
    class    (computationalDomainSpherical), intent(inout) :: self
    integer  (c_size_t                    )                :: i     , countFail
    integer                                                :: status
#ifdef USEMPI
    integer                                                :: p
#endif
    type     (timer                       )                :: timer_
    type     (varying_string              )                :: message
    character(len=12                      )                :: label
    
    ! Establish a timer.
    timer_=timer()
#ifdef USEMPI
    if (mpiSelf%isMaster()) call displayIndent  ('accumulating absorptions across processes',verbosityLevelStandard)
    call timer_%start()
    ! Reduce accumulated properties across all MPI processes.
    do i=1,self%countCells
       call self%radiativeTransferMatter_%accumulationReduction(self%properties(i))            
    end do
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayUnindent('done ['//timer_%reportText()//']'        ,verbosityLevelStandard)
#endif
    ! Solve for state in domain cells local to this process.
    if (mpiSelf%isMaster()) call displayIndent  ('solving for domain state',verbosityLevelStandard)
    call timer_%start()
    countFail=0_c_size_t
    do i=self%sliceMinimum(mpiSelf%rank()),self%sliceMaximum(mpiSelf%rank())
       if (mpiSelf%isMaster()) call displayCounter(int(100.0d0*dble(i-self%sliceMinimum(mpiSelf%rank()))/dble(self%sliceMaximum(mpiSelf%rank())-self%sliceMinimum(mpiSelf%rank())+1)),isNew=i==self%sliceMinimum(mpiSelf%rank()),verbosity=verbosityLevelWorking)
       call self%radiativeTransferMatter_%stateSolve(self%properties(i),status)
       if (status /= errorStatusSuccess) countFail=countFail+1_c_size_t
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
       do i=self%sliceMinimum(p),self%sliceMaximum(p)
          call self%radiativeTransferMatter_%broadcastState(p,self%properties(i))
          call mpiBarrier()
       end do
    end do
    call mpiBarrier     ()
    call timer_    %stop()
    if (mpiSelf%isMaster()) call displayCounterClear(                                   verbosityLevelWorking )
    if (mpiSelf%isMaster()) call displayUnindent     ('done ['//timer_%reportText()//']',verbosityLevelStandard)
#endif
    return
  end subroutine sphericalStateSolve

  logical function sphericalConverged(self)
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
    class           (computationalDomainSpherical), intent(inout)                         :: self
    double precision                              , dimension(self%countCellsConvergence) :: convergenceMeasures
    integer         (c_size_t                    )                                        :: i                  ,     l                  , &
         &                                                                                   m
    double precision                                                                      :: convergenceMeasure , convergenceMeasureRatio
    character       (len=128                     )                                        :: message
    logical                                                                               :: convergedCriteria
    type            (enumerationStatusCellType   )                                        :: statusCell
    type            (timer                       )                                        :: timer_

    ! Establish a timer.
    timer_=timer()
    if (mpiSelf%isMaster()) call displayIndent('computing convergence state',verbosityLevelStandard)
    call timer_%start()
    ! Find the top k(=self%countCellsConvergence) convergence measures across all cells.
    do i=1,self%countCellsConvergence
       convergenceMeasures(i)=-dble(self%countCellsConvergence+1_c_size_t-i)
    end do
    statusCell=statusCellFirst
    do i=1,self%countCells
       convergenceMeasure=self%radiativeTransferMatter_%convergenceMeasure(self%properties(i))
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
       if (i == self%countCells) statusCell=statusCellLast
       call self%radiativeTransferConvergence_%testConvergence(self%radiativeTransferMatter_,self%properties(i),statusCell,convergedCriteria)
       statusCell=statusCellOther
    end do
    ! Determine if the domain is converged.
    convergenceMeasureRatio=Disparity_Ratio(convergenceMeasures(1),self%convergenceMeasurePrevious)
    sphericalConverged     = convergenceMeasures(1)  < self%convergenceThreshold      &
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
  end function sphericalConverged

  subroutine sphericalOutput(self,outputGroup)
    !!{
    Output the computational domain.
    !!}
    !$ use :: HDF5_Access                     , only : hdf5Access
    use    :: ISO_Varying_String              , only : char
    use    :: Numerical_Constants_Astronomical, only : megaparsec
    implicit none
    class           (computationalDomainSpherical), intent(inout)               :: self
    type            (hdf5Object                  ), intent(inout)               :: outputGroup
    integer         (c_size_t                    )                              :: i             , countOutputs, &
         &                                                                         output
    double precision                              , allocatable  , dimension(:) :: propertyScalar
    type            (hdf5Object                  )                              :: dataset
    
    !$ call hdf5Access%set  ()
    call outputGroup%writeDataset  (self%boundariesCells%boundary                                                              ,'domainBoundariesRadial'  ,datasetReturned=dataset)
    call dataset    %writeAttribute(megaparsec                                                                                 ,'unitsInSI'                                       )
    call dataset    %writeAttribute('Mpc'                                                                                      ,'units'                                           )
    call dataset    %writeAttribute('boundaries of computational domain cells in the radial direction in spherical coordinates','description'                                     )
    !$ call hdf5Access%unset()
    countOutputs=self%radiativeTransferMatter_%countOutputs()
    allocate(propertyScalar(self%countCells))
    do output=1,countOutputs
       do i=1,self%countCells
          propertyScalar(i)=self%radiativeTransferMatter_%outputProperty(self%properties(i),output)
       end do
       !$ call hdf5Access%set  ()
       call outputGroup%writeDataset(propertyScalar,char(self%radiativeTransferMatter_%outputName(output)))
       !$ call hdf5Access%unset()
    end do
    return
  end subroutine sphericalOutput
