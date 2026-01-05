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

  use            :: Computational_Domains           , only : computationalDomainClass
  use, intrinsic :: ISO_C_Binding                   , only : c_size_t
  use            :: Numerical_Random_Numbers        , only : randomNumberGeneratorClass
  use            :: Radiative_Transfer_Convergences , only : radiativeTransferConvergenceClass
  use            :: Radiative_Transfer_Outputters   , only : radiativeTransferOutputterClass
  use            :: Radiative_Transfer_Photon_Packet, only : radiativeTransferPhotonPacketClass
  use            :: Radiative_Transfer_Sources      , only : radiativeTransferSourceClass
  
  !![
  <task name="taskRadiativeTransfer">
   <description>A task which performs radiative transfer.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskRadiativeTransfer
     !!{
     Implementation of a task which performs radiative transfer.
     !!}
     private
     type            (varying_string                    )                            :: outputGroupName
     class           (computationalDomainClass          ), pointer                   :: computationalDomain_           => null()
     class           (radiativeTransferConvergenceClass ), pointer                   :: radiativeTransferConvergence_  => null()
     class           (radiativeTransferPhotonPacketClass), pointer                   :: radiativeTransferPhotonPacket_ => null()
     class           (radiativeTransferSourceClass      ), pointer                   :: radiativeTransferSource_       => null()
     class           (radiativeTransferOutputterClass   ), pointer                   :: radiativeTransferOutputter_    => null()
     class           (randomNumberGeneratorClass        ), pointer                   :: randomNumberGenerator_         => null()
     integer                                                                         :: wavelengthCountPerDecade
     integer         (c_size_t                          )                            :: countPhotonsPerWavelength               , countIterationsMaximum                         , &
          &                                                                             countIterationsMinimum                  , countPhotonsPerWavelengthFinalIteration
     double precision                                                                :: wavelengthMinimum                       , wavelengthMaximum
     double precision                                    , allocatable, dimension(:) :: wavelengthsMinimum                      , wavelengthsMaximum                             , &
          &                                                                             wavelengths
     logical                                                                         :: outputIterations                        , nodeComponentsInitialized              =.false.
   contains
     final     ::            radiativeTransferDestructor
     procedure :: perform => radiativeTransferPerform
  end type taskRadiativeTransfer

  interface taskRadiativeTransfer
     !!{
     Constructors for the \refClass{taskRadiativeTransfer} task.
     !!}
     module procedure radiativeTransferConstructorParameters
     module procedure radiativeTransferConstructorInternal
  end interface taskRadiativeTransfer
  
contains

  function radiativeTransferConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskRadiativeTransfer} task class which takes a parameter set as input.
    !!}
    use :: ISO_Varying_String, only : var_str
    use :: Input_Parameters  , only : inputParameter              , inputParameters
    use :: Galacticus_Nodes  , only : nodeClassHierarchyInitialize
    use :: Node_Components   , only : Node_Components_Initialize
    implicit none
    type            (taskRadiativeTransfer             )                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (computationalDomainClass          ), pointer       :: computationalDomain_
    class           (radiativeTransferConvergenceClass ), pointer       :: radiativeTransferConvergence_
    class           (radiativeTransferPhotonPacketClass), pointer       :: radiativeTransferPhotonPacket_
    class           (radiativeTransferSourceClass      ), pointer       :: radiativeTransferSource_
    class           (radiativeTransferOutputterClass   ), pointer       :: radiativeTransferOutputter_
    class           (randomNumberGeneratorClass        ), pointer       :: randomNumberGenerator_
    type            (inputParameters                   ), pointer       :: parametersRoot
    integer                                                             :: wavelengthCountPerDecade
    integer         (c_size_t                          )                :: countPhotonsPerWavelength     , countIterationsMinimum                 , &
         &                                                                 countIterationsMaximum        , countPhotonsPerWavelengthFinalIteration
    double precision                                                    :: wavelengthMinimum             , wavelengthMaximum
    type            (varying_string                    )                :: outputGroupName
    logical                                                             :: outputIterations

    ! Ensure the nodes objects are initialized. This is necessary to ensure that the abundances class is initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    self%nodeComponentsInitialized=.true.
    !![
    <inputParameter>
      <name>wavelengthMinimum</name>
      <defaultValue>0.3d4</defaultValue>
      <description>The minimum wavelength at which to sample photon packets.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <defaultValue>10.0d4</defaultValue>
      <description>The maximum wavelength at which to sample photon packets.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthCountPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of wavelengths per decade at which to sample photon packets.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countPhotonsPerWavelength</name>
      <defaultValue>10_c_size_t</defaultValue>
      <description>The number of photon packets to generate at each wavelength.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countPhotonsPerWavelengthFinalIteration</name>
      <defaultValue>countPhotonsPerWavelength</defaultValue>
      <description>The number of photon packets to generate at each wavelength on the final iteration.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countIterationsMinimum</name>
      <defaultValue>1_c_size_t</defaultValue>
      <description>The minimum number of iterations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countIterationsMaximum</name>
      <defaultValue>10_c_size_t</defaultValue>
      <description>The maximum number of iterations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroupName</name>
      <defaultValue>var_str('radiativeTransferModel')</defaultValue>
      <description>The name of the group to which results should be output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputIterations</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, output data for all iterations, not just the final iteration.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="computationalDomain"           name="computationalDomain_"           source="parameters"/>
    <objectBuilder class="radiativeTransferConvergence"  name="radiativeTransferConvergence_"  source="parameters"/>
    <objectBuilder class="radiativeTransferPhotonPacket" name="radiativeTransferPhotonPacket_" source="parameters"/>
    <objectBuilder class="radiativeTransferSource"       name="radiativeTransferSource_"       source="parameters"/>
    <objectBuilder class="radiativeTransferOutputter"    name="radiativeTransferOutputter_"    source="parameters"/>
    <objectBuilder class="randomNumberGenerator"         name="randomNumberGenerator_"         source="parameters"/>
    !!]
    self=taskRadiativeTransfer(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,countPhotonsPerWavelength,countPhotonsPerWavelengthFinalIteration,countIterationsMinimum,countIterationsMaximum,outputGroupName,outputIterations,computationalDomain_,radiativeTransferPhotonPacket_,radiativeTransferSource_,radiativeTransferOutputter_,radiativeTransferConvergence_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="computationalDomain_"          />
    <objectDestructor name="radiativeTransferConvergence_" />
    <objectDestructor name="radiativeTransferPhotonPacket_"/>
    <objectDestructor name="radiativeTransferSource_"      />
    <objectDestructor name="radiativeTransferOutputter_"   />
    <objectDestructor name="randomNumberGenerator_"        />
    !!]
   return
  end function radiativeTransferConstructorParameters

  function radiativeTransferConstructorInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,countPhotonsPerWavelength,countPhotonsPerWavelengthFinalIteration,countIterationsMinimum,countIterationsMaximum,outputGroupName,outputIterations,computationalDomain_,radiativeTransferPhotonPacket_,radiativeTransferSource_,radiativeTransferOutputter_,radiativeTransferConvergence_,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{taskRadiativeTransfer} task class which takes a parameter set as input.
    !!}
    use :: ISO_Varying_String, only : char
    use :: Numerical_Ranges  , only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (taskRadiativeTransfer             )                        :: self
    integer                                             , intent(in   )         :: wavelengthCountPerDecade
    integer         (c_size_t                          ), intent(in   )         :: countPhotonsPerWavelength     , countIterationsMinimum                 , &
         &                                                                         countIterationsMaximum        , countPhotonsPerWavelengthFinalIteration
    double precision                                    , intent(in   )         :: wavelengthMinimum             , wavelengthMaximum
    type            (varying_string                    ), intent(in   )         :: outputGroupName
    logical                                             , intent(in   )         :: outputIterations
    class           (computationalDomainClass          ), intent(in   ), target :: computationalDomain_
    class           (radiativeTransferConvergenceClass ), intent(in   ), target :: radiativeTransferConvergence_
    class           (radiativeTransferPhotonPacketClass), intent(in   ), target :: radiativeTransferPhotonPacket_
    class           (radiativeTransferSourceClass      ), intent(in   ), target :: radiativeTransferSource_
    class           (radiativeTransferOutputterClass   ), intent(in   ), target :: radiativeTransferOutputter_
    class           (randomNumberGeneratorClass        ), intent(in   ), target :: randomNumberGenerator_
    integer                                                                     :: wavelengthCount
    double precision                                                            :: wavelengthFactor
    !![
    <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, wavelengthCountPerDecade, countPhotonsPerWavelength, countPhotonsPerWavelengthFinalIteration, countIterationsMinimum, countIterationsMaximum, outputGroupName, outputIterations, *computationalDomain_, *radiativeTransferPhotonPacket_, *radiativeTransferSource_, *radiativeTransferOutputter_, *radiativeTransferConvergence_, *randomNumberGenerator_"/>
    !!]

    wavelengthCount =int(log10(wavelengthMaximum/wavelengthMinimum)*dble(wavelengthCountPerDecade)+1)
    wavelengthFactor=exp(log(wavelengthMaximum/wavelengthMinimum)/dble(wavelengthCount))
    allocate(self%wavelengths       (wavelengthCount))
    allocate(self%wavelengthsMinimum(wavelengthCount))
    allocate(self%wavelengthsMaximum(wavelengthCount))
    self%wavelengths       =Make_Range(wavelengthMinimum,wavelengthMaximum,wavelengthCount,rangeTypeLogarithmic,rangeBinned=.true.)
    self%wavelengthsMinimum=self%wavelengths/sqrt(wavelengthFactor)
    self%wavelengthsMaximum=self%wavelengths*sqrt(wavelengthFactor)
    return
  end function radiativeTransferConstructorInternal

  subroutine radiativeTransferDestructor(self)
    !!{
    Destructor for the \refClass{taskRadiativeTransfer} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskRadiativeTransfer), intent(inout) :: self

    !![
    <objectDestructor name="self%computationalDomain_"          />
    <objectDestructor name="self%radiativeTransferConvergence_" />
    <objectDestructor name="self%radiativeTransferPhotonPacket_"/>
    <objectDestructor name="self%radiativeTransferSource_"      />
    <objectDestructor name="self%radiativeTransferOutputter_"   />
    <objectDestructor name="self%randomNumberGenerator_"        />
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine radiativeTransferDestructor

  subroutine radiativeTransferPerform(self,status)
    !!{
    Perform radiative transfer and output results.
    !!}
    use :: Display                 , only : displayIndent                            , displayMessage    , displayUnindent, verbosityLevelStandard
    use :: Error                   , only : Error_Report                             , errorStatusSuccess
    use :: Output_HDF5             , only : outputFile
    use :: IO_HDF5                 , only : hdf5Object
    use :: MPI_Utilities           , only : mpiBarrier                               , mpiSelf
    use :: Statistics_Distributions, only : distributionFunction1DNegativeExponential
    use :: Timers                  , only : timer
    implicit none
    class           (taskRadiativeTransfer                    ), intent(inout), target       :: self
    integer                                                    , intent(  out), optional     :: status
    integer         (c_size_t                                 ), allocatable  , dimension(:) :: domainIndices            , domainIndicesNeighbor
    double precision                                                          , dimension(3) :: positionBoundary
    integer                                                                                  :: iWavelength
    integer         (c_size_t                                 )                              :: iPhoton                  , countIterations           , &
         &                                                                                      countPhotonsPerWavelength
    type            (distributionFunction1DNegativeExponential)                              :: opticalDepthDistribution
    logical                                                                                  :: photonPacketInDomain     , photonPacketAlive         , &
         &                                                                                      updateDomainIndices      , converged                 , &
         &                                                                                      photonPacketInteracts    , finalIteration
    double precision                                                                         :: opticalDepthToInteraction, opticalDepthToCellBoundary, &
         &                                                                                      absorptionCoefficient    , lengthToCellBoundary      , &
         &                                                                                      lengthTraversed
    type            (hdf5Object                               )                              :: outputGroup              , iterationOutputGroup
    character       (len=128                                  )                              :: message                  , label
    type            (timer                                    )                              :: timer_                   , timerTotal_               , &
         &                                                                                      timerIteration_

    if (mpiSelf%isMaster()) call displayIndent('Begin task: radiative transfer')
    ! Establish timers.
    timer_         =timer()
    timerTotal_    =timer()
    timerIteration_=timer()
    ! Open group for output of our model data.
    if (mpiSelf%isMaster()) outputGroup=outputFile%openGroup(char(self%outputGroupName),'Radiative transfer model.')
    call timerTotal_%start()
    ! Initialize the computational domain.
    call self%computationalDomain_       %initialize      (                                         )
    ! Compute and output properties of the sources.
    call self%radiativeTransferOutputter_%sourceProperties(self%radiativeTransferSource_,outputGroup)
    ! Construct a negative exponential distribution from which to sample optical depths.
    opticalDepthDistribution=distributionFunction1DNegativeExponential(1.0d0)
    ! Iterate until convergence.
    converged      =.false.
    finalIteration =.false.
    countIterations=0_c_size_t
    iterations : do while ((.not.converged .or. finalIteration .or. countIterations < self%countIterationsMinimum) .and. countIterations < self%countIterationsMaximum)
       countIterations=countIterations+1
       if (mpiSelf%isMaster()) then
          write (message,'(a,i6)') 'begin iteration ',countIterations
          call displayIndent(trim(message),verbosityLevelStandard)
       end if
       call mpiBarrier()
       call timerIteration_%start()
       ! Reset the computational domain and outputter.
       call self%computationalDomain_       %reset()
       call self%radiativeTransferOutputter_%reset()
       ! Iterate over photon wavelengths.
       call timer_%start()
       if (mpiSelf%isMaster()) call displayIndent('cast photon packets',verbosityLevelStandard)
       wavelengths : do iWavelength=1,size(self%wavelengths)
          ! Iterate over photon packets.
          countPhotonsPerWavelength=self%countPhotonsPerWavelength
          if (finalIteration) countPhotonsPerWavelength=self%countPhotonsPerWavelengthFinalIteration
          photonPackets : do iPhoton=1,countPhotonsPerWavelength
#ifdef USEMPI
             ! Skip photons which do not belong to this MPI process.
             if (mod(iPhoton,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Initialize the photon packet. Wavelength is sampled uniformly in the log of wavelength within this wavelength bin.
             call self%radiativeTransferPhotonPacket_%wavelengthSet         (                                                                  &
                  &                                                          +exp(                                                             &
                  &                                                               +log(                                                        &
                  &                                                                    +self%wavelengthsMaximum                  (iWavelength) &
                  &                                                                    /self%wavelengthsMinimum                  (iWavelength) &
                  &                                                                   )                                                        &
                  &                                                               *     self%randomNumberGenerator_%uniformSample(           ) &
                  &                                                               +log(                                                        &
                  &                                                                    +self%wavelengthsMinimum                  (iWavelength) &
                  &                                                                   )                                                        &
                  &                                                              )                                                             &
                  &                                                         )
             call self%radiativeTransferPhotonPacket_%wavelengthMinimumSet  (self%wavelengthsMinimum                              (iWavelength              ))
             call self%radiativeTransferPhotonPacket_%wavelengthMaximumSet  (self%wavelengthsMaximum                              (iWavelength              ))
             call self%radiativeTransferSource_      %initializePhotonPacket(self%radiativeTransferPhotonPacket_                                             )
             call self%radiativeTransferPhotonPacket_%luminositySet         (self%radiativeTransferPhotonPacket_%luminosity()/dble(countPhotonsPerWavelength))
             ! If the photon packet has no luminosity, we can simply ignore it.
             if (self%radiativeTransferPhotonPacket_%luminosity() <= 0.0d0) cycle
             ! Find the indices of the domain cell in which the photon is located.
             call self%computationalDomain_%indicesFromPosition(self%radiativeTransferPhotonPacket_%position(),domainIndices)
             ! Iterate photon packet steps until the photon is absorbed or reaches the domain boundary.
             photonPacketAlive   =.true.
             photonPacketInDomain=all(domainIndices > -huge(0_c_size_t))
             propagate : do while (photonPacketAlive .and. photonPacketInDomain)
                ! Draw an optical depth from a negative exponential distribution. Note that we choose a new optical depth on each step
                ! even if the photon packet did not interact on the previous step. This is allowable due to the memorylessness property
                ! of the negative exponential distribution.
                opticalDepthToInteraction=opticalDepthDistribution%sample(randomNumberGenerator_=self%randomNumberGenerator_)
                ! Get absorption coefficient in domain cell.
                absorptionCoefficient=self%computationalDomain_%absorptionCoefficient(self%radiativeTransferPhotonPacket_,domainIndices                                       )
                ! Find length to domain cell boundary, and coordinates of the neighboring cell (or if the domain edge is reached).
                lengthToCellBoundary =self%computationalDomain_%lengthToCellBoundary (self%radiativeTransferPhotonPacket_,domainIndices,domainIndicesNeighbor,positionBoundary)
                ! Determine if the photon packet interacts in this domain cell.
                opticalDepthToCellBoundary=+absorptionCoefficient &
                     &                     *lengthToCellBoundary
                ! Propagate the packet.
                updateDomainIndices  =.false.
                photonPacketInteracts=.false.
                if (opticalDepthToCellBoundary > opticalDepthToInteraction) then
                   ! Interaction occurs.
                   lengthTraversed      =+opticalDepthToInteraction &
                        &                /absorptionCoefficient
                   photonPacketInteracts=.true.
                else
                   ! No interaction occurs - check if the photon packet remains in the computational domain, and update domain indices if so.
                   lengthTraversed=lengthToCellBoundary
                   call self%radiativeTransferPhotonPacket_%positionSet(positionBoundary)
                   photonPacketInDomain=all(domainIndicesNeighbor > -huge(0_c_size_t))
                   if (photonPacketInDomain) updateDomainIndices=.true.
                end if
                ! Allow matter to absorb the photon packet according to the Lucy (1999; A&A; 344; 282;
                ! https://ui.adsabs.harvard.edu/abs/1999A%26A...344..282L; see section 3.4) methodology.
                call self%computationalDomain_%accumulatePhotonPacket(self%radiativeTransferPhotonPacket_,domainIndices,absorptionCoefficient,lengthTraversed)
                ! Allow interaction with matter.
                if (photonPacketInteracts) photonPacketAlive=self%computationalDomain_%interactWithPhotonPacket(self%radiativeTransferPhotonPacket_,domainIndices)
                ! Move to a new domain cell if necessary.
                if (updateDomainIndices) domainIndices=domainIndicesNeighbor
             end do propagate
             ! Photon packet propagation is done - determine how to handle it. If the packet is no longer alive, there's nothing to do.
             if (photonPacketAlive) then
                ! Photon packet is still alive - it should have left the domain.
                if (photonPacketInDomain) call Error_Report('photon packet propagation stopped while still in computational domain'//{introspection:location})
                call self%radiativeTransferOutputter_  %photonPacketEscapes(self%radiativeTransferPhotonPacket_)
                call self%radiativeTransferConvergence_%photonPacketEscapes(self%radiativeTransferPhotonPacket_)
             end if
          end do photonPackets
       end do wavelengths
       call mpiBarrier()
       call timer_%stop()
       if (mpiSelf%isMaster()) call displayUnindent('done ['//timer_%reportText()//']',verbosityLevelStandard)
       ! Skip state solving and convergence on any final iteration.
       if (.not.finalIteration) then
          ! Solve for state of matter in the computational domain.
          if (mpiSelf%isMaster()) call displayIndent('solve for matter state',verbosityLevelStandard)
          call timer_%start()
          call self%computationalDomain_%stateSolve()
          call mpiBarrier()
          call timer_%stop()
          if (mpiSelf%isMaster()) call displayUnindent('done ['//timer_%reportText()//']',verbosityLevelStandard)       
          ! Test convergence on the computational domain.
          if (mpiSelf%isMaster()) call displayIndent('test for convergence',verbosityLevelStandard)
          call timer_%start()
          converged=self%computationalDomain_%converged()
          call mpiBarrier()
          call timer_%stop()
          if (mpiSelf%isMaster()) then
             if (converged) then
                call displayMessage('converged'    ,verbosityLevelStandard)
             else
                call displayMessage('not converged',verbosityLevelStandard)
             end if
             call displayUnindent('done ['//timer_%reportText()//']',verbosityLevelStandard)
          end if
          call mpiBarrier()
          call timerIteration_%stop()
          if (mpiSelf%isMaster()) call displayUnindent('done ['//timerIteration_%reportText()//']',verbosityLevelStandard)
       end if
       ! Output the computational domain.
       if (mpiSelf%isMaster().and.self%outputIterations) then
          write (label,'(i6)') countIterations
          iterationOutputGroup=outputGroup%openGroup('iteration'//trim(adjustl(label)),'Data for iteration '//trim(adjustl(label)))
          call self                %computationalDomain_%output        (iterationOutputGroup            )
          call iterationOutputGroup                     %writeAttribute(converged           ,'converged')
       end if
       ! If converged and we need to use a different number of photons in the final iteration, do that now.       
       finalIteration=.not.finalIteration .and. countIterations >= self%countIterationsMinimum .and. converged .and. self%countPhotonsPerWavelengthFinalIteration /= self%countPhotonsPerWavelength
    end do iterations
    ! Output convergence status.
    if (mpiSelf%isMaster()) call outputGroup%writeAttribute(converged,'converged')
    ! Output the computational domain.
    if (mpiSelf%isMaster()) call self%computationalDomain_       %output(outputGroup)
    ! Output outputter results.
    call self%radiativeTransferOutputter_%finalize()
    if (mpiSelf%isMaster()) call self%radiativeTransferOutputter_%output(outputGroup)
    ! Done.
    call timerTotal_%stop()
    if (present(status)) status=errorStatusSuccess
    if (mpiSelf%isMaster()) call displayUnindent('Done task: radiative transfer ['//timerTotal_%reportText()//']')
    return
  end subroutine radiativeTransferPerform
