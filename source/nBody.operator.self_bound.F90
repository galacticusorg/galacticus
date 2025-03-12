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
Contains a module which implements an N-body data operator which determines the subset of particles that are self-bound.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorSelfBound">
   <description>An N-body data operator which determines the subset of particles that are self-bound.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSelfBound
     !!{
     An N-body data operator which determines the subset of particles that are self-bound.
     !!}
     private
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     integer         (c_size_t                  )          :: bootstrapSampleCount            , representativeMinimumCount
     double precision                                      :: tolerance                       , bootstrapSampleRate       , &
          &                                                   representativeFraction
     logical                                               :: analyzeAllParticles             , useVelocityMostBound
   contains
     final     ::            selfBoundDestructor
     procedure :: operate => selfBoundOperate
  end type nbodyOperatorSelfBound

  interface nbodyOperatorSelfBound
     !!{
     Constructors for the ``selfBound'' N-body operator class.
     !!}
     module procedure selfBoundConstructorParameters
     module procedure selfBoundConstructorInternal
  end interface nbodyOperatorSelfBound

contains

  function selfBoundConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``selfBound'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSelfBound    )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    integer         (c_size_t                  )                :: bootstrapSampleCount  , representativeMinimumCount
    double precision                                            :: tolerance             , bootstrapSampleRate       , &
         &                                                         representativeFraction
    logical                                                     :: analyzeAllParticles   , useVelocityMostBound

    !![
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <inputParameter>
      <name>representativeMinimumCount</name>
      <source>parameters</source>
      <defaultValue>10_c_size_t</defaultValue>
      <description>Minimum number of representative particles used to compute the center of a halo.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerance</name>
      <source>parameters</source>
      <defaultValue>1.0d-2</defaultValue>
      <description>The tolerance in the summed weight of bound particles which must be attained to declare convergence.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleRate</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The sampling rate for particles.</description>
    </inputParameter>
    <inputParameter>
      <name>representativeFraction</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>Fraction of bound particles used to compute the center of a halo.</description>
    </inputParameter>
    <inputParameter>
      <name>analyzeAllParticles</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, all particles are assumed to be self-bound at the beginning of the analysis. Unbound particles at previous times are allowed to become bound in the current snapshot. If false and the self-bound information from the previous snapshot is available, only the particles that are self-bound at the previous snapshot are assumed to be bound at the beginning of the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>useVelocityMostBound</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the velocity of the most bound particle in velocity space is used as the representative velocity of the satellite. If false, use the mass weighted mean velocity (center-of-mass velocity) of self-bound particles instead.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorSelfBound(tolerance,bootstrapSampleCount,bootstrapSampleRate,representativeMinimumCount,representativeFraction,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function selfBoundConstructorParameters

  function selfBoundConstructorInternal(tolerance,bootstrapSampleCount,bootstrapSampleRate,representativeMinimumCount,representativeFraction,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the ``selfBound'' N-body operator class
    !!}
    implicit none
    type            (nbodyOperatorSelfBound    )                        :: self
    double precision                            , intent(in   )         :: tolerance             , bootstrapSampleRate       , &
         &                                                                 representativeFraction
    integer         (c_size_t                  ), intent(in   )         :: bootstrapSampleCount  , representativeMinimumCount
    logical                                     , intent(in   )         :: analyzeAllParticles   , useVelocityMostBound
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="tolerance, bootstrapSampleCount, bootstrapSampleRate, representativeMinimumCount, representativeFraction, analyzeAllParticles, useVelocityMostBound, *randomNumberGenerator_"/>
    !!]
    if (representativeMinimumCount < 1) then
       call Error_Report('"representativeMinimumCount" must be larger than or equal to 1.'//{introspection:location})
    end if

    return
  end function selfBoundConstructorInternal

  subroutine selfBoundDestructor(self)
    !!{
    Destructor for the ``selfBound'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSelfBound), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine selfBoundDestructor
  
  subroutine selfBoundOperate(self,simulations)
    !!{
    Determine the subset of N-body particles which are self-bound.
    !!}
    use :: Display                         , only : displayIndent                 , displayUnindent  , displayMessage
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : var_str
    use :: String_Handling                 , only : operator(//)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Sorting                         , only : sortIndex                     , sortSmallestIndex
    implicit none
    class           (nbodyOperatorSelfBound), intent(inout)                          :: self
    type            (nBodyData             ), intent(inout), dimension(:  )          :: simulations
    integer                                 , parameter                              :: countIterationMaximum      =30
    logical                                 , allocatable  , dimension(:,:)          :: isBound                       , isBoundNew            , &
         &                                                                              isBoundCompute                , isBoundComputeActual
    logical                                 , allocatable  , dimension(:  )          :: compute                       , computeActual
    integer         (c_size_t              ), pointer      , dimension(:,:)          :: boundStatus                   , boundStatusPrevious
    double precision                        , pointer      , dimension(:,:)          :: position                      , velocity              , &
         &                                                                              sampleWeight                  , sampleWeightPrevious
    integer         (c_size_t              ), pointer      , dimension(:  )          :: particleIDs                   , particleIDsPrevious
    double precision                        , allocatable  , dimension(:,:)          :: positionRelative              , positionOffset
    double precision                        , allocatable  , dimension(:  )          :: separation                    , potential             , &
         &                                                                              separationSquared             , potentialActual
    double precision                        , allocatable  , dimension(:,:)          :: energyPotential               , velocityPotential     , &
         &                                                                              energyKinetic                 , energyPotentialChange , &
         &                                                                              velocityPotentialChange
    double precision                        , allocatable  , dimension(:,:)          :: velocityCenterOfMass
    double precision                                       , dimension(3  )          :: velocityRepresentative
    integer         (c_size_t              ), pointer      , dimension(:,:)          :: indexMostBound                , indexVelocityMostBound
    integer         (c_size_t              ), pointer      , dimension(:  )          :: indexSorted                   , indexSortedPrevious
    integer         (c_size_t              )                                         :: particleCount                 , i                     , &
         &                                                                              k                             , iSample               , &
         &                                                                              current                       , previous
    integer         (c_size_t              ), allocatable  , dimension(:  )          :: countBound                    , countBoundPrevious
    double precision                        , allocatable  , dimension(:  )          :: weightBound                   , weightBoundPrevious
    logical                                 , allocatable  , dimension(:  )          :: isConverged
    integer         (c_size_t              )                                         :: representativeParticleCount
    integer                                                                          :: addSubtract                   , countIteration
    double precision                                                                 :: lengthSoftening               , velocitySoftening     , &
         &                                                                              massParticle                  , convergenceFactor     , &
         &                                                                              weightRepresentative
    type            (varying_string        )                                         :: message
    character       (len=12                )                                         :: label

    ! Determine current and previous simulations.
    if      (size(simulations) == 1) then
       ! Single simulation.
       current =+1
       previous=-1
    else if (size(simulations) == 2) then
       ! Two simulations - one should be active, the other previous. Previous simulation must provide the "isBound" property.
       current =-1
       previous=-1
       do i=1,2
          if (simulations(i)%label == "active"  ) current =i
          if (simulations(i)%label == "previous") previous=i
       end do
       if (current  == -1) call Error_Report('no "active" simulation found'  //{introspection:location})
       if (previous == -1) call Error_Report('no "previous" simulation found'//{introspection:location})
       if (.not.simulations(previous)%propertiesIntegerRank1%exists('isBound')) &
            & call Error_Report('"previous" simulation must provide the "isBound" property'//{introspection:location})
    else
       current =-1
       previous=-1
       call Error_Report('either 1 or 2 simulations (labelled "active" and "previous" in the case of 2 simulations) should be provided'//{introspection:location})
    end if
    ! Get simulation attributes.
    lengthSoftening  =simulations(current)%attributesReal%value('lengthSoftening')
    massParticle     =simulations(current)%attributesReal%value('massParticle'   )
    velocitySoftening=sqrt(gravitationalConstant_internal*massParticle/lengthSoftening)
    ! Get particle data.
    position    => simulations(current)%propertiesRealRank1%value('position'  )
    velocity    => simulations(current)%propertiesRealRank1%value('velocity'  )
    particleIDs => simulations(current)%propertiesInteger  %value('particleID')
    ! Allocate workspaces.
    particleCount=size(position,dim=2)
    allocate(isBound                (           particleCount,self%bootstrapSampleCount))
    allocate(isBoundNew             (           particleCount,self%bootstrapSampleCount))
    allocate(isBoundCompute         (           particleCount,self%bootstrapSampleCount))
    allocate(energyKinetic          (           particleCount,self%bootstrapSampleCount))
    allocate(energyPotential        (           particleCount,self%bootstrapSampleCount))
    allocate(velocityPotential      (           particleCount,self%bootstrapSampleCount))
    allocate(compute                (           particleCount                          ))
    allocate(energyPotentialChange  (           particleCount,self%bootstrapSampleCount))
    allocate(velocityPotentialChange(           particleCount,self%bootstrapSampleCount))
    allocate(sampleWeight           (           particleCount,self%bootstrapSampleCount))
    allocate(velocityCenterOfMass   (3_c_size_t              ,self%bootstrapSampleCount))
    allocate(positionOffset         (3_c_size_t,particleCount                          ))
    allocate(boundStatus            (           particleCount,self%bootstrapSampleCount))
    allocate(indexSorted            (           particleCount                          ))
    allocate(countBound             (                         self%bootstrapSampleCount))
    allocate(countBoundPrevious     (                         self%bootstrapSampleCount))
    allocate(weightBound            (                         self%bootstrapSampleCount))
    allocate(weightBoundPrevious    (                         self%bootstrapSampleCount))
    allocate(isConverged            (                         self%bootstrapSampleCount))
    ! Iterate over bootstrap samples.
    call displayIndent('Performing self-bound analysis on bootstrap samples')
    ! If previous bound status is available, read in the self-bound status and sampling weights. If not, generate new values.
    if (previous > 0) then
       ! Get the self-bound status from the previous snapshot.
       particleIDsPrevious  => simulations(previous)%propertiesInteger     %value('particleID'  )
       boundStatusPrevious  => simulations(previous)%propertiesIntegerRank1%value('isBound'     )
       sampleWeightPrevious => simulations(previous)%propertiesRealRank1   %value('sampleWeight')
       if (self%bootstrapSampleCount /= size(boundStatusPrevious,dim=2)) &
            & call Error_Report('The number of bootstrap samples is not consistent with the previous snapshot.'//{introspection:location})
       ! Sort particles according to their particle IDs.
       if (simulations(previous)%propertiesInteger%exists('particleOrder')) then
          indexSortedPrevious => simulations(previous)%propertiesInteger%value('particleOrder')
       else
          allocate(indexSortedPrevious(particleCount))
          indexSortedPrevious=sortIndex(particleIDsPrevious)
       end if
       indexSorted           =sortIndex(particleIDs        )
       !$omp parallel do private(i)
       do i=1,particleCount
          sampleWeight(indexSorted(i),:)=dble(sampleWeightPrevious(indexSortedPrevious(i),:))
          isBound     (indexSorted(i),:)=      boundStatusPrevious(indexSortedPrevious(i),:) > 0
       end do
       !$omp end parallel do
       nullify(particleIDsPrevious )
       nullify(boundStatusPrevious )
       nullify(sampleWeightPrevious)
       if (simulations(previous)%propertiesInteger%exists('particleOrder')) then
          nullify   (indexSortedPrevious)
       else
          deallocate(indexSortedPrevious)
       end if
    else
       ! Generate new sampling weights.
       do iSample=1,self%bootstrapSampleCount
          ! Determine weights for particles.
          do i=1,particleCount
             sampleWeight(i,iSample)=dble(self%randomNumberGenerator_%poissonSample(self%bootstrapSampleRate))
          end do
          isBound(:,iSample)=sampleWeight(:,iSample) > 0.0d0
       end do
    end if
    compute=.false.
    do iSample=1,self%bootstrapSampleCount
       ! Initialize count of bound particles.
       countBoundPrevious (iSample)=count(                             isBound(:,iSample))
       weightBoundPrevious(iSample)=sum  (sampleWeight(:,iSample),mask=isBound(:,iSample))
       ! Compute the center-of-mass velocity.
       forall(k=1:3)
          velocityCenterOfMass(k,iSample)=sum(velocity(k,:)*sampleWeight(:,iSample),mask=isBound(:,iSample)) &
               &                          /weightBoundPrevious(iSample)
       end forall
       ! If the self-bound status from the previous snapshot is used, reset the self-bound status. All particles in the sample are
       ! assumed to be self-bound at the beginning of the first iteration.
       if (previous > 0 .and. self%analyzeAllParticles) then
         isBound        (:,iSample)      =sampleWeight(:,iSample) > 0.0d0
       end if
       ! Initialize potentials.
       energyPotential  (:,iSample)      =0.0d0
       velocityPotential(:,iSample)      =0.0d0
       isBoundCompute   (:,iSample)      =             isBound(:,iSample)
       compute                           =compute .or. isBound(:,iSample)
    end do
    ! Determine the number of representative particles.
    representativeParticleCount=int(self%representativeFraction*minval(countBoundPrevious))
    representativeParticleCount=max(representativeParticleCount,self%representativeMinimumCount)
    allocate(indexMostBound        (representativeParticleCount,self%bootstrapSampleCount))
    allocate(indexVelocityMostBound(representativeParticleCount,self%bootstrapSampleCount))
    addSubtract =+1
    isConverged = .false.
    ! Begin iterations.
    countIteration=0
    do while (.true.)
       countIteration         =countIteration+1
       velocityPotentialChange=0.0d0
       energyPotentialChange  =0.0d0
       call displayIndent(var_str('iteration ')//countIteration)
       !$omp parallel private(i,k,positionRelative,separationSquared,separation,potential,potentialActual,computeActual,isBoundComputeActual,velocityRepresentative,weightRepresentative)
       allocate(positionRelative    (3_c_size_t,particleCount                          ))
       allocate(separation          (           particleCount                          ))
       allocate(separationSquared   (           particleCount                          ))
       allocate(potential           (           particleCount                          ))
       allocate(potentialActual     (           particleCount                          ))
       allocate(computeActual       (           particleCount                          ))
       allocate(isBoundComputeActual(           particleCount,self%bootstrapSampleCount))
       separation       =0.0d0
       separationSquared=0.0d0
       positionRelative =0.0d0
       potential        =0.0d0
       potentialActual  =0.0d0
       !$omp do schedule(dynamic) reduction(+: energyPotentialChange, velocityPotentialChange)
       do i=1,particleCount-1
          ! Skip particles we do not need to compute.
          if (.not.compute(i)) cycle
          ! Check how to accumulate changes to potentials.
          if (addSubtract == +1) then
             ! Recompute potentials for all self-bound particles.
             computeActual        = compute
             isBoundComputeActual = isBoundCompute
          else
             ! Subtract potentials due to particles which are marked as unbound in the previous iteration.
             computeActual = .false.
             do iSample=1,self%bootstrapSampleCount
                where(isBoundCompute(:,iSample))
                   isBoundComputeActual(:,iSample) = isBound(:,iSample) .neqv. isBound(i,iSample)
                else where
                   isBoundComputeActual(:,iSample) = .false.
                end where
                computeActual = computeActual .or. isBoundComputeActual(:,iSample)
             end do
          end if
          ! Compute "potentials" in velocity space in order to find a particle which best represents the velocity of the particles.
          ! Find particle velocity separations.
          forall(k=1:3)
             where(computeActual(i+1:particleCount))
                positionRelative(k,i+1:particleCount)=+velocity(k,i+1:particleCount) &
                     &                                -velocity(k,i                )
             end where
          end forall
          where(computeActual(i+1:particleCount))
             separationSquared(i+1:particleCount)=+sum (positionRelative (:,i+1:particleCount)**2,dim=1) &
                  &                               /velocitySoftening**2
             separation       (i+1:particleCount)=+sqrt(separationSquared(  i+1:particleCount)         )
             ! Compute potentials.
             potential        (i+1:particleCount)=+selfBoundPotential(                                      &
                  &                                                   separation       (i+1:particleCount), &
                  &                                                   separationSquared(i+1:particleCount)  &
                  &                                                  )
          elsewhere
             ! For particles not participating in this computation, set potential to zero.
             potential        (i+1:particleCount)=0.0d0
          end where
          do iSample=1,self%bootstrapSampleCount
             ! Skip bootstrap samples for which we have already got converged results.
             if (isConverged(iSample)) cycle
             if (isBoundCompute(i,iSample)) then
                ! Compute "potentials" needed in the analysis of a particular sample.
                where(isBoundComputeActual(i+1:particleCount,iSample))
                   potentialActual(i+1:particleCount) = potential(i+1:particleCount)
                else where
                   potentialActual(i+1:particleCount) = 0.0d0
                end where
                ! Accumulate potential energy.
                velocityPotentialChange(i                ,iSample)=+velocityPotentialChange(i                ,iSample)                             &
                     &                                             +sum(potentialActual(i+1:particleCount)*sampleWeight(i+1:particleCount,iSample))
                velocityPotentialChange(i+1:particleCount,iSample)=+velocityPotentialChange(i+1:particleCount,iSample)                             &
                     &                                             +    potentialActual(i+1:particleCount)*sampleWeight(i                ,iSample)
             end if
          end do
          ! Compute gravitational potential energies.
          ! Find particle separations.
          forall(k=1:3)
             where(computeActual(i+1:particleCount))
                positionRelative(k,i+1:particleCount)=+position(k,i+1:particleCount) &
                     &                                -position(k,i                )
             end where
          end forall
          where(computeActual(i+1:particleCount))
             separationSquared(i+1:particleCount)=+sum (positionRelative (:,i+1:particleCount)**2,dim=1) &
                  &                               /lengthSoftening**2
             separation       (i+1:particleCount)=+sqrt(separationSquared(  i+1:particleCount)         )
             ! Compute potentials.
             potential        (i+1:particleCount)=+selfBoundPotential(                                      &
                  &                                                   separation       (i+1:particleCount), &
                  &                                                   separationSquared(i+1:particleCount)  &
                  &                                                  )
          elsewhere
             ! For particles not participating in this computation, set potential to zero.
             potential        (i+1:particleCount)=0.0d0
          end where
          do iSample=1,self%bootstrapSampleCount
             ! Skip bootstrap samples for which we have already got converged results.
             if (isConverged(iSample)) cycle
             if (isBoundCompute(i,iSample)) then
                ! Compute potentials needed in the analysis of a particular sample.
                where(isBoundComputeActual(i+1:particleCount,iSample))
                   potentialActual(i+1:particleCount) = potential(i+1:particleCount)
                else where
                   potentialActual(i+1:particleCount) = 0.0d0
                end where
                ! Accumulate potential energy.
                energyPotentialChange(i                ,iSample)=+energyPotentialChange(i                ,iSample)                               &
                     &                                           +sum(potentialActual(i+1:particleCount)*sampleWeight(i+1:particleCount,iSample))
                energyPotentialChange(i+1:particleCount,iSample)=+energyPotentialChange(i+1:particleCount,iSample)                               &
                     &                                           +    potentialActual(i+1:particleCount)*sampleWeight(i                ,iSample)
             end if
          end do
       end do
       !$omp end do
       !$omp workshare
       ! Apply constant multipliers to potential energy.
       where(isBoundCompute)
          energyPotentialChange=+energyPotentialChange                                     &
               &                *gravitationalConstant_internal                            &
               &                *massParticle                                              &
               &                /self%bootstrapSampleRate                                  &
               &                /lengthSoftening                                           &
               &                * dble(particleCount)                                      &
               &                /(dble(particleCount-1_c_size_t)-self%bootstrapSampleRate)
       end where
       !$omp end workshare
       ! Accumulate changes to potentials.
       if (addSubtract == +1) then
          !$omp workshare
          where(isBoundCompute)
             energyPotential  =  energyPotential+  energyPotentialChange
             velocityPotential=velocityPotential+velocityPotentialChange
          end where
          !$omp end workshare
       else
          !$omp workshare
          where(isBoundCompute)
             energyPotential  =  energyPotential-  energyPotentialChange
             velocityPotential=velocityPotential-velocityPotentialChange
          end where
          !$omp end workshare
       end if
       do iSample=1,self%bootstrapSampleCount
          ! Skip bootstrap samples for which we have already got converged results.
          if (isConverged(iSample)) cycle
          !$omp single
          ! Find the indices of the k most bound particles.
          if (countBoundPrevious(iSample) >= representativeParticleCount) then
             indexMostBound        (:,iSample)=sortSmallestIndex(energyPotential  (:,iSample),representativeParticleCount,mask=isBound(:,iSample))
             indexVelocityMostBound(:,iSample)=sortSmallestIndex(velocityPotential(:,iSample),representativeParticleCount,mask=isBound(:,iSample))
          else
             call Error_Report('There are not enough bound particles to estimate the center of the halo.'//{introspection:location})
          end if
          !$omp end single
          ! Check whether we should use the velocity of the most bound particle in velocity space as the
          ! representative velocity of the satellite. If not, use the center-of-mass velocity instead.
          if (self%useVelocityMostBound) then
             velocityRepresentative=0.0d0
             weightRepresentative  =0.0d0
             do k=1,representativeParticleCount
                velocityRepresentative=+velocityRepresentative                                              &
                     &                 +velocity              (:,indexVelocityMostBound(k,iSample)        ) &
                     &                 *sampleWeight          (  indexVelocityMostBound(k,iSample),iSample)
                weightRepresentative  =+weightRepresentative                                                &
                     &                 +sampleWeight          (  indexVelocityMostBound(k,iSample),iSample)
             end do
             velocityRepresentative=velocityRepresentative/weightRepresentative
          else
             velocityRepresentative=velocityCenterOfMass(:,iSample)
          end if
          !$omp workshare
          ! Compute kinetic energies.
          forall(k=1:3)
             where(isBound(:,iSample))
                positionOffset(k,:)=+velocity(k,:)-velocityRepresentative(k)
             end where
          end forall
          where(isBound(:,iSample))
             energyKinetic(:,iSample)=0.5d0*sum(positionOffset**2,dim=1)
          end where
          ! Determine which particles are bound.
          where (isBound(:,iSample))
             isBoundNew(:,iSample)=(energyKinetic(:,iSample)+energyPotential(:,iSample) < 0.0d0)
          elsewhere
             isBoundNew(:,iSample)=.false.
          end where
          ! Count bound particles.
          countBound (iSample)=count(                             isBoundNew(:,iSample))
          weightBound(iSample)=sum  (sampleWeight(:,iSample),mask=isBoundNew(:,iSample))
          ! Compute the center-of-mass velocity.
          forall(k=1:3)
             velocityCenterOfMass(k,iSample)=sum(velocity(k,:)*sampleWeight(:,iSample),mask=isBoundNew(:,iSample)) &
                  &                          /weightBound(iSample)
          end forall
          !$omp end workshare
       end do
       ! Free workspaces.
       deallocate(positionRelative    )
       deallocate(separation          )
       deallocate(separationSquared   )
       deallocate(potential           )
       deallocate(potentialActual     )
       deallocate(computeActual       )
       deallocate(isBoundComputeActual)
       !$omp end parallel
       ! Decide if it is faster to recompute potentials fully for the new bound set, or to subtract potential due
       ! to particles which are now unbound.
       if (2*sum(countBoundPrevious-countBound) < sum(countBound+1)) then
          addSubtract   =-1
       else
          addSubtract   =+1
       end if
       ! Test for convergence.
       compute = .false.
       do iSample=1,self%bootstrapSampleCount
          convergenceFactor=+2.0d0                             &
               &            *abs(                              &
               &                 +weightBound        (iSample) &
               &                 -weightBoundPrevious(iSample) &
               &                )                              &
               &            /   (                              &
               &                 +weightBound        (iSample) &
               &                 +weightBoundPrevious(iSample) &
               &                )
          if (convergenceFactor < self%tolerance) then
             ! Converged.
             isBound       (:,iSample)=isBoundNew(:,iSample)
             isBoundCompute(:,iSample)=.false.
             isConverged   (  iSample)=.true.
          else
             ! Not converged.
             if (addSubtract==-1) then
                ! Number of particles that have become unbound is sufficiently small that it is faster to simply subtract off their
                ! contribution to the potential.
                compute                     =compute .or. isBound(:,iSample)
                isBoundCompute(:,iSample)   =             isBound(:,iSample)
             else
                ! Number of particles that have become unbound is sufficiently large that it will be faster to simply recompute
                ! potentials completely.
                compute                     =compute .or. isBoundNew(:,iSample)
                isBoundCompute(:,iSample)   =             isBoundNew(:,iSample)
                ! Reset potentials.
                energyPotential  (:,iSample)=0.0d0
                velocityPotential(:,iSample)=0.0d0
             end if
             ! Update bound status.
             isBound            (:,iSample)=isBoundNew (:,iSample)
             countBoundPrevious (  iSample)=countBound (  iSample)
             weightBoundPrevious(  iSample)=weightBound(  iSample)
          end if
          write (label,'(e12.6)') convergenceFactor
          message=var_str('sample ')//iSample//' convergence factor = '//trim(adjustl(label))
          call displayMessage(message)
          ! Check for excess iterations.
          if (countIteration > countIterationMaximum) call Error_Report('maximum iterations exceeded'//{introspection:location})
       end do
       call displayUnindent('done')
       if (count(isConverged)==self%bootstrapSampleCount) exit
    end do
    ! Write bound status to file.
    !$omp parallel workshare
    where (isBound)
       boundStatus=nint(sampleWeight)
    elsewhere
       boundStatus=0
    end where
    !$omp end parallel workshare
    call displayUnindent('done')
    ! Store the self bound status, the particle order and the index of the most bound particle.
    call simulations(current)%propertiesIntegerRank1%set         ('isBound'                 , boundStatus            )
    call simulations(current)%propertiesRealRank1   %set         ('sampleWeight'            , sampleWeight           )
    call simulations(current)%propertiesInteger     %set         ('particleOrder'           , indexSorted            )
    call simulations(current)%propertiesIntegerRank1%set         ('indexMostBound'          , indexMostBound         )
    call simulations(current)%propertiesIntegerRank1%set         ('indexVelocityMostBound'  , indexVelocityMostBound )
    ! Write indices of most bound particles to file.
    call simulations(current)%analysis              %writeDataset( indexMostBound           ,'indexMostBound'        )
    call simulations(current)%analysis              %writeDataset( indexVelocityMostBound   ,'indexVelocityMostBound')
    ! Write bound status to file.
    call simulations(current)%analysis              %writeDataset( boundStatus              ,'selfBoundStatus'       )
    call simulations(current)%analysis              %writeDataset( nint(sampleWeight)       ,'weight'                )
    call simulations(current)%analysis              %writeDataset([self%bootstrapSampleRate],'bootstrapSampleRate'   )
    ! Free workspaces.
    nullify   (boundStatus            )
    nullify   (sampleWeight           )
    nullify   (indexSorted            )
    deallocate(compute                )
    deallocate(isBound                )
    deallocate(isBoundNew             )
    deallocate(isBoundCompute         )
    deallocate(energyKinetic          )
    deallocate(energyPotential        )
    deallocate(velocityPotential      )
    deallocate(energyPotentialChange  )
    deallocate(velocityPotentialChange)
    deallocate(velocityCenterOfMass   )
    nullify   (indexMostBound         )
    nullify   (indexVelocityMostBound )
    deallocate(positionOffset         )
    deallocate(countBound             )
    deallocate(countBoundPrevious     )
    deallocate(weightBound            )
    deallocate(weightBoundPrevious    )
    deallocate(isConverged            )
    return
  end subroutine selfBoundOperate

  pure function selfBoundPotential(separation,separationSquared)
    !!{
    Compute the potential for an array of particle separations. Currently assumes the functional form of the softening used by
    Gadget.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:               ) :: separation        , separationSquared
    double precision               , dimension(size(separation)) :: selfBoundPotential

    where     (separation == 0.0d0)
       selfBoundPotential=0.0d0 ! No self-energy.
    elsewhere (separation <= 0.5d0)
       selfBoundPotential=-14.0d0              &
            &             /5.0d0               &
            &             +separationSquared   &
            &             *(                   &
            &               +16.0d0            &
            &               /3.0d0             &
            &               +separationSquared &
            &               *(                 &
            &                 -48.0d0          &
            &                 / 5.0d0          &
            &                 +32.0d0          &
            &                 *separation      &
            &                 /5.0d0           &
            &                )                 &
            &              )
    elsewhere (separation <= 1.0d0)
       selfBoundPotential=-1.0d0                 &
            &             +(                     &
            &               +1.0d0               &
            &               +separation          &
            &               *(                   &
            &                 -33.0d0            &
            &                 +separationSquared &
            &                 *(                 &
            &                   +160.0d0         &
            &                   +separation      &
            &                   *(               &
            &                     -240.0d0       &
            &                     +separation    &
            &                     *(             &
            &                       +144.0d0     &
            &                       -32.0d0      &
            &                       *separation  &
            &                      )             &
            &                    )               &
            &                  )                 &
            &                )                   &
            &              )                     &
            &             /15.0d0                &
            &             /separation
    elsewhere
       selfBoundPotential=-1.0d0                &
            &             /separation
    end where
    return
  end function selfBoundPotential

