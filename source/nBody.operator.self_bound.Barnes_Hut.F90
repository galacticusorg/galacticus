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
Contains a module which implements an N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorSelfBoundBarnesHut">
   <description>An N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorSelfBound) :: nbodyOperatorSelfBoundBarnesHut
     !!{
     An N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.
     !!}
     private
     double precision :: thetaTolerance
   contains
     procedure :: operate => selfBoundBarnesHutOperate
  end type nbodyOperatorSelfBoundBarnesHut

  interface nbodyOperatorSelfBoundBarnesHut
     !!{
     Constructors for the ``selfBoundBarnesHut'' N-body operator class.
     !!}
     module procedure selfBoundBarnesHutConstructorParameters
     module procedure selfBoundBarnesHutConstructorInternal
  end interface nbodyOperatorSelfBoundBarnesHut

contains

  function selfBoundBarnesHutConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``selfBoundBarnesHut'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSelfBoundBarnesHut)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: thetaTolerance

    !![
    <inputParameter>
      <name>thetaTolerance</name>
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <description>The criterion for the opening angle.</description>
    </inputParameter>
    !!]
    self%nbodyOperatorSelfBound=nbodyOperatorSelfBound(parameters)
    self%thetaTolerance        =thetaTolerance
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function selfBoundBarnesHutConstructorParameters

  function selfBoundBarnesHutConstructorInternal(thetaTolerance,tolerance,bootstrapSampleCount,bootstrapSampleRate,representativeMinimumCount,representativeFraction,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the ``selfBoundBarnesHut'' N-body operator class
    !!}
    implicit none
    type            (nbodyOperatorSelfBoundBarnesHut)                        :: self
    double precision                                 , intent(in   )         :: thetaTolerance        , tolerance                 , &
         &                                                                      bootstrapSampleRate   , representativeFraction
    integer         (c_size_t                       ), intent(in   )         :: bootstrapSampleCount  , representativeMinimumCount
    logical                                          , intent(in   )         :: analyzeAllParticles   , useVelocityMostBound
    class           (randomNumberGeneratorClass     ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="thetaTolerance"/>
    !!]

    self%nbodyOperatorSelfBound=nbodyOperatorSelfBound(tolerance,bootstrapSampleCount,bootstrapSampleRate,representativeMinimumCount,representativeFraction,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_)
    return
  end function selfBoundBarnesHutConstructorInternal
  
  subroutine selfBoundBarnesHutOperate(self,simulations)
    !!{
    Determine the subset of N-body particles which are self-bound.
    !!}
    use :: Display                         , only : displayIndent                 , displayUnindent  , displayMessage
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : var_str
    use :: String_Handling                 , only : operator(//)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Octree_Data_Structure           , only : octreeData
    use :: Sorting                         , only : sortIndex                     , sortSmallestIndex
    implicit none
    class           (nbodyOperatorSelfBoundBarnesHut), intent(inout)                 :: self
    type            (nBodyData                      ), intent(inout), dimension(:  ) :: simulations
    type            (octreeData                     )                                :: octreePosition                , octreeVelocity
    integer                                          , parameter                     :: countIterationMaximum      =30
    logical                                          , allocatable  , dimension(:,:) :: isBound                       , isBoundNew            , &
         &                                                                              isRemoved
    integer         (c_size_t                       ), pointer      , dimension(:,:) :: boundStatus                   , boundStatusPrevious
    double precision                                 , pointer      , dimension(:,:) :: position                      , velocity              , &
         &                                                                              sampleWeight                  , sampleWeightPrevious
    integer         (c_size_t                       ), pointer      , dimension(:  ) :: particleIDs                   , particleIDsPrevious
    double precision                                 , allocatable  , dimension(:,:) :: positionOffset                , positionRescaled      , &
         &                                                                              velocityRescaled
    double precision                                 , allocatable  , dimension(:,:) :: energyPotential               , velocityPotential     , &
         &                                                                              energyKinetic                 , sampleWeightActual
    double precision                                 , allocatable  , dimension(:,:) :: velocityCenterOfMass
    double precision                                                , dimension(3  ) :: velocityRepresentative
    integer         (c_size_t                       ), pointer      , dimension(:,:) :: indexMostBound                , indexVelocityMostBound
    integer         (c_size_t                       ), pointer      , dimension(:  ) :: indexSorted                   , indexSortedPrevious
    integer         (c_size_t                       )                                :: particleCount                 , i                     , &
         &                                                                              k                             , iSample               , &
         &                                                                              current                       , previous
    integer         (c_size_t                       ), allocatable  , dimension(:  ) :: countBound                    , countBoundPrevious
    double precision                                 , allocatable  , dimension(:  ) :: weightBound                   , weightBoundPrevious
    logical                                          , allocatable  , dimension(:  ) :: isConverged
    integer         (c_size_t                       )                                :: representativeParticleCount
    integer                                                                          :: countIteration
    double precision                                                                 :: lengthSoftening               , velocitySoftening     , &
         &                                                                              massParticle                  , convergenceFactor     , &
         &                                                                              weightRepresentative
    type            (varying_string                 )                                :: message
    character       (len=12                         )                                :: label

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
    allocate(isBound               (           particleCount,self%bootstrapSampleCount))
    allocate(isBoundNew            (           particleCount,self%bootstrapSampleCount))
    allocate(isRemoved             (           particleCount,self%bootstrapSampleCount))
    allocate(energyKinetic         (           particleCount,self%bootstrapSampleCount))
    allocate(energyPotential       (           particleCount,self%bootstrapSampleCount))
    allocate(velocityPotential     (           particleCount,self%bootstrapSampleCount))
    allocate(sampleWeight          (           particleCount,self%bootstrapSampleCount))
    allocate(sampleWeightActual    (           particleCount,self%bootstrapSampleCount))
    allocate(velocityCenterOfMass  (3_c_size_t              ,self%bootstrapSampleCount))
    allocate(positionOffset        (3_c_size_t,particleCount                          ))
    allocate(positionRescaled      (3_c_size_t,particleCount                          ))
    allocate(velocityRescaled      (3_c_size_t,particleCount                          ))
    allocate(boundStatus           (           particleCount,self%bootstrapSampleCount))
    allocate(indexSorted           (           particleCount                          ))
    allocate(countBound            (                         self%bootstrapSampleCount))
    allocate(countBoundPrevious    (                         self%bootstrapSampleCount))
    allocate(weightBound           (                         self%bootstrapSampleCount))
    allocate(weightBoundPrevious   (                         self%bootstrapSampleCount))
    allocate(isConverged           (                         self%bootstrapSampleCount))
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
    isConverged = .false.
    ! Rescale the particle coordinates by the softening length.
    positionRescaled=position/lengthSoftening
    velocityRescaled=velocity/velocitySoftening
    call displayIndent('Performing self-bound analysis on bootstrap samples')
    do iSample=1,self%bootstrapSampleCount
       ! Initialize count of bound particles.
       countBoundPrevious (iSample)=count(                             isBound(:,iSample))
       weightBoundPrevious(iSample)=sum  (sampleWeight(:,iSample),mask=isBound(:,iSample))
    end do
    ! Determine the number of representative particles.
    representativeParticleCount=int(self%representativeFraction*minval(countBoundPrevious))
    representativeParticleCount=max(representativeParticleCount,self%representativeMinimumCount)
    allocate(indexMostBound        (representativeParticleCount,self%bootstrapSampleCount))
    allocate(indexVelocityMostBound(representativeParticleCount,self%bootstrapSampleCount))
    ! Iterate over bootstrap samples.
    do iSample=1,self%bootstrapSampleCount
       call displayIndent(var_str('sample ')//iSample//' of '//self%bootstrapSampleCount)
       ! Compute the center-of-mass velocity.
       forall(k=1:3)
          velocityCenterOfMass(k,iSample)=sum(velocity(k,:)*sampleWeight(:,iSample),mask=isBound(:,iSample)) &
               &                          /weightBoundPrevious(iSample)
       end forall
       ! If the self-bound status from the previous snapshot is used, reset the self-bound status. All particles in the
       ! sample are assumed to be self-bound at the beginning of the first iteration.
       if (previous > 0 .and. self%analyzeAllParticles) then
         isBound        (:,iSample)      =sampleWeight(:,iSample) > 0.0d0
       end if
       ! Initialize potentials.
       energyPotential  (:,iSample)      =0.0d0
       velocityPotential(:,iSample)      =0.0d0
       where(isBound(:,iSample))
          sampleWeightActual(:,iSample)  =sampleWeight(:,iSample)
       else where
          sampleWeightActual(:,iSample)  =0.0d0
       end where
       ! Build octrees.
       call displayIndent('build octrees')
       call octreePosition%build(positionRescaled,sampleWeightActual(:,iSample))
       call octreeVelocity%build(velocityRescaled,sampleWeightActual(:,iSample))
       call displayUnindent('done')
       ! Begin iterations.
       countIteration=0
       do while (.true.)
          countIteration=countIteration+1
          !$omp parallel private(i,k,velocityRepresentative,weightRepresentative)
          !$omp do schedule(dynamic)
          do i=1,particleCount
             ! Skip particles we do not need to analyze.
             if (.not.isBound(i,iSample)) cycle
             call octreePosition%traverseCompute(positionRescaled(:,i),sampleWeightActual(i,iSample),self%thetaTolerance,energyPotential  (i,iSample),selfBoundBarnesHutPotential)
             call octreeVelocity%traverseCompute(velocityRescaled(:,i),sampleWeightActual(i,iSample),self%thetaTolerance,velocityPotential(i,iSample),selfBoundBarnesHutPotential)
          end do
          !$omp end do
          !$omp workshare
          ! Apply constant multipliers to potential energy.
          where(isBound(:,iSample))
             energyPotential(:,iSample)=+energyPotential(:,iSample)                                &
                  &                     *gravitationalConstant_internal                            &
                  &                     *massParticle                                              &
                  &                     /self%bootstrapSampleRate                                  &
                  &                     /lengthSoftening                                           &
                  &                     * dble(particleCount)                                      &
                  &                     /(dble(particleCount-1_c_size_t)-self%bootstrapSampleRate)
          end where
          !$omp end workshare
          !$omp single
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
          !$omp end parallel
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
             isConverged   (  iSample)=.true.
          else
             ! Not converged.
             ! Reset potentials.
             energyPotential  (:,iSample)=0.0d0
             velocityPotential(:,iSample)=0.0d0
             ! Remove unbound particles from the octree.
             isRemoved        (:,iSample)=(isBound(:,iSample) .and. .not.isBoundNew(:,iSample))
             do i=1,particleCount
                if (isRemoved(i,iSample)) then
                   call octreePosition%removeParticle(positionRescaled(:,i),sampleWeight(i,iSample))
                   call octreeVelocity%removeParticle(velocityRescaled(:,i),sampleWeight(i,iSample))
                end if
             end do
             ! Update bound status.
             isBound            (:,iSample)=isBoundNew (:,iSample)
             countBoundPrevious (  iSample)=countBound (  iSample)
             weightBoundPrevious(  iSample)=weightBound(  iSample)
             where(isBound(:,iSample))
                sampleWeightActual(:,iSample)=sampleWeight(:,iSample)
             else where
                sampleWeightActual(:,iSample)=0.0d0
             end where
          end if
          write (label,'(e12.6)') convergenceFactor
          message=var_str('iteration ')//countIteration//' convergence factor = '//trim(adjustl(label))
          call displayMessage(message)
          ! Check for excess iterations.
          if (countIteration > countIterationMaximum) call Error_Report('maximum iterations exceeded'//{introspection:location})
          if (isConverged(iSample)) exit
       end do
       call displayUnindent('done')
       ! Destroy the octrees.
       call octreePosition%destroy()
       call octreeVelocity%destroy()
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
    deallocate(isBound                )
    deallocate(isBoundNew             )
    deallocate(isRemoved              )
    deallocate(energyKinetic          )
    deallocate(energyPotential        )
    deallocate(velocityPotential      )
    deallocate(sampleWeightActual     )
    deallocate(velocityCenterOfMass   )
    nullify   (indexMostBound         )
    nullify   (indexVelocityMostBound )
    deallocate(positionOffset         )
    deallocate(positionRescaled       )
    deallocate(velocityRescaled       )
    deallocate(countBound             )
    deallocate(countBoundPrevious     )
    deallocate(weightBound            )
    deallocate(weightBoundPrevious    )
    deallocate(isConverged            )
    return
  end subroutine selfBoundBarnesHutOperate

  subroutine selfBoundBarnesHutPotential(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquared)
    !!{
    Compute the potential given the separation between a particle and a node in the octree. Currently assumes the functional form of the softening used by
    Gadget.
    !!}
    implicit none
    double precision              , intent(inout) :: value
    double precision, dimension(3), intent(in   ) :: centerOfMass     , relativePosition
    double precision,               intent(in   ) :: nodeWeight       , separation      , &
         &                                           separationSquared
    double precision                              :: potential
    !$GLC attributes unused :: centerOfMass, relativePosition

    if     (separation == 0.0d0) then
       potential=0.0d0 ! No self-energy.
    else if (separation <= 0.5d0) then
       potential=-14.0d0              &
            &    /5.0d0               &
            &    +separationSquared   &
            &    *(                   &
            &      +16.0d0            &
            &      /3.0d0             &
            &      +separationSquared &
            &      *(                 &
            &        -48.0d0          &
            &        / 5.0d0          &
            &        +32.0d0          &
            &        *separation      &
            &        /5.0d0           &
            &       )                 &
            &      )
    else if (separation <= 1.0d0) then
       potential=-1.0d0                 &
            &    +(                     &
            &      +1.0d0               &
            &      +separation          &
            &      *(                   &
            &        -33.0d0            &
            &        +separationSquared &
            &        *(                 &
            &          +160.0d0         &
            &          +separation      &
            &          *(               &
            &            -240.0d0       &
            &            +separation    &
            &            *(             &
            &              +144.0d0     &
            &              -32.0d0      &
            &              *separation  &
            &             )             &
            &           )               &
            &         )                 &
            &       )                   &
            &     )                     &
            &    /15.0d0                &
            &    /separation
    else
       potential=-1.0d0                &
            &     /separation
    end if
    value=value+nodeWeight*potential
    return
  end subroutine selfBoundBarnesHutPotential
