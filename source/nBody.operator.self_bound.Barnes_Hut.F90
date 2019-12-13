!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements an N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !# <nbodyOperator name="nbodyOperatorSelfBoundBarnesHut">
  !#  <description>An N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorSelfBound) :: nbodyOperatorSelfBoundBarnesHut
     !% An N-body data operator which determines the subset of particles that are self-bound. The potential is computed using a tree method following \cite{barnes_hierarchical_1986}.
     private
     double precision :: thetaTolerance
   contains
     procedure :: operate => selfBoundBarnesHutOperate
  end type nbodyOperatorSelfBoundBarnesHut

  interface nbodyOperatorSelfBoundBarnesHut
     !% Constructors for the ``selfBoundBarnesHut'' N-body operator class.
     module procedure selfBoundBarnesHutConstructorParameters
     module procedure selfBoundBarnesHutConstructorInternal
  end interface nbodyOperatorSelfBoundBarnesHut

contains

  function selfBoundBarnesHutConstructorParameters(parameters) result (self)
    !% Constructor for the ``selfBoundBarnesHut'' N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSelfBoundBarnesHut)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: thetaTolerance

    !# <inputParameter>
    !#   <name>thetaTolerance</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.5d0</defaultValue>
    !#   <description>The criterion for the opening angle.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self%nbodyOperatorSelfBound=nbodyOperatorSelfBound(parameters)
    self%thetaTolerance        =thetaTolerance
    !# <inputParametersValidate source="parameters"/>
    return
  end function selfBoundBarnesHutConstructorParameters

  function selfBoundBarnesHutConstructorInternal(thetaTolerance,tolerance,bootstrapSampleCount,bootstrapSampleRate,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_) result (self)
    !% Internal constructor for the ``selfBoundBarnesHut'' N-body operator class
    implicit none
    type            (nbodyOperatorSelfBoundBarnesHut)                        :: self
    double precision                                 , intent(in   )         :: thetaTolerance        , tolerance           , &
         &                                                                      bootstrapSampleRate
    integer         (c_size_t                       ), intent(in   )         :: bootstrapSampleCount
    logical                                          , intent(in   )         :: analyzeAllParticles   , useVelocityMostBound
    class           (randomNumberGeneratorClass     ), intent(in   ), target :: randomNumberGenerator_
    !# <constructorAssign variables="thetaTolerance"/>

    self%nbodyOperatorSelfBound=nbodyOperatorSelfBound(tolerance,bootstrapSampleCount,bootstrapSampleRate,analyzeAllParticles,useVelocityMostBound,randomNumberGenerator_)
    return
  end function selfBoundBarnesHutConstructorInternal
  
  subroutine selfBoundBarnesHutOperate(self,simulation)
    !% Determine the subset of N-body particles which are self-bound.
    use :: Galacticus_Display          , only : Galacticus_Display_Message
    use :: Galacticus_Error            , only : Galacticus_Error_Report
    use :: ISO_Varying_String          , only : varying_string                 , assignment(=)  , operator(//)
    use :: String_Handling             , only : operator(//)
    use :: Memory_Management           , only : allocateArray                  , deallocateArray
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    use :: Octree_Data_Structure       , only : octreeData
    use :: Sort                        , only : Sort_Index_Do
    implicit none
    class           (nbodyOperatorSelfBoundBarnesHut), intent(inout)                          :: self
    type            (nBodyData                      ), intent(inout)                          :: simulation
    type            (octreeData                     )                                         :: octreePosition          , octreeVelocity
    integer                                          , parameter                              :: countIterationMaximum=30
    logical                                          , allocatable  , dimension(:,:)          :: isBound                 , isBoundNew             , &
         &                                                                                       isRemoved
    integer                                          , allocatable  , dimension(:,:)          :: boundStatus
    double precision                                 , allocatable  , dimension(:,:)          :: positionOffset          , positionRescaled
    double precision                                 , allocatable  , dimension(:,:)          :: energyPotential         , velocityPotential      , &
         &                                                                                       energyKinetic           , sampleWeight           , &
         &                                                                                       sampleWeightActual
    double precision                                 , allocatable  , dimension(:,:)          :: velocityCenterOfMass
    double precision                                                , dimension(3  )          :: velocityRepresentative
    integer         (c_size_t                       ), allocatable  , dimension(:  )          :: indexMostBound          , indexVelocityMostBound , &
         &                                                                                       indexSorted             , indexSortedPrevious
    integer         (c_size_t                       )                                         :: particleCount           , i                      , &
         &                                                                                       k                       , iSample
    integer         (c_size_t                       ), allocatable  , dimension(:  )          :: countBound              , countBoundPrevious
    double precision                                 , allocatable  , dimension(:  )          :: weightBound             , weightBoundPrevious
    logical                                          , allocatable  , dimension(:  )          :: isConverged
    logical                                                                                   :: isPreviousSnapshotAvail
    integer                                                                                   :: countIteration
    type            (varying_string                 )                                         :: message

    ! Allocate workspaces.
    particleCount=size(simulation%position,dim=2)
    call allocateArray(isBound                ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(isBoundNew             ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(isRemoved              ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(energyKinetic          ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(energyPotential        ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(velocityPotential      ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(sampleWeight           ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(sampleWeightActual     ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(velocityCenterOfMass   ,[3_c_size_t              ,self%bootstrapSampleCount])
    call allocateArray(positionOffset         ,[3_c_size_t,particleCount                          ])
    call allocateArray(positionRescaled       ,[3_c_size_t,particleCount                          ])
    call allocateArray(boundStatus            ,[           particleCount,self%bootstrapSampleCount])
    call allocateArray(indexMostBound         ,[                         self%bootstrapSampleCount])
    call allocateArray(indexVelocityMostBound ,[                         self%bootstrapSampleCount])
    call allocateArray(indexSorted            ,[           particleCount                          ])
    call allocateArray(indexSortedPrevious    ,[           particleCount                          ])
    call allocateArray(countBound             ,[                         self%bootstrapSampleCount])
    call allocateArray(countBoundPrevious     ,[                         self%bootstrapSampleCount])
    call allocateArray(weightBound            ,[                         self%bootstrapSampleCount])
    call allocateArray(weightBoundPrevious    ,[                         self%bootstrapSampleCount])
    call allocateArray(isConverged            ,[                         self%bootstrapSampleCount])
    ! Check whether the self-bound status from the previous snapshot is available. If it is, read in the self-bound status
    ! and sanpling weights. If not, generate new values.
    isPreviousSnapshotAvail=.false.
    if (allocated(simulation%boundStatusPrevious)) then
       ! Read in the self-bound status from the previous snapshot.
       isPreviousSnapshotAvail=.true.
       if (self%bootstrapSampleCount /= size(simulation%boundStatusPrevious,dim=2)) then
          call Galacticus_Error_Report('The number of bootstrap samples is not consistent with the previous snapshot.'//{introspection:location})
       end if
       ! Sort particles according to their particle IDs.
       indexSortedPrevious=Sort_Index_Do(simulation%ParticleIDsPrevious)
       indexSorted        =Sort_Index_Do(simulation%ParticleIDs        )
       !$omp parallel do private(i)
       do i=1,particleCount
          sampleWeight(indexSorted(i),:)=dble(simulation%sampleWeightPrevious(indexSortedPrevious(i),:))
          isBound     (indexSorted(i),:)=     simulation% boundStatusPrevious(indexSortedPrevious(i),:) > 0
       end do
       !$omp end parallel do
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
    ! Rescle the particle coordinates by the softening length.
    positionRescaled=simulation%position/simulation%lengthSoftening
    ! Iterate over bootstrap samplings.
    do iSample=1,self%bootstrapSampleCount
       message='Performing self-bound analysis on bootstrap sample '
       message=message//iSample//' of '//self%bootstrapSampleCount
       call Galacticus_Display_Message(message)
       ! Initialize count of bound particles.
       countBoundPrevious (iSample)=count(                             isBound(:,iSample))
       weightBoundPrevious(iSample)=sum  (sampleWeight(:,iSample),mask=isBound(:,iSample))
       ! Compute the center-of-mass velocity.
       forall(k=1:3)
          velocityCenterOfMass(k,iSample)=sum(simulation%velocity(k,:)*sampleWeight(:,iSample),mask=isBound(:,iSample)) &
               &                          /weightBoundPrevious(iSample)
       end forall
       ! If the self-bound status from the previous snapshot is used, reset the self-boud status. All particles in the
       ! sample are assumed to be self-bound at the beginning of the first iteration.
       if (isPreviousSnapshotAvail .and. self%analyzeAllParticles) then
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
       message='Build octrees.'
       call Galacticus_Display_Message(message)
       call octreePosition%build(positionRescaled   ,sampleWeightActual(:,iSample))
       call octreeVelocity%build(simulation%velocity,sampleWeightActual(:,iSample))
       ! Begin iterations.
       countIteration=0
       do while (.true.)
          countIteration=countIteration+1
          !$omp parallel private(i,k)
          !$omp do schedule(dynamic)
          do i=1,particleCount
             ! Skip particles we do not need to analyze.
             if (.not.isBound(i,iSample)) cycle
             call octreePosition%traverseCompute(positionRescaled   (:,i),sampleWeightActual(i,iSample),self%thetaTolerance,energyPotential  (i,iSample),selfBoundBarnesHutPotential)
             call octreeVelocity%traverseCompute(simulation%velocity(:,i),sampleWeightActual(i,iSample),self%thetaTolerance,velocityPotential(i,iSample),selfBoundBarnesHutPotential)
          end do
          !$omp end do
          !$omp workshare
          ! Apply constant multipliers to potential energy.
          where(isBound(:,iSample))
             energyPotential(:,iSample)=+energyPotential(:,iSample)                                &
                  &                     *gravitationalConstantGalacticus                           &
                  &                     *simulation%massParticle                                   &
                  &                     /self%bootstrapSampleRate                                  &
                  &                     /simulation%lengthSoftening                                &
                  &                     * dble(particleCount)                                      &
                  &                     /(dble(particleCount-1_c_size_t)-self%bootstrapSampleRate)
          end where
          ! Find the index of the most bound particle.
          indexMostBound        (iSample) =minloc(energyPotential  (:,iSample),dim=1,mask=isBound(:,iSample))
          ! Find the index of the most bound particle in velocity space.
          indexVelocityMostBound(iSample) =minloc(velocityPotential(:,iSample),dim=1,mask=isBound(:,iSample))
          !$omp end workshare
          ! Check whether we should use the velocity of the most bound particle in velocity space as the
          ! representative velocity of the satellite. If not, use the center-of-mass velocity instead.
          if (self%useVelocityMostBound) then
             velocityRepresentative=simulation%velocity (:,indexVelocityMostBound(iSample))
          else
             velocityRepresentative=velocityCenterOfMass(:,iSample)
          end if
          !$omp workshare
          ! Compute kinetic energies.
          forall(k=1:3)
             where(isBound(:,iSample))
                positionOffset(k,:)=+simulation%velocity(k,:)-velocityRepresentative(k)
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
             velocityCenterOfMass(k,iSample)=sum(simulation%velocity(k,:)*sampleWeight(:,iSample),mask=isBoundNew(:,iSample)) &
                  &                          /weightBound(iSample)
          end forall
          !$omp end workshare
          !$omp end parallel
          if (abs(weightBound(iSample)-weightBoundPrevious(iSample)) < 0.5d0*self%tolerance                 &
               &                                                            *(                              &
               &                                                              +weightBound        (iSample) &
               &                                                              +weightBoundPrevious(iSample) &
               &                                                             )) then
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
                   call octreePosition%removeParticle(positionRescaled   (:,i),sampleWeight(i,iSample))
                   call octreeVelocity%removeParticle(simulation%velocity(:,i),sampleWeight(i,iSample))
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
          ! Check for excess iterations.
          if (countIteration > countIterationMaximum) call Galacticus_Error_Report('maximum iterations exceeded'//{introspection:location})
          if (isConverged(iSample)) exit
       end do
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
    ! Write indices of most bound particles to file.
    call simulation%analysis%writeDataset(indexMostBound        ,'indexMostBound'        )
    call simulation%analysis%writeDataset(indexVelocityMostBound,'indexVelocityMostBound')
    ! Write bound status to file.
    call simulation%analysis%writeDataset(boundStatus           ,'selfBoundStatus'       )
    call simulation%analysis%writeDataset(nint(sampleWeight)    ,'weight'                )
    ! Free workspaces.
    call deallocateArray(isBound                )
    call deallocateArray(isBoundNew             )
    call deallocateArray(isRemoved              )
    call deallocateArray(energyKinetic          )
    call deallocateArray(energyPotential        )
    call deallocateArray(velocityPotential      )
    call deallocateArray(boundStatus            )
    call deallocateArray(sampleWeight           )
    call deallocateArray(sampleWeightActual     )
    call deallocateArray(velocityCenterOfMass   )
    call deallocateArray(indexMostBound         )
    call deallocateArray(indexVelocityMostBound )
    call deallocateArray(indexSorted            )
    call deallocateArray(indexSortedPrevious    )
    call deallocateArray(positionOffset         )
    call deallocateArray(positionRescaled       )
    call deallocateArray(countBound             )
    call deallocateArray(countBoundPrevious     )
    call deallocateArray(weightBound            )
    call deallocateArray(weightBoundPrevious    )
    call deallocateArray(isConverged            )
    return
  end subroutine selfBoundBarnesHutOperate

  subroutine selfBoundBarnesHutPotential(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquared)
    !% Compute the potential given the separation between a particle and a node in the octree. Currently assumes the functional form of the softening used by
    !% Gadget.
    implicit none
    double precision              , intent(inout) :: value
    double precision, dimension(3), intent(in   ) :: centerOfMass     , relativePosition
    double precision,               intent(in   ) :: nodeWeight       , separation      , &
         &                                           separationSquared
    double precision                              :: potential
    !GCC$ attributes unused :: centerOfMass, relativePosition

    if     (separation == 0.0d0) then
       potential=0.0d0 ! No self-energy.
    else if (separation <= 0.5d0) then
       potential=-14.0d0              &
            &                      /5.0d0               &
            &                      +separationSquared   &
            &                      *(                   &
            &                        +16.0d0            &
            &                        /3.0d0             &
            &                        +separationSquared &
            &                        *(                 &
            &                          -48.0d0          &
            &                          / 5.0d0          &
            &                          +32.0d0          &
            &                          *separation      &
            &                          /5.0d0           &
            &                         )                 &
            &                       )
    else if (separation <= 1.0d0) then
       potential=-1.0d0                 &
            &                      +(                     &
            &                        +1.0d0               &
            &                        +separation          &
            &                        *(                   &
            &                          -33.0d0            &
            &                          +separationSquared &
            &                          *(                 &
            &                            +160.0d0         &
            &                            +separation      &
            &                            *(               &
            &                              -240.0d0       &
            &                              +separation    &
            &                              *(             &
            &                                +144.0d0     &
            &                                -32.0d0      &
            &                                *separation  &
            &                               )             &
            &                             )               &
            &                           )                 &
            &                         )                   &
            &                       )                     &
            &                      /15.0d0                &
            &                      /separation
    else
       potential=-1.0d0                &
            &                      /separation
    end if
    value=value+nodeWeight*potential
    return
  end subroutine selfBoundBarnesHutPotential
