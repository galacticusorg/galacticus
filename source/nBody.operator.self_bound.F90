!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either versteeion 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements an N-body data operator which determines the subset of particles that are self-bound.

  use, intrinsic :: ISO_C_Binding

  !# <nbodyOperator name="nbodyOperatorSelfBound">
  !#  <description>An N-body data operator which determines the subset of particles that are self-bound.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorSelfBound
     !% An N-body data operator which determines the subset of particles that are self-bound.
     private
     integer         (c_size_t) :: bootstrapSampleCount
     double precision           :: tolerance
   contains
     procedure :: operate => selfBoundOperate
  end type nbodyOperatorSelfBound

  interface nbodyOperatorSelfBound
     !% Constructors for the ``selfBound'' N-body operator class.
     module procedure selfBoundConstructorParameters
     module procedure selfBoundConstructorInternal
  end interface nbodyOperatorSelfBound

contains

  function selfBoundConstructorParameters(parameters) result (self)
    !% Constructor for the ``selfBound'' N-body operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (nbodyOperatorSelfBound)                :: self
    type            (inputParameters       ), intent(inout) :: parameters
    integer         (c_size_t              )                :: bootstrapSampleCount
    double precision                                        :: tolerance
    
    !# <inputParameter>
    !#   <name>bootstrapSampleCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>30_c_size_t</defaultValue>
    !#   <description>The number of bootstrap resamples of the particles that should be used.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>tolerance</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d-2</defaultValue>
    !#   <description>The tolerance in the summed weight of bound particles which must be attained to declare convergence.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyOperatorSelfBound(tolerance,bootstrapSampleCount)
    return
  end function selfBoundConstructorParameters

  function selfBoundConstructorInternal(tolerance,bootstrapSampleCount) result (self)
    !% Internal constructor for the ``selfBound'' N-body operator class
    use Input_Parameters2
    implicit none
    type            (nbodyOperatorSelfBound)                :: self
    double precision                        , intent(in   ) :: tolerance
    integer         (c_size_t              ), intent(in   ) :: bootstrapSampleCount
    !# <constructorAssign variables="tolerance, bootstrapSampleCount"/>

    return
  end function selfBoundConstructorInternal

  subroutine selfBoundOperate(self,simulation)
    !% Determine the subset of N-body particles which are self-bound.
    use FGSL
    use Memory_Management
    use Numerical_Constants_Physical
    use Galacticus_Error
    use Poisson_Random
    implicit none
    class           (nbodyOperatorSelfBound), intent(inout)                          :: self
    type            (nBodyData             ), intent(inout)                          :: simulation
    integer                                 , parameter                              :: countIterationMaximum=30
    double precision                        , parameter                              :: sampleRate           = 1.0d0
    logical                                 , allocatable  , dimension(:  )          :: isBound                      , compute                , &
         &                                                                              isBoundNew                   , isBoundCompute
    integer                                 , allocatable  , dimension(:,:), target  :: boundStatus
    integer                                                , dimension(:  ), pointer :: boundStatusSingle
    double precision                        , allocatable  , dimension(:,:)          :: positionRelative
    double precision                        , allocatable  , dimension(:  )          :: separation                   , potential              , &
         &                                                                              energyPotential              , separationSquared      , &
         &                                                                              velocityPotential            , energyKinetic          , &
         &                                                                              energyPotentialChange        , velocityPotentialChange, &
         &                                                                              sampleWeight
    integer         (c_size_t              ), allocatable  , dimension(:  )          :: indexMostBound               , indexVelocityMostBound
    integer         (c_size_t              )                                         :: particleCount                , i                      , &
         &                                                                              k                            , iSample                , &
         &                                                                              countBound                   , countBoundPrevious
    double precision                                                                 :: weightBound                  , weightBoundPrevious
    integer                                                                          :: addSubtract                  , countIteration
    type            (fgsl_rng              )                                         :: pseudoSequenceObject
    logical                                                                          :: pseudoSequenceReset   =.true.

    ! Allocate workspaces.
    particleCount=size(simulation%position,dim=2)
    call allocateArray(isBound                ,[particleCount                          ])    
    call allocateArray(isBoundNew             ,[particleCount                          ])    
    call allocateArray(isBoundCompute         ,[particleCount                          ])    
    call allocateArray(energyKinetic          ,[particleCount                          ])    
    call allocateArray(energyPotential        ,[particleCount                          ])    
    call allocateArray(velocityPotential      ,[particleCount                          ])
    call allocateArray(compute                ,[particleCount                          ])
    call allocateArray(energyPotentialChange  ,[particleCount                          ])
    call allocateArray(velocityPotentialChange,[particleCount                          ])
    call allocateArray(sampleWeight           ,[particleCount                          ])
    call allocateArray(boundStatus            ,[particleCount,self%bootstrapSampleCount])
    call allocateArray(indexMostBound         ,[              self%bootstrapSampleCount])
    call allocateArray(indexVelocityMostBound ,[              self%bootstrapSampleCount])
    ! Iterate over bootstrap samplings.
    do iSample=1,self%bootstrapSampleCount
       ! Determine weights for particles.
       do i=1,particleCount
          sampleWeight(i)=dble(Poisson_Random_Get(pseudoSequenceObject,sampleRate,pseudoSequenceReset))
       end do
       ! Initialize count of bound particles.
       countBoundPrevious =     particleCount
       weightBoundPrevious=dble(particleCount)
       ! Initialize potentials.
       energyPotential  =0.0d0
       velocityPotential=0.0d0
       isBound          =sampleWeight > 0.0d0
       isBoundCompute   =isBound
       compute          =isBound
       addSubtract      =+1
       ! Begin iterations.
       countIteration=0
       do while (.true.)
          countIteration=countIteration+1
          !$omp parallel private(i,k,positionRelative,separationSquared,separation,potential)
          call allocateArray(positionRelative ,[3_c_size_t,particleCount])
          call allocateArray(separation       ,[           particleCount])
          call allocateArray(separationSquared,[           particleCount])
          call allocateArray(potential        ,[           particleCount])
          separation       =0.0d0
          separationSquared=0.0d0
          positionRelative =0.0d0
          potential        =0.0d0
           !$omp do schedule(dynamic) reduction(+: energyPotentialChange, velocityPotentialChange)
          do i=1,particleCount-1
             ! Skip particles we do not need to compute.
             if (.not.compute(i)) cycle
             ! Compute "potentials" in velocity space in order to find a particle which best represents the velocity of the particles.
             ! Find particle velocity separations.
             forall(k=1:3)
                where(isBoundCompute(i+1:particleCount))
                   positionRelative(k,i+1:particleCount)=+simulation%velocity(k,i+1:particleCount) &
                        &                                -simulation%velocity(k,i                )
                end where
             end forall
             where(isBoundCompute(i+1:particleCount))
                separationSquared(i+1:particleCount)=+sum (positionRelative (:,i+1:particleCount)**2,dim=1)
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
             ! Accumulate potential energy.
             velocityPotentialChange(i                )=velocityPotentialChange(i                )+sum(potential(i+1:particleCount)*sampleWeight(i+1:particleCount))
             velocityPotentialChange(i+1:particleCount)=velocityPotentialChange(i+1:particleCount)+    potential(i+1:particleCount)*sampleWeight(i                )
             ! Compute gravitational potential energies.
             ! Find particle separations.
             forall(k=1:3)
                where(isBoundCompute(i+1:particleCount))
                   positionRelative(k,i+1:particleCount)=+simulation%position(k,i+1:particleCount) &
                        &                                -simulation%position(k,i                )
                end where
             end forall
             where(isBoundCompute(i+1:particleCount))
                separationSquared(i+1:particleCount)=+sum (positionRelative (:,i+1:particleCount)**2,dim=1) &
                     &                               /simulation%lengthSoftening
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
             ! Accumulate potential energy.
             energyPotentialChange(i                )=energyPotentialChange(i                )+sum(potential(i+1:particleCount)*sampleWeight(i+1:particleCount))
             energyPotentialChange(i+1:particleCount)=energyPotentialChange(i+1:particleCount)+    potential(i+1:particleCount)*sampleWeight(i                )
          end do
          !$omp end do
          !$omp workshare
          ! Apply constant multipliers to potential energy.
          where(isBoundCompute)
             energyPotentialChange=+energyPotentialChange                       &
                  &                *gravitationalConstantGalacticus             &
                  &                *simulation%massParticle                     &
                  &                /simulation%lengthSoftening                  &
                  &                * dble(particleCount)                        &
                  &                /(dble(particleCount-1_c_size_t)-sampleRate)
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
          !$omp workshare
          ! Find the index of the most bound particle.
          indexMostBound        (iSample)=minloc(energyPotential  ,dim=1,mask=isBound)
          ! Find the index of the most bound particle in velocity space.
          indexVelocityMostBound(iSample)=minloc(velocityPotential,dim=1,mask=isBound)
          ! Compute kinetic energies.
          forall(k=1:3)
             where(isBound)
                positionRelative(k,:)=+simulation%velocity(k,                             : ) &
                     &                -simulation%velocity(k,indexVelocityMostBound(iSample))
             end where
          end forall
          where(isBound)
             energyKinetic=0.5d0*sum(positionRelative**2,dim=1)
          end where
          ! Determine which particles are bound.
          where (isBound)
             isBoundNew=(energyKinetic+energyPotential < 0.0d0)
          elsewhere
             isBoundNew=.false.
          end where
          ! Count bound particles.
          countBound =count(                  isBoundNew)
          weightBound=sum  (sampleWeight,mask=isBoundNew)
          !$omp end workshare
          ! Free workspaces.
          call deallocateArray(positionRelative       )
          call deallocateArray(separation             )
          call deallocateArray(separationSquared      )
          call deallocateArray(potential              )
          !$omp end parallel
          ! Test for convergence.
          if (abs(weightBound-weightBoundPrevious) < 0.5d0*self%tolerance*(weightBound+weightBoundPrevious)) then
             ! Converged.
             isBound=isBoundNew
             exit
          else
             ! Not converged. Decide if it is faster to recompute potentials fully for the new bound set, or to subtract potential due
             ! to particles which are now unbound.
             if (2*(countBoundPrevious-countBound) < countBound+1) then
                ! Number of particles that have become unbound is sufficiently small that it is faster to simply subtract off their
                ! contribution to the potential.
                compute       =isBound.and..not.isBoundNew
                isBoundCompute=isBound
                addSubtract   =-1
             else
                ! Number of particles that have become unbound is sufficiently large that it will be faster to simply recompute
                ! potentials completely.          
                compute       =isBoundNew
                isBoundCompute=isBoundNew
                addSubtract   =+1
                ! Reset potentials.
                energyPotential  =0.0d0
                velocityPotential=0.0d0
             end if
             ! Update bound status.
             isBound            =isBoundNew
             countBoundPrevious =countBound
             weightBoundPrevious=weightBound
          end if
          ! Check for excess iterations.
          if (countIteration > countIterationMaximum) call Galacticus_Error_Report('selfBoundOperate','maximum iterations exceeded')
       end do
       ! Write bound status to file.
       boundStatusSingle => boundStatus(:,iSample)
       where (isBound)
          boundStatusSingle=nint(sampleWeight)
       elsewhere
          boundStatusSingle=0
       end where
    end do
    ! Write indices of most bound particles to file.
    call simulation%analysis%writeDataset(indexMostBound        ,'indexMostBound'        )
    call simulation%analysis%writeDataset(indexVelocityMostBound,'indexVelocityMostBound')
    ! Write bound status to file.
    call simulation%analysis%writeDataset(boundStatus           ,'selfBoundStatus'       )
    ! Free workspaces.
    call deallocateArray(compute                )
    call deallocateArray(isBound                )
    call deallocateArray(isBoundNew             )
    call deallocateArray(isBoundCompute         )
    call deallocateArray(energyKinetic          )
    call deallocateArray(energyPotential        )
    call deallocateArray(velocityPotential      )
    call deallocateArray(boundStatus            )
    call deallocateArray(energyPotentialChange  )
    call deallocateArray(velocityPotentialChange)
    call deallocateArray(indexMostBound         )
    call deallocateArray(indexVelocityMostBound )
    return
  end subroutine selfBoundOperate

  pure function selfBoundPotential(separation,separationSquared)
    !% Compute the potential for an array of particle separations. Currently assumes the functional form of the softening used by
    !% Gadget.
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

