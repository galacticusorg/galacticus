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
Implements an N-body data operator which determines the kinetic and Chandrasekhar potential energy tensors.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <nbodyOperator name="nbodyOperatorEnergyTensors">
   <description>An N-body data operator which determines the kinetic and Chandrasekhar potential energy tensors.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorEnergyTensors
     !!{
     An N-body data operator which determines the kinetic and Chandrasekhar potential energy tensors.
     !!}
     private
     logical                                      :: selfBoundParticlesOnly
     integer(c_size_t                  )          :: bootstrapSampleCount
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            energyTensorsDestructor
     procedure :: operate => energyTensorsOperate
  end type nbodyOperatorEnergyTensors

  interface nbodyOperatorEnergyTensors
     !!{
     Constructors for the \refClass{nbodyOperatorEnergyTensors} N-body operator class.
     !!}
     module procedure energyTensorsConstructorParameters
     module procedure energyTensorsConstructorInternal
  end interface nbodyOperatorEnergyTensors

contains

  function energyTensorsConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorEnergyTensors} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyOperatorEnergyTensors)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    class  (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    logical                                            :: selfBoundParticlesOnly
    integer(c_size_t                  )                :: bootstrapSampleCount

    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the mean position and velocity are computed only for self-bound particles</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorEnergyTensors(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function energyTensorsConstructorParameters

  function energyTensorsConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorEnergyTensors} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorEnergyTensors)                        :: self
    logical                            , intent(in   )         :: selfBoundParticlesOnly
    integer(c_size_t                  ), intent(in   )         :: bootstrapSampleCount
    class  (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, *randomNumberGenerator_"/>
    !!]

    return
  end function energyTensorsConstructorInternal

  subroutine energyTensorsDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorEnergyTensors} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorEnergyTensors), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine energyTensorsDestructor
  
  subroutine energyTensorsOperate(self,simulations)
    !!{
    Determine the kinetic and potential energy tensors of N-body particles.
    !!}
    use    :: Display                         , only : displayCounter                , displayCounterClear   , displayIndent, displayMessage, &
         &                                             displayUnindent               , verbosityLevelStandard
    use    :: Error                           , only : Error_Report
    use    :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use    :: Linear_Algebra                  , only : matrix                        , vector                , operator(*)  , assignment(=)
#ifdef USEMPI
    use    :: MPI_Utilities                   , only : mpiSelf
#endif
    !$ use :: OMP_Lib                         , only : OMP_Get_Thread_Num
    implicit none
    class          (nbodyOperatorEnergyTensors), intent(inout)                   :: self
    type           (nBodyData                 ), intent(inout), dimension(:  )   :: simulations
    double precision                           , parameter                       :: sampleRate                  =1.0d0
    double precision                           , allocatable  , dimension(:,:,:) :: energyTensorKinetic               , energyTensorPotential          , &
         &                                                                          energyTensorKineticPrincipal      , energyTensorPotentialPrincipal
    double precision                           , pointer      , dimension(:,:  ) :: position                          , velocity                       , &
         &                                                                          separations                       , angularVelocityVector          , &
         &                                                                          positionMean                      , velocityMean
    integer         (c_size_t                 ), pointer      , dimension(:,:  ) :: selfBoundStatus
    double precision                           , allocatable  , dimension(:    ) :: separationsCubed                  , separationsCubedInverseWeighted
    double precision                           , allocatable  , dimension(:,:  ) :: positionOffset                    , velocityOffset
    integer         (c_size_t                 )                                  :: i                                 , j                              , &
         &                                                                          k                                 , m                              , &
         &                                                                          iSimulation
    double precision                                                             :: massParticle                      , lengthSoftening
    type            (matrix                   )                                  :: energyTensor                      , eigenVectors                   , &
         &                                                                          energyTensorKinetic_              , energyTensorPotential_
    type            (vector                   )                                  :: eigenValues
    
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
    call displayIndent('compute energy tensors',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    do iSimulation=1,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayMessage(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Determine the particle mask to use.
       if (self%selfBoundParticlesOnly) then
          if (simulations(iSimulation)%propertiesIntegerRank1%exists('isBound')) then
             selfBoundStatus => simulations(iSimulation)%propertiesIntegerRank1%value('isBound')
             if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of selfBoundStatus samples must equal number of requested bootstrap samples'//{introspection:location})
          else
             call Error_Report('self-bound status not available - apply a self-bound operator first'//{introspection:location})
          end if
       else
          position => simulations(iSimulation)%propertiesRealRank1%value('position')
          allocate(selfBoundStatus(size(position,dim=2,kind=c_size_t),self%bootstrapSampleCount))
          do i=1,self%bootstrapSampleCount
             do j=1,size(position,dim=2)
                selfBoundStatus(j,i)=self%randomNumberGenerator_%poissonSample(sampleRate)
             end do
          end do
       end if
       ! Get simulation attributes.
       lengthSoftening=simulations(iSimulation)%attributesReal%value('lengthSoftening')
       massParticle   =simulations(iSimulation)%attributesReal%value('massParticle'   )
       ! Get mean position and velocity, and the angular velocity vector of the rotating frame in which the angular momentum of the system is zero.
       if (simulations(iSimulation)%propertiesRealRank1%exists('positionMean')) then
          positionMean => simulations(iSimulation)%propertiesRealRank1%value('positionMean')          
          if (size(positionMean,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of positionMean samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('mean position not available - apply a "positionMean" operator first'//{introspection:location})
       end if
       if (simulations(iSimulation)%propertiesRealRank1%exists('velocityMean')) then
          velocityMean => simulations(iSimulation)%propertiesRealRank1%value('velocityMean')          
          if (size(velocityMean,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of velocityMean samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('mean velocity not available - apply a "positionMean" operator first'//{introspection:location})
       end if
       if (simulations(iSimulation)%propertiesRealRank1%exists('angularVelocityVector')) then
          angularVelocityVector => simulations(iSimulation)%propertiesRealRank1%value('angularVelocityVector')          
          if (size(angularVelocityVector,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of angularVelocityVector samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('mean angularVelocityVector not available - apply a "angularMomentum" operator first'//{introspection:location})
       end if
       ! Get position and velocity.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       velocity => simulations(iSimulation)%propertiesRealRank1%value('velocity')
       ! Compute the kinetic and potential energy tensors.
       allocate(positionOffset                (3_c_size_t           ,size(position,dim=2)     ))
       allocate(velocityOffset                (3_c_size_t           ,size(position,dim=2)     ))
       allocate(energyTensorKinetic           (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(energyTensorPotential         (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(energyTensorKineticPrincipal  (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(energyTensorPotentialPrincipal(3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       energyTensorPotential=0.0d0
       energyTensorKinetic  =0.0d0
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       do i=1,self%bootstrapSampleCount
          !$omp parallel private(j,k,m,separations,separationsCubed,separationsCubedInverseWeighted) reduction(+:energyTensorKinetic,energyTensorPotential)
          allocate(separations                    (3,size(position,dim=2)))
          allocate(separationsCubed               (  size(position,dim=2)))
          allocate(separationsCubedInverseWeighted(  size(position,dim=2)))
          !$omp workshare
          forall(j=1:3)
             positionOffset(j,:)=position(j,:)-positionMean(j,i)
             velocityOffset(j,:)=velocity(j,:)-velocityMean(j,i)
          end forall
          velocityOffset(1,:)=velocityOffset(1,:)+(+angularVelocityVector(2,i)*positionOffset(3,:)-angularVelocityVector(3,i)*positionOffset(2,:))
          velocityOffset(2,:)=velocityOffset(2,:)+(-angularVelocityVector(1,i)*positionOffset(3,:)+angularVelocityVector(3,i)*positionOffset(1,:))
          velocityOffset(3,:)=velocityOffset(3,:)+(+angularVelocityVector(1,i)*positionOffset(2,:)-angularVelocityVector(2,i)*positionOffset(1,:))
          !$omp end workshare
          do j=1,3
             do k=j,3
                !$omp workshare
                energyTensorKinetic(j,k,i)=+sum(                            &
                     &                          +dble(selfBoundStatus(:,i)) &
                     &                          *     velocityOffset (j,:)  &
                     &                          *     velocityOffset (k,:)  &
                     &                         )
                !$omp end workshare
             end do
          end do
          !$omp do schedule(dynamic)
          do m=2,size(position,dim=2)
             ! Skip particles with zero weight.
             if (selfBoundStatus(m,i) == 0) cycle
             ! Sum over all prior particles. This reduces calculations by a factor of 2, and avoids the double-counting of
             ! potential energy so that we can simply drop the factor of ½ in the final expression for the potential energy
             ! tensor.
             !! Compute separations and related quantities.
             forall(j=1:3)
                separations(j,1:m-1)=+position(j,  m  ) &
                     &               -position(j,1:m-1)
             end forall
             separationsCubed               (1:m-1)=+sum (separations     (:,1:m-1  )**2,dim=1)**(3.0d0/2.0d0)
             separationsCubedInverseWeighted(1:m-1)=+dble(selfBoundStatus (  1:m-1,i)         )                &
                  &                                 /     separationsCubed(  1:m-1  )
             do j=1,3
                do k=j,3
                   energyTensorPotential(j,k,i)=+      energyTensorPotential          (j,k      ,i)  &
                        &                       +dble(selfBoundStatus                 (    m    ,i)) &
                        &                       *sum (                                               &
                        &                             +separations                    (j  ,1:m-1  )  &
                        &                             *separations                    (  k,1:m-1  )  &
                        &                             *separationsCubedInverseWeighted(    1:m-1  )  &
                        &                            )
                end do
             end do             
             ! Update progress.
             !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
             if (mpiSelf%isMaster()) then
#endif
                call displayCounter(                                                     &
                     &                          int(                                     &
                     &                              +100.0d0                             &
                     &                              *(                                   &
                     &                                +float(i                   -1   )  &
                     &                                *float(size(position,dim=2)-1   )  &
                     &                                +float(m                   -1   )  &
                     &                               )                                   &
                     &                              /  float(size(position,dim=2)-1   )  &
                     &                              /  float(self%bootstrapSampleCount)  &
                     &                             )                                   , &
                     &                          .false.                                  &
                     &                         )
#ifdef USEMPI
             end if
#endif
             !$ end if
          end do
          !$omp end do
          deallocate(separations                    )
          deallocate(separationsCubed               )
          deallocate(separationsCubedInverseWeighted)
          !$omp end parallel
       end do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
       ! Apply constant factors to the tensors.
       energyTensorPotential=-gravitationalConstant_internal    &
            &                *massParticle                  **2 &
            &                /sampleRate                    **2 &
            &                *  dble(size(position,dim=2))      &
            &                /(                                 &
            &                  +dble(size(position,dim=2))      &
            &                  -1.0d0                           &
            &                  -sampleRate                      &
            &                 )                                 &
            &                *energyTensorPotential
       energyTensorKinetic  =+0.5d0                             &
            &                *massParticle                      &
            &                /sampleRate                        &
            &                *energyTensorKinetic
       ! Symmetrize the tensors.
       do j=2,3
          do k=1,j-1
             energyTensorKinetic  (j,k,:)=energyTensorKinetic  (k,j,:)
             energyTensorPotential(j,k,:)=energyTensorPotential(k,j,:)
          end do
       end do
       ! Rotate tensors into the frame of the principal axes of the total energy tensor.
       do i=1,self%bootstrapSampleCount
          energyTensor=matrix          (                              &
               &                        +energyTensorKinetic  (:,:,i) &
               &                        +energyTensorPotential(:,:,i) &
               &                       )
          energyTensorKinetic_  =matrix(                              &
               &                        +energyTensorKinetic  (:,:,i) &
               &                       )
          energyTensorPotential_=matrix(                              &
               &                        +energyTensorPotential(:,:,i) &
               &                       )
          call energyTensor%eigenSystem(eigenVectors,eigenValues)
          energyTensorKineticPrincipal  (:,:,i)=eigenVectors%transpose()*energyTensorKinetic_  *eigenVectors
          energyTensorPotentialPrincipal(:,:,i)=eigenVectors%transpose()*energyTensorPotential_*eigenVectors
       end do
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(energyTensorKinetic           ,'energyTensorKinetic'           )
       call simulations(iSimulation)%analysis%writeDataset(energyTensorKineticPrincipal  ,'energyTensorKineticPrincipal'  )
       call simulations(iSimulation)%analysis%writeDataset(energyTensorPotential         ,'energyTensorPotential'         )
       call simulations(iSimulation)%analysis%writeDataset(energyTensorPotentialPrincipal,'energyTensorPotentialPrincipal')
       ! Deallocate workspace.
       if (self%selfBoundParticlesOnly) then
          nullify   (selfBoundStatus      )
       else
          deallocate(selfBoundStatus      )
       end if
       nullify      (positionMean         )
       nullify      (velocityMean         )
       nullify      (angularVelocityVector)
       deallocate   (energyTensorKinetic  )
       deallocate   (energyTensorPotential)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine energyTensorsOperate

