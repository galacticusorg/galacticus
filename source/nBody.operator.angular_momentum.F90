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
Implements an N-body data operator which determines the mean angular momentum of particles.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <nbodyOperator name="nbodyOperatorAngularMomentum">
   <description>An N-body data operator which determines the mean angular momentum of particles. Also finds the angular velocity vector for the rotating frame in which the angular momentum would be zero.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorAngularMomentum
     !!{
     An N-body data operator which determines the mean angular momentum of particles.
     !!}
     private
     logical                                      :: selfBoundParticlesOnly
     integer(c_size_t                  )          :: bootstrapSampleCount
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            angularMomentumDestructor
     procedure :: operate => angularMomentumOperate
  end type nbodyOperatorAngularMomentum

  interface nbodyOperatorAngularMomentum
     !!{
     Constructors for the ``angularMomentum'' N-body operator class.
     !!}
     module procedure angularMomentumConstructorParameters
     module procedure angularMomentumConstructorInternal
  end interface nbodyOperatorAngularMomentum

contains

  function angularMomentumConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``angularMomentum'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyOperatorAngularMomentum)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    logical                                              :: selfBoundParticlesOnly
    integer(c_size_t                    )                :: bootstrapSampleCount

    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the mean angular momentum is computed only for self-bound particles</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorAngularMomentum(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function angularMomentumConstructorParameters

  function angularMomentumConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the ``angularMomentum'' N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorAngularMomentum)                        :: self
    logical                              , intent(in   )         :: selfBoundParticlesOnly
    integer(c_size_t                    ), intent(in   )         :: bootstrapSampleCount
    class  (randomNumberGeneratorClass  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, *randomNumberGenerator_"/>
    !!]

    return
  end function angularMomentumConstructorInternal

  subroutine angularMomentumDestructor(self)
    !!{
    Destructor for the ``angularMomentum'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorAngularMomentum), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine angularMomentumDestructor
  
  subroutine angularMomentumOperate(self,simulations)
    !!{
    Determine the mean position and velocity of N-body particles.
    !!}
    use :: Error            , only : Error_Report
    use :: Linear_Algebra   , only : vector       , matrix         , assignment(=), operator(*)
    implicit none
    class          (nbodyOperatorAngularMomentum), intent(inout)                   :: self
    type           (nBodyData                   ), intent(inout), dimension(:    ) :: simulations
    integer         (c_size_t                   ), pointer      , dimension(:,:  ) :: selfBoundStatus
    double precision                             , parameter                       :: sampleRate          =1.0d0
    double precision                             , pointer      , dimension(:,:  ) :: positionMean              , velocityMean         , &
         &                                                                            position                  , velocity             , &
         &                                                                            angularMomentumMean       , angularVelocityVector
    double precision                             , allocatable  , dimension(:,:  ) :: positionOffset            , velocityOffset
    double precision                             , allocatable  , dimension(:,:,:) :: inertiaTensor
    type            (vector                     )                                  :: angularMomentum_
    type            (matrix                     )                                  :: inertiaTensor_
    double precision                                                               :: weight
    integer         (c_size_t                   )                                  :: i                         , j                    , &
         &                                                                            iSimulation
    
    do iSimulation=1,size(simulations)
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
       ! Get the mean position/velocity.
       if (simulations(iSimulation)%propertiesRealRank1%exists('positionMean')) then
          positionMean => simulations(iSimulation)%propertiesRealRank1%value('positionMean')
          if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of selfBoundStatus samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('mean position not available - apply a meanPosition operator first'//{introspection:location})
       end if
       if (simulations(iSimulation)%propertiesRealRank1%exists('velocityMean')) then
          velocityMean => simulations(iSimulation)%propertiesRealRank1%value('velocityMean')
          if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of selfBoundStatus samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('mean velocity not available - apply a meanVelocity operator first'//{introspection:location})
       end if
       ! Compute mean angular momentum.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       velocity => simulations(iSimulation)%propertiesRealRank1%value('velocity')
       allocate(angularMomentumMean  (3_c_size_t           ,self%bootstrapSampleCount))
       allocate(angularVelocityVector(3_c_size_t           ,self%bootstrapSampleCount))
       allocate(inertiaTensor        (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(positionOffset       (3_c_size_t           ,size(position,dim=2)     ))
       allocate(velocityOffset       (3_c_size_t           ,size(position,dim=2)     ))
       do i=1,self%bootstrapSampleCount
          !$omp parallel workshare
          forall(j=1:3)
             positionOffset(j,:)=position(j,:)-positionMean(j,i)
             velocityOffset(j,:)=velocity(j,:)-velocityMean(j,i)
          end forall
          weight                    =dble(sum(selfBoundStatus(:,i)))
          angularMomentumMean(1  ,i)=+sum(                                           &
               &                          +(                                         &
               &                            +velocityOffset(2,:)*positionOffset(3,:) &
               &                            -velocityOffset(3,:)*positionOffset(2,:) &
               &                           )                                         &
               &                          *dble(selfBoundStatus(:,i))                &
               &                         )                                           &
               &                     /weight
          angularMomentumMean(2  ,i)=+sum(                                           &
               &                          +(                                         &
               &                            -velocityOffset(1,:)*positionOffset(3,:) &
               &                            +velocityOffset(3,:)*positionOffset(1,:) &
               &                           )                                         &
               &                          *dble(selfBoundStatus(:,i))                &
               &                         )                                           &
               &                     /weight
          angularMomentumMean(3  ,i)=+sum(                                           &
               &                          +(                                         &
               &                            +velocityOffset(1,:)*positionOffset(2,:) &
               &                            -velocityOffset(2,:)*positionOffset(1,:) &
               &                           )                                         &
               &                          *dble(selfBoundStatus(:,i))                &
               &                         )                                           &
               &                     /weight
          inertiaTensor      (1,1,i)=+sum(                                                  &
               &                          +(+positionOffset(2,:)**2+positionOffset(3,:)**2) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          inertiaTensor      (2,2,i)=+sum(                                                  &
               &                          +(+positionOffset(3,:)**2+positionOffset(1,:)**2) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          inertiaTensor      (3,3,i)=+sum(                                                  &
               &                          +(+positionOffset(1,:)**2+positionOffset(2,:)**2) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          inertiaTensor      (1,2,i)=-sum(                                                  &
               &                          +(+positionOffset(1,:)   *positionOffset(2,:)   ) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          inertiaTensor      (1,3,i)=-sum(                                                  &
               &                          +(+positionOffset(1,:)   *positionOffset(3,:)   ) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          inertiaTensor      (2,3,i)=-sum(                                                  &
               &                          +(+positionOffset(2,:)   *positionOffset(3,:)   ) &
               &                          *dble(selfBoundStatus(:,i))                       &
               &                         )                                                  &
               &                     /weight
          !$omp end parallel workshare
          inertiaTensor      (2,1,i)=inertiaTensor(1,2,i)
          inertiaTensor      (3,1,i)=inertiaTensor(1,3,i)
          inertiaTensor      (3,2,i)=inertiaTensor(2,3,i)
          ! Find the angular velocity vector corresponding to the zero angular momentum frame. To find this we make use of the
          ! definition:
          !
          !  J = ℐ ω
          !
          ! where J is the angular momentum, ℐ is the inertia tensor, and ω the angular velocity.
          inertiaTensor_            = matrix                (inertiaTensor      (:,:,i))
          angularMomentum_          = vector                (angularMomentumMean(:  ,i))
          angularVelocityVector(:,i)= inertiaTensor_%inverse()                           &
               &                     *angularMomentum_          
       end do
       ! Store the results.
       call simulations(iSimulation)%propertiesRealRank1%set         ('angularMomentumMean'  ,angularMomentumMean    )
       call simulations(iSimulation)%propertiesRealRank1%set         ('angularVelocityVector',angularVelocityVector  )
       ! Store results to file.
       call simulations(iSimulation)%analysis           %writeDataset(angularMomentumMean    ,'angularMomentumMean'  )
       call simulations(iSimulation)%analysis           %writeDataset(angularVelocityVector  ,'angularVelocityVector')
       ! Deallocate workspace.
       if (self%selfBoundParticlesOnly) then
          nullify   (      selfBoundStatus)
       else
          deallocate(      selfBoundStatus)
       end if
       nullify      (       positionMean  )
       nullify      (       velocityMean  )
       deallocate   (       positionOffset)
       deallocate   (       velocityOffset)
       nullify      (angularMomentumMean  )
       nullify      (angularVelocityVector)
       deallocate   (        inertiaTensor)
    end do
    return
  end subroutine angularMomentumOperate

