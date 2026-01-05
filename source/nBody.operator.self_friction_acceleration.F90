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
Implements an N-body data operator which determines the acceleration of self-bound particles from unbound ones. The interaction between particles is computed using a tree method following \cite{barnes_hierarchical_1986}.
!!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorSelfFrictionAcceleration">
   <description>An N-body data operator which determines the acceleration of self-bound particles from unbound ones. The interaction between particles is computed using a tree method following \cite{barnes_hierarchical_1986}.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSelfFrictionAcceleration
     !!{
     An N-body data operator which determines the acceleration of self-bound particles from unbound ones. The interaction is computed using a tree method following \cite{barnes_hierarchical_1986}.
     !!}
     private
     integer         (c_size_t) :: bootstrapSampleCount
     double precision           :: thetaTolerance
   contains
     final     ::            selfFrictionAccelerationDestructor
     procedure :: operate => selfFrictionAccelerationOperate
  end type nbodyOperatorSelfFrictionAcceleration

  interface nbodyOperatorSelfFrictionAcceleration
     !!{
     Constructors for the \refClass{nbodyOperatorSelfFrictionAcceleration} N-body operator class.
     !!}
     module procedure selfFrictionAccelerationConstructorParameters
     module procedure selfFrictionAccelerationConstructorInternal
  end interface nbodyOperatorSelfFrictionAcceleration

contains

  function selfFrictionAccelerationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSelfFrictionAcceleration} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSelfFrictionAcceleration)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    integer         (c_size_t                             )                :: bootstrapSampleCount
    double precision                                                       :: thetaTolerance

    !![
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <inputParameter>
      <name>thetaTolerance</name>
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <description>The criterion for the opening angle.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorSelfFrictionAcceleration(bootstrapSampleCount,thetaTolerance)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function selfFrictionAccelerationConstructorParameters

  function selfFrictionAccelerationConstructorInternal(bootstrapSampleCount,thetaTolerance) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSelfFrictionAcceleration} N-body operator class
    !!}
    implicit none
    type            (nbodyOperatorSelfFrictionAcceleration)                :: self
    integer         (c_size_t                             ), intent(in   ) :: bootstrapSampleCount
    double precision                                       , intent(in   ) :: thetaTolerance
    !![
    <constructorAssign variables="bootstrapSampleCount,thetaTolerance"/>
    !!]

    return
  end function selfFrictionAccelerationConstructorInternal

  subroutine selfFrictionAccelerationDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorSelfFrictionAcceleration} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSelfFrictionAcceleration), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine selfFrictionAccelerationDestructor
 
  subroutine selfFrictionAccelerationOperate(self,simulations)
    !!{
    Determine the acceleration of bound particles from unbound ones.
    !!}
    use :: Display                         , only : displayIndent                 , displayUnindent
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : var_str
    use :: String_Handling                 , only : operator(//)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, gigaYear       , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Octree_Data_Structure           , only : octreeData
    implicit none
    class           (nbodyOperatorSelfFrictionAcceleration), intent(inout)                 :: self
    type            (nBodyData                            ), intent(inout), dimension(:  ) :: simulations
    type            (octreeData                           )                                :: octreePosition
    integer         (c_size_t                             ), pointer      , dimension(:,:) :: selfBoundStatus
    double precision                                       , pointer      , dimension(:,:) :: position                , sampleWeight
    double precision                                       , allocatable  , dimension(:,:) :: unBoundStatus
    double precision                                       , allocatable  , dimension(:,:) :: positionRescaled
    double precision                                       , allocatable  , dimension(:,:) :: accelerationSelfFriction
    double precision                                       , allocatable  , dimension(:  ) :: sampleRate
    double precision                                                      , dimension(3  ) :: acceleration            , accelerationAccumulate
    double precision                                                                       :: lengthSoftening         , massParticle 
    integer         (c_size_t                             )                                :: particleCount           , iSimulation           , &
         &                                                                                    i                       , j

    do iSimulation=1,size(simulations)
       ! Get simulation attributes.
       lengthSoftening=simulations(iSimulation)%attributesReal%value('lengthSoftening')
       massParticle   =simulations(iSimulation)%attributesReal%value('massParticle'   )
       ! Get particle data.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       ! Total number of particles.
       particleCount=size(position,dim=2,kind=c_size_t)
       if (simulations(iSimulation)%propertiesIntegerRank1%exists('isBound')) then
          selfBoundStatus => simulations(iSimulation)%propertiesIntegerRank1%value('isBound'     )
          sampleWeight    => simulations(iSimulation)%propertiesRealRank1   %value('sampleWeight')
          if (simulations(iSimulation)%analysis%hasDataset('bootstrapSampleRate')) then
             call simulations(iSimulation)%analysis%readDataset('bootstrapSampleRate',sampleRate)
          else
             allocate(sampleRate(1))
             sampleRate=1.0d0
          end if
          if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of selfBoundStatus samples must equal number of requested bootstrap samples'//{introspection:location})
       else
          call Error_Report('self-bound status not available - apply a self-bound operator first'//{introspection:location})
       end if
       allocate(unBoundStatus(particleCount,self%bootstrapSampleCount))
       do i=1,self%bootstrapSampleCount
          do j=1,particleCount
             if (selfBoundStatus(j,i) <= 0) then
                unBoundStatus(j,i)=sampleWeight(j,i)
             else
                unBoundStatus(j,i)=0.0d0
             end if
          end do
       end do
       allocate(accelerationSelfFriction(3_c_size_t,             self%bootstrapSampleCount))
       allocate(positionRescaled        (3_c_size_t,particleCount                         ))
       positionRescaled=position/lengthSoftening
       ! Iterate over bootstrap samples.
       call displayIndent('Performing self-friction analysis on bootstrap samples')
       do i=1,self%bootstrapSampleCount
          call displayIndent(var_str('sample ')//i//' of '//self%bootstrapSampleCount)
          if (sum(unBoundStatus(:,i)) > 0.0d0) then
             ! Build octree.
             call octreePosition%build(positionRescaled,unBoundStatus(:,i))
             accelerationAccumulate=0.0d0
             !$omp parallel private(j,acceleration)
             !$omp do reduction(+:accelerationAccumulate)
             do j=1,particleCount
                if (selfBoundStatus(j,i) > 0) then
                   call octreePosition%traverseCompute(positionRescaled(:,j),sampleWeight(j,i),self%thetaTolerance,acceleration,gravitationalAcceleration,isExternalParticle=.true.)
                   accelerationAccumulate=accelerationAccumulate+acceleration*sampleWeight(j,i)
                end if
             end do
             !$omp end do
             !$omp end parallel
             accelerationSelfFriction(:,i)=+accelerationAccumulate          &
                  &                        /sum(dble(selfBoundStatus(:,i))) &
                  &                        *gravitationalConstant_internal  &
                  &                        *massParticle                    &
                  &                        /sampleRate(1)                   &
                  &                        /lengthSoftening**2              &
                  &                        *kilo                            &
                  &                        *gigaYear                        &
                  &                        /megaParsec
             ! Destroy the octree.
             call octreePosition%destroy()
          else
             accelerationSelfFriction(:,i)=0.0d0
          end if
          call displayUnindent('done')
       end do
       call displayUnindent('done')
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(accelerationSelfFriction,'accelerationSelfFriction')
       ! Deallocate workspace.
       nullify   (selfBoundStatus         )
       nullify   (sampleWeight            )
       deallocate(unBoundStatus           )
       deallocate(positionRescaled        )
       deallocate(accelerationSelfFriction)
       deallocate(sampleRate              )
    end do
    return
  end subroutine selfFrictionAccelerationOperate

  subroutine gravitationalAcceleration(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquared)
    !!{
    Compute the interaction between a particle and a node in the octree. Currently assumes the functional form of the softening used by
    Gadget.
    !!}:
    implicit none
    double precision, dimension(:), intent(inout) :: value
    double precision, dimension(3), intent(in   ) :: centerOfMass     , relativePosition
    double precision,               intent(in   ) :: nodeWeight       , separation      , &
         &                                           separationSquared
    double precision, dimension(3)                :: acceleration
    !$GLC attributes unused :: centerOfMass, separationSquared

    if      (separation == 0.0d0) then
       acceleration=0.0d0
    else if (separation <= 0.5d0) then
       acceleration=-(                   &
            &          32.0d0/3.0d0      &
            &         +32.0d0            &
            &         *separationSquared &
            &         *(                 &
            &           -6.0d0/5.0d0     &
            &           +separation      &
            &          )                 &
            &        )                   &
            &       *relativePosition
    else if (separation <= 1.0d0) then
       acceleration=-(                             &
            &         -  1.0d0/15.0d0              &
            &         + 64.0d0*separation**3/3.0d0 &
            &         - 48.0d0*separation**4       &
            &         +192.0d0*separation**5/5.0d0 &
            &         - 32.0d0*separation**6/3.0d0 &
            &        )                             &
            &       *relativePosition              &
            &       /separation**3
    else
       acceleration=-relativePosition  &
            &       /separation**3
    end if
    value=value+nodeWeight*acceleration
    return
  end subroutine gravitationalAcceleration
