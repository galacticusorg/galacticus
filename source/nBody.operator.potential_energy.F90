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
Implements an N-body data operator which determines the potential energy of each particle.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <nbodyOperator name="nbodyOperatorPotentialEnergy">
   <description>An N-body data operator which determines the potential energy of each particle.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorPotentialEnergy
     !!{
     An N-body data operator which determines the potential energy of each particle.
     !!}
     private
     logical                                               :: selfBoundParticlesOnly
     integer         (c_size_t                  )          :: bootstrapSampleCount
     double precision                                      :: bootstrapSampleRate             , thetaTolerance
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            potentialEnergyDestructor
     procedure :: operate => potentialEnergyOperate
  end type nbodyOperatorPotentialEnergy

  interface nbodyOperatorPotentialEnergy
     !!{
     Constructors for the ``potentialEnergy'' N-body operator class.
     !!}
     module procedure potentialEnergyConstructorParameters
     module procedure potentialEnergyConstructorInternal
  end interface nbodyOperatorPotentialEnergy

contains

  function potentialEnergyConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``potentialEnergy'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorPotentialEnergy)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    logical                                                       :: selfBoundParticlesOnly
    integer         (c_size_t                    )                :: bootstrapSampleCount
    double precision                                              :: bootstrapSampleRate   , thetaTolerance 

    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the gravitational potential is computed only from self-bound particles.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleRate</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The sampling rate for particles.</description>
    </inputParameter>
    <inputParameter>
      <name>thetaTolerance</name>
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <description>The criterion for the opening angle.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorPotentialEnergy(selfBoundParticlesOnly,bootstrapSampleCount,bootstrapSampleRate,thetaTolerance,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function potentialEnergyConstructorParameters

  function potentialEnergyConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,bootstrapSampleRate,thetaTolerance,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the ``potentialEnergy'' N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorPotentialEnergy)                        :: self
    logical                                       , intent(in   )         :: selfBoundParticlesOnly
    integer         (c_size_t                    ), intent(in   )         :: bootstrapSampleCount
    double precision                              , intent(in   )         :: bootstrapSampleRate   , thetaTolerance
    class           (randomNumberGeneratorClass  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, bootstrapSampleRate, thetaTolerance, *randomNumberGenerator_"/>
    !!]

    return
  end function potentialEnergyConstructorInternal

  subroutine potentialEnergyDestructor(self)
    !!{
    Destructor for the ``potentialEnergy'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorPotentialEnergy), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine potentialEnergyDestructor
  
  subroutine potentialEnergyOperate(self,simulations)
    !!{
    Determine the acceleration of bound particles from stripped ones.
    !!}
    use :: Error                           , only : Error_Report
    use :: Octree_Data_Structure           , only : octreeData
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nbodyOperatorPotentialEnergy), intent(inout)                 :: self
    type            (nBodyData                   ), intent(inout), dimension(:  ) :: simulations
    type            (octreeData                  )                                :: octreePosition
    integer         (c_size_t                    ), pointer      , dimension(:,:) :: selfBoundStatus
    double precision                              , pointer      , dimension(:,:) :: position
    double precision                              , allocatable  , dimension(:,:) :: positionRescaled
    double precision                              , allocatable  , dimension(:,:) :: potentialEnergy
    double precision                              , allocatable  , dimension(:  ) :: sampleRate 
    double precision                                                              :: lengthSoftening  , massParticle
    integer         (c_size_t                    )                                :: particleCount
    integer         (c_size_t                    )                                :: i                , j           , &
         &                                                                           iSimulation

    do iSimulation=1,size(simulations)
       ! Get simulation attributes.
       lengthSoftening=simulations(iSimulation)%attributesReal%value('lengthSoftening')
       massParticle   =simulations(iSimulation)%attributesReal%value('massParticle'   )
       ! Get particle data.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       ! Total number of particles.
       particleCount=size(position,dim=2,kind=c_size_t)
       ! Determine the particle mask to use.
       if (self%selfBoundParticlesOnly) then
          if (simulations(iSimulation)%propertiesIntegerRank1%exists('isBound')) then
             selfBoundStatus => simulations(iSimulation)%propertiesIntegerRank1%value('isBound')
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
       else
          allocate(selfBoundStatus(particleCount,self%bootstrapSampleCount))
          allocate(sampleRate     (1                                      ))
          sampleRate=self%bootstrapSampleRate
          do i=1,self%bootstrapSampleCount
             do j=1,particleCount
                selfBoundStatus(j,i)=self%randomNumberGenerator_%poissonSample(sampleRate(1))
             end do
          end do
       end if
       ! Compute the potential energy.
       allocate(potentialEnergy (           particleCount,self%bootstrapSampleCount))
       allocate(positionRescaled(3_c_size_t,particleCount                          ))
       positionRescaled=position/lengthSoftening
       potentialEnergy =0.0d0
       do i=1,self%bootstrapSampleCount
          if (sum(selfBoundStatus(:,i)) > 0) then
             ! Build octree.
             call octreePosition%build(positionRescaled,dble(selfBoundStatus(:,i)))
             !$omp parallel private(j)
             !$omp do schedule(dynamic)
             do j=1,particleCount
                if (selfBoundStatus(j,i) > 0 ) then
                   call octreePosition%traverseCompute(positionRescaled(:,j),dble(selfBoundStatus(j,i)),self%thetaTolerance,potentialEnergy(j,i),potentialEnergyPotential)
                   potentialEnergy(j,i)=+potentialEnergy(j,i)                                      &
                        &               *gravitationalConstant_internal                            &
                        &               *massParticle                                              &
                        &               /sampleRate(1)                                             &
                        &               /lengthSoftening                                           &
                        &               * dble(particleCount)                                      &
                        &               /(dble(particleCount-1_c_size_t)-self%bootstrapSampleRate)
                end if
             end do
             !$omp end do
             !$omp end parallel
             ! Destroy the octree.
             call octreePosition%destroy()
          end if
       end do
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(potentialEnergy,'energyPotential')
       ! Deallocate workspace.
       if (self%selfBoundParticlesOnly) then
          nullify   (selfBoundStatus )
       else
          deallocate(selfBoundStatus )
       end if
       deallocate   (positionRescaled)
       deallocate   (potentialEnergy )
       deallocate   (sampleRate      )
    end do
    return
  end subroutine potentialEnergyOperate

  subroutine potentialEnergyPotential(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquared)
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
  end subroutine potentialEnergyPotential

