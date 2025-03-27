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
Implements an N-body data operator which determines the mean position and velocity of particles.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <nbodyOperator name="nbodyOperatorMeanPosition">
   <description>An N-body data operator which determines the mean position and velocity of particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorMeanPosition
     !!{
     An N-body data operator which determines the mean position and velocity of particles.
     !!}
     private
     logical                                      :: selfBoundParticlesOnly
     integer(c_size_t                  )          :: bootstrapSampleCount
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            meanPositionDestructor
     procedure :: operate => meanPositionOperate
  end type nbodyOperatorMeanPosition

  interface nbodyOperatorMeanPosition
     !!{
     Constructors for the {\normalfont \ttfamily meanPosition} N-body operator class.
     !!}
     module procedure meanPositionConstructorParameters
     module procedure meanPositionConstructorInternal
  end interface nbodyOperatorMeanPosition

contains

  function meanPositionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily meanPosition} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyOperatorMeanPosition )                :: self
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
    self=nbodyOperatorMeanPosition(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function meanPositionConstructorParameters

  function meanPositionConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily meanPosition} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorMeanPosition)                         :: self
    logical                            , intent(in   )         :: selfBoundParticlesOnly
    integer(c_size_t                  ), intent(in   )         :: bootstrapSampleCount
    class  (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, *randomNumberGenerator_"/>
    !!]

    return
  end function meanPositionConstructorInternal

  subroutine meanPositionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily meanPosition} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorMeanPosition), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine meanPositionDestructor
  
  subroutine meanPositionOperate(self,simulations)
    !!{
    Determine the mean position and velocity of N-body particles.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class          (nbodyOperatorMeanPosition), intent(inout)                 :: self
    type           (nBodyData                ), intent(inout), dimension(:  ) :: simulations
    integer         (c_size_t                ), pointer      , dimension(:,:) :: selfBoundStatus
    double precision                          , parameter                     :: sampleRate     =1.0d0
    double precision                          , pointer      , dimension(:,:) :: positionMean          , velocityMean
    double precision                          , pointer      , dimension(:,:) :: position              , velocity
    double precision                                                          :: weight
    integer         (c_size_t                )                                :: i                     , j           , &
         &                                                                       iSimulation
    
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
       ! Compute mean position and velocity.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       velocity => simulations(iSimulation)%propertiesRealRank1%value('velocity')
       allocate(positionMean(3_c_size_t,self%bootstrapSampleCount))
       allocate(velocityMean(3_c_size_t,self%bootstrapSampleCount))
       do i=1,self%bootstrapSampleCount
          !$omp parallel workshare
          weight=dble(sum(selfBoundStatus(:,i)))
          forall(j=1:3)
             positionMean(j,i)=+sum(position(j,:)*dble(selfBoundStatus(:,i))) &
                  &            /weight
             velocityMean(j,i)=+sum(velocity(j,:)*dble(selfBoundStatus(:,i))) &
                  &            /weight
          end forall
          !$omp end parallel workshare
       end do
       ! Store the results.
       call simulations(iSimulation)%propertiesRealRank1%set('positionMean',positionMean)
       call simulations(iSimulation)%propertiesRealRank1%set('velocityMean',velocityMean)
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(positionMean,'positionMean')
       call simulations(iSimulation)%analysis%writeDataset(velocityMean,'velocityMean')
       ! Deallocate workspace.
       if (self%selfBoundParticlesOnly) then
          nullify   (selfBoundStatus)
       else
          deallocate(selfBoundStatus)
       end if
       nullify      (positionMean   )
       nullify      (velocityMean   )
    end do
    return
  end subroutine meanPositionOperate

