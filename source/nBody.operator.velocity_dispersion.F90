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
Implements an N-body data operator which computes the velocity dispersion in a set of given spherical shells.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorVelocityDispersion">
   <description>An N-body data operator which computes the rotation curve at a set of given radii.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorVelocityDispersion
     !!{
     An N-body data operator which computes the rotation curve at a set of given radii.
     !!}
     private
     logical                                                                 :: selfBoundParticlesOnly
     integer         (c_size_t                  )                            :: bootstrapSampleCount
     double precision                            , allocatable, dimension(:) :: radiusInner                     , radiusOuter
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
   contains
     final     ::            velocityDispersionDestructor
     procedure :: operate => velocityDispersionOperate
  end type nbodyOperatorVelocityDispersion

  interface nbodyOperatorVelocityDispersion
     !!{
     Constructors for the {\normalfont \ttfamily velocityDispersion} N-body operator class.
     !!}
     module procedure velocityDispersionConstructorParameters
     module procedure velocityDispersionConstructorInternal
  end interface nbodyOperatorVelocityDispersion

contains

  function velocityDispersionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily velocityDispersion} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorVelocityDispersion)                              :: self
    type            (inputParameters                ), intent(inout)               :: parameters
    class           (randomNumberGeneratorClass     ), pointer                     :: randomNumberGenerator_
    double precision                                 , allocatable  , dimension(:) :: radiusInner           , radiusOuter
    logical                                                                        :: selfBoundParticlesOnly
    integer         (c_size_t                       )                              :: bootstrapSampleCount

    allocate(radiusInner(parameters%count('radiusInner')))
    allocate(radiusOuter(parameters%count('radiusOuter')))
    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the velocity dispersion is computed only for self-bound particles.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusInner</name>
      <source>parameters</source>
      <description>Inner radii of spherical shells within which the velocity dispersion should be computed.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusOuter</name>
      <source>parameters</source>
      <description>Outer radii of spherical shells within which the velocity dispersion should be computed.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorVelocityDispersion(selfBoundParticlesOnly,bootstrapSampleCount,radiusInner,radiusOuter,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function velocityDispersionConstructorParameters

  function velocityDispersionConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,radiusInner,radiusOuter,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily velocityDispersion} N-body operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyOperatorVelocityDispersion)                              :: self
    logical                                          , intent(in   )               :: selfBoundParticlesOnly
    integer         (c_size_t                       ), intent(in   )               :: bootstrapSampleCount
    double precision                                 , intent(in   ), dimension(:) :: radiusInner            , radiusOuter
    class           (randomNumberGeneratorClass     ), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, radiusInner, radiusOuter, *randomNumberGenerator_"/>
    !!]

    if (size(self%radiusInner) /= size(self%radiusOuter)) call Error_Report('number of inner and outer radii should be equal'//{introspection:location})
    return
  end function velocityDispersionConstructorInternal

  subroutine velocityDispersionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily velocityDispersion} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorVelocityDispersion), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine velocityDispersionDestructor

  subroutine velocityDispersionOperate(self,simulations)
    !!{
    Determine the mean position and velocity of N-body particles.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (nbodyOperatorVelocityDispersion), intent(inout)                 :: self
    type            (nBodyData                      ), intent(inout), dimension(:  ) :: simulations
    integer                                          , allocatable  , dimension(:,:) :: selfBoundStatus
    double precision                                 , parameter                     :: sampleRate           =1.0d0
    double precision                                 , allocatable  , dimension(:,:) :: positionMean                , velocityDispersion
    double precision                                 , pointer      , dimension(:,:) :: position                    , velocity
    double precision                                 , allocatable  , dimension(:  ) :: distanceRadialSquared
    double precision                                 , allocatable  , dimension(:,:) :: positionRelative
    logical                                          , allocatable  , dimension(:  ) :: mask
    double precision                                                , dimension(3  ) :: velocityMean                , velocityMeanSquared
    integer                                                                          :: k
    integer         (c_size_t                       )                                :: i                           , j                  , &
         &                                                                              iSimulation

    do iSimulation=1,size(simulations)
       ! Get particle data.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       velocity => simulations(iSimulation)%propertiesRealRank1%value('velocity')
       ! Allocate workspace.
       allocate(distanceRadialSquared(  size(     position   ,dim =2)                          ))
       allocate(positionRelative     (3,size(     position   ,dim =2)                          ))
       allocate(velocityDispersion   (  size(self%radiusInner       ),self%bootstrapSampleCount))
       allocate(mask                 (  size(     position   ,dim =2)                          ))
       ! Determine the particle mask to use.
       if (self%selfBoundParticlesOnly) then
          if (simulations(iSimulation)%analysis%hasDataset('selfBoundStatus')) then
             call simulations(iSimulation)%analysis%readDataset('selfBoundStatus',selfBoundStatus)
             if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of selfBoundStatus samples must equal number of requested bootstrap samples'//{introspection:location})
          else
             call Error_Report('self-bound status not available - apply a self-bound operator first'//{introspection:location})
          end if
       else
          allocate(selfBoundStatus(size(position,dim=2),self%bootstrapSampleCount))
          do i=1,self%bootstrapSampleCount
             do j=1,size(position,dim=2)
                selfBoundStatus(j,i)=self%randomNumberGenerator_%poissonSample(sampleRate)
             end do
          end do
       end if
       ! Get mean position.
       if (.not.simulations(iSimulation)%analysis%hasDataset('positionMean')) call Error_Report('mean position not available - apply the mean position operator first'//{introspection:location})
       call simulations(iSimulation)%analysis%readDataset('positionMean',positionMean)
       if (size(positionMean,dim=2) /= self%bootstrapSampleCount) call Error_Report('number of positionMean samples must equal number of requested bootstrap samples'//{introspection:location})
       do i=1,self%bootstrapSampleCount
          !$omp parallel workshare
          ! Compute radial distance from the mean position.
          forall(k=1:3)
             positionRelative(k,:)=+position    (k,:) &
                  &                -positionMean(k,i)
          end forall
          distanceRadialSquared=sum(positionRelative**2,dim=1)
          !$omp end parallel workshare
          ! Compute velocity dispersion within each shell.
          do k=1,size(self%radiusInner)
             !$omp parallel workshare
             mask   = distanceRadialSquared >= self%radiusInner(k)**2 &
                  &  .and.                                            &
                  &   distanceRadialSquared <  self%radiusOuter(k)**2
             forall(j=1:3)
                velocityMean       (j)=sum(velocity(j,:)   *dble(selfBoundStatus(:,i)),mask=mask)/dble(sum(selfBoundStatus(:,i),mask))
                velocityMeanSquared(j)=sum(velocity(j,:)**2*dble(selfBoundStatus(:,i)),mask=mask)/dble(sum(selfBoundStatus(:,i),mask))
             end forall
             velocityDispersion(k,i)=sqrt(sum(velocityMeanSquared-velocityMean**2)/3.0d0)
             !$omp end parallel workshare
          end do
       end do
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(self%radiusInner  ,'velocityDispersionRadiusInner')
       call simulations(iSimulation)%analysis%writeDataset(self%radiusInner  ,'velocityDispersionRadiusOuter')
       call simulations(iSimulation)%analysis%writeDataset(velocityDispersion,'velocityDispersion'           )
       ! Deallocate workspace.
       deallocate(selfBoundStatus      )
       deallocate(distanceRadialSquared)
       deallocate(positionRelative     )
       deallocate(velocityDispersion   )
       deallocate(mask                 )
    end do
    return
  end subroutine velocityDispersionOperate

