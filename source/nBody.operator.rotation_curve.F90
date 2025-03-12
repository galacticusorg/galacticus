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
Contains a module which implements an N-body data operator which computes the rotation curve at a set of given radii.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorRotationCurve">
   <description>An N-body data operator which computes the rotation curve at a set of given radii.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorRotationCurve
     !!{
     An N-body data operator which computes the rotation curve at a set of given radii.
     !!}
     private
     logical                                                                 :: selfBoundParticlesOnly
     integer         (c_size_t                  )                            :: bootstrapSampleCount
     double precision                            , allocatable, dimension(:) :: radius
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
   contains
     final     ::            rotationCurveDestructor
     procedure :: operate => rotationCurveOperate
  end type nbodyOperatorRotationCurve

  interface nbodyOperatorRotationCurve
     !!{
     Constructors for the ``rotationCurve'' N-body operator class.
     !!}
     module procedure rotationCurveConstructorParameters
     module procedure rotationCurveConstructorInternal
  end interface nbodyOperatorRotationCurve

contains

  function rotationCurveConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``rotationCurve'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorRotationCurve)                              :: self
    type            (inputParameters           ), intent(inout)               :: parameters
    class           (randomNumberGeneratorClass), pointer                     :: randomNumberGenerator_
    double precision                            , allocatable  , dimension(:) :: radius
    logical                                                                   :: selfBoundParticlesOnly
    integer         (c_size_t                  )                              :: bootstrapSampleCount

    allocate(radius(parameters%count('radius')))
    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the mean position and velocity are computed only for self-bound particles.</description>
    </inputParameter>
    <inputParameter>
      <name>radius</name>
      <source>parameters</source>
      <description>Radii at which the rotation curve should be computed.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorRotationCurve(selfBoundParticlesOnly,bootstrapSampleCount,radius,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function rotationCurveConstructorParameters

  function rotationCurveConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,radius,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the ``rotationCurve'' N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorRotationCurve)                              :: self
    logical                                     , intent(in   )               :: selfBoundParticlesOnly
    integer         (c_size_t                  ), intent(in   )               :: bootstrapSampleCount
    double precision                            , intent(in   ), dimension(:) :: radius
    class           (randomNumberGeneratorClass), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, radius, *randomNumberGenerator_"/>
    !!]

    return
  end function rotationCurveConstructorInternal

  subroutine rotationCurveDestructor(self)
    !!{
    Destructor for the ``meanPosition'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorRotationCurve), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine rotationCurveDestructor

  subroutine rotationCurveOperate(self,simulations)
    !!{
    Determine the mean position and velocity of N-body particles.
    !!}
    use :: Error                           , only : Error_Report
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nbodyOperatorRotationCurve), intent(inout)                 :: self
    type            (nBodyData                 ), intent(inout), dimension(:  ) :: simulations
    integer                                     , allocatable  , dimension(:,:) :: selfBoundStatus
    double precision                            , parameter                     :: sampleRate           =1.0d0
    double precision                            , pointer      , dimension(:,:) :: position
    double precision                            , allocatable  , dimension(:,:) :: positionMean                , rotationCurve
    double precision                            , allocatable  , dimension(:  ) :: distanceRadialSquared
    double precision                            , allocatable  , dimension(:,:) :: positionRelative
    integer                                                                     :: k
    type            (hdf5Object                )                                :: rotationCurveGroup
    integer         (c_size_t                  )                                :: i                           , j            , &
         &                                                                         iSimulation
    double precision                                                            :: massParticle
    
    do iSimulation=1,size(simulations)
       ! Allocate workspace.
       position     => simulations(iSimulation)%propertiesRealRank1%value('position'    )
       massparticle =  simulations(iSimulation)%attributesReal     %value('massParticle')
       allocate(distanceRadialSquared(  size(     position,dim =2)                          ))
       allocate(positionRelative     (3,size(     position,dim =2)                          ))
       allocate(rotationCurve        (  size(self%radius         ),self%bootstrapSampleCount))
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
          ! Count particles within each radius.
          forall(k=1:size(self%radius))
             rotationCurve(k,i)=dble(sum(selfBoundStatus(:,i),mask=distanceRadialSquared <= self%radius(k)**2))
          end forall
          !$omp end parallel workshare
          ! Compute corresponding rotation curve.
          rotationCurve(:,i)=sqrt(gravitationalConstant_internal*massParticle*rotationCurve(:,i)/self%radius)
       end do
       ! Store results to file.
       rotationCurveGroup=simulations(iSimulation)%analysis%openGroup('rotationCurve')
       call rotationCurveGroup%writeDataset(self%radius  ,'rotationCurveRadius'  )
       call rotationCurveGroup%writeDataset(rotationCurve,'rotationCurveVelocity')
       ! Deallocate workspace.
       deallocate(selfBoundStatus      )
       deallocate(distanceRadialSquared)
       deallocate(positionRelative     )
       deallocate(rotationCurve        )
    end do
    return
  end subroutine rotationCurveOperate

