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

!% Contains a module which implements an N-body data operator which computes the rotation curve at a set of given radii.
  
  use, intrinsic :: ISO_C_Binding

  !# <nbodyOperator name="nbodyOperatorRotationCurve">
  !#  <description>An N-body data operator which computes the rotation curve at a set of given radii.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorRotationCurve
     !% An N-body data operator which computes the rotation curve at a set of given radii.
     private
     logical                                               :: selfBoundParticlesOnly
     integer         (c_size_t)                            :: bootstrapSampleCount 
     double precision          , allocatable, dimension(:) :: radius
   contains
     procedure :: operate => rotationCurveOperate
  end type nbodyOperatorRotationCurve

  interface nbodyOperatorRotationCurve
     !% Constructors for the ``rotationCurve'' N-body operator class.
     module procedure rotationCurveConstructorParameters
     module procedure rotationCurveConstructorInternal
  end interface nbodyOperatorRotationCurve

contains

  function rotationCurveConstructorParameters(parameters) result (self)
    !% Constructor for the ``rotationCurve'' N-body operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (nbodyOperatorRotationCurve)                              :: self
    type            (inputParameters           ), intent(inout)               :: parameters
    double precision                            , allocatable  , dimension(:) :: radius 
    logical                                                                   :: selfBoundParticlesOnly
    integer         (c_size_t                  )                              :: bootstrapSampleCount

    allocate(radius(parameters%count('radius')))
    !# <inputParameter>
    !#   <name>selfBoundParticlesOnly</name>
    !#   <source>parameters</source>
    !#   <description>If true, the mean position and velocity are computed only for self-bound particles.</description>
    !#   <type>logical</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radius</name>
    !#   <source>parameters</source>
    !#   <description>Radii at which the rotation curve should be computed.</description>
    !#   <type>float</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bootstrapSampleCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>30_c_size_t</defaultValue>
    !#   <description>The number of bootstrap resamples of the particles that should be used.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyOperatorRotationCurve(selfBoundParticlesOnly,bootstrapSampleCount,radius)
    !# <inputParametersValidate source="parameters"/>
    return
  end function rotationCurveConstructorParameters

  function rotationCurveConstructorInternal(selfBoundParticlesOnly,bootstrapSampleCount,radius) result (self)
    !% Internal constructor for the ``rotationCurve'' N-body operator class.
    implicit none
    type            (nbodyOperatorRotationCurve)                              :: self
    logical                                     , intent(in   )               :: selfBoundParticlesOnly
    integer         (c_size_t                  ), intent(in   )               :: bootstrapSampleCount
    double precision                            , intent(in   ), dimension(:) :: radius 
    !# <constructorAssign variables="selfBoundParticlesOnly, bootstrapSampleCount, radius"/>

    return
  end function rotationCurveConstructorInternal

  subroutine rotationCurveOperate(self,simulation)
    !% Determine the mean position and velocity of N-body particles.
    use Memory_Management
    use Galacticus_Error
    use Numerical_Constants_Physical
    use IO_HDF5
    use FGSL
    use Poisson_Random
    implicit none
    class           (nbodyOperatorRotationCurve), intent(inout)                 :: self
    type            (nBodyData                 ), intent(inout)                 :: simulation
    integer                                     , allocatable  , dimension(:,:) :: selfBoundStatus
    double precision                            , parameter                     :: sampleRate           =1.0d0
    double precision                            , allocatable  , dimension(:,:) :: positionMean                , rotationCurve
    double precision                            , allocatable  , dimension(:  ) :: distanceRadialSquared
    double precision                            , allocatable  , dimension(:,:) :: positionRelative
    integer                                                                     :: k
    type            (hdf5Object                )                                :: rotationCurveGroup
    integer         (c_size_t                  )                                :: i                           , j
    type            (fgsl_rng                  )                                :: pseudoSequenceObject
    logical                                                                     :: pseudoSequenceReset  =.true.

    ! Allocate workspace.
    call allocateArray(distanceRadialSquared,[  size(simulation%position,dim =2       )                          ])
    call allocateArray(positionRelative     ,[3,size(simulation%position,dim =2       )                          ])
    call allocateArray(rotationCurve        ,[  size(self      %radius  ,kind=c_size_t),self%bootstrapSampleCount])
    ! Determine the particle mask to use.
    if (self%selfBoundParticlesOnly) then
       if (simulation%analysis%hasDataset('selfBoundStatus')) then
          call simulation%analysis%readDataset('selfBoundStatus',selfBoundStatus)
          if (size(selfBoundStatus,dim=2) /= self%bootstrapSampleCount) call Galacticus_Error_Report('rotationCurveOperate','number of selfBoundStatus samples must equal number of requested bootstrap samples')
       else
          call Galacticus_Error_Report('rotationCurveOperate','self-bound status not available - apply a self-bound operator first')
       end if
    else
       call allocateArray(selfBoundStatus,[size(simulation%position,dim=2,kind=c_size_t),self%bootstrapSampleCount])
       do i=1,self%bootstrapSampleCount
          do j=1,size(simulation%position,dim=2)
             selfBoundStatus(j,i)=Poisson_Random_Get(pseudoSequenceObject,sampleRate,pseudoSequenceReset)
          end do
       end do
    end if
    ! Get mean position.
    if (.not.simulation%analysis%hasDataset('positionMean')) call Galacticus_Error_Report('rotationCurveOperate','mean position not available - apply the mean position operator first')
    call simulation%analysis%readDataset('positionMean',positionMean)
    if (size(positionMean,dim=2) /= self%bootstrapSampleCount) call Galacticus_Error_Report('rotationCurveOperate','number of positionMean samples must equal number of requested bootstrap samples')
    do i=1,self%bootstrapSampleCount
       !$omp parallel workshare
       ! Compute radial distance from the mean position.
       forall(k=1:3)
          positionRelative(k,:)=+simulation  %position(k,:) &
               &                -positionMean         (k,i)
       end forall
       distanceRadialSquared=sum(positionRelative**2,dim=1)
       ! Count particles within each radius.
       forall(k=1:size(self%radius))
          rotationCurve(k,i)=dble(sum(selfBoundStatus(:,i),mask=distanceRadialSquared <= self%radius(k)**2))
       end forall
       !$omp end parallel workshare
       ! Compute corresponding rotation curve.
       rotationCurve(:,i)=sqrt(gravitationalConstantGalacticus*simulation%massParticle*rotationCurve(:,i)/self%radius)
    end do
    ! Store results to file.
    rotationCurveGroup=simulation%analysis%openGroup('rotationCurve')
    call rotationCurveGroup%writeDataset(self%radius  ,'rotationCurveRadius'  )
    call rotationCurveGroup%writeDataset(rotationCurve,'rotationCurveVelocity')
    call rotationCurveGroup%close()
    ! Deallocate workspace.
    call deallocateArray(selfBoundStatus      )
    call deallocateArray(distanceRadialSquared)
    call deallocateArray(positionRelative     )
    call deallocateArray(rotationCurve        )
     return
  end subroutine rotationCurveOperate

