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
Implements an N-body data operator which computes the inertia tensor eigenvalues and eigenvectors, along with axis ratios.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <nbodyOperator name="nbodyOperatorInertiaTensor">
   <description>An N-body data operator which computes the inertia tensor eigenvalues and eigenvectors, along with axis ratios.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorInertiaTensor
     !!{
     An N-body data operator which computes the inertia tensor eigenvalues and eigenvectors, along with axis ratios.
     !!}
     private
     logical                                               :: selfBoundParticlesOnly
     integer         (c_size_t                  )          :: bootstrapSampleCount
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     double precision                                      :: radiusMaximum
   contains
     final     ::            inertiaTensorDestructor
     procedure :: operate => inertiaTensorOperate
  end type nbodyOperatorInertiaTensor

  interface nbodyOperatorInertiaTensor
     !!{
     Constructors for the \refClass{nbodyOperatorInertiaTensor} N-body operator class.
     !!}
     module procedure inertiaTensorConstructorParameters
     module procedure inertiaTensorConstructorInternal
  end interface nbodyOperatorInertiaTensor

contains

  function inertiaTensorConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorInertiaTensor} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter  , inputParameters
    implicit none
    type            (nbodyOperatorInertiaTensor)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    logical                                                     :: selfBoundParticlesOnly
    integer         (c_size_t                  )                :: bootstrapSampleCount
    double precision                                            :: radiusMaximum
    
    !![
    <inputParameter>
      <name>selfBoundParticlesOnly</name>
      <source>parameters</source>
      <description>If true, the maximum velocity is computed only for self-bound particles</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusMaximum</name>
      <source>parameters</source>
      <description>The maximum radius (relative to the position of the most-bound particle) to include when computing the inertia tensor.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorInertiaTensor(radiusMaximum,selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function inertiaTensorConstructorParameters

  function inertiaTensorConstructorInternal(radiusMaximum,selfBoundParticlesOnly,bootstrapSampleCount,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorInertiaTensor} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorInertiaTensor)                        :: self
    double precision                            , intent(in   )         :: radiusMaximum
    logical                                     , intent(in   )         :: selfBoundParticlesOnly
    integer         (c_size_t                  ), intent(in   )         :: bootstrapSampleCount
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="radiusMaximum, selfBoundParticlesOnly, bootstrapSampleCount, *randomNumberGenerator_"/>
    !!]

    return
  end function inertiaTensorConstructorInternal

  subroutine inertiaTensorDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorInertiaTensor} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorInertiaTensor), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine inertiaTensorDestructor
  
  subroutine inertiaTensorOperate(self,simulations)
    !!{
    Determine the inertia tensor of the particle distribution, along with its eigenvectors, eigenvalues, and axis ratios.
    !!}
    use    :: Display        , only : displayCounter    , displayCounterClear   , displayIndent, displayMessage, &
         &                            displayUnindent   , verbosityLevelStandard
    use    :: Error          , only : Error_Report
    use    :: Linear_Algebra , only : matrix            , vector                , operator(*)  , assignment(=)
    use    :: Sorting        , only : sort
    use    :: Array_Utilities, only : Array_Reverse
#ifdef USEMPI
    use    :: MPI_Utilities , only : mpiSelf
#endif
    !$ use :: OMP_Lib       , only : OMP_Get_Thread_Num
    implicit none
    class          (nbodyOperatorInertiaTensor), intent(inout)                   :: self
    type           (nBodyData                 ), intent(inout), dimension(:  )   :: simulations
    double precision                           , parameter                       :: sampleRate      =1.0d0
    double precision                           , allocatable  , dimension(:,:,:) :: inertiaTensor         , eigenVectors
    double precision                           , pointer      , dimension(:,:  ) :: position
    integer         (c_size_t                 ), pointer      , dimension(:,:  ) :: selfBoundStatus
    integer         (c_size_t                 ), pointer      , dimension(:,:  ) :: indexMostBound
    double precision                           , allocatable  , dimension(:,:  ) :: positionOffset        , axisRatios   , &
         &                                                                          eigenValues
    double precision                           , allocatable  , dimension(:    ) :: radiusOffset
    logical                                    , allocatable  , dimension(:    ) :: mask
    integer         (c_size_t                 )                                  :: i                     , j            , &
         &                                                                          k                     , iSimulation
    double precision                                                             :: massParticle
    type            (matrix                   )                                  :: inertiaTensor_        , eigenVectors_
    type            (vector                   )                                  :: eigenValues_
    double precision                                          , dimension(3    ) :: axisLengths

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
    call displayIndent('compute inertia tensor',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Iterate over simulations.
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
       massParticle=simulations(iSimulation)%attributesReal%value('massParticle')
       ! Get index of the most-bound particle.
       if (simulations(iSimulation)%propertiesIntegerRank1%exists('indexMostBound')) then
          indexMostBound => simulations(iSimulation)%propertiesIntegerRank1%value('indexMostBound')
       else
          call Error_Report('index of most bound particle not available - apply a self-bound operator first'//{introspection:location})
       end if
       ! Get position.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       ! Compute the inertia tensor.
       allocate(positionOffset(3_c_size_t           ,size(position,dim=2)     ))
       allocate(radiusOffset  (                      size(position,dim=2)     ))
       allocate(mask          (                      size(position,dim=2)     ))
       allocate(inertiaTensor (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(eigenVectors  (3_c_size_t,3_c_size_t,self%bootstrapSampleCount))
       allocate(eigenValues   (3_c_size_t           ,self%bootstrapSampleCount))
       allocate(axisRatios    (3_c_size_t           ,self%bootstrapSampleCount))
       inertiaTensor=0.0d0
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       do i=1,self%bootstrapSampleCount
          !$omp parallel private(j,k) reduction(+:inertiaTensor)
          ! First find the position (and radius) relative to the most-bound particle.
          !$omp workshare
          forall(j=1:3)
             positionOffset(j,:)=+position(j,               :   ) &
                  &              -position(j,indexMostBound(1,i))
          end forall
          radiusOffset=sqrt(sum(positionOffset**2,dim=1))
          ! Construct a mask of particles to include in the calculation.
          mask=radiusOffset < self%radiusMaximum
          ! Compute the inertia tensor, defined as (https://scienceworld.wolfram.com/physics/MomentofInertia.html):
          !
          !     ⌠                   ⎡ +y²+z² -xy    -xz    ⎤
          ! I = ⎮ dx dy dz ρ(x,y,z) ⎢ -xy    +z²+x² -yz    ⎢
          !     ⌡                   ⎣ -xz    -yz    +x²+y² ⎦          
          !
          ! which we estimate as a discrete sum over particles.
          inertiaTensor(1,1,i)=+sum(                                   &
               &                         +dble(selfBoundStatus(:,i))   &
               &                         *(                            &
               &                           +positionOffset (2,:)**2    &
               &                           +positionOffset (3,:)**2    &
               &                          )                          , &
               &                    mask= mask                         &
               &                   )
          inertiaTensor(2,2,i)=+sum(                                   &
               &                         +dble(selfBoundStatus(:,i))   &
               &                         *(                            &
               &                           +positionOffset (1,:)**2    &
               &                           +positionOffset (3,:)**2    &
               &                          )                          , &
               &                    mask= mask                         &
               &                   )
          inertiaTensor(3,3,i)=+sum(                                   &
               &                         +dble(selfBoundStatus(:,i))   &
               &                         *(                            &
               &                           +positionOffset (1,:)**2    &
               &                           +positionOffset (2,:)**2    &
               &                          )                          , &
               &                    mask= mask                         &
               &                   )
          !$omp end workshare
          do j=1,2
             do k=j+1,3
                !$omp workshare
                inertiaTensor(j,k,i)=+sum(                                  &
                     &                         -dble(selfBoundStatus(:,i))  &
                     &                         *     positionOffset (j,:)   &
                     &                         *     positionOffset (k,:) , &
                     &                    mask= mask                        &
                     &                   )
                !$omp end workshare
             end do
          end do
          !$omp end parallel
       end do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
       ! Symmetrize the tensor and scale to particle mass.
       do j=2,3
          do k=1,j-1
             inertiaTensor(j,k,:)=inertiaTensor(k,j,:)
          end do
       end do
       inertiaTensor=+inertiaTensor &
            &        *massParticle
       ! Find the eigenvectors, eigenvalues, and axis ratios.
       do i=1,self%bootstrapSampleCount
          inertiaTensor_=matrix(inertiaTensor(:,:,i))
          call inertiaTensor_%eigenSystem(eigenVectors_,eigenValues_)
          eigenVectors(:,:,i)=eigenVectors_
          eigenValues (  :,i)=eigenValues_
          ! For an inhomogeneous ellipsoid with axis lengths (a,b,c), the eigenvalues of the inertia tensor, (A,B,C), are
          ! (https://scienceworld.wolfram.com/physics/MomentofInertiaEllipsoid.html):          
          !
          !  A ∝ (b²+c²),
          !  B ∝ (a²+c²),
          !  C ∝ (a²+b²).
          !
          ! Therefore, the axis lengths are related to these eigenvalues via:
          !
          !  a² ∝ (B+C-A)/2,
          !  b² ∝ (A+C-B)/2,
          !  c² ∝ (A+B-C)/2,
          !
          ! from which the axis ratios are trivially found.
          axisLengths(1  )=sqrt(0.5d0*(-eigenValues(1,i)+eigenValues(2,i)+eigenValues(3,i)))
          axisLengths(2  )=sqrt(0.5d0*(+eigenValues(1,i)-eigenValues(2,i)+eigenValues(3,i)))
          axisLengths(3  )=sqrt(0.5d0*(+eigenValues(1,i)+eigenValues(2,i)-eigenValues(3,i)))
          call sort(axisLengths)
          axisLengths     =Array_Reverse(axisLengths)
          axisRatios (:,i)=axisLengths/axisLengths(1)
       end do
       ! Store results to file.
       call simulations(iSimulation)%analysis%writeDataset(inertiaTensor,'inertiaTensor'            )
       call simulations(iSimulation)%analysis%writeDataset(eigenVectors ,'inertiaTensorEigenVectors')
       call simulations(iSimulation)%analysis%writeDataset(eigenValues  ,'inertiaTensorEigenValues' )
       call simulations(iSimulation)%analysis%writeDataset(axisRatios   ,'ellipsoidAxisRatios'      )
       ! Deallocate workspace.
       if (self%selfBoundParticlesOnly) then
          nullify   (selfBoundStatus)
       else
          deallocate(selfBoundStatus)
       end if
       nullify      (indexMostBound )
       deallocate   (inertiaTensor  )
       deallocate   (eigenVectors   )
       deallocate   (eigenValues    )
       deallocate   (positionOffset )
       deallocate   (radiusOffset   )
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine inertiaTensorOperate

