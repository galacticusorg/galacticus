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
Implements an N-body data operator which determines the environmental overdensity around particles.
!!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorEnvironmentalOverdensity">
   <description>An N-body data operator which determines the environmental overdensity around particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorEnvironmentalOverdensity
     !!{
     An N-body data operator which determines the environmental overdensity around particles.
     !!}
     private
     double precision           :: radiusSphere     , densityParticleMean, &
          &                        particleCountMean, lengthBox
     integer         (c_size_t) :: sampleRate
     logical                    :: periodic
   contains
     procedure :: operate => environmentalOverdensityOperate
  end type nbodyOperatorEnvironmentalOverdensity

  interface nbodyOperatorEnvironmentalOverdensity
     !!{
     Constructors for the \refClass{nbodyOperatorEnvironmentalOverdensity} N-body operator class.
     !!}
     module procedure environmentalOverdensityConstructorParameters
     module procedure environmentalOverdensityConstructorInternal
  end interface nbodyOperatorEnvironmentalOverdensity

contains

  function environmentalOverdensityConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorEnvironmentalOverdensity} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorEnvironmentalOverdensity)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: radiusSphere, densityParticleMean, &
         &                                                                    lengthBox
    integer         (c_size_t                             )                :: sampleRate
    logical                                                                :: periodic

    !![
    <inputParameter>
      <name>radiusSphere</name>
      <source>parameters</source>
      <description>The radius of the sphere within which to measure environmental overdensity.</description>
    </inputParameter>
    <inputParameter>
      <name>densityParticleMean</name>
      <source>parameters</source>
      <description>The mean density of particles in the simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>sampleRate</name>
      <source>parameters</source>
      <description>One in {\normalfont \ttfamily [sampleRate]} particles will be sampled when computed environmental overdensities.</description>
      <defaultValue>1_c_size_t</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>lengthBox</name>
      <source>parameters</source>
      <description>The length of the periodic box.</description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>periodic</name>
      <source>parameters</source>
      <description>If true, periodic boundary conditions will be used.</description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    !!]
    self=nbodyOperatorEnvironmentalOverdensity(radiusSphere,densityParticleMean,sampleRate,lengthBox,periodic)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function environmentalOverdensityConstructorParameters

  function environmentalOverdensityConstructorInternal(radiusSphere,densityParticleMean,sampleRate,lengthBox,periodic) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorEnvironmentalOverdensity} N-body operator class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (nbodyOperatorEnvironmentalOverdensity)                :: self
    double precision                                       , intent(in   ) :: radiusSphere, densityParticleMean, &
         &                                                                    lengthBox
    integer         (c_size_t                             ), intent(in   ) :: sampleRate
    logical                                                , intent(in   ) :: periodic
    !![
    <constructorAssign variables="radiusSphere, densityParticleMean, sampleRate, lengthBox, periodic"/>
    !!]

    ! Evaluate mean number of particles in the sphere.
    self%particleCountMean=+self%densityParticleMean    &
         &                 *self%radiusSphere       **3 &
         &                 *4.0d0                       &
         &                 *Pi                          &
         &                 /3.0d0
    return
  end function environmentalOverdensityConstructorInternal

  subroutine environmentalOverdensityOperate(self,simulations)
    !!{
    Determine the mean position and velocity of N-body particles.
    !!}
    use    :: Display          , only : displayCounter    , displayCounterClear
    use    :: Nearest_Neighbors, only : nearestNeighbors
    !$ use :: OMP_Lib          , only : omp_get_thread_num
    implicit none
    class           (nbodyOperatorEnvironmentalOverdensity), intent(inout)                 :: self
    type            (nBodyData                            ), intent(inout), dimension(:  ) :: simulations
    double precision                                       , parameter                     :: toleranceZero          =0.0d0
    integer                                                , allocatable  , dimension(:  ) :: neighborIndex
    double precision                                       , allocatable  , dimension(:  ) :: neighborDistance             , overdensity
    double precision                                       , allocatable  , dimension(:,:) :: position
    double precision                                       , pointer      , dimension(:,:) :: position_
    logical                                                , allocatable  , dimension(:  ) :: particleMask
    type            (nearestNeighbors                     )                                :: neighborFinder
    integer                                                                                :: neighborCount
    integer                                                               , dimension(3  ) :: i
    integer                                                                                :: l                            , m                     , &
         &                                                                                    n                            , iSimulation
    integer         (c_size_t                             )                                :: particleCount                , particleCountReplicant, &
         &                                                                                    j                            , k

    do iSimulation=1,size(simulations)
       position_ => simulations(iSimulation)%propertiesRealRank1%value('position')
       if (self%periodic) then
          ! Periodic boundary conditions - add buffers of periodically replicated particles.
          ! First count how many particles total we have including the replicated buffers.
          particleCount=0_c_size_t
          do l=-1,+1
             i(1)=l
             do m=-1,+1
                i(2)=m
                do n=-1,+1
                   i(3)=n
                   particleCount=+particleCount                              &
                        &        +count(                                     &
                        &                position_(1,:) >= boundLower(i(1))  &
                        &               .and.                                &
                        &                position_(1,:) <  boundUpper(i(1))  &
                        &               .and.                                &
                        &                position_(2,:) >= boundLower(i(2))  &
                        &               .and.                                &
                        &                position_(2,:) <  boundUpper(i(2))  &
                        &               .and.                                &
                        &                position_(3,:) >= boundLower(i(3))  &
                        &               .and.                                &
                        &                position_(3,:) <  boundUpper(i(3)), &
                        &               kind=c_size_t                        &
                        &              )
                end do
             end do
          end do
          allocate(position    (3_c_size_t,particleCount        ))
          allocate(particleMask(           size(position_,dim=2)))
          particleCount=0_c_size_t
          do l=-1,+1
             i(1)=l
             do m=-1,+1
                i(2)=m
                do n=-1,+1
                   i(3)=n
                   particleMask= position_(1,:) >= boundLower(i(1)) &
                        &       .and.                               &
                        &        position_(1,:) <  boundUpper(i(1)) &
                        &       .and.                               &
                        &        position_(2,:) >= boundLower(i(2)) &
                        &       .and.                               &
                        &        position_(2,:) <  boundUpper(i(2)) &
                        &       .and.                               &
                        &        position_(3,:) >= boundLower(i(3)) &
                        &       .and.                               &
                        &        position_(3,:) <  boundUpper(i(3))
                   particleCountReplicant=count(particleMask,kind=c_size_t)
                   do j=1,3
                      position(j,particleCount+1_c_size_t:particleCount+particleCountReplicant)=pack(position_(j,:),particleMask)+float(i(j))*self%lengthBox
                   end do
                   particleCount=particleCount+particleCountReplicant
                end do
             end do
          end do
          deallocate(particleMask)
       else
          ! Non-periodic boundaries - no need to replicate particles.
          allocate(position,mold=position_)
          position=position_
       end if
       neighborFinder=nearestNeighbors(transpose(position))
       allocate(overdensity(size(position_,dim=2)))
       overdensity=-2.0d0
       ! Iterate over particles.
       call displayCounter(0,.true.)
       j=0_c_size_t
       !$omp parallel do private(neighborCount,neighborIndex,neighborDistance) schedule(dynamic)
       do k=1_c_size_t,size(position_,dim=2,kind=c_size_t),self%sampleRate
          ! Locate particles nearby.
          call neighborFinder%searchFixedRadius(position_(:,k),self%radiusSphere,toleranceZero,neighborCount,neighborIndex,neighborDistance)
          overdensity(k)=float(neighborCount)/self%particleCountMean-1.0d0
          !$omp atomic
          j=j+self%sampleRate
          !$ if (omp_get_thread_num() == 0) then
          call displayCounter(int(100.0d0*float(j)/float(size(position_,dim=2,kind=c_size_t))),.false.)
          !$ end if
       end do
       !$omp end parallel do
       call displayCounterClear()
       deallocate(position)
       call simulations(iSimulation)%analysis%writeDataset(overdensity,'overdensityEnvironmental')
    end do
    return

  contains

    double precision function boundLower(l)
      !!{
      Compute lower bounds for particle inclusion in periodically replicated volumes.
      !!}
      use :: Error, only : Error_Report
      implicit none
      integer, intent(in   ) :: l

      select case (l)
      case (-1)
         boundLower=+self%lengthBox-self%radiusSphere
      case ( 0)
         boundLower=+0.0d0
      case (+1)
         boundLower=+0.0d0
      case default
         boundLower=+0.0d0
         call Error_Report('replicant index must be -1, 0, or +1'//{introspection:location})
      end select
      return
    end function boundLower

    double precision function boundUpper(l)
      !!{
      Compute upper bounds for particle inclusion in periodically replicated volumes.
      !!}
      use :: Error, only : Error_Report
      implicit none
      integer, intent(in   ) :: l

      select case (l)
      case (-1)
         boundUpper=+self%lengthBox
      case ( 0)
         boundUpper=+self%lengthBox
      case (+1)
         boundUpper=+self%radiusSphere
      case default
         boundUpper=+0.0d0
         call Error_Report('replicant index must be -1, 0, or +1'//{introspection:location})
      end select
      return
    end function boundUpper

  end subroutine environmentalOverdensityOperate

