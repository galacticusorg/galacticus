!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements an N-body data operator which determines the environmental overoverdensity around particles.
  
  use, intrinsic :: ISO_C_Binding

  !# <nbodyOperator name="nbodyOperatorEnvironmentalOverdensity">
  !#  <description>An N-body data operator which determines the environmental overoverdensity around particles.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorEnvironmentalOverdensity
     !% An N-body data operator which determines the environmental overoverdensity around particles.
     private
     double precision           :: radiusSphere     , densityParticleMean, &
          &                        particleCountMean, lengthBox
     integer         (c_size_t) :: sampleRate
     logical                    :: periodic
   contains
     procedure :: operate => environmentalOverdensityOperate
  end type nbodyOperatorEnvironmentalOverdensity

  interface nbodyOperatorEnvironmentalOverdensity
     !% Constructors for the ``environmentalOverdensity'' N-body operator class.
     module procedure environmentalOverdensityConstructorParameters
     module procedure environmentalOverdensityConstructorInternal
  end interface nbodyOperatorEnvironmentalOverdensity

contains

  function environmentalOverdensityConstructorParameters(parameters) result (self)
    !% Constructor for the ``environmentalOverdensity'' N-body operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (nbodyOperatorEnvironmentalOverdensity)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: radiusSphere, densityParticleMean, &
         &                                                                    lengthBox
    integer         (c_size_t                             )                :: sampleRate
    logical                                                                   periodic

    !# <inputParameter>
    !#   <name>radiusSphere</name>
    !#   <source>parameters</source>
    !#   <description>The radius of the sphere within which to measure environmental overdensity.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>densityParticleMean</name>
    !#   <source>parameters</source>
    !#   <description>The mean density of particles in the simulation.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sampleRate</name>
    !#   <source>parameters</source>
    !#   <description>One in {\normalfont \ttfamily [sampleRate]} particles will be sampled when computed environmental overdensities.</description>
    !#   <defaultValue>1_c_size_t</defaultValue>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>lengthBox</name>
    !#   <source>parameters</source>
    !#   <description>The length of the periodic box.</description>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>periodic</name>
    !#   <source>parameters</source>
    !#   <description>If true, periodic boundary conditions will be used.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyOperatorEnvironmentalOverdensity(radiusSphere,densityParticleMean,sampleRate,lengthBox,periodic)
    !# <inputParametersValidate source="parameters"/>
    return
  end function environmentalOverdensityConstructorParameters

  function environmentalOverdensityConstructorInternal(radiusSphere,densityParticleMean,sampleRate,lengthBox,periodic) result (self)
    !% Internal constructor for the ``environmentalOverdensity'' N-body operator class.
    use Numerical_Constants_Math
    implicit none
    type            (nbodyOperatorEnvironmentalOverdensity)                :: self
    double precision                                       , intent(in   ) :: radiusSphere, densityParticleMean, &
         &                                                                    lengthBox
    integer         (c_size_t                             ), intent(in   ) :: sampleRate
    logical                                                , intent(in   ) :: periodic
    !# <constructorAssign variables="radiusSphere, densityParticleMean, sampleRate, lengthBox, periodic"/>

    ! Evaluate mean number of particles in the sphere.
    self%particleCountMean=+self%densityParticleMean    &
         &                 *self%radiusSphere       **3 &
         &                 *4.0d0                       &
         &                 *Pi                          &
         &                 /3.0d0
    return
  end function environmentalOverdensityConstructorInternal

  subroutine environmentalOverdensityOperate(self,simulation)
    !% Determine the mean position and velocity of N-body particles.
    use Nearest_Neighbors
    use Memory_Management
    use Galacticus_Display
    use OMP_Lib
    implicit none
    class           (nbodyOperatorEnvironmentalOverdensity), intent(inout)                 :: self
    type            (nBodyData                            ), intent(inout)                 :: simulation
    double precision                                       , parameter                     :: toleranceZero          =0.0d0
    integer                                                , allocatable  , dimension(:  ) :: neighborIndex
    double precision                                       , allocatable  , dimension(:  ) :: neighborDistance             , overdensity
    double precision                                       , allocatable  , dimension(:,:) :: position
    logical                                                , allocatable  , dimension(:  ) :: particleMask
    type            (nearestNeighbors                     )                                :: neighborFinder
    integer                                                                                :: neighborCount
    integer                                                               , dimension(3  ) :: i
    integer                                                                                :: l                            , m                     , &
         &                                                                                    n
    integer         (c_size_t                             )                                :: particleCount                , particleCountReplicant, &
         &                                                                                    j                            , k

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
                particleCount=+particleCount                                        &
                     &        +count(                                               &
                     &                simulation%position(1,:) >= boundLower(i(1))  &
                     &               .and.                                          &
                     &                simulation%position(1,:) <  boundUpper(i(1))  &
                     &               .and.                                          &
                     &                simulation%position(2,:) >= boundLower(i(2))  &
                     &               .and.                                          &
                     &                simulation%position(2,:) <  boundUpper(i(2))  &
                     &               .and.                                          &
                     &                simulation%position(3,:) >= boundLower(i(3))  &
                     &               .and.                                          &
                     &                simulation%position(3,:) <  boundUpper(i(3)), &
                     &               kind=c_size_t                                  &
                     &              )
             end do
          end do
       end do
       call allocateArray(position    ,[3_c_size_t,particleCount                                ])
       call allocateArray(particleMask,[           size(simulation%position,dim=2,kind=c_size_t)])
       particleCount=0_c_size_t
       do l=-1,+1          
          i(1)=l
          do m=-1,+1
             i(2)=m
             do n=-1,+1
                i(3)=n
                particleMask= simulation%position(1,:) >= boundLower(i(1)) &
                     &       .and.                                         &
                     &        simulation%position(1,:) <  boundUpper(i(1)) &
                     &       .and.                                         &
                     &        simulation%position(2,:) >= boundLower(i(2)) &
                     &       .and.                                         &
                     &        simulation%position(2,:) <  boundUpper(i(2)) &
                     &       .and.                                         &
                     &        simulation%position(3,:) >= boundLower(i(3)) &
                     &       .and.                                         &
                     &        simulation%position(3,:) <  boundUpper(i(3))
                particleCountReplicant=count(particleMask,kind=c_size_t)
                do j=1,3
                   position(j,particleCount+1_c_size_t:particleCount+particleCountReplicant)=pack(simulation%position(j,:),particleMask)+float(i(j))*self%lengthBox
                end do
                particleCount=particleCount+particleCountReplicant
             end do
          end do
       end do
      call deallocateArray(particleMask)
     else
       ! Non-periodic boundaries - no need to replicate particles.
       call allocateArray(position,shape(simulation%position))
       position=simulation%position
    end if
    neighborFinder=nearestNeighbors(transpose(position))
    call allocateArray(overdensity,[size(simulation%position,dim=2,kind=c_size_t)])
    overdensity=-2.0d0
    ! Iterate over particles.
    call Galacticus_Display_Counter(0,.true.)
    j=0_c_size_t
    !$omp parallel do private(neighborCount,neighborIndex,neighborDistance) schedule(dynamic)
    do k=1_c_size_t,size(simulation%position,dim=2,kind=c_size_t),self%sampleRate
       ! Locate particles nearby.
       call neighborFinder%searchFixedRadius(simulation%position(:,k),self%radiusSphere,toleranceZero,neighborCount,neighborIndex,neighborDistance)
       overdensity(k)=float(neighborCount)/self%particleCountMean-1.0d0
       !$omp atomic
       j=j+self%sampleRate
       if (omp_get_thread_num() == 0) call Galacticus_Display_Counter(int(100.0d0*float(j)/float(size(simulation%position,dim=2,kind=c_size_t))),.false.)
    end do
    !$omp end parallel do
    call Galacticus_Display_Counter_Clear()
    call deallocateArray(position)
    call simulation%analysis%writeDataset(overdensity,'overdensityEnvironmental')
    return

  contains

    double precision function boundLower(l)
      !% Compute lower bounds for particle inclusion in periodically replicated volumes.
      use Galacticus_Error
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
         call Galacticus_Error_Report('replicant index must be -1, 0, or +1'//{introspection:location})
      end select
      return
    end function boundLower
    
    double precision function boundUpper(l)
      !% Compute upper bounds for particle inclusion in periodically replicated volumes.
      use Galacticus_Error
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
         call Galacticus_Error_Report('replicant index must be -1, 0, or +1'//{introspection:location})
      end select
      return
    end function boundUpper
    
  end subroutine environmentalOverdensityOperate

