!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data operator which computes pair counts in bins of separation.

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !# <nbodyOperator name="nbodyOperatorPairCounts">
  !#  <description>An N-body data operator which computes pair counts in bins of separation.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorPairCounts
     !% An N-body data operator which determines the environmental overoverdensity around particles.
     private
     double precision                                       :: separationMinimum              , separationMaximum   , &
          &                                                   bootstrapSampleRate
     integer         (c_size_t                  )          :: separationCount                 , bootstrapSampleCount
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            pairCountsDestructor
     procedure :: operate => pairCountsOperate
  end type nbodyOperatorPairCounts

  interface nbodyOperatorPairCounts
     !% Constructors for the ``pairCounts'' N-body operator class.
     module procedure pairCountsConstructorParameters
     module procedure pairCountsConstructorInternal
  end interface nbodyOperatorPairCounts

contains

  function pairCountsConstructorParameters(parameters) result (self)
    !% Constructor for the ``pairCounts'' N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorPairCounts   )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    double precision                                            :: separationMinimum     , separationMaximum   , &
         &                                                         bootstrapSampleRate
    integer         (c_size_t                  )                :: separationCount       , bootstrapSampleCount

    !# <inputParameter>
    !#   <name>bootstrapSampleCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>30_c_size_t</defaultValue>
    !#   <description>The number of bootstrap resamples of the particles that should be used.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bootstrapSampleRate</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The sampling rate for particles.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMinimum</name>
    !#   <source>parameters</source>
    !#   <description>The minimum separation to consider for pair counts.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMaximum</name>
    !#   <source>parameters</source>
    !#   <description>The maximum separation to consider for pair counts.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationCount</name>
    !#   <source>parameters</source>
    !#   <description>The number of bins in separation for pair counts.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    self=nbodyOperatorPairCounts(separationMinimum,separationMaximum,separationCount,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="randomNumberGenerator_"/>
    return
  end function pairCountsConstructorParameters

  function pairCountsConstructorInternal(separationMinimum,separationMaximum,separationCount,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_) result (self)
    !% Internal constructor for the ``pairCounts'' N-body operator class.
    implicit none
    type            (nbodyOperatorPairCounts)                           :: self
    double precision                            , intent(in   )         :: separationMinimum     , separationMaximum   , &
         &                                                                 bootstrapSampleRate
    integer         (c_size_t                  ), intent(in   )         :: separationCount       , bootstrapSampleCount
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !# <constructorAssign variables="separationMinimum, separationMaximum, separationCount, bootstrapSampleCount, bootstrapSampleRate, *randomNumberGenerator_"/>

    return
  end function pairCountsConstructorInternal

  subroutine pairCountsDestructor(self)
    !% Destructor for the ``meanPosition'' N-body operator class.
    implicit none
    type(nbodyOperatorPairCounts), intent(inout) :: self

    !# <objectDestructor name="self%randomNumberGenerator_"/>
    return
  end subroutine pairCountsDestructor
  
  subroutine pairCountsOperate(self,simulation)
    !% Determine the mean position and velocity of N-body particles.
    use :: Galacticus_Display, only : Galacticus_Display_Counter, Galacticus_Display_Counter_Clear
    use :: Nearest_Neighbors , only : nearestNeighbors
    use :: Numerical_Ranges  , only : Make_Range                , rangeTypeLogarithmic
    implicit none
    class           (nbodyOperatorPairCounts), intent(inout)                 :: self
    type            (nBodyData              ), intent(inout)                 :: simulation
    double precision                         , parameter                     :: toleranceZero       =0.0d0
    integer                                  , allocatable  , dimension(:  ) :: neighborIndex
    double precision                         , allocatable  , dimension(:  ) :: neighborDistance          , separationCentralBin, &
         &                                                                      separationMinimumBin      , separationMaximumBin
    integer                                  , allocatable  , dimension(:  ) :: weight
    integer         (c_size_t               ), allocatable  , dimension(:,:) :: pairCountBin
    integer         (c_size_t               )                                :: i                         , j                   , &
         &                                                                      iSample
    type            (nearestNeighbors       )                                :: neighborFinder
    integer                                                                  :: neighborCount

    ! Construct bins of separation.
    allocate(separationCentralBin(self%separationCount                                     ))
    allocate(separationMinimumBin(self%separationCount                                     ))
    allocate(separationMaximumBin(self%separationCount                                     ))
    allocate(pairCountBin        (self%separationCount           ,self%bootstrapSampleCount))
    allocate(weight              (size(simulation%position,dim=2)                          ))
    separationCentralBin=Make_Range(self%separationMinimum,self%separationMaximum,int(self%separationCount),rangeTypeLogarithmic,rangeBinned=.true.)
    separationMinimumBin=separationCentralBin/sqrt(separationCentralBin(2)/separationCentralBin(1))
    separationMaximumBin=separationCentralBin*sqrt(separationCentralBin(2)/separationCentralBin(1))
    pairCountBin=0_c_size_t
    ! Iterate over bootstrap samplings.
    call Galacticus_Display_Counter(0,.true.)
    do iSample=1,self%bootstrapSampleCount
       ! Determine weights for particles.
       do i=1,size(simulation%position,dim=2,kind=c_size_t)
          weight(i)=self%randomNumberGenerator_%poissonSample(self%bootstrapSampleRate)
       end do
       ! Iterate over particles.
       !$omp parallel private(neighborCount,neighborIndex,neighborDistance,neighborFinder)
       ! Construct the nearest neighbor finder object.
       neighborFinder=nearestNeighbors(transpose(simulation%position))
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(simulation%position,dim=2,kind=c_size_t)
          if (weight(i) <= 0.0d0) cycle
          ! Locate particles nearby.
          call neighborFinder%searchFixedRadius(simulation%position(:,i),self%separationMaximum,toleranceZero,neighborCount,neighborIndex,neighborDistance)
          ! Accumulate particles into bins.
          do j=1,self%separationCount
             if (weight(j) <= 0.0d0) cycle
             !$omp atomic
             pairCountBin(j,iSample)=pairCountBin(j,iSample)+int(weight(i),c_size_t)*int(weight(j),c_size_t)*count(neighborDistance >= separationMinimumBin(j) .and. neighborDistance < separationMaximumBin(j))
          end do
          ! Update progress.
          call Galacticus_Display_Counter(int(100.0d0*float(i+(iSample-1_c_size_t)*size(simulation%position,dim=2,kind=c_size_t))/float(size(simulation%position,dim=2,kind=c_size_t))/float(self%bootstrapSampleCount)),.false.)
       end do
       !$omp end do
       !$omp end parallel
    end do
    call Galacticus_Display_Counter_Clear()
    call simulation%analysis%writeDataset(pairCountBin        ,'pairCountCount'     )
    call simulation%analysis%writeDataset(separationCentralBin,'pairCountSeparation')
    return
  end subroutine pairCountsOperate

