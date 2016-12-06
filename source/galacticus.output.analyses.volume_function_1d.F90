!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a generic 1D volume function (i.e. number density of objects binned by some property, e.g. a
!% mass function) output analysis class.

  !$ use            OMP_Lib
  use, intrinsic :: ISO_C_Binding
  use               Output_Analysis_Property_Extractions
  use               Output_Analysis_Property_Operators
  use               Output_Analysis_Distribution_Operators
  use               Output_Analysis_Distribution_Normalizers
  use               Galactic_Filters

  !# <outputAnalysis name="outputAnalysisVolumeFunction1D" defaultThreadPrivate="yes">
  !#  <description>A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisClass) :: outputAnalysisVolumeFunction1D
     !% A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output
     !% analysis class.
     private
     class           (outputAnalysisPropertyExtractorClass     ), pointer                     :: outputAnalysisPropertyExtractor_      => null()
     class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_       => null()
     class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_ => null()
     class           (galacticFilterClass                      ), pointer                     :: galacticFilter_                       => null()
     double precision                                           , dimension(:,:), allocatable :: outputWeight
     double precision                                           , dimension(:  ), allocatable :: binCenter                                      , functionValue
     double precision                                           , dimension(:  ), allocatable :: binMinimum                                     , binMaximum
     integer         (c_size_t                                 )                              :: binCount
     !$ integer      (omp_lock_kind                            )                              :: accumulateLock
   contains
     final     ::             volumeFunction1DDestructor
     procedure :: analyze  => volumeFunction1DAnalyze
     procedure :: finalize => volumeFunction1DFinalize
  end type outputAnalysisVolumeFunction1D

  interface outputAnalysisVolumeFunction1D
     !% Constructors for the ``volumeFunction1D'' output analysis class.
     module procedure volumeFunction1DConstructorParameters
     module procedure volumeFunction1DConstructorInternal
  end interface outputAnalysisVolumeFunction1D

contains

  function volumeFunction1DConstructorParameters(parameters) result(self)
    !% Constructor for the ``volumeFunction1D'' output analysis class which takes a parameter set as input.
    use Input_Parameters2
    use Galacticus_Output_Times
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (outputAnalysisPropertyExtractorClass     ), pointer                     :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    double precision                                           , dimension(:  ), allocatable :: binCenter                            , outputWeight
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    call allocateArray(binCenter   ,[int(parameters%count('binCenter'),kind=c_size_t)                               ])
    call allocateArray(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t)*Galacticus_Output_Time_Count()])
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*Galacticus_Output_Time_Count()) &
         & call Galacticus_Error_Report('volumeFunction1DConstructorParameters','incorrect number of output weights provided')
    !# <inputParameter>
    !#   <name>binCenter</name>
    !#   <source>parameters</source>
    !#   <variable>binCenter</variable>
    !#   <description>The value of the property at the center of each bin.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputWeight</name>
    !#   <source>parameters</source>
    !#   <variable>outputWeight</variable>
    !#   <description>The weight to assign to each bin at each output.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisPropertyExtractor_"      source="parameters"/>
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"/>
    !# <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"/>
    !# <objectBuilder class="outputAnalysisDistributionNormalizer" name="outputAnalysisDistributionNormalizer_" source="parameters"/>
    !# <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"/>
    ! Build the object.
    self=outputAnalysisVolumeFunction1D(                                                                                                         &
         &                              binCenter                                                                                              , &
         &                              reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),Galacticus_Output_Time_Count()]), &
         &                              outputAnalysisPropertyExtractor_                                                                       , &
         &                              outputAnalysisPropertyOperator_                                                                        , &
         &                              outputAnalysisDistributionOperator_                                                                    , &
         &                              outputAnalysisDistributionNormalizer_                                                                  , &
         &                              galacticFilter_                                                                                          &
         &                             )
    return
  end function volumeFunction1DConstructorParameters

  function volumeFunction1DConstructorInternal(binCenter,outputWeight,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisDistributionNormalizer_,galacticFilter_) result (self)
    !% Constructor for the ``volumeFunction1D'' output analysis class for internal use.
    use Memory_Management
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                                :: self
    double precision                                           , intent(in   ), dimension(:  ) :: binCenter
    double precision                                           , intent(in   ), dimension(:,:) :: outputWeight
    class           (outputAnalysisPropertyExtractorClass     ), intent(in   ), target         :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), intent(in   ), target         :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass  ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), intent(in   ), target         :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), intent(in   ), target         :: galacticFilter_
    integer         (c_size_t                                 )                                :: i
    !# <constructorAssign variables="binCenter, outputWeight, *outputAnalysisPropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisDistributionOperator_, *outputAnalysisDistributionNormalizer_, *galacticFilter_"/>

    ! Count bins.
    self%binCount=size(binCenter,kind=c_size_t)
    ! Determine bin minima and maxima.
    call allocateArray(self%binMinimum,[self%binCount])
    call allocateArray(self%binMaximum,[self%binCount])
    do i=1,self%binCount
       if (i == 1) then
          self%binMinimum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i+1))
       else
          self%binMinimum(i)=             +0.5d0*(self%binCenter(i-1)+self%binCenter(i  ))
       end if
       if (i == self%binCount) then
          self%binMaximum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i-1))
       else
          self%binMaximum(i)=             +0.5d0*(self%binCenter(i  )+self%binCenter(i+1))
       end if
    end do
    ! Allocate and initialize function values.
    call allocateArray(self%functionValue,[self%binCount])
    self%functionValue=0.0d0
    ! Initialize OpenMP accumulation lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
   return
  end function volumeFunction1DConstructorInternal

  subroutine volumeFunction1DDestructor(self)
    !% Destructor for  the ``volumeFunction1D'' output analysis class.
    type(outputAnalysisVolumeFunction1D), intent(inout) :: self
    
    !# <objectDestructor name="self%outputAnalysisPropertyExtractor_"     />
    !# <objectDestructor name="self%outputAnalysisPropertyOperator_"      />
    !# <objectDestructor name="self%outputAnalysisDistributionOperator_"  />
    !# <objectDestructor name="self%outputAnalysisDistributionNormalizer_"/>
    !# <objectDestructor name="self%galacticFilter_"                      />
    ! Destroy OpenMP lock.
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine volumeFunction1DDestructor
  
  subroutine volumeFunction1DAnalyze(self,node,iOutput)
    !% Implement a volumeFunction1D output analysis.
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout)            :: self
    type            (treeNode                      ), intent(inout)            :: node
    integer         (c_size_t                      ), intent(in   )            :: iOutput
    double precision                                , dimension(self%binCount) :: distribution
    double precision                                                           :: propertyValue
    
    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    !! AJB HACK
    write (0,*) " -> post weight"
 
    ! Extract the property from the node.
    propertyValue=self%outputAnalysisPropertyExtractor_%extract(node)
    !! AJB HACK
    write (0,*) " -> raw value: ",propertyValue

    ! Filter on the property.
    if (.not.self%galacticFilter_%passes(node)) return
    !! AJB HACK
    write (0,*) " -> post filter raw value: ",propertyValue
    
    ! Apply property operators.
    propertyValue=self%outputAnalysisPropertyOperator_%operate(propertyValue)
    !! AJB HACK
    write (0,*) " -> post operator value: ",propertyValue
    !! Instances:
    !!!  Cosmology shift.
    !!!  Systematic shift.
    
    ! Apply distribution operators.
    distribution=self%outputAnalysisDistributionOperator_%operateScalar(propertyValue,self%binMinimum,self%binMaximum)
    !! Need to add in some buffer to the edges of the distribution to allow operators to be unaffected by edge effects.
    !! Instances:
    !!!  Gravitational lensing.
    !!!  Cosmology (volume) shift.
    !! AJB HACK
    write (0,*) " -> post operator distribution: ",distribution
 
    ! Accumulate the property, including weights from both the host tree and the output.
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%functionValue=+self%functionValue           &
         &             +distribution                 &
         &             *self%outputWeight(:,iOutput) &
         &             *node%hostTree%volumeWeight
    !$ call OMP_Unset_Lock(self%accumulateLock)
    !! AJB HACK
    write (0,*) " -> post accumulation function: ",self%functionValue
 
    ! Accumulate covariance.
    
    return
  end subroutine volumeFunction1DAnalyze

  subroutine volumeFunction1DFinalize(self)
    !% Implement a volumeFunction1D output analysis finalization.
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout) :: self

    ! Apply final distribution operators.
    self%functionValue=self%outputAnalysisDistributionNormalizer_%normalize(self%functionValue,self%binMinimum,self%binMaximum)
    !! Bin width
    !! dlog10->dlog
    !! AJB HACK
    write (0,*) " -> post normalize ",self%functionValue
    
    ! Output.
    
    return
  end subroutine volumeFunction1DFinalize
