!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson, Christoph Behrens, Xiaolong Du.

!% Contains a module which implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.

  use, intrinsic :: ISO_C_Binding
  use            :: FGSL                   , only : fgsl_interp_accel
  use            :: Excursion_Sets_Barriers
  use            :: File_Utilities
  use            :: Cosmology_Functions

  !# <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahi">
  !#  <description>An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.</description>
  !# </excursionSetFirstCrossing>
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingFarahi
     !% An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.
     private
     class           (cosmologyFunctionsClass ), pointer                       :: cosmologyFunctions_                     => null()
     class           (excursionSetBarrierClass), pointer                       :: excursionSetBarrier_                    => null()
     ! Variables used in tabulation the first crossing function.
     double precision                                                          :: timeMaximum                                      , timeMinimum                         , &
          &                                                                       varianceMaximum                         
     integer                                                                   :: timeTableCount                                   , varianceTableCount
     double precision                          , allocatable, dimension(:,:)   :: firstCrossingProbabilityTable
     double precision                          , allocatable, dimension(:  )   :: timeTable                                        , varianceTable
     double precision                                                          :: varianceTableStep
     logical                                                                   :: tableInitialized                        
     type            (fgsl_interp_accel       )                                :: interpolationAcceleratorTime                     , interpolationAcceleratorVariance
     logical                                                                   :: interpolationResetTime                           , interpolationResetVariance          
     ! Variables used in tabulation the first crossing rate function.
     double precision                                                          :: timeMaximumRate                                  , timeMinimumRate                     , &
          &                                                                       varianceMaximumRate                     
     integer                                                                   :: timeTableCountRate                               , varianceTableCountRate              , &
          &                                                                       varianceTableCountRateBase
     double precision                          , allocatable, dimension(:,:,:) :: firstCrossingTableRate
     double precision                          , allocatable, dimension(:,:  ) :: nonCrossingTableRate
     double precision                          , allocatable, dimension(:    ) :: timeTableRate                                    , varianceTableRate                   , &
          &                                                                       varianceTableRateBase
     logical                                                                   :: tableInitializedRate                    
     type            (fgsl_interp_accel       )                                :: interpolationAcceleratorTimeRate                 , interpolationAcceleratorVarianceRate, &
          &                                                                       interpolationAcceleratorVarianceRateBase
     logical                                                                   :: interpolationResetTimeRate                       , interpolationResetVarianceRate      , &
          &                                                                       interpolationResetVarianceRateBase      
     ! File name used to store tabulations.
     type            (varying_string          )                                :: fileName
     logical                                                                   :: useFile
     ! The fractional step in time used to compute barrier crossing rates.
     double precision                                                          :: timeStepFractional
     ! Record of variance and time in previous call to rate functions.
     double precision                                                          :: timeRatePrevious                                 , varianceRatePrevious
     double precision                                       , dimension(0:1)   :: hTimeRate                                        , hVarianceRate
     integer         (c_size_t                )                                :: iTimeRate                                        , iVarianceRate
   contains
     !@ <objectMethods>
     !@   <object>excursionSetFirstCrossingFarahi</object>
     !@   <objectMethod>
     !@     <method>rateTabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ varianceProgenitor\argin, \doublezero\ time\argin</arguments>
     !@     <description>Tabulate excursion set barrier crossing rates ensuring that they span the given progenitor variance and time.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>varianceRange</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\doublezero\ rangeMinimum\argin, \doublezero\ rangeMaximum\argin, \intzero\ rangeNumber, \doublezero\ ratioAtMaximum</arguments>
     !@     <description>Build a range of variances at which to tabulate the excursion set solutions.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fileRead</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Read excursion set solutions from file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fileWrite</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Write excursion set solutions to file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                    farahiDestructor
     procedure :: probability     => farahiProbability
     procedure :: rate            => farahiRate
     procedure :: rateNonCrossing => farahiRateNonCrossing
     procedure :: rateTabulate    => farahiRateTabulate
     procedure :: fileRead        => farahiFileRead
     procedure :: fileWrite       => farahiFileWrite
     procedure :: varianceRange   => farahiVarianceRange
  end type excursionSetFirstCrossingFarahi

  interface excursionSetFirstCrossingFarahi
     !% Constructors for the Farahi excursion set barrier class.
     module procedure farahiConstructorParameters
     module procedure farahiConstructorInternal
  end interface excursionSetFirstCrossingFarahi

  ! Parameters controlling tabulation range and granularity.
  integer                         , parameter :: farahiVarianceNumberPerUnitProbability=1000
  integer                         , parameter :: farahiTimeNumberPerDecade             =  10    , farahiVarianceNumberPerDecade=400    , &
       &                                         farahiVarianceNumberPerUnit           =  40
  double precision                , parameter :: farahiRateRedshiftMaximum             =  30.0d0, farahiRateRedshiftMinimum    =  0.0d0

  ! Lock used for file access.
  type            (lockDescriptor)            :: farahiFileLock
  logical                                     :: farahiFileLockInitialized              =.false.

contains

  function farahiConstructorParameters(parameters) result(self)
    !% Constructor for the Farahi excursion set class first crossing class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (excursionSetFirstCrossingFarahi)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (excursionSetBarrierClass       ), pointer       :: excursionSetBarrier_
    double precision                                                 :: timeStepFractional
    type            (varying_string                 )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <defaultValue>var_str('none')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !#   <description>The name of the file to/from which tabulations of barrier first crossing probabilities should be written/read. If set to ``{\normalfont \ttfamily none}'' tables will not be stored.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeStepFractional</name>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !#   <description>The fractional time step used when computing barrier crossing rates (i.e. the step used in finite difference calculations).</description>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="excursionSetBarrier" name="excursionSetBarrier_" source="parameters"/>
    self=excursionSetFirstCrossingFarahi(timeStepFractional,fileName,cosmologyFunctions_,excursionSetBarrier_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="excursionSetBarrier_"/>
    return
  end function farahiConstructorParameters

  function farahiConstructorInternal(timeStepFractional,fileName,cosmologyFunctions_,excursionSetBarrier_) result(self)
    !% Internal constructor for the Farahi excursion set class first crossing class.
    use Input_Parameters
    use Galacticus_Paths
    use Galacticus_Display
    use System_Command
    use File_Utilities
    implicit none
    type            (excursionSetFirstCrossingFarahi)                        :: self
    double precision                                 , intent(in   )         :: timeStepFractional
    type            (varying_string                 ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass       ), intent(in   ), target :: excursionSetBarrier_
    !# <constructorAssign variables="timeStepFractional, fileName, *cosmologyFunctions_, *excursionSetBarrier_"/>

    self%tableInitialized                  =.false.
    self%interpolationResetTime            =.true.
    self%interpolationResetVariance        =.true.
    self%tableInitializedRate              =.false.
    self%interpolationResetTimeRate        =.true.
    self%interpolationResetVarianceRate    =.true.
    self%interpolationResetVarianceRateBase=.true.
    self%timeMaximum                       =-huge(0.0d0)
    self%timeMinimum                       =+huge(0.0d0)
    self%varianceMaximum                   =-huge(0.0d0)
    self%timeMaximumRate                   =-huge(0.0d0)
    self%timeMinimumRate                   =+huge(0.0d0)
    self%varianceMaximumRate               =-huge(0.0d0)
    self%timeRatePrevious                  =-huge(0.0d0)
    self%varianceRatePrevious              =-huge(0.0d0)
    self%useFile                           =(self%fileName /= 'none')
    ! Build an automatic file name based on the descriptor for this object.
    if (self%fileName == "auto") then
       self%fileName=galacticusPath(pathTypeDataDynamic)//'largeScaleStructure/excursionSets/firstCrossDistributionFarahi_'//self%hashedDescriptor(includeSourceDigest=.true.)//'.hdf5'
       call Galacticus_Display_Message('excursion set data will be read from/written to "'//char(self%fileName)//'"',verbosityWorking)
    end if
    ! Expand file name.
    self%fileName=File_Name_Expand(char(self%fileName))
    ! Ensure directory exists.
    call System_Command_Do("mkdir -p `dirname "//char(self%fileName)//"`")
    ! Initialize file lock.    
    if (.not.farahiFileLockInitialized) then
       !$omp critical(farahiFileLockInitialize)
       if (.not.farahiFileLockInitialized) then
          call File_Lock_Initialize(farahiFileLock)
          farahiFileLockInitialized=.true.
       end if
       !$omp end critical(farahiFileLockInitialize)
    end if
    return
  end function farahiConstructorInternal

  subroutine farahiDestructor(self)
    !% Destructor for the Farahi excursion set first crossing class.    
    implicit none
    type(excursionSetFirstCrossingFarahi), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"  />
    !# <objectDestructor name="self%excursionSetBarrier_" />
    return
  end subroutine farahiDestructor

  double precision function farahiProbability(self,variance,time,node)
    !% Return the excursion set barrier at the given variance and time.
    use Numerical_Ranges
    use Numerical_Interpolation
    use Memory_Management
    use Galacticus_Display
    use Kind_Numbers
    use Error_Functions
    use MPI_Utilities
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)  :: self
    double precision                                 , intent(in   )  :: variance                     , time
    type            (treeNode                       ), intent(inout)  :: node
    double precision                                 , dimension(0:1) :: hTime                        , hVariance
    double precision                                 , parameter      :: varianceTableTolerance=1.0d-6
    class           (excursionSetBarrierClass       ), pointer        :: excursionSetBarrier_
    logical                                                           :: makeTable
    integer         (c_size_t                       )                 :: iTime                        , iVariance     , &
         &                                                               loopCount                    , loopCountTotal, &
         &                                                               i                            , j             , &
         &                                                               jTime                        , jVariance
    double precision                                                  :: sigma1f
    character       (len=6                          )                 :: label
    type            (varying_string                 )                 :: message
    logical                                                           :: locked

    ! Read tables from file if possible.
    locked=.false.
    if (self%useFile.and..not.self%tableInitialized) then
       call File_Lock(char(self%fileName),farahiFileLock)
       locked=.true.
       call self%fileRead()
    end if
    ! Construct the table if necessary.
    makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
#ifdef USEMPI
    if (self%coordinatedMPI_) then
       if (locked) then
          call File_Unlock(farahiFileLock)
          locked=.false.
       end if
       call mpiBarrier()
    end if
#endif
    if (makeTable) then
#ifdef USEMPI
       ! If coordinating under MPI then only the rank-0 process locks the file.
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          if (self%useFile.and..not.locked) then
             call File_Lock(char(self%fileName),farahiFileLock)
             locked=.true.
          end if
#ifdef USEMPI
       end if
#endif
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
       if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
       self%varianceMaximum   =max(self%varianceMaximum,variance)
       self%varianceTableCount=int(self%varianceMaximum*dble(farahiVarianceNumberPerUnitProbability))
       if (self%tableInitialized) then
          self%timeMinimum=min(      self%timeMinimum                                          ,0.5d0*time)
          self%timeMaximum=max(      self%timeMaximum                                          ,2.0d0*time)
       else
          self%timeMinimum=                                                                     0.5d0*time
          self%timeMaximum=max(2.0d0*self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0),2.0d0*time)
       end if
       self%timeTableCount=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(farahiTimeNumberPerDecade))+1)
       call allocateArray(self%varianceTable                ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
       call allocateArray(self%timeTable                    ,[                          self%timeTableCount]                  )
       call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount,self%timeTableCount],lowerBounds=[0,1])
       self%timeTable        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%timeTableCount      ,rangeType=rangeTypeLogarithmic)
       self%varianceTable    =Make_Range(0.0d0           ,self%varianceMaximum,self%varianceTableCount+1,rangeType=rangeTypeLinear     )
       self%varianceTableStep=self%varianceTable(1)-self%varianceTable(0)
       ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
          message="    time: "
          write (label,'(f6.3)') self%timeMinimum
          message=message//label//" to "
          write (label,'(f6.3)') self%timeMaximum
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
          message="variance: "
          write (label,'(f6.3)') self%varianceMaximum
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
#ifdef USEMPI
       end if
#endif
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
          loopCountTotal=(int(self%timeTableCount,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t)*(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
       else
#endif
          loopCountTotal= int(self%timeTableCount,kind=c_size_t)                                               *(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
#ifdef USEMPI
       end if
#endif
       loopCount=0
#ifdef USEMPI
       if (self%coordinatedMPI_) self%firstCrossingProbabilityTable=0.0d0
#endif
       !$omp parallel private(iTime,i,j,sigma1f,excursionSetBarrier_) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
       allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
       !# <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
       !$omp do schedule(dynamic)
       do iTime=1,self%timeTableCount
#ifdef USEMPI
          if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          self%firstCrossingProbabilityTable(0,iTime)=0.0d0
          self%firstCrossingProbabilityTable(1,iTime)=                                                                                                                                          &
               &                                 real(                                                                                                                                          &
               &                                        +2.0_kind_quad                                                                                                                          &
               &                                      *(                                                                                                                                        &
               &                                        +1.0_kind_quad                                                                                                                          &
               &                                        -erfApproximate(                                                                                                                        &
               &                                                        +excursionSetBarrier_%barrier(                   self%varianceTable(1),self%timeTable(iTime),node,rateCompute=.false.)  &
               &                                                        /                             sqrt(2.0_kind_quad*self%varianceTable(1)                                               )  &
               &                                                       )                                                                                                                        &
               &                                       )                                                                                                                                        &
               &                                      /self%varianceTableStep                                                                                                                 , &
               &                                      kind=kind_dble                                                                                                                            &
               &                                     )
          do i=2,self%varianceTableCount
#ifdef USEMPI
             if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityWorking)
#ifdef USEMPI
             end if
#endif
             !$omp atomic
             loopCount=loopCount+(i-1)
             sigma1f  =0.0d0
             do j=1,i-1
                sigma1f=+sigma1f                                                                                                                     &
                     &  +self%firstCrossingProbabilityTable(j,iTime)                                                                                 &
                     &  *real(                                                                                                                       &
                     &        +1.0_kind_quad                                                                                                         &
                     &        -erfApproximate(                                                                                                       &
                     &                        +(                                                                                                     &
                     &                          +excursionSetBarrier_%barrier(self%varianceTable(i),self%timeTable(iTime),node,rateCompute=.false.)  &
                     &                          -excursionSetBarrier_%barrier(self%varianceTable(j),self%timeTable(iTime),node,rateCompute=.false.)  &
                     &                         )                                                                                                     &
                     &                        /sqrt(2.0_kind_quad*(self%varianceTable(i)-self%varianceTable(j)))                                     &
                     &                       )                                                                                                     , &
                     &        kind=kind_dble                                                                                                         &
                     &       )
             end do
             self%firstCrossingProbabilityTable(i,iTime)=+real(                                                                                                                                               &
                  &                                            +max(                                                                                                                                          &
                  &                                                 +0.0_kind_quad,                                                                                                                           &
                  &                                                 +2.0_kind_quad                                                                                                                            &
                  &                                                 *(                                                                                                                                        &
                  &                                                   +1.0_kind_quad                                                                                                                          &
                  &                                                   -erfApproximate(                                                                                                                        &
                  &                                                                   +excursionSetBarrier_%barrier(                   self%varianceTable(i),self%timeTable(iTime),node,rateCompute=.false.)  &
                  &                                                                   /                             sqrt(2.0_kind_quad*self%varianceTable(i)                                               )  &
                  &                                                                  )                                                                                                                        &
                  &                                                  )                                                                                                                                        &
                  &                                                 /self%varianceTableStep                                                                                                                   &
                  &                                                 -2.0_kind_quad                                                                                                                            &
                  &                                                 *sigma1f                                                                                                                                  &
                  &                                                )                                                                                                                                        , &
                  &                                            kind=kind_dble                                                                                                                                 &
                  &                                           )
          end do
          ! Force the probability at maximum variance to zero.
          self%firstCrossingProbabilityTable(self%varianceTableCount,iTime)=0.0d0
       end do
       !$omp end do
       !# <objectDestructor name="excursionSetBarrier_"/>
       !$omp end parallel
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          call Galacticus_Display_Counter_Clear(verbosityWorking)
          call Galacticus_Display_Unindent("done",verbosityWorking)
#ifdef USEMPI
       end if
       if (self%coordinatedMPI_) then
          call mpiBarrier()
          self%firstCrossingProbabilityTable=mpiSelf%sum(self%firstCrossingProbabilityTable)
       end if
#endif
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVariance,reset=self%interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorTime    ,reset=self%interpolationResetTime    )
       self%interpolationResetVariance=.true.
       self%interpolationResetTime    =.true.
       ! Record that the table is now built.
       self%tableInitialized=.true.
       ! Write the table to file if possible.
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          if (self%useFile) call self%fileWrite()
#ifdef USEMPI
       end if
#endif
    end if
    if (locked) call File_Unlock(farahiFileLock)
    ! Get interpolation in time.
    iTime    =Interpolate_Locate                 (self%timeTable    ,self%interpolationAcceleratorTime    ,time    ,reset=self%interpolationResetTime    )
    hTime    =Interpolate_Linear_Generate_Factors(self%timeTable    ,iTime    ,time    )
    ! Get interpolation in variance.
    iVariance=Interpolate_Locate                 (self%varianceTable,self%interpolationAcceleratorVariance,variance,reset=self%interpolationResetVariance)
    hVariance=Interpolate_Linear_Generate_Factors(self%varianceTable,iVariance,variance)
    ! Compute first crossing probability by interpolating.
    farahiProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiProbability=+farahiProbability                                                     &
               &            +hTime                             (                            jTime) &
               &            *hVariance                         (            jVariance            ) &
               &            *self%firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function farahiProbability

  double precision function farahiRate(self,variance,varianceProgenitor,time,node)
    !% Return the excursion set barrier at the given variance and time.
    use Numerical_Interpolation
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)  :: self
    double precision                                 , intent(in   )  :: variance           , varianceProgenitor, &
         &                                                               time
    type            (treeNode                       ), intent(inout)  :: node
    double precision                                 , dimension(0:1) :: hVarianceProgenitor
    integer                                                           :: jVarianceProgenitor, jTime             , &
         &                                                               jVariance
    integer         (c_size_t                       )                 :: iVarianceProgenitor

    ! For progenitor variances less than or equal to the original variance, return zero.
    if (varianceProgenitor <= variance) then
       farahiRate=0.0d0
       return
    end if
    ! Ensure that the rate is tabulated.
    call self%rateTabulate(varianceProgenitor,time,node)
    ! For progenitor variances greater than the maximum allowed variance, return zero.
    if (varianceProgenitor > self%varianceMaximumRate) then
       farahiRate=0.0d0
       return
    end if
    ! Get interpolation in time.
    if (time /= self%timeRatePrevious) then
       self%timeRatePrevious    =time
       self%iTimeRate           =Interpolate_Locate                 (self%timeTableRate        ,self%interpolationAcceleratorTimeRate        ,time                       ,reset=self%interpolationResetTimeRate        )
       self%hTimeRate           =Interpolate_Linear_Generate_Factors(self%timeTableRate        ,self%iTimeRate     ,time    )
    end if
    ! Get interpolation in variance.
    if (variance /= self%varianceRatePrevious) then
       self%varianceRatePrevious=variance
       self%iVarianceRate       =Interpolate_Locate                 (self%varianceTableRateBase,self%interpolationAcceleratorVarianceRateBase,                   variance,reset=self%interpolationResetVarianceRateBase)
       self%hVarianceRate       =Interpolate_Linear_Generate_Factors(self%varianceTableRateBase,self%iVarianceRate,variance)
    end if

    ! Get interpolation in progenitor variance.
    iVarianceProgenitor         =Interpolate_Locate                 (self%varianceTableRate    ,self%interpolationAcceleratorVarianceRate    ,varianceProgenitor-variance,reset=self%interpolationResetVarianceRate    )
    ! Catch cases where the maximum variance is approached.
    if (self%varianceTableRate(iVarianceProgenitor)+variance > self%varianceMaximumRate) then
       ! Force the rate to drop to zero at the maximum variance. (Necessary because we will not have a tabulated point precisely
       ! at the maximum variance.)
       hVarianceProgenitor=[                                                                                      &
            &               +1.0d0                                                                                &
            &               -((     varianceProgenitor -variance)-self%varianceTableRate(iVarianceProgenitor-1))  &
            &               /((self%varianceMaximumRate-variance)-self%varianceTableRate(iVarianceProgenitor-1)), &
            &               +0.0d0                                                                                &
            &              ]
    else
       hVarianceProgenitor=Interpolate_Linear_Generate_Factors(self%varianceTableRate,iVarianceProgenitor,varianceProgenitor-variance)
    end if
    ! Compute first crossing probability by interpolating.
    farahiRate=0.0d0
    do jTime=0,1
       do jVariance=0,1
          do jVarianceProgenitor=0,1
             farahiRate=+farahiRate                                                                                                                 &
                  &     +self%hTimeRate             (                                                                                        jTime) &
                  &     *self%hVarianceRate         (                                                               jVariance                     ) &
                  &     *hVarianceProgenitor        (                      jVarianceProgenitor                                                    ) &
                  &     *self%firstCrossingTableRate(iVarianceProgenitor-1+jVarianceProgenitor,self%iVarianceRate-1+jVariance,self%iTimeRate+jTime)
          end do
       end do
    end do
    return
  end function farahiRate

  double precision function farahiRateNonCrossing(self,variance,time,node)
   !% Return the rate for excursion set non-crossing.
    use Numerical_Interpolation
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: time , variance
    type            (treeNode                       ), intent(inout) :: node
    integer                                                          :: jTime, jVariance

    ! Ensure that the rate is tabulated.
    call self%rateTabulate(variance,time,node)
    ! Get interpolation in time.
    if (time /= self%timeRatePrevious) then
       self%timeRatePrevious    =time
       self%iTimeRate           =Interpolate_Locate                 (self%timeTableRate        ,self%interpolationAcceleratorTimeRate        ,time    ,reset=self%interpolationResetTimeRate        )
       self%hTimeRate           =Interpolate_Linear_Generate_Factors(self%timeTableRate        ,self%iTimeRate    ,time    )
    end if
    ! Get interpolation in variance.
    if (variance /= self%varianceRatePrevious) then
       self%varianceRatePrevious=variance
       self%iVarianceRate       =Interpolate_Locate                 (self%varianceTableRateBase,self%interpolationAcceleratorVarianceRateBase,variance,reset=self%interpolationResetVarianceRateBase)
       self%hVarianceRate       =Interpolate_Linear_Generate_Factors(self%varianceTableRateBase,self%iVarianceRate,variance)
    end if
    ! Compute non-crossing probability by interpolating.
    farahiRateNonCrossing=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiRateNonCrossing=+farahiRateNonCrossing                                                          &
               &                +self%hTimeRate           (                                              jTime) &
               &                *self%hVarianceRate       (                     jVariance                     ) &
               &                *self%nonCrossingTableRate(self%iVarianceRate-1+jVariance,self%iTimeRate+jTime)
       end do
    end do
    return
  end function farahiRateNonCrossing

  subroutine farahiRateTabulate(self,varianceProgenitor,time,node)
    !% Tabulate the excursion set crossing rate.
    use Numerical_Ranges
    use Numerical_Interpolation
    use Memory_Management
    use Galacticus_Display
    use Kind_Numbers
    use Error_Functions
    use MPI_Utilities
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)               :: self
    double precision                                 , intent(in   )               :: time                             , varianceProgenitor
    type            (treeNode                       ), intent(inout)               :: node
    double precision                                 , parameter                   :: varianceMinimumDefault    =1.0d-2
    double precision                                 , parameter                   :: varianceTolerance         =1.0d-6
    real            (kind=kind_quad                 ), allocatable  , dimension(:) :: firstCrossingTableRateQuad       , varianceTableRateBaseQuad, &
         &                                                                            varianceTableRateQuad
    class           (excursionSetBarrierClass       ), pointer                     :: excursionSetBarrier_
#ifdef USEMPI
    integer                                                                        :: taskCount
#endif
    logical                                                                        :: makeTable
    integer         (c_size_t                       )                              :: loopCount                        , loopCountTotal
    integer                                                                        :: i                                , iTime                    , &
         &                                                                            iVariance                        , j
    double precision                                                               :: timeProgenitor                   , varianceMinimumRate
    character       (len=6                          )                              :: label
    type            (varying_string                 )                              :: message
    real            (kind=kind_quad                 )                              :: crossingFraction                 , effectiveBarrierInitial  , &
         &                                                                            sigma1f                          , varianceTableStepRate    , &
         &                                                                            barrier
    logical                                                                        :: locked
    
    ! Determine if we need to make the table.
    ! Read tables from file if possible.
    locked=.false.
    if (self%useFile.and..not.self%tableInitializedRate) then
       call File_Lock(char(self%fileName),farahiFileLock)
       locked=.true.
       call self%fileRead()
    end if
    makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
#ifdef USEMPI
    if (self%coordinatedMPI_) then
       if (locked) then
          call File_Unlock(farahiFileLock)
          locked=.false.
       end if
       call mpiBarrier()
    end if
#endif
    if (makeTable) then
#ifdef USEMPI
       ! If coordinating under MPI then only the rank-0 process locks the file.
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (self%useFile.and..not.locked) then
             call File_Lock(char(self%fileName),farahiFileLock)
             locked=.true.
          end if
#ifdef USEMPI
       end if
#endif
       if (allocated(self%varianceTableRate     )) call deallocateArray(self%varianceTableRate     )
       if (allocated(self%varianceTableRateBase )) call deallocateArray(self%varianceTableRateBase )
       if (allocated(self%timeTableRate         )) call deallocateArray(self%timeTableRate         )
       if (allocated(self%firstCrossingTableRate)) call deallocateArray(self%firstCrossingTableRate)
       if (allocated(self%nonCrossingTableRate  )) call deallocateArray(self%nonCrossingTableRate  )
       if (self%tableInitializedRate) then
          self%timeMinimumRate   =min(self%timeMinimumRate,0.5d0*time)
          self%timeMaximumRate   =max(self%timeMaximumRate,2.0d0*time)
          self%timeTableCountRate=int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(farahiTimeNumberPerDecade))+1
       else
          self%timeMinimumRate   =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(farahiRateRedshiftMaximum))
          self%timeMaximumRate   =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(farahiRateRedshiftMinimum))
          self%timeMinimumRate   =min(self%timeMinimumRate,0.5d0*time)
          self%timeMaximumRate   =max(self%timeMaximumRate,2.0d0*time)
          self%timeTableCountRate=max(int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(farahiTimeNumberPerDecade))+1,2)
       end if
       ! Set the default minimum variance.
       varianceMinimumRate       =varianceMinimumDefault
       ! Next reduce the variance if necessary such that the typical amplitude of fluctuations is less (by a factor of sqrt[10])
       ! than the effective barrier height at zero variance for the minimum and maximum times that we must consider.
       allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
       !# <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
       varianceMinimumRate            =min(                                                                                                                      &
            &                              +varianceMinimumRate                                                                                                , &
            &                              +1.0d-2                                                                                                               &
            &                              *(                                                                                                                    &
            &                                +excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate*(1.0d0-self%timeStepFractional),node,rateCompute=.true.)  &
            &                                -excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate                                ,node,rateCompute=.true.)  &
            &                               )**2                                                                                                                 &
            &                             )
       !# <objectDestructor name="excursionSetBarrier_"/>
       self%varianceMaximumRate       =max(self%varianceMaximumRate,varianceProgenitor)
       self%varianceTableCountRate    =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(farahiVarianceNumberPerDecade))+1
       self%varianceTableCountRateBase=int(self%varianceMaximumRate*dble(farahiVarianceNumberPerUnit))
       call allocateArray(self%varianceTableRate     ,[1+self%varianceTableCountRate                                                          ],lowerBounds=[0    ])
       call allocateArray(self%varianceTableRateBase ,[                              1+self%varianceTableCountRateBase                        ],lowerBounds=[0    ])
       call allocateArray(self%timeTableRate         ,[                                                                self%timeTableCountRate]                    )
       call allocateArray(self%firstCrossingTableRate,[1+self%varianceTableCountRate,1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[0,0,1])
       call allocateArray(self%nonCrossingTableRate  ,[                              1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[  0,1])
       ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
       self%varianceTableRate    (0                                )=0.0d0
       self%varianceTableRate    (1:self%varianceTableCountRate    )=self%varianceRange(varianceMinimumRate,self%varianceMaximumRate,self%varianceTableCountRate      ,ratioAtMaximum=10.0d0    )
       self%varianceTableRateBase(0:self%varianceTableCountRateBase)=Make_Range        (0.0d0              ,self%varianceMaximumRate,self%varianceTableCountRateBase+1,rangeType=rangeTypeLinear)
       ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
       allocate(varianceTableRateQuad     (0:self%varianceTableCountRate    ))
       varianceTableRateQuad    =self%varianceTableRate
       allocate(varianceTableRateBaseQuad (0:self%varianceTableCountRateBase))
       varianceTableRateBaseQuad=self%varianceTableRateBase
       ! The time table is logarithmically distributed in time.
       self%timeTableRate=Make_Range(self%timeMinimumRate,self%timeMaximumRate,self%timeTableCountRate,rangeType=rangeTypeLogarithmic)
       ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          call Galacticus_Display_Indent("solving for excursion set barrier crossing rates",verbosityWorking)
          message="    time: "
          write (label,'(f6.3)') self%timeMinimumRate
          message=message//label//" to "
          write (label,'(f6.3)') self%timeMaximumRate
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
          message="variance: "
          write (label,'(f6.3)') self%varianceMaximumRate
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
#ifdef USEMPI
       end if
#endif
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
          loopCountTotal=int(self%timeTableCountRate,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t
       else
#endif
          loopCountTotal=int(self%timeTableCountRate,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)
#ifdef USEMPI
       end if
#endif
       loopCount=0
#ifdef USEMPI
       if (self%coordinatedMPI_) then
          self%firstCrossingTableRate=0.0d0
          self%nonCrossingTableRate  =0.0d0
       end if
       taskCount=-1
#endif
       !$omp parallel private(iTime,timeProgenitor,iVariance,varianceTableStepRate,i,j,sigma1f,crossingFraction,barrier,effectiveBarrierInitial,firstCrossingTableRateQuad,excursionSetBarrier_) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
       allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
       !# <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
       !$omp do schedule(dynamic)
       do iTime=1,self%timeTableCountRate
          if (.not.allocated(firstCrossingTableRateQuad)) allocate(firstCrossingTableRateQuad(0:self%varianceTableCountRate))
          ! Compute a suitable progenitor time.
          timeProgenitor=self%timeTableRate(iTime)*(1.0d0-self%timeStepFractional)
          ! Loop through the starting variances.
          do iVariance=0,self%varianceTableCountRateBase
#ifdef USEMPI
             taskCount=taskCount+1
             if (self%coordinatedMPI_ .and. mod(taskCount,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
#ifdef USEMPI
             if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityWorking)
#ifdef USEMPI
             end if
#endif
             !$omp atomic
             loopCount=loopCount+1_c_size_t
             ! For zero variance, the rate is initialized to zero.
             firstCrossingTableRateQuad(0)=0.0d0
             ! Compute the step in variance across this first grid cell.
             varianceTableStepRate=varianceTableRateQuad(1)-varianceTableRateQuad(0)
             ! Compute the barrier for the descendent.
             barrier=real(excursionSetBarrier_%barrier(real(varianceTableRateBaseQuad(iVariance),kind=8),self%timeTableRate(iTime),node,rateCompute=.true.),kind=kind_quad)
             ! Compute the first crossing distribution at the first grid cell.
             if (varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                firstCrossingTableRateQuad(1)= 0.0d0
             else
                firstCrossingTableRateQuad(1)=+2.0_kind_quad                                                                                                                                                                             &
                     &                        *(                                                                                                                                                                                         &
                     &                          +1.0_kind_quad                                                                                                                                                                           &
                     &                          -erfApproximate(                                                                                                                                                                         &
                     &                                          +(                                                                                                                                                                       &
                     &                                            +real(excursionSetBarrier_%barrier(real(+varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance),kind=8),timeProgenitor,node,rateCompute=.true.),kind=kind_quad) &
                     &                                            -barrier                                                                                                                                                               &
                     &                                           )                                                                                                                                                                       &
                     &                                          /sqrt(2.0_kind_quad*varianceTableRateQuad(1))                                                                                                                            &
                     &                                         )                                                                                                                                                                         &
                     &                         )                                                                                                                                                                                         &
                     &                        /varianceTableStepRate
             end if
             do i=2,self%varianceTableCountRate
                if (varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                   firstCrossingTableRateQuad(i)=0.0d0
                else
                   effectiveBarrierInitial=+real(excursionSetBarrier_%barrier(real(+varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance),kind=8),timeProgenitor,node,rateCompute=.true.),kind=kind_quad) &
                        &                  -barrier
                   sigma1f                =+0.0d0
                   do j=1,i-1
                      varianceTableStepRate=(varianceTableRateQuad(j+1)-varianceTableRateQuad(j-1))/2.0_kind_quad
                      sigma1f=+sigma1f                                                                                                                                                                                   &
                           &  +firstCrossingTableRateQuad(j)                                                                                                                                                             &
                           &  *varianceTableStepRate                                                                                                                                                                     &
                           &  *(                                                                                                                                                                                         &
                           &    +1.0_kind_quad                                                                                                                                                                           &
                           &    -erfApproximate(                                                                                                                                                                         &
                           &                    +(                                                                                                                                                                       &
                           &                      +effectiveBarrierInitial                                                                                                                                               &
                           &                      -real(excursionSetBarrier_%barrier(real(+varianceTableRateQuad(j)+varianceTableRateBaseQuad(iVariance),kind=8),timeProgenitor,node,rateCompute=.true.),kind=kind_quad) &
                           &                      +barrier                                                                                                                                                               &
                           &                     )                                                                                                                                                                       &
                           &                    /sqrt(2.0_kind_quad*(varianceTableRateQuad(i)-varianceTableRateQuad(j)))                                                                                                 &
                           &                   )                                                                                                                                                                         &
                           &   )
                   end do
                   varianceTableStepRate=varianceTableRateQuad(i)-varianceTableRateQuad(i-1)
                   firstCrossingTableRateQuad(i)=max(                                                                  &
                        &                            +0.0_kind_quad,                                                   &
                        &                            +(                                                                &
                        &                              +2.0_kind_quad                                                  &
                        &                              *(                                                              &
                        &                                +1.0_kind_quad                                                &
                        &                                -erfApproximate(                                              &
                        &                                                +effectiveBarrierInitial                      &
                        &                                                /sqrt(2.0_kind_quad*varianceTableRateQuad(i)) &
                        &                                               )                                              &
                        &                               )                                                              &
                        &                              -2.0_kind_quad*sigma1f                                          &
                        &                             )                                                                &
                        &                            /varianceTableStepRate                                            &
                        &                           )
                end if
             end do
             ! Compute the fraction of trajectories which never cross the barrier.
             crossingFraction=0.0_kind_quad
             do j=0,self%varianceTableCountRate-1
                varianceTableStepRate=varianceTableRateQuad(j+1)-varianceTableRateQuad(j)
                crossingFraction=+crossingFraction                  &
                     &           +0.5_kind_quad                     &
                     &           *(                                 &
                     &              firstCrossingTableRateQuad(j  ) &
                     &             +firstCrossingTableRateQuad(j+1) &
                     &            )                                 &
                     &           *varianceTableStepRate
             end do
             ! Compute the rate for trajectories which never cross the barrier.
             self%nonCrossingTableRate(iVariance,iTime)=real(                                   &
                  &                                          +(1.0_kind_quad-crossingFraction)  &
                  &                                          /self%timeTableRate(iTime)         &
                  &                                          /self%timeStepFractional         , &
                  &                                          kind=kind_dble                     &
                  &                                         )
             ! Store the compute crossing rate in our table.
             self%firstCrossingTableRate(:,iVariance,iTime)=real(firstCrossingTableRateQuad,kind=kind_dble)

          end do
          ! Divide through by the time step to get the rate of barrier crossing.
          self%firstCrossingTableRate(:,:,iTime)=+self%firstCrossingTableRate(:,:,iTime) &
               &                                 /self%timeTableRate         (    iTime) &
               &                                 /self%timeStepFractional
       end do
       !$omp end do
       !# <objectDestructor name="excursionSetBarrier_"/>
       !$omp end parallel
       ! Deallocate work arrays.
       deallocate(varianceTableRateBaseQuad )
       deallocate(varianceTableRateQuad     )
       if (allocated(firstCrossingTableRateQuad)) deallocate(firstCrossingTableRateQuad)
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          call Galacticus_Display_Counter_Clear(       verbosityWorking)
          call Galacticus_Display_Unindent     ("done",verbosityWorking)
#ifdef USEMPI
       end if
       if (self%coordinatedMPI_) then
          call mpiBarrier()
          self%firstCrossingTableRate=mpiSelf%sum(self%firstCrossingTableRate)
          self%  nonCrossingTableRate=mpiSelf%sum(self%  nonCrossingTableRate)
       end if
#endif
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVarianceRate    ,reset=self%interpolationResetVarianceRate    )
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVarianceRateBase,reset=self%interpolationResetVarianceRateBase)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorTimeRate        ,reset=self%interpolationResetTimeRate        )
       self%interpolationResetVarianceRate    =.true.
       self%interpolationResetVarianceRateBase=.true.
       self%interpolationResetTimeRate        =.true.
       ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
       self%varianceRatePrevious=-1.0d0
       self%timeRatePrevious    =-1.0d0
       ! Record that the table is now built.
       self%tableInitializedRate=.true.
       ! Write the table to file if possible.
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          if (self%useFile) call self%fileWrite()
#ifdef USEMPI
       end if
#endif
    end if
    if (locked) call File_Unlock(farahiFileLock)
    return
  end subroutine farahiRateTabulate

  subroutine farahiFileRead(self)
    !% Read tabulated data on excursion set first crossing probabilities from file.
    use IO_HDF5
    use File_Utilities
    use Memory_Management
    use Numerical_Interpolation
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                   :: self
    type            (hdf5Object                     )                                  :: dataFile                   , dataGroup
    double precision                                 , allocatable  , dimension(:    ) :: varianceTableBaseTmp       , varianceTableTmp
    double precision                                 , allocatable  , dimension(:,:  ) :: firstCrossingProbabilityTmp, nonCrossingTableRate
    double precision                                 , allocatable  , dimension(:,:,:) :: firstCrossingRateTmp
    type            (varying_string                 )                                  :: message
    character       (len=32                         )                                  :: label
    
    ! Return immediately if the file does not exist.
    if (.not.File_Exists(File_Name_Expand(char(self%fileName)))) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(char(self%fileName))
    ! Check if the standard table is populated.
    if (dataFile%hasGroup('probability')) then
       ! Deallocate arrays if necessary.
       if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
       if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
       ! Read the datasets.
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%readDataset('variance'                ,varianceTableTmp           )
       call dataGroup%readDataset('time'                    ,self%timeTable             )
       call dataGroup%readDataset('firstCrossingProbability',firstCrossingProbabilityTmp)
       call dataGroup%close()
       ! Set table sizes and limits.
       self%varianceTableCount=size(varianceTableTmp)-1
       self%timeTableCount    =size(self%timeTable  )
       ! Transfer to tables.
       call allocateArray(self%varianceTable                ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
       call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount,self%timeTableCount],lowerBounds=[0,1])
       self%varianceTable                (0:self%varianceTableCount  )=varianceTableTmp           (1:self%varianceTableCount+1  )
       self%firstCrossingProbabilityTable(0:self%varianceTableCount,:)=firstCrossingProbabilityTmp(1:self%varianceTableCount+1,:)
       call deallocateArray(varianceTableTmp           )
       call deallocateArray(firstCrossingProbabilityTmp)
       ! Set table limits.
       self%timeMinimum      =self%timeTable    (                      1)
       self%timeMaximum      =self%timeTable    (self%    timeTableCount)
       self%varianceMaximum  =self%varianceTable(self%varianceTableCount)
       self%varianceTableStep=self%varianceTable(1)-self%varianceTable(0)
       self%tableInitialized =.true.
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVariance,reset=self%interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorTime    ,reset=self%interpolationResetTime    )
       self%interpolationResetVariance=.true.
       self%interpolationResetTime    =.true.
       ! Report.
       message=var_str('read excursion set first crossing probability from: ')//File_Name_Expand(char(self%fileName))
       call Galacticus_Display_Indent  (message,verbosityWorking)
       write (label,'(e12.6)') self%timeMinimum
       message=var_str('    time minimum: ')//label//' Gyr'
       call Galacticus_Display_Message (message,verbosityWorking)
       write (label,'(e12.6)') self%timeMaximum
       message=var_str('    time maximum: ')//label//' Gyr'
       write (label,'(e12.6)') self%varianceMaximum
       message=var_str('variance maximum: ')//label
       call Galacticus_Display_Message (message,verbosityWorking)
       call Galacticus_Display_Unindent(''     ,verbosityWorking)
    end if
    ! Check if the rate table is populated.
    if (dataFile%hasGroup('rate')) then
       ! Deallocate arrays if necessary.
       if (allocated(self%varianceTableRate     )) call deallocateArray(self%varianceTableRate     )
       if (allocated(self%varianceTableRateBase )) call deallocateArray(self%varianceTableRateBase )
       if (allocated(self%timeTableRate         )) call deallocateArray(self%timeTableRate         )
       if (allocated(self%firstCrossingTableRate)) call deallocateArray(self%firstCrossingTableRate)
       if (allocated(self%nonCrossingTableRate  )) call deallocateArray(self%nonCrossingTableRate  )
       ! Read the datasets.
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%readDataset('variance'         ,varianceTableTmp    )
       call dataGroup%readDataset('varianceBase'     ,varianceTableBaseTmp)
       call dataGroup%readDataset('time'             ,self%timeTableRate  )
       call dataGroup%readDataset('firstCrossingRate',firstCrossingRateTmp)
       call dataGroup%readDataset('nonCrossingRate'  ,nonCrossingTableRate)
       call dataGroup%close()
       ! Set table sizes and limits.
       self%varianceTableCountRate    =size(varianceTableTmp    )-1
       self%varianceTableCountRateBase=size(varianceTableBaseTmp)-1
       self%timeTableCountRate        =size(self%timeTableRate  )
       ! Transfer to tablse.
       call allocateArray(self%varianceTableRate     ,[1+self%varianceTableCountRate                                                          ],lowerBounds=[0    ])
       call allocateArray(self%varianceTableRateBase ,[                              1+self%varianceTableCountRateBase                        ],lowerBounds=[  0  ])
       call allocateArray(self%firstCrossingTableRate,[1+self%varianceTableCountRate,1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[0,0,1])
       call allocateArray(self%nonCrossingTableRate  ,[                              1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[  0,1])
       self%varianceTableRate     (0:self%varianceTableCountRate                                    )=varianceTableTmp    (1:self%varianceTableCountRate+1                                      )
       self%varianceTableRateBase (                              0:self%varianceTableCountRateBase  )=varianceTableBaseTmp(                                1:self%varianceTableCountRateBase+1  )
       self%firstCrossingTableRate(0:self%varianceTableCountRate,0:self%varianceTableCountRateBase,:)=firstCrossingRateTmp(1:self%varianceTableCountRate+1,1:self%varianceTableCountRateBase+1,:)
       self%nonCrossingTableRate  (                              0:self%varianceTableCountRateBase,:)=nonCrossingTableRate(                                1:self%varianceTableCountRateBase+1,:)
       call deallocateArray(varianceTableTmp    )
       call deallocateArray(varianceTableBaseTmp)
       ! Set table limits.
       self%varianceMaximumRate =self%varianceTableRate(self%varianceTableCountRate)
       self%timeMinimumRate     =self%timeTableRate    (                          1)
       self%timeMaximumRate     =self%timeTableRate    (    self%timeTableCountRate)
       self%tableInitializedRate=.true.
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVarianceRate    ,reset=self%interpolationResetVarianceRate    )
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVarianceRateBase,reset=self%interpolationResetVarianceRateBase)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorTimeRate        ,reset=self%interpolationResetTimeRate        )
       self%interpolationResetVarianceRate    =.true.
       self%interpolationResetVarianceRateBase=.true.
       self%interpolationResetTimeRate        =.true.
       ! Report.
       message=var_str('read excursion set first crossing rates from: ')//File_Name_Expand(char(self%fileName))
       call Galacticus_Display_Indent  (message,verbosityWorking)
       write (label,'(e12.6)') self%timeMinimumRate
       message=var_str('    time minimum: ')//label//' Gyr'
       call Galacticus_Display_Message (message,verbosityWorking)
       write (label,'(e12.6)') self%timeMaximumRate
       message=var_str('    time maximum: ')//label//' Gyr'
       write (label,'(e12.6)') self%varianceMaximumRate
       message=var_str('variance minimum: ')//label
       call Galacticus_Display_Message (message,verbosityWorking)
       call Galacticus_Display_Message (message,verbosityWorking)
       call Galacticus_Display_Unindent(''     ,verbosityWorking)
    end if
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileRead

  subroutine farahiFileWrite(self)
    !% Write tabulated data on excursion set first crossing probabilities to file.
    use HDF5
    use IO_HDF5
    implicit none
    class(excursionSetFirstCrossingFarahi), intent(inout) :: self
    type (hdf5Object                     )                :: dataFile, dataGroup

    ! Don't write anything if neither table is initialized.
    if (.not.(self%tableInitialized.or.self%tableInitializedRate)) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(char(self%fileName),overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    ! Check if the standard table is populated.
    if (self%tableInitialized) then
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%writeDataset(self%varianceTable                ,'variance'                ,'The variance at which results are tabulated.'                         )
       call dataGroup%writeDataset(self%timeTable                    ,'time'                    ,'The cosmic times at which results are tabulated.'                     )
       call dataGroup%writeDataset(self%firstCrossingProbabilityTable,'firstCrossingProbability','The probability of first crossing as a function of variance and time.')
       call dataGroup%close()
    end if
    ! Check if the rate table is populated.
    if (self%tableInitializedRate) then
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%writeDataset(self%varianceTableRate     ,'variance'         ,'The variance at which results are tabulated.'                               )
       call dataGroup%writeDataset(self%varianceTableRateBase ,'varianceBase'     ,'The variance of the base halo at which results are tabulated.'              )
       call dataGroup%writeDataset(self%timeTableRate         ,'time'             ,'The cosmic times at which results are tabulated.'                           )
       call dataGroup%writeDataset(self%firstCrossingTableRate,'firstCrossingRate','The probability rate of first crossing as a function of variances and time.')
       call dataGroup%writeDataset(self%nonCrossingTableRate  ,'nonCrossingRate'  ,'The probability rate of non crossing as a function of variance and time.')
       call dataGroup%close()
    end if
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileWrite

  function farahiVarianceRange(self,rangeMinimum,rangeMaximum,rangeNumber,ratioAtMaximum) result (rangeValues)
    !% Builds a numerical range between {\normalfont \ttfamily rangeMinimum} and {\normalfont \ttfamily rangeMaximum} using {\normalfont \ttfamily rangeNumber} points with spacing that
    !% varies from logarithmic to linear spacing with the transition point controlled by {\normalfont \ttfamily ratioAtMaximum}.
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)          :: self
    double precision                                 , intent(in   )          :: rangeMaximum  , rangeMinimum    , &
         &                                                                       ratioAtMaximum
    integer                                          , intent(in   )          :: rangeNumber
    double precision                                 , dimension(rangeNumber) :: rangeValues
    integer                                                                   :: iRange
    double precision                                                          :: rangeLinear   , rangeLogarithmic
    !GCC$ attributes unused :: self

    do iRange=1,rangeNumber
       rangeLinear        =        rangeMinimum +   (rangeMaximum-rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1)
       rangeLogarithmic   =exp(log(rangeMinimum)+log(rangeMaximum/rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1))
       rangeValues(iRange)=(1.0d0+1.0d0/ratioAtMaximum)/(1.0d0/rangeLinear+1.0d0/rangeLogarithmic/ratioAtMaximum)
    end do
    return
  end function farahiVarianceRange
