!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson.

!% Contains a module which implements a fast and accurate method to
!% solve the excursion set barrier crossing problem for generic
!% barriers.

module Excursion_Sets_First_Crossing_Farahi
  !% Implements a fast and accurate method to solve the excursion set
  !% barrier crossing problem for generic barriers.
  use FGSL
  use ISO_Varying_String
  private
  public :: Excursion_Sets_First_Crossing_Farahi_Initialize

  ! Variables used in tabulation the first crossing function.
  double precision                                                              :: timeMaximum                                      =0.0d0  , timeMinimum                         =0.0d0 , &
       &                                                                           varianceMaximum                                  =0.0d0
  integer                                                                       :: timeTableCount                                           , varianceTableCount
  integer                                                           , parameter :: varianceTableNumberPerUnitProbability            =1000
  integer                                                           , parameter :: timeTableNumberPerDecade                         =10     , varianceTableNumberPerDecade        =400   , &
       &                                                                           varianceTableNumberPerUnit                       =40
  double precision                   , allocatable, dimension(:,:)              :: firstCrossingProbabilityTable
  double precision                   , allocatable, dimension(:  )              :: timeTable                                                , varianceTable
  double precision                                                              :: varianceTableStep
  logical                                                                       :: tableInitialized                                 =.false.
  type            (fgsl_interp_accel)                                           :: interpolationAcceleratorTime                             , interpolationAcceleratorVariance
  logical                                                                       :: interpolationResetTime                           =.true. , interpolationResetVariance          =.true.

  ! Variables used in tabulation the first crossing rate function.
  double precision                                                  , parameter :: redshiftMaximumRate                              =30.0d0 , redshiftMinimumRate                 =0.0d0
  double precision                                                              :: timeMaximumRate                                  =0.0d0  , timeMinimumRate                     =0.0d0 , &
       &                                                                           varianceMaximumRate                              =0.0d0
  integer                                                                       :: timeTableCountRate                                       , varianceTableCountRate                     , &
       &                                                                           varianceTableCountRateBase
  double precision                   , allocatable, dimension(:,:,:)            :: firstCrossingTableRate
  double precision                   , allocatable, dimension(:,:  )            :: nonCrossingTableRate
  double precision                   , allocatable, dimension(:    )            :: timeTableRate                                            , varianceTableRate                          , &
       &                                                                           varianceTableRateBase
  logical                                                                       :: tableInitializedRate                             =.false.
  type            (fgsl_interp_accel)                                           :: interpolationAcceleratorTimeRate                         , interpolationAcceleratorVarianceRate       , &
       &                                                                           interpolationAcceleratorVarianceRateBase
  logical                                                                       :: interpolationResetTimeRate                       =.true. , interpolationResetVarianceRate      =.true., &
       &                                                                           interpolationResetVarianceRateBase               =.true.

  ! File name used to store tabulations.
  type            (varying_string   )                                           :: excursionSetFirstCrossingFarahiFileName

  ! The fractional step in time used to compute barrier crossing rates.
  double precision                                                              :: excursionSetFirstCrossingFarahiFractionalTimeStep

  ! Record of variance and time in previous call to rate functions.
  double precision                                                              :: timeRatePrevious                                         , varianceRatePrevious
  double precision                                , dimension(0:1)              :: hTimeRate                                                , hVarianceRate
  integer                                                                       :: iTimeRate                                                , iVarianceRate

contains

  !# <excursionSetFirstCrossingMethod>
  !#  <unitName>Excursion_Sets_First_Crossing_Farahi_Initialize</unitName>
  !# </excursionSetFirstCrossingMethod>
  subroutine Excursion_Sets_First_Crossing_Farahi_Initialize(excursionSetFirstCrossingMethod &
       &,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get&
       &,Excursion_Sets_Non_Crossing_Rate_Get)
    !% Initialize the ``Farahi'' first crossing distribution method for excursion sets module.
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: excursionSetFirstCrossingMethod
    procedure(double precision), intent(inout), pointer :: Excursion_Sets_First_Crossing_Probability_Get, Excursion_Sets_First_Crossing_Rate_Get, &
         &                                                 Excursion_Sets_Non_Crossing_Rate_Get

    if (excursionSetFirstCrossingMethod == 'Farahi') then
       Excursion_Sets_First_Crossing_Probability_Get => Excursion_Sets_First_Crossing_Probability_Farahi
       Excursion_Sets_First_Crossing_Rate_Get        => Excursion_Sets_First_Crossing_Rate_Farahi
       Excursion_Sets_Non_Crossing_Rate_Get          => Excursion_Sets_Non_Crossing_Rate_Farahi
       ! Get the name of any file to which results should be stored.
       !@ <inputParameter>
       !@   <name>excursionSetFirstCrossingFarahiFileName</name>
       !@   <defaultValue>none</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <type>string</type>
       !@   <cardinality>0..1</cardinality>
       !@   <description>
       !@     The name of the file to/from which tabulations of barrier first crossing probabilities should be written/read. If set to ``{\tt none}'' tables will not be stored.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('excursionSetFirstCrossingFarahiFileName',excursionSetFirstCrossingFarahiFileName,defaultValue='none')
       ! Get the size of the fractional time step to use when computing barrier crossing rates.
       !@ <inputParameter>
       !@   <name>excursionSetFirstCrossingFarahiFractionalTimeStep</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <type>real</type>
       !@   <cardinality>0..1</cardinality>
       !@   <description>
       !@     The fractional time step used when computing barrier crossing rates in the Farahi excursion set solver
       !@     (i.e. the step used in finite difference calculations).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('excursionSetFirstCrossingFarahiFractionalTimeStep',excursionSetFirstCrossingFarahiFractionalTimeStep,defaultValue=1.0d-2)
    end if
    return
  end subroutine Excursion_Sets_First_Crossing_Farahi_Initialize

  double precision function Excursion_Sets_First_Crossing_Probability_Farahi(variance,time)
    !% Return the probability for excursion set first crossing using the methodology of Farahi.
    use Numerical_Ranges
    use Numerical_Interpolation
    use Excursion_Sets_Barriers
    use Memory_Management
    use Galacticus_Display
    use Kind_Numbers
    !$ use OMP_Lib
    implicit none
    double precision                , intent(in   )  :: time                         , variance
    double precision                , dimension(0:1) :: hTime                        , hVariance
    double precision                , parameter      :: varianceTableTolerance=1.0d-6
    logical                                          :: makeTable
    integer                                          :: i                            , iTime         , &
         &                                              iVariance                    , j             , &
         &                                              jTime                        , jVariance     , &
         &                                              loopCount                    , loopCountTotal
    double precision                                 :: sigma1f
    character       (len=6         )                 :: label
    type            (varying_string)                 :: message

    ! Determine if we need to make the table.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    ! Read tables from file if possible.
    if (.not.tableInitialized.and.excursionSetFirstCrossingFarahiFileName /= 'none') call Excursion_Sets_First_Crossing_Farahi_Read_File()
    ! Construct the table if necessary.
    makeTable=.not.tableInitialized.or.(variance > varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < timeMinimum).or.(time > timeMaximum)
    if (makeTable) then
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(varianceTable                )) call Dealloc_Array(varianceTable                )
       if (allocated(timeTable                    )) call Dealloc_Array(timeTable                    )
       if (allocated(firstCrossingProbabilityTable)) call Dealloc_Array(firstCrossingProbabilityTable)
       varianceMaximum   =max(varianceMaximum,variance)
       varianceTableCount=int(varianceMaximum*dble(varianceTableNumberPerUnitProbability))
       if (tableInitialized) then
          timeMinimum=min(timeMinimum,0.5d0*time)
          timeMaximum=max(timeMaximum,2.0d0*time)
       else
          timeMinimum=0.5d0*time
          timeMaximum=2.0d0*time
       end if
       timeTableCount=max(2,int(log10(timeMaximum/timeMinimum)*dble(timeTableNumberPerDecade))+1)
       call Alloc_Array(varianceTable                ,[1+varianceTableCount               ],lowerBounds=[0  ])
       call Alloc_Array(timeTable                    ,[                     timeTableCount]                  )
       call Alloc_Array(firstCrossingProbabilityTable,[1+varianceTableCount,timeTableCount],lowerBounds=[0,1])
       varianceTable    =Make_Range(0.0d0,varianceMaximum,varianceTableCount+1,rangeType=rangeTypeLinear)
       varianceTableStep=varianceTable(1)-varianceTable(0)
       timeTable        =Make_Range(timeMinimum,timeMaximum,timeTableCount,rangeType=rangeTypeLogarithmic)

       ! Loop through the table and solve for the first crossing distribution.
       call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
       message="    time: "
       write (label,'(f6.3)') timeMinimum
       message=message//label//" to "
       write (label,'(f6.3)') timeMaximum
       message=message//label
       call Galacticus_Display_Message(message,verbosityWorking)
       message="variance: "
       write (label,'(f6.3)') varianceMaximum
       message=message//label
       call Galacticus_Display_Message(message,verbosityWorking)
       loopCountTotal=timeTableCount*(varianceTableCount-1)
       loopCount=0
       !$omp parallel do private(iTime,i,j,sigma1f)
       do iTime=1,timeTableCount
          firstCrossingProbabilityTable(0,iTime)=0.0d0
          firstCrossingProbabilityTable(1,iTime)=                                                                                                   &
               &                                 real(                                                                                              &
               &                                         2.0_kind_quad                                                                              &
               &                                      *(                                                                                            &
               &                                         1.0_kind_quad                                                                              &
               &                                        -erfApproximation(                                                                          &
               &                                                           Excursion_Sets_Barrier(              varianceTable(1),timeTable(iTime))  &
               &                                                          /                  sqrt(2.0_kind_quad*varianceTable(1)                 )  &
               &                                                         )                                                                          &
               &                                       )                                                                                            &
               &                                      /varianceTableStep                                                                          , &
               &                                      kind=kind_dble                                                                                &
               &                                     )
          do i=2,varianceTableCount
             call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityWorking)
             !$omp atomic
             loopCount=loopCount+1
             sigma1f=0.0d0
             do j=1,i-1
                sigma1f= sigma1f                                                                             &
                     &  +firstCrossingProbabilityTable(j,iTime)                                              &
                     &  *real(                                                                               &
                     &         1.0_kind_quad                                                                 &
                     &        -erfApproximation(                                                             &
                     &                          (                                                            &
                     &                            Excursion_Sets_Barrier(varianceTable(i),timeTable(iTime))  &
                     &                           -Excursion_Sets_Barrier(varianceTable(j),timeTable(iTime))  &
                     &                          )                                                            &
                     &                          /sqrt(2.0_kind_quad*(varianceTable(i)-varianceTable(j)))     &
                     &                         )                                                           , &
                     &        kind=kind_dble                                                                 &
                     &       )
             end do
             firstCrossingProbabilityTable(i,iTime)=                                                                                                       &
                  &                                 real(                                                                                                  &
                  &                                      max(                                                                                              &
                  &                                           0.0_kind_quad,                                                                               &
                  &                                           2.0_kind_quad                                                                                &
                  &                                          *(                                                                                            &
                  &                                             1.0_kind_quad                                                                              &
                  &                                            -erfApproximation(                                                                          &
                  &                                                               Excursion_Sets_Barrier(              varianceTable(i),timeTable(iTime))  &
                  &                                                              /                  sqrt(2.0_kind_quad*varianceTable(i)                 )  &
                  &                                                             )                                                                          &
                  &                                           )                                                                                            &
                  &                                          /varianceTableStep                                                                            &
                  &                                          -2.0_kind_quad*sigma1f                                                                        &
                  &                                         )                                                                                            , &
                  &                                      kind=kind_dble                                                                                    &
                  &                                     )
          end do
          ! Force the probability at maximum variance to zero.
          firstCrossingProbabilityTable(varianceTableCount,iTime)=0.0d0
       end do
       !$omp end parallel do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       call Galacticus_Display_Unindent("done",verbosityWorking)
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVariance,reset=interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorTime    ,reset=interpolationResetTime    )
       interpolationResetVariance=.true.
       interpolationResetTime    =.true.
       ! Record that the table is now built.
       tableInitialized=.true.
       ! Write the table to file if possible.
       if (excursionSetFirstCrossingFarahiFileName /= 'none') call Excursion_Sets_First_Crossing_Farahi_Write_File()
    end if

    ! Get interpolation in time.
    iTime    =Interpolate_Locate                 (timeTableCount      ,timeTable    ,interpolationAcceleratorTime    ,time    ,reset=interpolationResetTime    )
    hTime    =Interpolate_Linear_Generate_Factors(timeTableCount      ,timeTable    ,iTime    ,time    )

    ! Get interpolation in variance.
    iVariance=Interpolate_Locate                 (varianceTableCount+1,varianceTable,interpolationAcceleratorVariance,variance,reset=interpolationResetVariance)
    hVariance=Interpolate_Linear_Generate_Factors(varianceTableCount+1,varianceTable,iVariance,variance)

    ! Compute first crossing probability by interpolating.
    Excursion_Sets_First_Crossing_Probability_Farahi=0.0d0
    do jTime=0,1
       do jVariance=0,1
          Excursion_Sets_First_Crossing_Probability_Farahi=                                                                  &
               &                                            Excursion_Sets_First_Crossing_Probability_Farahi                 &
               &                                           +hTime                        (                            jTime) &
               &                                           *hVariance                    (            jVariance            ) &
               &                                           *firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    return
  end function Excursion_Sets_First_Crossing_Probability_Farahi

  subroutine Excursion_Sets_First_Crossing_Rate_Tabulate_Farahi(varianceProgenitor,time)
    !% Tabulate the excursion set crossing rate.
    use Numerical_Ranges
    use Numerical_Interpolation
    use Memory_Management
    use Galacticus_Display
    use Cosmology_Functions
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )               :: time                             , varianceProgenitor
    double precision                , parameter                   :: varianceMinimumDefault    =1.0d-2
    double precision                , parameter                   :: varianceTolerance         =1.0d-6
    real            (kind=kind_quad), allocatable  , dimension(:) :: firstCrossingTableRateQuad       , varianceTableRateBaseQuad, &
         &                                                           varianceTableRateQuad
    logical                                                       :: makeTable
    integer                                                       :: i                                , iTime                    , &
         &                                                           iVariance                        , j                        , &
         &                                                           loopCount                        , loopCountTotal
    double precision                                              :: timeProgenitor                   , varianceMinimumRate
    character       (len=6         )                              :: label
    type            (varying_string)                              :: message
    real            (kind=kind_quad)                              :: crossingFraction                 , effectiveBarrierInitial  , &
         &                                                           sigma1f                          , varianceTableStepRate

    ! Determine if we need to make the table.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    ! Read tables from file if possible.
    if (.not.tableInitializedRate.and.excursionSetFirstCrossingFarahiFileName /= 'none') call Excursion_Sets_First_Crossing_Farahi_Read_File()
    makeTable=.not.tableInitializedRate.or.(varianceProgenitor > varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < timeMinimumRate).or.(time > timeMaximumRate)
    if (makeTable) then
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(varianceTableRate     )) call Dealloc_Array(varianceTableRate     )
       if (allocated(varianceTableRateBase )) call Dealloc_Array(varianceTableRateBase )
       if (allocated(timeTableRate         )) call Dealloc_Array(timeTableRate         )
       if (allocated(firstCrossingTableRate)) call Dealloc_Array(firstCrossingTableRate)
       if (allocated(nonCrossingTableRate  )) call Dealloc_Array(nonCrossingTableRate  )
       if (tableInitializedRate) then
          timeMinimumRate   =min(timeMinimumRate,0.5d0*time)
          timeMaximumRate   =max(timeMaximumRate,2.0d0*time)
          timeTableCountRate=int(log10(timeMaximumRate/timeMinimumRate)*dble(timeTableNumberPerDecade))+1
       else
          timeMinimumRate   =Cosmology_Age(Expansion_Factor_From_Redshift(redshiftMaximumRate))
          timeMaximumRate   =Cosmology_Age(Expansion_Factor_From_Redshift(redshiftMinimumRate))
          timeMinimumRate   =min(timeMinimumRate,0.5d0*time)
          timeMaximumRate   =max(timeMaximumRate,2.0d0*time)
          timeTableCountRate=max(int(log10(timeMaximumRate/timeMinimumRate)*dble(timeTableNumberPerDecade))+1,2)
       end if
       ! Set the default minimum variance.
       varianceMinimumRate       =varianceMinimumDefault
       ! Next reduce the variance if necessary such that the typical amplitude of fluctuations is less (by a factor of sqrt[10])
       ! than the effective barrier height at zero variance for the minimum and maximum times that we must consider.
       varianceMinimumRate       =min(varianceMinimumRate,1.0d-2*real(Excursion_Sets_Barrier_Effective(0.0_kind_quad,timeMinimumRate,0.0_kind_quad,timeMinimumRate*(1.0d0-excursionSetFirstCrossingFarahiFractionalTimeStep)),kind=8)**2)
       varianceMinimumRate       =min(varianceMinimumRate,1.0d-2*real(Excursion_Sets_Barrier_Effective(0.0_kind_quad,timeMaximumRate,0.0_kind_quad,timeMaximumRate*(1.0d0-excursionSetFirstCrossingFarahiFractionalTimeStep)),kind=8)**2)
       varianceMaximumRate       =max(varianceMaximumRate,varianceProgenitor)
       varianceTableCountRate    =int(log10(varianceMaximumRate/varianceMinimumRate)*dble(varianceTableNumberPerDecade))+1
       varianceTableCountRateBase=int(varianceMaximumRate*dble(varianceTableNumberPerUnit))
       call Alloc_Array(varianceTableRate     ,[1+varianceTableCountRate                                                ],lowerBounds=[0    ])
       call Alloc_Array(varianceTableRateBase ,[                         1+varianceTableCountRateBase                   ],lowerBounds=[0    ])
       call Alloc_Array(timeTableRate         ,[                                                      timeTableCountRate]                    )
       call Alloc_Array(firstCrossingTableRate,[1+varianceTableCountRate,1+varianceTableCountRateBase,timeTableCountRate],lowerBounds&
            &=[0,0,1])
       call Alloc_Array(nonCrossingTableRate  ,[                         1+varianceTableCountRateBase,timeTableCountRate],lowerBounds&
            &=[  0,1])
       ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
       varianceTableRate    (0                           )=0.0d0
       varianceTableRate    (1:varianceTableCountRate    )=Make_Variance_Range(varianceMinimumRate,varianceMaximumRate,varianceTableCountRate      ,ratioAtMaximum=10.0d0    )
       varianceTableRateBase(0:varianceTableCountRateBase)=Make_Range         (0.0d0              ,varianceMaximumRate,varianceTableCountRateBase+1,rangeType=rangeTypeLinear)
       ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
       allocate(varianceTableRateQuad     (0:varianceTableCountRate    ))
       varianceTableRateQuad    =varianceTableRate
       allocate(varianceTableRateBaseQuad (0:varianceTableCountRateBase))
       varianceTableRateBaseQuad=varianceTableRateBase
       allocate(firstCrossingTableRateQuad(0:varianceTableCountRate    ))
       ! The time table is logarithmically distributed in time.
       timeTableRate=Make_Range(timeMinimumRate,timeMaximumRate,timeTableCountRate,rangeType=rangeTypeLogarithmic)
       ! Loop through the table and solve for the first crossing distribution.
       call Galacticus_Display_Indent("solving for excursion set barrier crossing rates",verbosityWorking)
       message="    time: "
       write (label,'(f6.3)') timeMinimumRate
       message=message//label//" to "
       write (label,'(f6.3)') timeMaximumRate
       message=message//label
       call Galacticus_Display_Message(message,verbosityWorking)
       message="variance: "
       write (label,'(f6.3)') varianceMaximumRate
       message=message//label
       call Galacticus_Display_Message(message,verbosityWorking)
       loopCountTotal=timeTableCountRate*(1+varianceTableCountRateBase)
       loopCount=0
       !$omp parallel do private(iTime,timeProgenitor,iVariance,varianceTableStepRate,i,j,sigma1f,crossingFraction,effectiveBarrierInitial,firstCrossingTableRateQuad)
       do iTime=1,timeTableCountRate
          ! Compute a suitable progenitor time.
          timeProgenitor=timeTableRate(iTime)*(1.0d0-excursionSetFirstCrossingFarahiFractionalTimeStep)

          ! Loop through the starting variances.
          do iVariance=0,varianceTableCountRateBase
             call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityWorking)
             !$omp atomic
             loopCount=loopCount+1
             ! For zero variance, the rate is initialized to zero.
             firstCrossingTableRateQuad(0)=0.0d0
             ! Compute the step in variance across this first grid cell.
             varianceTableStepRate=varianceTableRateQuad(1)-varianceTableRateQuad(0)
             ! Compute the first crossing distribution at the first grid cell.
             if (varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance) > varianceMaximumRate) then
                firstCrossingTableRateQuad(1)= 0.0d0
             else
                firstCrossingTableRateQuad(1)=                                                                                                                                          &
                     &                         2.0_kind_quad                                                                                                                            &
                     &                        *(                                                                                                                                        &
                     &                           1.0_kind_quad                                                                                                                          &
                     &                          -erfApproximation(                                                                                                                      &
                     &                                             Excursion_Sets_Barrier_Effective(                                                                                    &
                     &                                                                                                       varianceTableRateBaseQuad(iVariance),timeTableRate(iTime), &
                     &                                                                              varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance),timeProgenitor        &
                     &                                                                             )                                                                                    &
                     &                                            /sqrt(2.0_kind_quad*varianceTableRateQuad(1))                                                                         &
                     &                                           )                                                                                                                      &
                     &                         )                                                                                                                                        &
                     &                        /varianceTableStepRate
             end if
             do i=2,varianceTableCountRate
                if (varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance) > varianceMaximumRate) then
                   firstCrossingTableRateQuad(i)=0.0d0
                else
                   effectiveBarrierInitial=Excursion_Sets_Barrier_Effective(                                                                                    &
                        &                                                                            varianceTableRateBaseQuad(iVariance),timeTableRate(iTime), &
                        &                                                   varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance),timeProgenitor        &
                        &                                                  )
                   sigma1f=0.0d0
                   do j=1,i-1
                      varianceTableStepRate=(varianceTableRateQuad(j+1)-varianceTableRateQuad(j-1))/2.0_kind_quad
                      sigma1f= sigma1f                                                                                                                                   &
                           &  +firstCrossingTableRateQuad(j)                                                                                                             &
                           &  *varianceTableStepRate                                                                                                                     &
                           &  *(                                                                                                                                         &
                           &     1.0_kind_quad                                                                                                                           &
                           &    -erfApproximation(                                                                                                                       &
                           &                      (                                                                                                                      &
                           &                       +effectiveBarrierInitial                                                                                              &
                           &                       -Excursion_Sets_Barrier_Effective(                                                                                    &
                           &                                                                                  varianceTableRateBaseQuad(iVariance),timeTableRate(iTime), &
                           &                                                         varianceTableRateQuad(j)+varianceTableRateBaseQuad(iVariance),timeProgenitor        &
                           &                                                        )                                                                                    &
                           &                      )                                                                                                                      &
                           &                      /sqrt(2.0_kind_quad*(varianceTableRateQuad(i)-varianceTableRateQuad(j)))                                               &
                           &                     )                                                                                                                       &
                           &   )
                   end do
                   varianceTableStepRate=varianceTableRateQuad(i)-varianceTableRateQuad(i-1)
                   firstCrossingTableRateQuad(i)=                                                                        &
                        &                         max(                                                                   &
                        &                             0.0_kind_quad,                                                     &
                        &                             (                                                                  &
                        &                               2.0_kind_quad                                                    &
                        &                              *(                                                                &
                        &                                 1.0_kind_quad                                                  &
                        &                                -erfApproximation(                                              &
                        &                                                  +effectiveBarrierInitial                      &
                        &                                                  /sqrt(2.0_kind_quad*varianceTableRateQuad(i)) &
                        &                                                 )                                              &
                        &                               )                                                                &
                        &                               -2.0_kind_quad*sigma1f                                           &
                        &                             )                                                                  &
                        &                             /varianceTableStepRate                                             &
                        &                            )
                end if
             end do
             ! Compute the fraction of trajectories which never cross the barrier.
             crossingFraction=0.0_kind_quad
             do j=0,varianceTableCountRate-1
                varianceTableStepRate=varianceTableRateQuad(j+1)-varianceTableRateQuad(j)
                crossingFraction=                                   &
                     &            crossingFraction                  &
                     &           +0.5_kind_quad                     &
                     &           *(                                 &
                     &              firstCrossingTableRateQuad(j  ) &
                     &             +firstCrossingTableRateQuad(j+1) &
                     &            )                                 &
                     &           *varianceTableStepRate
             end do
             ! Compute the rate for trajectories which never cross the barrier.
             nonCrossingTableRate(iVariance,iTime)=                                                         &
                  &                                real(                                                    &
                  &                                      (1.0_kind_quad-crossingFraction)                   &
                  &                                     /timeTableRate(iTime)                               &
                  &                                     /excursionSetFirstCrossingFarahiFractionalTimeStep, &
                  &                                     kind=kind_dble                                      &
                  &                                    )
             ! Store the compute crossing rate in our table.
             firstCrossingTableRate(:,iVariance,iTime)=real(firstCrossingTableRateQuad,kind=kind_dble)
           end do
           ! Divide through by the time step to get the rate of barrier crossing.
           firstCrossingTableRate(:,:,iTime)= firstCrossingTableRate(:,:,iTime)                 &
                &                            /timeTableRate         (    iTime)                 &
                &                            /excursionSetFirstCrossingFarahiFractionalTimeStep
        end do
       !$omp end parallel do
       ! Deallocate work arrays.
       deallocate(varianceTableRateBaseQuad )
       deallocate(varianceTableRateQuad     )
       deallocate(firstCrossingTableRateQuad)
       call Galacticus_Display_Counter_Clear(       verbosityWorking)
       call Galacticus_Display_Unindent     ("done",verbosityWorking)
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVarianceRate    ,reset=interpolationResetVarianceRate    )
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVarianceRateBase,reset=interpolationResetVarianceRateBase)
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorTimeRate        ,reset=interpolationResetTimeRate        )
       interpolationResetVarianceRate    =.true.
       interpolationResetVarianceRateBase=.true.
       interpolationResetTimeRate        =.true.
       ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
       varianceRatePrevious=-1.0d0
       timeRatePrevious    =-1.0d0
       ! Record that the table is now built.
       tableInitializedRate=.true.
       ! Write the table to file if possible.
       if (excursionSetFirstCrossingFarahiFileName /= 'none') call Excursion_Sets_First_Crossing_Farahi_Write_File()
    end if
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    return
  end subroutine Excursion_Sets_First_Crossing_Rate_Tabulate_Farahi

  double precision function Excursion_Sets_First_Crossing_Rate_Farahi(variance,varianceProgenitor,time)
    !% Return the rate for excursion set first crossing.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   )  :: time               , variance, varianceProgenitor
    double precision, dimension(0:1) :: hVarianceProgenitor
    integer                          :: iVarianceProgenitor, jTime   , jVariance         , &
         &                              jVarianceProgenitor

    ! For progenitor variances less than or equal to the original variance, return zero.
    if (varianceProgenitor <= variance) then
       Excursion_Sets_First_Crossing_Rate_Farahi=0.0d0
       return
    end if

    ! Ensure that the rate is tabulated.
    call Excursion_Sets_First_Crossing_Rate_Tabulate_Farahi(varianceProgenitor,time)

    ! For progenitor variances greater than the maximum allowed variance, return zero.
    if (varianceProgenitor > varianceMaximumRate) then
       Excursion_Sets_First_Crossing_Rate_Farahi=0.0d0
       return
    end if

    ! Get interpolation in time.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    if (time /= timeRatePrevious) then
       timeRatePrevious    =time
       iTimeRate       =Interpolate_Locate                 (timeTableCountRate      ,timeTableRate    ,interpolationAcceleratorTimeRate    ,time    ,reset=interpolationResetTimeRate    )
       hTimeRate       =Interpolate_Linear_Generate_Factors(timeTableCountRate      ,timeTableRate    ,iTimeRate ,time    )
    end if

    ! Get interpolation in variance.
    if (variance /= varianceRatePrevious) then
       varianceRatePrevious=variance
       iVarianceRate   =Interpolate_Locate                 (varianceTableCountRateBase+1,varianceTableRateBase,interpolationAcceleratorVarianceRateBase,variance,reset=interpolationResetVarianceRateBase)
       hVarianceRate   =Interpolate_Linear_Generate_Factors(varianceTableCountRateBase+1,varianceTableRateBase,iVarianceRate,variance)
    end if

    ! Get interpolation in progenitor variance.
    iVarianceProgenitor=Interpolate_Locate                 (varianceTableCountRate+1,varianceTableRate,interpolationAcceleratorVarianceRate,varianceProgenitor-variance,reset=interpolationResetVarianceRate)

    ! Catch cases where the maximum variance is approached.
    if (varianceTableRate(iVarianceProgenitor)+variance > varianceMaximumRate) then
       ! Force the rate to drop to zero at the maximum variance. (Necessary because we will not have a tabulated point precisely
       ! at the maximum variance.)
       hVarianceProgenitor=[                                                                            &
            &               +1.0d0                                                                      &
            &               -((varianceProgenitor -variance)-varianceTableRate(iVarianceProgenitor-1))  &
            &               /((varianceMaximumRate-variance)-varianceTableRate(iVarianceProgenitor-1)), &
            &               +0.0d0                                                                      &
            &              ]
    else
       hVarianceProgenitor=Interpolate_Linear_Generate_Factors(varianceTableCountRate+1,varianceTableRate,iVarianceProgenitor,varianceProgenitor-variance)
    end if

    ! Compute first crossing probability by interpolating.
    Excursion_Sets_First_Crossing_Rate_Farahi=0.0d0
    do jTime=0,1
       do jVariance=0,1
          do jVarianceProgenitor=0,1
             Excursion_Sets_First_Crossing_Rate_Farahi=                                                                                                             &
                  &                                     Excursion_Sets_First_Crossing_Rate_Farahi                                                                   &
                  &                                    +hTimeRate             (                                                                              jTime) &
                  &                                    *hVarianceRate         (                                                          jVariance                ) &
                  &                                    *hVarianceProgenitor   (                      jVarianceProgenitor                                          ) &
                  &                                    *firstCrossingTableRate(iVarianceProgenitor-1+jVarianceProgenitor,iVarianceRate-1+jVariance,iTimeRate+jTime)
          end do
       end do
    end do
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    return
  end function Excursion_Sets_First_Crossing_Rate_Farahi

  double precision function Excursion_Sets_Non_Crossing_Rate_Farahi(variance,time)
    !% Return the rate for excursion set non-crossing.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: time , variance
    integer                         :: jTime, jVariance

    ! Ensure that the rate is tabulated.
    call Excursion_Sets_First_Crossing_Rate_Tabulate_Farahi(variance,time)

    ! Get interpolation in time.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    if (time /= timeRatePrevious) then
       timeRatePrevious    =time
       iTimeRate    =Interpolate_Locate                 (timeTableCountRate      ,timeTableRate    ,interpolationAcceleratorTimeRate    ,time    ,reset=interpolationResetTimeRate    )
       hTimeRate    =Interpolate_Linear_Generate_Factors(timeTableCountRate      ,timeTableRate    ,iTimeRate    ,time    )
    end if

    ! Get interpolation in variance.
    if (variance /= varianceRatePrevious) then
       varianceRatePrevious=variance
       iVarianceRate=Interpolate_Locate                 (varianceTableCountRateBase+1,varianceTableRateBase,interpolationAcceleratorVarianceRateBase,variance,reset=interpolationResetVarianceRateBase)
       hVarianceRate=Interpolate_Linear_Generate_Factors(varianceTableCountRateBase+1,varianceTableRateBase,iVarianceRate,variance)
    end if

    ! Compute non-crossing probability by interpolating.
    Excursion_Sets_Non_Crossing_Rate_Farahi=0.0d0
    do jTime=0,1
       do jVariance=0,1
          Excursion_Sets_Non_Crossing_Rate_Farahi=                                                                 &
               &                                   Excursion_Sets_Non_Crossing_Rate_Farahi                         &
               &                                  +hTimeRate           (                                    jTime) &
               &                                  *hVarianceRate       (                jVariance                ) &
               &                                  *nonCrossingTableRate(iVarianceRate-1+jVariance,iTimeRate+jTime)
       end do
    end do
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Farahi_Init)
    return
  end function Excursion_Sets_Non_Crossing_Rate_Farahi

 function Excursion_Sets_Barrier_Effective(variance0,time0,variance,time)
    !% The effective barrier for conditional excursion sets.
    use Kind_Numbers
    use Excursion_Sets_Barriers
    implicit none
    real            (kind=kind_quad)                :: Excursion_Sets_Barrier_Effective
    real            (kind=kind_quad), intent(in   ) :: variance                                       , variance0
    double precision                , intent(in   ) :: time                                           , time0
    real            (kind=kind_quad), save          :: variance0Previous               =-1.0_kind_quad, time0Previous=-1.0_kind_quad, &
         &                                             barrier0Previous
    !$omp threadprivate(variance0Previous,time0Previous,barrier0Previous)

    if (variance0 /= variance0Previous .or. time0 /= time0Previous) then
       variance0Previous=variance0
       time0Previous    =time0
       barrier0Previous =Excursion_Sets_Barrier(real(variance0,kind=8),time0,ratesCalculation=.true.)
    end if
    Excursion_Sets_Barrier_Effective=                                                                            &
         &                            Excursion_Sets_Barrier(real(variance,kind=8),time,ratesCalculation=.true.) &
         &                           -barrier0Previous
    return
  end function Excursion_Sets_Barrier_Effective

 function erfApproximation(x)
    !% An \href{http://sites.google.com/site/winitzki/sergei-winitzkis-files/erf-approx.pdf}{approximation to the error function}
    !% that is designed to be very accurate in the vicinity of zero and infinity.
    use Numerical_Constants_Math
    implicit none
    real(kind=kind_quad)                :: erfApproximation
    real(kind=kind_quad), intent(in   ) :: x
    real(kind=kind_quad), parameter     :: a               =8.0_kind_quad*(PiQuadPrecision-3.0_kind_quad)/3.0_kind_quad/PiQuadPrecision/(4.0_kind_quad-PiQuadPrecision)

    erfApproximation=sqrt(1.0_kind_quad-exp(-x**2*(4.0_kind_quad/PiQuadPrecision+a*x**2)/(1.0_kind_quad+a*x**2)))
    return
  end function erfApproximation

  subroutine Excursion_Sets_First_Crossing_Farahi_Read_File()
    !% Read tabulated data on excursion set first crossing probabilities from file.
    use IO_HDF5
    use File_Utilities
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    type            (hdf5Object)                            :: dataFile                  , dataGroup
    double precision            , allocatable, dimension(:) :: varianceTableBaseTemporary, varianceTableTemporary

    ! Return immediately if the file does not exist.
    if (.not.File_Exists(excursionSetFirstCrossingFarahiFileName)) return
    ! Open the data file.
    call dataFile%openFile(char(excursionSetFirstCrossingFarahiFileName))
    ! Check if the standard table is populated.
    if (dataFile%hasGroup('probability')) then
       ! Deallocate arrays if necessary.
       if (allocated(varianceTable                )) call Dealloc_Array(varianceTable                )
       if (allocated(timeTable                    )) call Dealloc_Array(timeTable                    )
       if (allocated(firstCrossingProbabilityTable)) call Dealloc_Array(firstCrossingProbabilityTable)
       ! Read the datasets.
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%readDataset('variance'                ,varianceTableTemporary       )
       call dataGroup%readDataset('time'                    ,timeTable                    )
       call dataGroup%readDataset('firstCrossingProbability',firstCrossingProbabilityTable)
       call dataGroup%close()
       ! Set table sizes and limits.
       varianceTableCount=size(varianceTableTemporary)-1
       timeTableCount    =size(timeTable      )
       ! Transfer to variance table.
       call Alloc_Array(varianceTable,[1+varianceTableCount],lowerBounds=[0])
       varianceTable(0:varianceTableCount)=varianceTableTemporary(1:varianceTableCount+1)
       call Dealloc_Array(varianceTableTemporary)
       ! Set table limits.
       varianceMaximum   =varianceTable(varianceTableCount)
       varianceTableStep =varianceTable(1)-varianceTable(0)
       timeMinimum       =timeTable    (                 1)
       timeMaximum       =timeTable    (    timeTableCount)
       tableInitialized  =.true.
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVariance,reset=interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorTime    ,reset=interpolationResetTime    )
       interpolationResetVariance=.true.
       interpolationResetTime    =.true.
    end if
    ! Check if the rate table is populated.
    if (dataFile%hasGroup('rate')) then
       ! Deallocate arrays if necessary.
       if (allocated(varianceTableRate     )) call Dealloc_Array(varianceTableRate     )
       if (allocated(varianceTableRateBase )) call Dealloc_Array(varianceTableRateBase )
       if (allocated(timeTableRate         )) call Dealloc_Array(timeTableRate         )
       if (allocated(firstCrossingTableRate)) call Dealloc_Array(firstCrossingTableRate)
       if (allocated(nonCrossingTableRate  )) call Dealloc_Array(nonCrossingTableRate  )
       ! Read the datasets.
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%readDataset('variance'         ,varianceTableTemporary    )
       call dataGroup%readDataset('varianceBase'     ,varianceTableBaseTemporary)
       call dataGroup%readDataset('time'             ,timeTableRate             )
       call dataGroup%readDataset('firstCrossingRate',firstCrossingTableRate    )
       call dataGroup%readDataset('nonCrossingRate'  ,nonCrossingTableRate      )
       call dataGroup%close()
       ! Set table sizes and limits.
       varianceTableCountRate    =size(varianceTableTemporary    )-1
       varianceTableCountRateBase=size(varianceTableBaseTemporary)-1
       timeTableCountRate        =size(timeTableRate             )
       ! Transfer to variance table.
       call Alloc_Array(varianceTableRate    ,[1+varianceTableCountRate    ],lowerBounds=[0])
       call Alloc_Array(varianceTableRateBase,[1+varianceTableCountRateBase],lowerBounds=[0])
       varianceTableRate    (0:varianceTableCountRate    )=varianceTableTemporary    (1:varianceTableCountRate    +1)
       varianceTableRateBase(0:varianceTableCountRateBase)=varianceTableBaseTemporary(1:varianceTableCountRateBase+1)
       call Dealloc_Array(varianceTableTemporary    )
       call Dealloc_Array(varianceTableBaseTemporary)
       ! Set table limits.
       varianceMaximumRate =varianceTableRate(varianceTableCountRate)
       timeMinimumRate     =timeTableRate    (                     1)
       timeMaximumRate     =timeTableRate    (    timeTableCountRate)
       tableInitializedRate=.true.
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVarianceRate    ,reset=interpolationResetVarianceRate    )
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVarianceRateBase,reset=interpolationResetVarianceRateBase)
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorTimeRate        ,reset=interpolationResetTimeRate        )
       interpolationResetVarianceRate    =.true.
       interpolationResetVarianceRateBase=.true.
       interpolationResetTimeRate        =.true.
    end if
    ! Close the data file.
    call dataFile%close()
    return
  end subroutine Excursion_Sets_First_Crossing_Farahi_Read_File

  subroutine Excursion_Sets_First_Crossing_Farahi_Write_File()
    !% Write tabulated data on excursion set first crossing probabilities to file.
    use IO_HDF5
    implicit none
    type(hdf5Object) :: dataFile, dataGroup

    ! Don't write anything if neither table is initialized.
    if (.not.(tableInitialized.or.tableInitializedRate)) return
    ! Open the data file.
    call dataFile%openFile(char(excursionSetFirstCrossingFarahiFileName),overWrite=.true.,chunkSize=100,compressionLevel=9)
    ! Check if the standard table is populated.
    if (tableInitialized) then
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%writeDataset(varianceTable                ,'variance'                ,'The variance at which results are tabulated.'                         )
       call dataGroup%writeDataset(timeTable                    ,'time'                    ,'The cosmic times at which results are tabulated.'                     )
       call dataGroup%writeDataset(firstCrossingProbabilityTable,'firstCrossingProbability','The probability of first crossing as a function of variance and time.')
       call dataGroup%close()
    end if
    ! Check if the rate table is populated.
    if (tableInitializedRate) then
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%writeDataset(varianceTableRate     ,'variance'         ,'The variance at which results are tabulated.'                               )
       call dataGroup%writeDataset(varianceTableRateBase ,'varianceBase'     ,'The variance of the base halo at which results are tabulated.'              )
       call dataGroup%writeDataset(timeTableRate         ,'time'             ,'The cosmic times at which results are tabulated.'                           )
       call dataGroup%writeDataset(firstCrossingTableRate,'firstCrossingRate','The probability rate of first crossing as a function of variances and time.')
       call dataGroup%writeDataset(nonCrossingTableRate  ,'nonCrossingRate'  ,'The probability rate of non crossing as a function of variance and time.')
       call dataGroup%close()
    end if
    ! Close the data file.
    call dataFile%close()
    return
  end subroutine Excursion_Sets_First_Crossing_Farahi_Write_File

  function Make_Variance_Range(rangeMinimum,rangeMaximum,rangeNumber,ratioAtMaximum) result (rangeValues)
    !% Builds a numerical range between {\tt rangeMinimum} and {\tt rangeMaximum} using {\tt rangeNumber} points with spacing that
    !% varies from logarithmic to linear spacing with the transition point controlled by {\tt ratioAtMaximum}.
    implicit none
    double precision, intent(in   )          :: rangeMaximum               , rangeMinimum    , &
         &                                      ratioAtMaximum
    integer         , intent(in   )          :: rangeNumber
    double precision, dimension(rangeNumber) :: rangeValues   (rangeNumber)
    integer                                  :: iRange
    double precision                         :: rangeLinear                , rangeLogarithmic

    do iRange=1,rangeNumber
       rangeLinear        =        rangeMinimum +   (rangeMaximum-rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1)
       rangeLogarithmic   =exp(log(rangeMinimum)+log(rangeMaximum/rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1))
       rangeValues(iRange)=(1.0d0+1.0d0/ratioAtMaximum)/(1.0d0/rangeLinear+1.0d0/rangeLogarithmic/ratioAtMaximum)
    end do
    return
  end function Make_Variance_Range

end module Excursion_Sets_First_Crossing_Farahi
