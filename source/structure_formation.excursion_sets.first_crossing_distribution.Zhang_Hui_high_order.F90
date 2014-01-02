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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson.

!% Contains a module which implements the first crossing distribution for excursion set calculations of dark matter halo formation
!% using a high order modification of the methodology of \cite{zhang_random_2006}.

module Excursion_Sets_First_Crossing_Zhang_Hui_High
  !% Implements the first crossing distribution for excursion set calculations of dark matter halo formation using a high order
  !% modification of the methodology of \cite{zhang_random_2006}.
  use FGSL
  private
  public :: Excursion_Sets_First_Crossing_Zhang_Hui_High_Initialize

  double precision                                                 :: timeMaximum                  =0.0d0  , timeMinimum                     =0.0d0 , &
       &                                                              varianceMaximum              =0.0d0
  integer                                                          :: timeTableCount                       , varianceTableCount
  integer                            , parameter                   :: timeTableNumberPerDecade     =40     , varianceTableNumberPerUnit      =400
  double precision                   , allocatable, dimension(:,:) :: firstCrossingProbabilityTable
  double precision                   , allocatable, dimension(:  ) :: timeTable                            , varianceTable
  double precision                                                 :: varianceTableStep
  logical                                                          :: tableInitialized             =.false.
  type            (fgsl_interp_accel)                              :: interpolationAcceleratorTime         , interpolationAcceleratorVariance
  logical                                                          :: interpolationResetTime       =.true. , interpolationResetVariance      =.true.

contains

  !# <excursionSetFirstCrossingMethod>
  !#  <unitName>Excursion_Sets_First_Crossing_Zhang_Hui_High_Initialize</unitName>
  !# </excursionSetFirstCrossingMethod>
  subroutine Excursion_Sets_First_Crossing_Zhang_Hui_High_Initialize(excursionSetFirstCrossingMethod&
       &,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get&
       &,Excursion_Sets_Non_Crossing_Rate_Get)
    !% Initialize the ``ZhangHui2006 high order'' first crossing distribution for excursion sets module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: excursionSetFirstCrossingMethod
    procedure(Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High), intent(inout), pointer :: Excursion_Sets_First_Crossing_Probability_Get
    procedure(Excursion_Sets_First_Crossing_Rate_Zhang_Hui_High), intent(inout), pointer :: Excursion_Sets_First_Crossing_Rate_Get
    procedure(Excursion_Sets_Non_Crossing_Rate_Zhang_Hui_High), intent(inout), pointer :: Excursion_Sets_Non_Crossing_Rate_Get

    if (excursionSetFirstCrossingMethod == 'ZhangHui2006HighOrder') then
       Excursion_Sets_First_Crossing_Probability_Get => Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High
       Excursion_Sets_First_Crossing_Rate_Get        => Excursion_Sets_First_Crossing_Rate_Zhang_Hui_High
       Excursion_Sets_Non_Crossing_Rate_Get          => Excursion_Sets_Non_Crossing_Rate_Zhang_Hui_High
    end if
    return
  end subroutine Excursion_Sets_First_Crossing_Zhang_Hui_High_Initialize

  double precision function Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High(variance,time)
    !% Return the probability for excursion set first crossing using a high order modification of the methodology of \cite{zhang_random_2006}.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Display
    use Excursion_Sets_First_Crossing_Zhang_Hui_Utilities
    implicit none
    double precision, intent(in   )                 :: time                                 , variance
    double precision               , dimension(0:1) :: hTime                                , hVariance
    double precision, allocatable  , dimension(:,:) :: firstCrossingProbabilityTablePrevious
    logical                                         :: makeTable                            , tableIsExtendable
    integer                                         :: i                                    , iTime                     , &
         &                                             iVariance                            , j                         , &
         &                                             jTime                                , jVariance                 , &
         &                                             k                                    , varianceTableCountPrevious
    double precision                                :: g2Average                            , summation0                , &
         &                                             summation1                           , summation2                , &
         &                                             timeMaximumPrevious                  , timeMinimumPrevious

    ! Determine if we need to make the table. Only recompute if the variance is sufficiently larger than that tabulated that we
    ! will add at least one new point to the tabulation.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High_Init)
    makeTable=.not.tableInitialized.or.(variance > varianceMaximum+1.0d0/dble(varianceTableNumberPerUnit)).or.(time < timeMinimum).or.(time > timeMaximum)
    if (makeTable) then
       ! Compute the range of times to tabulate.
       timeMinimumPrevious=timeMinimum
       timeMaximumPrevious=timeMaximum
       if (tableInitialized) then
          timeMinimum=min(timeMinimum,time)
          timeMaximum=max(timeMaximum,time)
       else
          timeMinimum=0.5d0*time
          timeMaximum=2.0d0*time
       end if
       tableIsExtendable=(timeMinimum == timeMinimumPrevious .and. timeMaximum == timeMaximumPrevious .and. tableInitialized)
       ! Make a copy of the current table if possible.
       if (tableIsExtendable) then
          ! We have a pre-existing table with the same range of times. Therefore, we can simply extend this table. So make a copy
          ! of the existing table.
          call Move_Alloc(firstCrossingProbabilityTable,firstCrossingProbabilityTablePrevious)
          varianceTableCountPrevious=varianceTableCount
       else
          varianceTableCountPrevious=-1
       end if
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(varianceTable                )) call Dealloc_Array(varianceTable                )
       if (allocated(timeTable                    )) call Dealloc_Array(timeTable                    )
       if (allocated(firstCrossingProbabilityTable)) call Dealloc_Array(firstCrossingProbabilityTable)
       varianceMaximum   =max(varianceMaximum,variance)
       varianceTableCount=int(varianceMaximum*dble(varianceTableNumberPerUnit))
       varianceMaximum   =dble(varianceTableCount)/dble(varianceTableNumberPerUnit)
       timeTableCount    =int(log10(timeMaximum/timeMinimum)*dble(timeTableNumberPerDecade))+1
       call Alloc_Array(varianceTable                ,[1+varianceTableCount]                 ,lowerBounds=[0  ])
       call Alloc_Array(timeTable                                           ,[timeTableCount]                  )
       call Alloc_Array(firstCrossingProbabilityTable,[1+varianceTableCount , timeTableCount],lowerBounds=[0,1])
       varianceTable    =Make_Range(0.0d0,varianceMaximum,varianceTableCount+1,rangeType=rangeTypeLinear)
       varianceTableStep=varianceTable(4)-varianceTable(0)
       timeTable        =Make_Range(timeMinimum,timeMaximum,timeTableCount,rangeType=rangeTypeLogarithmic)

       if (tableIsExtendable) then
          firstCrossingProbabilityTable(0:varianceTableCountPrevious,:)=firstCrossingProbabilityTablePrevious
          call Dealloc_Array(firstCrossingProbabilityTablePrevious)
       end if

       ! Loop through the table and solve for the first crossing distribution.
       call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
       do iTime=1,timeTableCount
          do i=varianceTableCountPrevious+1,varianceTableCount
             call Galacticus_Display_Counter(                                                                            &
                  &                           int(                                                                       &
                  &                                100.0d0                                                               &
                  &                               *dble(                                                                 &
                  &                                                    (i                 -varianceTableCountPrevious-1) &
                  &                                     +(iTime-1)    *(varianceTableCount-varianceTableCountPrevious-1) &
                  &                                    )                                                                 &
                  &                               /dble(                                                                 &
                  &                                     timeTableCount*(varianceTableCount-varianceTableCountPrevious-1) &
                  &                                    )                                                                 &
                  &                              )                                                                       &
                  &                          ,i==varianceTableCountPrevious+1 .and. iTime==1,verbosityWorking            &
                  &                         )

             if (i  > 3) then
                summation0=0.0d0
                summation1=0.0d0
                summation2=0.0d0

                k=mod(i,4)
                if      (k == 3) then
                   summation0=(                                                                                          &
                        &       3.0d0*firstCrossingProbabilityTable(1,iTime)*g_2(varianceTable(i),varianceTable(1),timeTable(iTime)) &
                        &      +3.0d0*firstCrossingProbabilityTable(2,iTime)*g_2(varianceTable(i),varianceTable(2),timeTable(iTime)) &
                        &      +      firstCrossingProbabilityTable(3,iTime)*g_2(varianceTable(i),varianceTable(3),timeTable(iTime)) &
                        &     )                                                                                          &
                        &     *3.0d0*varianceTableStep/32.0d0
                else if (k == 2) then
                   summation0=(                                                                                          &
                        &       4.0d0*firstCrossingProbabilityTable(1,iTime)*g_2(varianceTable(i),varianceTable(1),timeTable(iTime)) &
                        &      +      firstCrossingProbabilityTable(2,iTime)*g_2(varianceTable(i),varianceTable(2),timeTable(iTime)) &
                        &     )                                                                                          &
                        &     *varianceTableStep/12.0d0
                else if (k == 1) then
                   summation0=        firstCrossingProbabilityTable(1,iTime)*g_2(varianceTable(i),varianceTable(1),timeTable(iTime)) &
                        &     *varianceTableStep/8.0d0
                else
                   summation0=0.0d0
                end if

                do j=k,i-8,4
                   summation1= summation1                                                                                    &
                        &     + 7.0d0*firstCrossingProbabilityTable(j  ,iTime)*g_2(varianceTable(i),varianceTable(j  ),timeTable(iTime)) &
                        &     +32.0d0*firstCrossingProbabilityTable(j+1,iTime)*g_2(varianceTable(i),varianceTable(j+1),timeTable(iTime)) &
                        &     +12.0d0*firstCrossingProbabilityTable(j+2,iTime)*g_2(varianceTable(i),varianceTable(j+2),timeTable(iTime)) &
                        &     +32.0d0*firstCrossingProbabilityTable(j+3,iTime)*g_2(varianceTable(i),varianceTable(j+3),timeTable(iTime)) &
                        &     + 7.0d0*firstCrossingProbabilityTable(j+4,iTime)*g_2(varianceTable(i),varianceTable(j+4),timeTable(iTime))
                end do
                summation1=summation1*varianceTableStep/90.0d0
                g2Average=g_2_Integrated(varianceTable(i),varianceTableStep,timeTable(iTime))
                summation2= g2Average                                         &
                     &     *(                                                 &
                     &         7.0d0*firstCrossingProbabilityTable(i-4,iTime) &
                     &       +32.0d0*firstCrossingProbabilityTable(i-3,iTime) &
                     &       +12.0d0*firstCrossingProbabilityTable(i-2,iTime) &
                     &       +32.0d0*firstCrossingProbabilityTable(i-1,iTime) &
                     &      )                                                 &
                     &     /90.0d0

                firstCrossingProbabilityTable(i,iTime)=(g_1(varianceTable(i),timeTable(iTime))+summation0+summation1+summation2)/(1.0d0-7.0d0*g2Average/90.0d0)
             else if (i == 3) then
                g2Average=g_2_Integrated(varianceTable(i),0.75d0*varianceTableStep,timeTable(iTime))
                firstCrossingProbabilityTable(i,iTime)=(g_1(varianceTable(i),timeTable(iTime))+g2Average*3.0d0*(firstCrossingProbabilityTable(1&
                     &,iTime)+firstCrossingProbabilityTable(2,iTime))/8.0d0)/(1.0d0-g2Average/8.0d0)
             else if (i == 2) then
                g2Average=g_2_Integrated(varianceTable(i),0.5d0*varianceTableStep,timeTable(iTime))
                firstCrossingProbabilityTable(i,iTime)=(g_1(varianceTable(i),timeTable(iTime))+g2Average*4.0d0*firstCrossingProbabilityTable(1&
                     &,iTime)/6.0d0)/(1.0d0-g2Average/6.0d0)
             else if (i == 1) then
                g2Average=g_2_Integrated(varianceTable(i),0.25d0*varianceTableStep,timeTable(iTime))
                firstCrossingProbabilityTable(i,iTime)= g_1(varianceTable(i),timeTable(iTime))/(1.0d0-g2Average/2.0d0)
             else if (i == 0) then
                firstCrossingProbabilityTable(i,iTime)= 0.0d0
             end if
          end do
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       call Galacticus_Display_Unindent("done",verbosityWorking)
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorVariance,reset=interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorTime    ,reset=interpolationResetTime    )
       interpolationResetVariance=.true.
       interpolationResetTime    =.true.
       ! Record that the table is now built.
       tableInitialized=.true.
    end if
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High_Init)

    ! Get interpolation in time.
    iTime    =Interpolate_Locate                 (timeTableCount      ,timeTable    ,interpolationAcceleratorTime    ,time    ,reset=interpolationResetTime    )
    hTime    =Interpolate_Linear_Generate_Factors(timeTableCount      ,timeTable    ,iTime    ,time    )

    ! Get interpolation in variance.
    iVariance=Interpolate_Locate                 (varianceTableCount+1,varianceTable,interpolationAcceleratorVariance,variance,reset=interpolationResetVariance)
    hVariance=Interpolate_Linear_Generate_Factors(varianceTableCount+1,varianceTable,iVariance,variance)

    ! Compute first crossing probability by interpolating.
    Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High=0.0d0
    do jTime=0,1
       do jVariance=0,1
          Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High=Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High+hTime(jTime)*hVariance(jVariance)*firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function Excursion_Sets_First_Crossing_Probability_Zhang_Hui_High

  double precision function Excursion_Sets_First_Crossing_Rate_Zhang_Hui_High(variance,varianceProgenitor,time)
    !% Return the rate for excursion set first crossing assuming a linear barrier.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: time, variance, varianceProgenitor

    call Galacticus_Error_Report('Excursion_Sets_First_Crossing_Rate_Zhang_Hui_High','barrier crossing rates are not implemented for this method [too slow]')
    return
  end function Excursion_Sets_First_Crossing_Rate_Zhang_Hui_High

  double precision function Excursion_Sets_Non_Crossing_Rate_Zhang_Hui_High(variance,time)
    !% Return the rate for excursion set non-crossing.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: time, variance

    call Galacticus_Error_Report('Excursion_Sets_Non_Crossing_Rate_Zhang_Hui_High','barrier non-crossing rates are not implemented for this method [too slow]')
    return
  end function Excursion_Sets_Non_Crossing_Rate_Zhang_Hui_High

end module Excursion_Sets_First_Crossing_Zhang_Hui_High
