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

!% Contains a module which implements the first crossing distribution for excursion set calculations
!% of dark matter halo formation using the methodology of \cite{zhang_random_2006}.

module Excursion_Sets_First_Crossing_Zhang_Hui
  !% Implements the first crossing distribution for excursion set calculations of dark matter halo formation using the methodology
  !% of \cite{zhang_random_2006}.
  use FGSL
  private
  public :: Excursion_Sets_First_Crossing_Zhang_Hui_Initialize

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
  !#  <unitName>Excursion_Sets_First_Crossing_Zhang_Hui_Initialize</unitName>
  !# </excursionSetFirstCrossingMethod>
  subroutine Excursion_Sets_First_Crossing_Zhang_Hui_Initialize(excursionSetFirstCrossingMethod&
       &,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get&
       &,Excursion_Sets_Non_Crossing_Rate_Get)
    !% Initialize the ``ZhangHui2006'' first crossing distribution for excursion sets module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: excursionSetFirstCrossingMethod
    procedure(Excursion_Sets_First_Crossing_Probability_Zhang_Hui), intent(inout), pointer :: Excursion_Sets_First_Crossing_Probability_Get
    procedure(Excursion_Sets_First_Crossing_Rate_Zhang_Hui), intent(inout), pointer :: Excursion_Sets_First_Crossing_Rate_Get
    procedure(Excursion_Sets_Non_Crossing_Rate_Zhang_Hui), intent(inout), pointer :: Excursion_Sets_Non_Crossing_Rate_Get

    if (excursionSetFirstCrossingMethod == 'ZhangHui2006') then
       Excursion_Sets_First_Crossing_Probability_Get => Excursion_Sets_First_Crossing_Probability_Zhang_Hui
       Excursion_Sets_First_Crossing_Rate_Get        => Excursion_Sets_First_Crossing_Rate_Zhang_Hui
       Excursion_Sets_Non_Crossing_Rate_Get          => Excursion_Sets_Non_Crossing_Rate_Zhang_Hui
    end if
    return
  end subroutine Excursion_Sets_First_Crossing_Zhang_Hui_Initialize

  double precision function Excursion_Sets_First_Crossing_Probability_Zhang_Hui(variance,time)
    !% Return the probability for excursion set first crossing using the methodology of \cite{zhang_random_2006}.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Display
    use Excursion_Sets_First_Crossing_Zhang_Hui_Utilities
    implicit none
    double precision, intent(in   )  :: time             , variance
    double precision, dimension(0:1) :: hTime            , hVariance
    logical                          :: makeTable
    integer                          :: i                , iTime    , iVariance, &
         &                              j                , jTime    , jVariance
    double precision                 :: summedProbability

    ! Determine if we need to make the table.
    !$omp critical (Excursion_Sets_First_Crossing_Probability_Zhang_Hui_Init)
    makeTable=.not.tableInitialized.or.(variance > varianceMaximum).or.(time < timeMinimum).or.(time > timeMaximum)
    if (makeTable) then
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(varianceTable                )) call Dealloc_Array(varianceTable                )
       if (allocated(timeTable                    )) call Dealloc_Array(timeTable                    )
       if (allocated(firstCrossingProbabilityTable)) call Dealloc_Array(firstCrossingProbabilityTable)
       varianceMaximum   =max(varianceMaximum,variance)
       varianceTableCount=int(varianceMaximum*dble(varianceTableNumberPerUnit))
       if (tableInitialized) then
          timeMinimum=min(timeMinimum,time)
          timeMaximum=max(timeMaximum,time)
       else
          timeMinimum=0.5d0*time
          timeMaximum=2.0d0*time
       end if
       timeTableCount=int(log10(timeMaximum/timeMinimum)*dble(timeTableNumberPerDecade))+1
       call Alloc_Array(varianceTable                ,[1+varianceTableCount]                 ,lowerBounds=[0  ])
       call Alloc_Array(timeTable                                           ,[timeTableCount]                  )
       call Alloc_Array(firstCrossingProbabilityTable,[1+varianceTableCount , timeTableCount],lowerBounds=[0,1])
       varianceTable    =Make_Range(0.0d0,varianceMaximum,varianceTableCount+1,rangeType=rangeTypeLinear)
       varianceTableStep=varianceTable(1)-varianceTable(0)
       timeTable        =Make_Range(timeMinimum,timeMaximum,timeTableCount,rangeType=rangeTypeLogarithmic)

       ! Loop through the table and solve for the first crossing distribution.
       call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
       do iTime=1,timeTableCount
          do i=0,varianceTableCount
             call Galacticus_Display_Counter(int(100.0d0*dble(i+(iTime-1)*varianceTableCount)/dble(varianceTableCount*timeTableCount)),i==0 .and. iTime==1,verbosityWorking)
             if      (i  > 2) then
                summedProbability=0.0d0
                do j=1,i-1
                   summedProbability=summedProbability+firstCrossingProbabilityTable(j,iTime)*(Delta(i,j,varianceTable(i),varianceTable(j),varianceTableStep,timeTable(iTime))+Delta(i,j+1,varianceTable(i),varianceTable(j+1),varianceTableStep,timeTable(iTime)))
                end do
                firstCrossingProbabilityTable(i,iTime)=(g_1(varianceTable(i),timeTable(iTime))+summedProbability)/(1.0d0-Delta(i,i,varianceTable(i),varianceTable(i),varianceTableStep,timeTable(iTime)))
             else if (i == 2) then
                firstCrossingProbabilityTable(i,iTime)=(g_1(varianceTable(i),timeTable(iTime))+firstCrossingProbabilityTable(i-1&
                     &,iTime)*(Delta(i,1,varianceTable(i),varianceTable(1),varianceTableStep,timeTable(iTime))+Delta(i,2,varianceTable(i),varianceTable(2),varianceTableStep,timeTable(iTime))))/(1.0d0-Delta(i,i,varianceTable(i),varianceTable(i),varianceTableStep,timeTable(iTime)))
             else if (i == 1) then
                firstCrossingProbabilityTable(i,iTime)= g_1(varianceTable(i),timeTable(iTime))/(1.0d0-Delta(1,1,varianceTable(1),varianceTable(1),varianceTableStep,timeTable(iTime)))
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
    !$omp end critical (Excursion_Sets_First_Crossing_Probability_Zhang_Hui_Init)

    ! Get interpolation in time.
    iTime    =Interpolate_Locate                 (timeTableCount      ,timeTable    ,interpolationAcceleratorTime    ,time    ,reset=interpolationResetTime    )
    hTime    =Interpolate_Linear_Generate_Factors(timeTableCount      ,timeTable    ,iTime    ,time    )

    ! Get interpolation in variance.
    iVariance=Interpolate_Locate                 (varianceTableCount+1,varianceTable,interpolationAcceleratorVariance,variance,reset=interpolationResetVariance)
    hVariance=Interpolate_Linear_Generate_Factors(varianceTableCount+1,varianceTable,iVariance,variance)

    ! Compute first crossing probability by interpolating.
    Excursion_Sets_First_Crossing_Probability_Zhang_Hui=0.0d0
    do jTime=0,1
       do jVariance=0,1
          Excursion_Sets_First_Crossing_Probability_Zhang_Hui=Excursion_Sets_First_Crossing_Probability_Zhang_Hui+hTime(jTime)*hVariance(jVariance)*firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function Excursion_Sets_First_Crossing_Probability_Zhang_Hui

  double precision function Excursion_Sets_First_Crossing_Rate_Zhang_Hui(variance,varianceProgenitor,time)
    !% Return the rate for excursion set first crossing.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: time, variance, varianceProgenitor

    call Galacticus_Error_Report('Excursion_Sets_First_Crossing_Rate_Zhang_Hui','barrier crossing rates are not implemented for this method [too slow]')
    return
  end function Excursion_Sets_First_Crossing_Rate_Zhang_Hui

  double precision function Excursion_Sets_Non_Crossing_Rate_Zhang_Hui(variance,time)
    !% Return the rate for excursion set non-crossing.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: time, variance

    call Galacticus_Error_Report('Excursion_Sets_Non_Crossing_Rate_Zhang_Hui','barrier non-crossing rates are not implemented for this method [too slow]')
    return
  end function Excursion_Sets_Non_Crossing_Rate_Zhang_Hui

end module Excursion_Sets_First_Crossing_Zhang_Hui
