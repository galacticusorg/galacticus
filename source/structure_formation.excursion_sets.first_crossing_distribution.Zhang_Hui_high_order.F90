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

  !+    Contributions to this file made by: Arya Farahi, Andrew Benson.

  !% Contains a module which implements a excursion set first crossing statistics class utilizing a higher order generalization of
  !% the algorithm of \cite{zhang_random_2006}.

  !# <excursionSetFirstCrossing name="excursionSetFirstCrossingZhangHuiHighOrder">
  !#  <description>An excursion set first crossing statistics class utilizing a higher order generalization of the algorithm of \cite{zhang_random_2006}.</description>
  !# </excursionSetFirstCrossing>
  type, extends(excursionSetFirstCrossingZhangHui) :: excursionSetFirstCrossingZhangHuiHighOrder
     !% An excursion set first crossing statistics class utilizing a higher order generalization of the algorithm of
     !% \cite{zhang_random_2006}.
     private
     double precision :: timeMinimumPrevious, timeMaximumPrevious
   contains
     !@ <objectMethods>
     !@   <object>excursionSetFirstCrossingZhangHuiHighOrder</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Initialize the high order \cite{zhang_random_2006} class.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize  => zhangHuiHighOrderInitialize
     procedure :: probability => zhangHuiHighOrderProbability
  end type excursionSetFirstCrossingZhangHuiHighOrder

  interface excursionSetFirstCrossingZhangHuiHighOrder
     !% Constructors for the \cite{zhang_random_2006} excursion set barrier class.
     module procedure zhangHuiHighOrderConstructorParameters
     module procedure zhangHuiHighOrderConstructorInternal
  end interface excursionSetFirstCrossingZhangHuiHighOrder

contains

  function zhangHuiHighOrderConstructorParameters(parameters) result(self)
    !% Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (excursionSetFirstCrossingZhangHuiHighOrder)                :: self
    type (inputParameters                           ), intent(inout) :: parameters

    self%excursionSetFirstCrossingZhangHui=excursionSetFirstCrossingZhangHui(parameters)
    call self%initialize()
    return
  end function zhangHuiHighOrderConstructorParameters

  function zhangHuiHighOrderConstructorInternal(excursionSetBarrier_) result(self)
    !% Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    implicit none
    type (excursionSetFirstCrossingZhangHuiHighOrder)                        :: self
    class(excursionSetBarrierClass                  ), intent(in   ), target :: excursionSetBarrier_

    self%excursionSetFirstCrossingZhangHui=excursionSetFirstCrossingZhangHui(excursionSetBarrier_)
    call self%initialize()
    return
  end function zhangHuiHighOrderConstructorInternal

  subroutine zhangHuiHighOrderInitialize(self)
    !% Initialize the high order \cite{zhang_random_2006} excursion set first crossing class.
    implicit none
    class(excursionSetFirstCrossingZhangHuiHighOrder), intent(inout) :: self

    self%timeMinimumPrevious=+huge(0.0d0)
    self%timeMaximumPrevious=-huge(0.0d0)
    return
  end subroutine zhangHuiHighOrderInitialize

  double precision function zhangHuiHighOrderProbability(self,variance,time,node)
    !% Return the excursion set barrier at the given variance and time.
    use            :: Galacticus_Display, only : Galacticus_Display_Counter, Galacticus_Display_Counter_Clear, Galacticus_Display_Indent, Galacticus_Display_Unindent, &
          &                                      verbosityWorking
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Memory_Management , only : allocateArray             , deallocateArray
    use            :: Numerical_Ranges  , only : Make_Range                , rangeTypeLinear                 , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingZhangHuiHighOrder), intent(inout)                 :: self
    double precision                                            , intent(in   )                 :: variance                             , time
    type            (treeNode                                  ), intent(inout)                 :: node
    double precision                                                           , dimension(0:1) :: hTime                                , hVariance
    double precision                                            , allocatable  , dimension(:,:) :: firstCrossingProbabilityTablePrevious
    logical                                                                                     :: makeTable                            , tableIsExtendable
    integer         (c_size_t                                  )                                :: iVariance                            , iTime
    integer                                                                                     :: i                                    , j                         , &
         &                                                                                         jTime                                , jVariance                 , &
         &                                                                                         k                                    , varianceTableCountPrevious
    double precision                                                                            :: g2Average                            , summation0                , &
         &                                                                                         summation1                           , summation2

    ! Determine if we need to make the table.
    makeTable=           .not. self%tableInitialized  &
         &    .or.                                    &
         &     (variance >     self%varianceMaximum ) &
         &    .or.                                    &
         &     (time     <     self%timeMinimum     ) &
         &    .or.                                    &
         &     (time     >     self%timeMaximum     )
    if (makeTable) then
       ! Compute the table extent.
       if (self%tableInitialized) then
          self%timeMinimum=min(self%timeMinimum,time)
          self%timeMaximum=max(self%timeMaximum,time)
       else
          self%timeMinimum=0.5d0*time
          self%timeMaximum=2.0d0*time
       end if
       self%varianceMaximum   =max(self%varianceMaximum,variance)
       ! Make a copy of the current table if possible.
       tableIsExtendable=(self%timeMinimum == self%timeMinimumPrevious .and. self%timeMaximum == self%timeMaximumPrevious .and. self%tableInitialized)
       if (tableIsExtendable) then
          ! We have a pre-existing table with the same range of times. Therefore, we can simply extend this table. So make a copy
          ! of the existing table.
          call Move_Alloc(self%firstCrossingProbabilityTable,firstCrossingProbabilityTablePrevious)
          varianceTableCountPrevious=self%varianceTableCount
       else
          varianceTableCountPrevious=-1
       end if
       self%timeMinimumPrevious=self%timeMinimum
       self%timeMaximumPrevious=self%timeMaximum
       self%timeTableCount     =int(log10(self%timeMaximum/self%timeMinimum)*dble(timeTableNumberPerDecade))+1
       self%varianceTableCount =int(self%varianceMaximum*dble(varianceTableNumberPerUnit))
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
       if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
       call allocateArray(self%varianceTable                ,[1+self%varianceTableCount]                      ,lowerBounds=[0  ])
       call allocateArray(self%timeTable                                                ,[self%timeTableCount]                  )
       call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount , self%timeTableCount],lowerBounds=[0,1])
       self%timeTable        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%timeTableCount      ,rangeType=rangeTypeLogarithmic)
       self%varianceTable    =Make_Range(0.0d0           ,self%varianceMaximum,self%varianceTableCount+1,rangeType=rangeTypeLinear     )
       self%varianceTableStep=+self%varianceTable(4) &
            &                 -self%varianceTable(0)
       if (tableIsExtendable) then
          self%firstCrossingProbabilityTable(0:varianceTableCountPrevious,:)=firstCrossingProbabilityTablePrevious
          call deallocateArray(firstCrossingProbabilityTablePrevious)
       end if
       ! Loop through the table and solve for the first crossing distribution.
       call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
       do iTime=1,self%timeTableCount
          do i=varianceTableCountPrevious+1,self%varianceTableCount
             call Galacticus_Display_Counter(                                                                                      &
                  &                           int(                                                                                 &
                  &                                100.0d0                                                                         &
                  &                               *dble(                                                                           &
                  &                                                         (i                      -varianceTableCountPrevious-1) &
                  &                                     +(iTime-1)         *(self%varianceTableCount-varianceTableCountPrevious-1) &
                  &                                    )                                                                           &
                  &                               /dble(                                                                           &
                  &                                     self%timeTableCount*(self%varianceTableCount-varianceTableCountPrevious-1) &
                  &                                    )                                                                           &
                  &                              )                                                                                 &
                  &                          ,i==varianceTableCountPrevious+1 .and. iTime==1,verbosityWorking                      &
                  &                         )

             if (i  > 3) then
                summation0=0.0d0
                summation1=0.0d0
                summation2=0.0d0
                k=mod(i,4)
                if      (k == 3) then
                   summation0=+(                                                                                                                                   &
                        &       +3.0d0*self%firstCrossingProbabilityTable(1,iTime)*self%g2(self%varianceTable(i),self%varianceTable(1),self%timeTable(iTime),node) &
                        &       +3.0d0*self%firstCrossingProbabilityTable(2,iTime)*self%g2(self%varianceTable(i),self%varianceTable(2),self%timeTable(iTime),node) &
                        &       +      self%firstCrossingProbabilityTable(3,iTime)*self%g2(self%varianceTable(i),self%varianceTable(3),self%timeTable(iTime),node) &
                        &      )                                                                                                                                   &
                        &     *3.0d0*self%varianceTableStep/32.0d0
                else if (k == 2) then
                   summation0=+(                                                                                                                                   &
                        &       +4.0d0*self%firstCrossingProbabilityTable(1,iTime)*self%g2(self%varianceTable(i),self%varianceTable(1),self%timeTable(iTime),node) &
                        &       +      self%firstCrossingProbabilityTable(2,iTime)*self%g2(self%varianceTable(i),self%varianceTable(2),self%timeTable(iTime),node) &
                        &      )                                                                                                                                   &
                        &     *self%varianceTableStep/12.0d0
                else if (k == 1) then
                   summation0=+       self%firstCrossingProbabilityTable(1,iTime)*self%g2(self%varianceTable(i),self%varianceTable(1),self%timeTable(iTime),node) &
                        &     *self%varianceTableStep/8.0d0
                else
                   summation0=0.0d0
                end if
                do j=k,i-8,4
                   summation1=+summation1                                                                                                                             &
                        &     + 7.0d0*self%firstCrossingProbabilityTable(j  ,iTime)*self%g2(self%varianceTable(i),self%varianceTable(j  ),self%timeTable(iTime),node) &
                        &     +32.0d0*self%firstCrossingProbabilityTable(j+1,iTime)*self%g2(self%varianceTable(i),self%varianceTable(j+1),self%timeTable(iTime),node) &
                        &     +12.0d0*self%firstCrossingProbabilityTable(j+2,iTime)*self%g2(self%varianceTable(i),self%varianceTable(j+2),self%timeTable(iTime),node) &
                        &     +32.0d0*self%firstCrossingProbabilityTable(j+3,iTime)*self%g2(self%varianceTable(i),self%varianceTable(j+3),self%timeTable(iTime),node) &
                        &     + 7.0d0*self%firstCrossingProbabilityTable(j+4,iTime)*self%g2(self%varianceTable(i),self%varianceTable(j+4),self%timeTable(iTime),node)
                end do
                summation1=+summation1             &
                     &     *self%varianceTableStep &
                     &     /90.0d0
                g2Average=self%g2Integrated(self%varianceTable(i),self%varianceTableStep,self%timeTable(iTime),node)
                summation2=+g2Average                                              &
                     &     *(                                                      &
                     &       + 7.0d0*self%firstCrossingProbabilityTable(i-4,iTime) &
                     &       +32.0d0*self%firstCrossingProbabilityTable(i-3,iTime) &
                     &       +12.0d0*self%firstCrossingProbabilityTable(i-2,iTime) &
                     &       +32.0d0*self%firstCrossingProbabilityTable(i-1,iTime) &
                     &      )                                                      &
                     &     /90.0d0
                self%firstCrossingProbabilityTable(i,iTime)=+(                                                           &
                     &                                        +self%g1(self%varianceTable(i),self%timeTable(iTime),node) &
                     &                                        +summation0                                                &
                     &                                        +summation1                                                &
                     &                                        +summation2                                                &
                     &                                       )                                                           &
                     &                                      /(                                                           &
                     &                                        +1.0d0                                                     &
                     &                                        -7.0d0                                                     &
                     &                                        *g2Average                                                 &
                     &                                        /90.0d0                                                    &
                     &                                       )
             else if (i == 3) then
                g2Average=self%g2Integrated(self%varianceTable(i),0.75d0*self%varianceTableStep,self%timeTable(iTime),node)
                self%firstCrossingProbabilityTable(i,iTime)=+(                                                           &
                     &                                        +self%g1(self%varianceTable(i),self%timeTable(iTime),node) &
                     &                                        +g2Average                                                 &
                     &                                        *3.0d0                                                     &
                     &                                        *(                                                         &
                     &                                          +self%firstCrossingProbabilityTable(1,iTime)             &
                     &                                          +self%firstCrossingProbabilityTable(2,iTime)             &
                     &                                         )                                                         &
                     &                                        /8.0d0                                                     &
                     &                                       )                                                           &
                     &                                      /(                                                           &
                     &                                        +1.0d0                                                     &
                     &                                        -g2Average                                                 &
                     &                                        /8.0d0                                                     &
                     &                                       )
             else if (i == 2) then
                g2Average=self%g2Integrated(self%varianceTable(i),0.5d0*self%varianceTableStep,self%timeTable(iTime),node)
                self%firstCrossingProbabilityTable(i,iTime)=+(                                                           &
                     &                                        +self%g1(self%varianceTable(i),self%timeTable(iTime),node) &
                     &                                        +g2Average                                                 &
                     &                                        *4.0d0                                                     &
                     &                                        *self%firstCrossingProbabilityTable(1,iTime)               &
                     &                                        /6.0d0                                                     &
                     &                                       )                                                           &
                     &                                      /(                                                           &
                     &                                        +1.0d0                                                     &
                     &                                        -g2Average                                                 &
                     &                                        /6.0d0                                                     &
                     &                                       )
             else if (i == 1) then
                g2Average=self%g2Integrated(self%varianceTable(i),0.25d0*self%varianceTableStep,self%timeTable(iTime),node)
                self%firstCrossingProbabilityTable(i,iTime)=+self%g1(self%varianceTable(i),self%timeTable(iTime),node) &
                     &                                      /(                                                         &
                     &                                        +1.0d0                                                   &
                     &                                        -g2Average                                               &
                     &                                        /2.0d0                                                   &
                     &                                       )
             else if (i == 0) then
                self%firstCrossingProbabilityTable(i,iTime)=+0.0d0
             end if
          end do
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       call Galacticus_Display_Unindent("done",verbosityWorking)
       ! Build the interpolators.
       if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
       if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
       allocate(self%interpolatorVariance)
       allocate(self%interpolatorTime    )
       self%interpolatorVariance=interpolator(self%varianceTable)
       self%interpolatorTime    =interpolator(self%timeTable    )
       ! Record that the table is now built.
       self%tableInitialized=.true.
    end if
    ! Get interpolating factors.
    call self%interpolatorTime    %linearFactors(time    ,iTime    ,hTime    )
    call self%interpolatorVariance%linearFactors(variance,iVariance,hVariance)
    ! Compute first crossing probability by interpolating.
    zhangHuiHighOrderProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          zhangHuiHighOrderProbability=zhangHuiHighOrderProbability+hTime(jTime)*hVariance(jVariance)*self%firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function zhangHuiHighOrderProbability
