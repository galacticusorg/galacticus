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

!% Defines methods for the ``GelmanRubin'' convergence class.

function convergenceGelmanRubinConstructor(thresholdHatR,burnCount,testCount,outlierCountMaximum,outlierSignificance,outlierLogLikelihoodOffset,reportCount)
  !% Constructor for ``GelmanRubin'' convergence class.
  use Memory_Management
  use MPI_Utilities
  type            (convergenceGelmanRubin)                :: convergenceGelmanRubinConstructor
  double precision                        , intent(in   ) :: thresholdHatR             , outlierSignificance, &
       &                                                     outlierLogLikelihoodOffset
  integer                                 , intent(in   ) :: burnCount                 , testCount          , &
       &                                                     outlierCountMaximum       , reportCount

  call Alloc_Array(convergenceGelmanRubinConstructor%chainMask,[mpiSelf%count()])
  convergenceGelmanRubinConstructor%thresholdHatR             =thresholdHatR
  convergenceGelmanRubinConstructor%burnCount                 =burnCount
  convergenceGelmanRubinConstructor%testCount                 =testCount
  convergenceGelmanRubinConstructor%stepCount                 =testCount
  convergenceGelmanRubinConstructor%outlierCountMaximum       =outlierCountMaximum
  convergenceGelmanRubinConstructor%outlierSignificance       =outlierSignificance
  convergenceGelmanRubinConstructor%outlierLogLikelihoodOffset=outlierLogLikelihoodOffset
  convergenceGelmanRubinConstructor%reportCount               =reportCount
  convergenceGelmanRubinConstructor%estimateCount             =0
  convergenceGelmanRubinConstructor%converged                 =.false.
  convergenceGelmanRubinConstructor%chainMask                 =.true.
 return
end function convergenceGelmanRubinConstructor

logical function convergenceGelmanRubinIsConverged(self,simulationState,logLikelihood)
  !% Return whether the simulation is converged.
  use Constraints_State
  use Memory_Management
  use MPI_Utilities
  use FGSL
  use ISO_Varying_String
  use Galacticus_Display
  use String_Handling
  implicit none
  class           (convergenceGelmanRubin), intent(inout)               :: self
  class           (state                 ), intent(in   ), optional     :: simulationState
  double precision                        , intent(in   ), optional     :: logLikelihood
  double precision                        , allocatable  , dimension(:) :: chainMean                    , chainVariance          , &
       &                                                                   interchainMean               , interchainMeanVariance , &
       &                                                                   B                            , W                      , &
       &                                                                   varianceW                    , covarianceVarianceMean , &
       &                                                                   covarianceVarianceMeanSquared, interchainMeanSquared  , &
       &                                                                   hatV                         , varianceHatV           , &
       &                                                                   degreesOfFreedom             , hatR                   , &
       &                                                                   currentStateVariance         , chainDeviation         , &
       &                                                                   chainDeviationMaximum        , currentStateMean       , &
       &                                                                   currentStateMeanSquared      , logLikelihoods
  integer                                 , allocatable  , dimension(:) :: chainDeviationIndex
  logical                                 , allocatable  , dimension(:) :: hasOutlier
  logical                                                               :: newOutliersFound
  integer                                                               :: outlierCount                 , activeChainCount       , &
       &                                                                   i                            , deviationMaximumChain
  double precision                                                      :: grubbsCriticalValue          , tStatisticCriticalValue, &
       &                                                                   logLikelihoodMaximum         , deviationMaximum
  type            (varying_string        )                              :: message
  character       (len=16                )                              :: label

  ! If no aguments were provided, return current convergence status without updating.
  if (.not.(present(simulationState).and.present(logLikelihood))) then
     convergenceGelmanRubinIsConverged=self%converged
     return
  end if
  ! Check if we have enough steps to check convergence.
  if (simulationState%count() < self%burnCount) then
     convergenceGelmanRubinIsConverged=.false.
     return
  end if
  ! Decide if we should check convergence.
  self%stepCount=self%stepCount+1
  if (self%stepCount < self%testCount) then
     convergenceGelmanRubinIsConverged=self%converged
     return
  end if
  self%stepCount=0
  ! Allocate chain mask and Rhat arrays.
  if (.not.allocated(self%correctedHatR)) call Alloc_Array(self%correctedHatR,[simulationState%dimension()])
  ! Find outlier chains using a Grubb's outlier test.
  ! Initialize chain mask to all chains accepted.
  self%chainMask=.true.
  ! Find the log likelihoods of all states.
  logLikelihoods=mpiSelf%gather(logLikelihood)
  ! Find the maximum likelihood over current states.
  logLikelihoodMaximum=mpiSelf%maxval(logLikelihood)
  ! Begin iterative removal of outliers.
  outlierCount    =0
  activeChainCount=mpiSelf%count()
  newOutliersFound=.true.
  do while (outlierCount < self%outlierCountMaximum .and. newOutliersFound)
     ! Find interchain mean and variance of current states.
     currentStateMean       =mpiSelf%average(simulationState%get()   ,self%chainMask)
     currentStateMeanSquared=mpiSelf%average(simulationState%get()**2,self%chainMask)
     currentStateVariance   = (currentStateMeanSquared-currentStateMean**2) &
          &                  *dble(activeChainCount  )                      &
          &                  /dble(activeChainCount-1)
     chainDeviation         = abs(simulationState%get()-currentStateMean)   &
          &                  /sqrt(currentStateVariance)
     ! Find the maximum value of the deviation and which chain has this deviation.
     chainDeviationMaximum=mpiSelf%maxval(chainDeviation,self%chainMask)
     chainDeviationIndex  =mpiSelf%maxloc(chainDeviation,self%chainMask)+1
     ! Evaluate the critical value for outlier rejection.
     tStatisticCriticalValue=                                 &
          &   FGSL_CDF_tDist_Qinv(                            &
          &                       (                           &
          &                        +1.0d0                     &
          &                        -self%outlierSignificance  &
          &                       )                           &
          &                       /2.0d0                      &
          &                       /dble(activeChainCount  ) , &
          &                        dble(activeChainCount-2)   &
          &                      )
     grubbsCriticalValue= dble(                              &
          &                         activeChainCount-1       &
          &                   )                              &
          &              /sqrt(                              &
          &                    dble(activeChainCount)        &
          &                   )                              &
          &              *sqrt(                              &
          &                       tStatisticCriticalValue**2 &
          &                    /(                            &
          &                      +  activeChainCount         &
          &                      -2.0d0                      &
          &                      +tStatisticCriticalValue**2 &
          &                     )                            &
          &                   )
     ! Test for outliers using the Grubb's statistic.
     hasOutlier=                                         &
          &        chainDeviationMaximum                 &
          &      >                                       &
          &        grubbsCriticalValue                   &
          &     .and.                                    &
          &        logLikelihoods(chainDeviationIndex+1) &
          &      <                                       &
          &        logLikelihoodMaximum                  &
          &       -                                      &
          &        self%outlierLogLikelihoodOffset
     ! Mask chains with outliers.
     newOutliersFound=any(hasOutlier)
     if (newOutliersFound) then
        deviationMaximum     =0.0d0
        deviationMaximumChain=-1
        do i=1,size(hasOutlier)
           if (hasOutlier(i).and.self%chainMask(chainDeviationIndex(i)).and.chainDeviationMaximum(i) > deviationMaximum) then
              deviationMaximum     =chainDeviationMaximum(i)
              deviationMaximumChain=chainDeviationIndex  (i)
           end if
        end do
        self%chainMask(deviationMaximumChain)=.false.
        outlierCount                         =outlierCount    +1
        activeChainCount                     =activeChainCount-1
     end if
     ! Proceed to next iteration.
     call mpiBarrier()
  end do
  ! Get the mean and variance of each parameter in our chain.
  chainMean    =simulationState%mean    ()
  chainVariance=simulationState%variance()
  ! Get the interchain mean.
  interchainMean       =mpiSelf%average(chainMean   ,self%chainMask)
  ! Get the interchain mean squared.
  interchainMeanSquared=mpiSelf%average(chainMean**2,self%chainMask)
  ! Get the interchain variance.
  interchainMeanVariance=(                          &
       &                  +interchainMeanSquared    &
       &                  -interchainMean       **2 &
       &                 )                          &
       &                 *dble(activeChainCount  )  &
       &                 /dble(activeChainCount-1)
  ! Compute B from Brooks & Gelman (section 1.2).
  B                     =dble(simulationState%count())*interchainMeanVariance
  ! Compute W from Brooks & Gelman (section 1.2).
  W                     = mpiSelf%average(chainVariance   ,self%chainMask)
  ! Check for zero chain variance.
  if (any(W <= 0.0d0)) then
     self%estimateCount=self%estimateCount+1
     if (mpiSelf%isMaster() .and. mod(self%estimateCount,self%reportCount) == 0) then
        write (label,'(i16)') simulationState%count()
        message="Gelman-Rubin statistic at "//trim(adjustl(label))
        write (label,'(f6.2)') minval(self%correctedHatR)
        message=message//" steps cannot be computed (zero variances)"
        call Galacticus_Display_Message(message)
     end if
     convergenceGelmanRubinIsConverged=.false.
  end if
  ! Compute variance of chain variances.
  varianceW             =(                                                  &
       &                  +mpiSelf%average(chainVariance**2,self%chainMask) &
       &                  -                            W**2                 &
       &                 )                                                  &
       &                 *dble(activeChainCount  )                          &
       &                 /dble(activeChainCount-1)
  ! Find the covariance of chain variances and means.
  covarianceVarianceMean       = (mpiSelf%average(chainVariance*chainMean   ,self%chainMask)-W*interchainMean       ) &
       &                        *dble(activeChainCount  )                                                        &
       &                        /dble(activeChainCount-2)
  covarianceVarianceMeanSquared= (mpiSelf%average(chainVariance*chainMean**2,self%chainMask)-W*interchainMeanSquared) &
       &                        *dble(activeChainCount  )                                                        &
       &                        /dble(activeChainCount-2)
  ! Estimate Vhat and variance in Vhat from Gelman & Rubin (eqn. 4).
  hatV=  +W*dble(simulationState%count()-1)/dble(simulationState%count()) &
       & +B/dble(simulationState%count()  )                               &
       & +B/dble(simulationState%count()  )/dble(mpiSelf        %count())
  varianceHatV=                                      &
       &        (                                    &
       &          dble(simulationState%count()-1)    &
       &         /dble(simulationState%count()  )    &
       &        )**2                                 &
       &       *varianceW                            &
       &       /  dble(activeChainCount         )    &
       &       +(                                    &
       &          dble(activeChainCount       +1)    &
       &         /dble(activeChainCount         )    &
       &         /dble(simulationState%count()  )    &
       &        )**2                                 &
       &       *(                                    &
       &         2.0d0                               &
       &         /dble(activeChainCount       -1)    &
       &        )                                    &
       &       *B**2                                 &
       &       +2.0d0                                &
       &       *  dble(activeChainCount       +1)    &
       &       *  dble(simulationState%count()-1)    &
       &       /  dble(activeChainCount         )    &
       &       /  dble(simulationState%count()  )**2 &
       &       *(                                    &
       &          dble(simulationState%count()  )    &
       &         /dble(activeChainCount         )    &
       &        )                                    &
       &       *(                                    &
       &         +covarianceVarianceMeanSquared      &
       &         -2.0d0                              &
       &         *interchainMean                     &
       &         *covarianceVarianceMean             &
       &        )
  ! Compute degrees of freedom.
  degreesOfFreedom=2.0d0*hatV**2/varianceHatV
  ! Compute R-hat statistic.
  hatR=                                               &
       & sqrt(                                        &
       &      max(                                    &
       &           0.0d0                            , &
       &           (                                  &
       &             dble(mpiSelf         %count()+1) &
       &            /dble(mpiSelf         %count()  ) &
       &           )                                  &
       &          *(                                  &
       &            (                                 &
       &              dble(simulationState%count()-1) &
       &             /dble(simulationState%count()  ) &
       &            )                                 &
       &            *W                                &
       &            +B                                &
       &            / dble(simulationState%count()  ) &
       &           )                                  &
       &           /W                                 &
       &           -  dble(simulationState%count()-1) &
       &           /  dble(simulationState%count()  ) &
       &           /  dble(activeChainCount         ) &
       &         )                                    &
       &     )
  ! Compute corrected R-hat statistic.
  self%correctedHatR=(degreesOfFreedom+3.0d0)*hatR/(degreesOfFreedom+1.0d0)
  ! Check for convergence.
  self%converged=all(self%correctedHatR < self%thresholdHatR)
  convergenceGelmanRubinIsConverged=self%converged
  ! Report.
  self%estimateCount=self%estimateCount+1
  if (mpiSelf%isMaster() .and. mod(self%estimateCount,self%reportCount) == 0 .and. .not.self%converged) then
     write (label,'(i16)') simulationState%count()
     message="Gelman-Rubin ̂R at "//trim(adjustl(label))
     write (label,'(f6.2)') minval(self%correctedHatR)
     message=message//" steps min/max="//trim(adjustl(label))//"/"
     write (label,'(f6.2)') maxval(self%correctedHatR)
     message=message//trim(adjustl(label))//")"
     call Galacticus_Display_Message(message)
     if (activeChainCount < mpiSelf%count()) then
        message='outlier chains:'
        label=''
        do i=1,mpiSelf%count()
           if (.not.self%chainMask(i)) then
              message=message//trim(label)//' '//i
              label=","
           end if
        end do
        call Galacticus_Display_Message(message)
     else
        call Galacticus_Display_Message('no outlier chains')
     end if
  end if
  return
end function convergenceGelmanRubinIsConverged

subroutine convergenceGelmanRubinReset(self)
  !% Reset the convergence object.
  implicit none
  class(convergenceGelmanRubin), intent(inout) :: self

  self%converged=.false.
  return
end subroutine convergenceGelmanRubinReset

subroutine convergenceGelmanRubinLogReport(self,fileUnit)
  !% Write a convergence report to the given {\tt fileUnit}.
  use MPI_Utilities
  implicit none
  class    (convergenceGelmanRubin), intent(inout) :: self
  integer                          , intent(in   ) :: fileUnit
  character(len=25                )                :: label
  
  write (label,'(a,i4.4,a)') '(a,f5.2,',size(self%correctedHatR)-1,'(", ",f5.2))'
  write (fileUnit,label) 'Gelman-Rubin covergence ̂R: ',self%correctedHatR
  write (label,'(a,i4.4,a)') '(a,l1,',size(self%correctedHatR)-1,'(", ",l1))'
  write (fileUnit,label) 'Gelman-Rubin chain mask  : ',self%chainMask
  return
end subroutine convergenceGelmanRubinLogReport

logical function convergenceGelmanRubinChainIsOutlier(self,chainIndex)
  !% Return true if the specified chain is deemed to be an outlier.
  implicit none
  class  (convergenceGelmanRubin), intent(inout) :: self
  integer                        , intent(in   ) :: chainIndex

  convergenceGelmanRubinChainIsOutlier=.not.self%chainMask(chainIndex+1)
  return
end function convergenceGelmanRubinChainIsOutlier
