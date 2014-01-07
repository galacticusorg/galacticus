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

function convergenceGelmanRubinConstructor(thresholdHatR,burnCount,testCount,outlierCountMaximum,outlierSignificance,outlierLogLikelihoodOffset)
  !% Constructor for ``GelmanRubin'' convergence class.
  implicit none
  type            (convergenceGelmanRubin)                :: convergenceGelmanRubinConstructor
  double precision                        , intent(in   ) :: thresholdHatR             , outlierSignificance, &
       &                                                     outlierLogLikelihoodOffset
  integer                                 , intent(in   ) :: burnCount                 , testCount          , &
       &                                                     outlierCountMaximum

  convergenceGelmanRubinConstructor%thresholdHatR             =thresholdHatR
  convergenceGelmanRubinConstructor%burnCount                 =burnCount
  convergenceGelmanRubinConstructor%testCount                 =testCount
  convergenceGelmanRubinConstructor%stepCount                 =testCount
  convergenceGelmanRubinConstructor%outlierCountMaximum       =outlierCountMaximum
  convergenceGelmanRubinConstructor%outlierSignificance       =outlierSignificance
  convergenceGelmanRubinConstructor%outlierLogLikelihoodOffset=outlierLogLikelihoodOffset
  convergenceGelmanRubinConstructor%converged                 =.false.
  return
end function convergenceGelmanRubinConstructor

logical function convergenceGelmanRubinIsConverged(self,simulationState,logLikelihood)
  !% Return whether the simulation is converged.
  use Constraints_State
  use MPI_Utilities
  use FGSL
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
       &                                                                   correctedHatR                , chainDeviation         , &
       &                                                                   chainDeviationMaximum        , currentStateMean       , &
       &                                                                   currentStateMeanSquared      , currentStateVariance
  integer                                 , allocatable  , dimension(:) :: chainDeviationIndex
  logical                                 , allocatable  , dimension(:) :: chainMask                    , hasOutlier
  logical                                                               :: newOutliersFound
  integer                                                               :: outlierCount                 , activeChainCount       , &
       &                                                                   i
  double precision                                                      :: grubbsCriticalValue          , tStatisticCriticalValue, &
       &                                                                   logLikelihoodMaximum

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
  ! Find outlier chains using a Grubb's outlier test.
  ! Initialize chain mask to all chains accepted.
  allocate(chainMask(mpiSelf%count()))
  chainMask=.true.
  ! Find the maximum likelihood over current states.
  logLikelihoodMaximum=mpiSelf%maxval(logLikelihood)
  ! Begin iterative removal of outliers.
  outlierCount    =0
  activeChainCount=mpiSelf%count()
  newOutliersFound=.true.
  do while (outlierCount < self%outlierCountMaximum .and. newOutliersFound)
     ! Find interchain mean and variance of current states.
     currentStateMean       =mpiSelf%average(simulationState%get()   ,chainMask)
     currentStateMeanSquared=mpiSelf%average(simulationState%get()**2,chainMask)
     currentStateVariance   = (currentStateMeanSquared-currentStateMean**2) &
          &                  *dble(activeChainCount  )                      &
          &                  /dble(activeChainCount-1)
     chainDeviation         = max(                                             &
          &                       abs(simulationState%get()-currentStateMean), &
          &                       abs(simulationState%get()-currentStateMean)  &
          &                      )                                             &
          &                  /sqrt(currentStateVariance)
     ! Find the maximum value of the deviation and which chain has this deviation.
     chainDeviationMaximum=mpiSelf%maxval(chainDeviation,chainMask)
     chainDeviationIndex  =mpiSelf%maxloc(chainDeviation,chainMask)+1
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
     hasOutlier=                                   &
          &       chainDeviationMaximum            &
          &      >                                 &
          &       grubbsCriticalValue              &
          &     .and.                              &
          &       logLikelihood                    &
          &      <                                 &
          &        logLikelihoodMaximum            &
          &       -                                &
          &        self%outlierLogLikelihoodOffset
     ! Mask chains with outliers.
     newOutliersFound=.false.
     do i=1,mpiSelf%count()
        if (hasOutlier(i).and.chainMask(chainDeviationIndex(i))) then
           chainMask(chainDeviationIndex(i))=.false.
           newOutliersFound                 =.true.
           outlierCount                     =outlierCount    +1
           activeChainCount                 =activeChainCount-1
        end if
     end do
     ! Proceed to next iteration.
     call mpiBarrier()
  end do
  ! Get the mean and variance of each parameter in our chain.
  chainMean    =simulationState%mean    ()
  chainVariance=simulationState%variance()
  ! Get the interchain mean.
  interchainMean       =mpiSelf%average(chainMean   ,chainMask)
  ! Get the interchain mean squared.
  interchainMeanSquared=mpiSelf%average(chainMean**2,chainMask)
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
  W                     = mpiSelf%average(chainVariance   ,chainMask)
  ! Compute variance of chain variances.
  varianceW             =(                                             &
       &                  +mpiSelf%average(chainVariance**2,chainMask) &
       &                  -                            W**2            &
       &                 )                                             &
       &                 *dble(activeChainCount  )                     &
       &                 /dble(activeChainCount-1)
  ! Find the covariance of chain variances and means.
  covarianceVarianceMean       = (mpiSelf%average(chainVariance*chainMean   ,chainMask)-W*interchainMean       ) &
       &                        *dble(activeChainCount  )                                                        &
       &                        /dble(activeChainCount-2)
  covarianceVarianceMeanSquared= (mpiSelf%average(chainVariance*chainMean**2,chainMask)-W*interchainMeanSquared) &
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
  correctedHatR=(degreesOfFreedom+3.0d0)*hatR/(degreesOfFreedom+1.0d0)
  ! Check for convergence.
  self%converged=all(correctedHatR < self%thresholdHatR)
  convergenceGelmanRubinIsConverged=self%converged
  return
end function convergenceGelmanRubinIsConverged
