!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which computes correlation statistics from point distributions.
!!}

module Statistics_Points_Correlations
  !!{
  Compute correlation statistics from point distributions.
  !!}
  private
  public :: Statistics_Points_Correlation

contains

  subroutine Statistics_Points_Correlation(dataPosition,randomPosition,separationMinimum,separationMaximum,separationCount,separation,correlation,projected,radialSeparationMaximum,halfIntegral)
    !!{
    Compute the correlation function from a set of points.
    !!}
    use :: Display          , only : displayCounter  , displayCounterClear , displayIndent, displayUnindent
    use :: Error            , only : Error_Report
    use :: Kind_Numbers     , only : kind_int8
    use :: Nearest_Neighbors, only : nearestNeighbors
    use :: Numerical_Ranges , only : Make_Range      , rangeTypeLogarithmic
    use :: Vectors          , only : Vector_Magnitude
    implicit none
    double precision                  , intent(in   ), dimension(:,:), target      :: dataPosition                    , randomPosition
    double precision                  , intent(in   )                              :: separationMinimum               , separationMaximum
    double precision                  , intent(  out), dimension(:  ), allocatable :: correlation                     , separation
    logical                           , intent(in   ), optional                    :: projected                       , halfIntegral
    double precision                  , intent(in   ), optional                    :: radialSeparationMaximum
    integer                                                                        :: separationCount
    double precision                  , allocatable  , dimension(:  )              :: separationSquaredMaximum        , neighborDistance            , &
         &                                                                            separationSquaredMinimum
    integer                           , allocatable  , dimension(:  )              :: neighborIndex
    integer         (kind=kind_int8  ), allocatable  , dimension(:,:), target      :: dataDataCount                   , randomRandomCount           , &
         &                                                                            dataRandomCount                 , pairCount
    double precision                  , allocatable  , dimension(:,:)              :: dataData                        , randomRandom                , &
         &                                                                            dataRandom                      , correlationTmp
    integer                           , parameter                                  :: radialBinsPerDecade       =30
    double precision                  , parameter                                  :: toleranceZero             =0.0d0
    double precision                  , parameter                                  :: radialMinimumRatio        =0.3d0
    class           (nearestNeighbors), allocatable  , save                        :: neighborFinder
    !$omp threadprivate(neighborFinder)
    double precision                                 , dimension(:,:), pointer     :: fromPoints                      , toPoints
    integer                                                                        :: iPoint                          , jPoint                      , &
         &                                                                            iBin                            , neighborCount               , &
         &                                                                            pointCount                      , iPass                       , &
         &                                                                            radialBinCount                  , jBin
    double precision                                                               :: separationLimit                 , separationLogarithmicMinimum, &
         &                                                                            separationLogarithmicStepInverse, radialPositionFrom          , &
         &                                                                            radialPositionTo                , separationRadial            , &
         &                                                                            separationProjected             , cosTheta                    , &
         &                                                                            radialLogarithmicStepInverse    , radialLogarithmicMinimum    , &
         &                                                                            radialMinimum
    logical                                                                        :: projectedActual

    ! Determine if we are to compute regular or projected correlation function.
    projectedActual=.false.
    if (present(projected)) projectedActual=projected
    ! Number of radial separation bins.
    if (projectedActual) then
       radialMinimum               =radialMinimumRatio*separationMinimum
       radialBinCount              =int(dble(radialBinsPerDecade)*log10(radialSeparationMaximum/radialMinimum))+1
       radialLogarithmicStepInverse=1.0d0/(log(radialSeparationMaximum/radialMinimum)/dble(radialBinCount-1))
       radialLogarithmicMinimum    =       log(                        radialMinimum)-0.5d0/radialLogarithmicStepInverse
    else
       radialMinimum               =-1.0d0
       radialBinCount              = 1
       radialLogarithmicStepInverse= 0.0d0
       radialLogarithmicMinimum    = 0.0d0
    end if
    ! Allocate arrays.
    if (allocated(separation )) deallocate(separation )
    if (allocated(correlation)) deallocate(correlation)
    allocate(separation              (separationCount               ))
    allocate(correlation             (separationCount               ))
    allocate(separationSquaredMinimum(separationCount               ))
    allocate(separationSquaredMaximum(separationCount               ))
    allocate(dataDataCount           (separationCount,radialBinCount))
    allocate(randomRandomCount       (separationCount,radialBinCount))
    allocate(dataRandomCount         (separationCount,radialBinCount))
    allocate(pairCount               (separationCount,radialBinCount))
    allocate(dataData                (separationCount,radialBinCount))
    allocate(randomRandom            (separationCount,radialBinCount))
    allocate(dataRandom              (separationCount,radialBinCount))
    allocate(correlationTmp          (separationCount,radialBinCount))
    ! Initialize counts.
    dataDataCount    =0
    randomRandomCount=0
    dataRandomCount  =0
    ! Generate the array of separations.
    separation=Make_Range(separationMinimum,separationMaximum,separationCount,rangeType=rangeTypeLogarithmic)
    do iBin=1,separationCount
       separationSquaredMinimum(iBin)=exp(2.0d0*(log(separation(iBin))-0.5d0*log(separation(2)/separation(1))))
       separationSquaredMaximum(iBin)=exp(2.0d0*(log(separation(iBin))+0.5d0*log(separation(2)/separation(1))))
    end do
    ! Find the minimum and maximum separations, and the logarithmic step size.
    separationLimit                 =separationMaximum*sqrt(separation(2)/separation(1))
    if (projectedActual) then
       if (.not.present(radialSeparationMaximum)) call Error_Report('maximum radial separation required for projected correlation functions'//{introspection:location})
       separationLimit=sqrt(separationLimit**2+radialSeparationMaximum**2)
    end if
    separationLogarithmicStepInverse=1.0d0/log(separation(2)/separation(1))
    separationLogarithmicMinimum    =log(separationMinimum)-0.5d0/separationLogarithmicStepInverse
    ! Iterate over points, counting neighbors.
    do iPass=1,3
       ! Set pointers to set of points to use in this pass.
       select case (iPass)
       case (1)
          ! Data-data pairs.
          call displayIndent('Counting data-data pairs')
          fromPoints => dataPosition
          toPoints   => dataPosition
      case (2)
          ! Data-random pairs.
          call displayIndent('Counting data-random pairs')
          fromPoints => dataPosition
          toPoints   => randomPosition
      case (3)
          ! Random-random pairs.
          call displayIndent('Counting random-random pairs')
          fromPoints => randomPosition
          toPoints   => randomPosition
       end select
       pointCount=0
       !$omp parallel private(iBin,jBin,iPoint,jPoint,neighborCount,neighborIndex,neighborDistance,separationProjected,radialPositionFrom,radialPositionTo,separationRadial,cosTheta)
       ! Construct nearest neighbor object.
       allocate(nearestNeighbors :: neighborFinder)
       select type (neighborFinder)
       type is (nearestNeighbors)
          neighborFinder=nearestNeighbors(transpose(toPoints))
       end select
       pairCount=0
       !$omp do reduction(+:pairCount)
       do iPoint=1,size(fromPoints,dim=2)
          !$omp atomic
          pointCount=pointCount+1
          call displayCounter(int(100.0d0*dble(pointCount-1)/dble(size(fromPoints,dim=2))),isNew=(iPoint==1))
          ! Get neighbors within the maximum separation.
          call neighborFinder%searchFixedRadius(fromPoints(:,iPoint),separationLimit,toleranceZero,neighborCount,neighborIndex,neighborDistance)
          select case (projectedActual)
          case (.false.)
             ! For regular correlation function, simply count pairs into bins based on separation. Points are assumed to be in
             ! order of increasing separation.
             iBin=1
             do jPoint=1,neighborCount
                do while (neighborDistance(jPoint) > separationSquaredMaximum(iBin))
                   iBin=iBin+1
                   if (iBin > separationCount) exit
                end do
                if (iBin > separationCount) exit
                if (neighborDistance(jPoint) > separationSquaredMinimum(iBin)) pairCount(iBin,1)=pairCount(iBin,1)+1
             end do
          case (.true.)
             ! For projected correlation function, find the radial and projected separations.
             do jPoint=1,neighborCount
                radialPositionFrom=Vector_Magnitude(fromPoints(:,              iPoint ))
                radialPositionTo  =Vector_Magnitude(toPoints  (:,neighborIndex(jPoint)))
                separationRadial=abs(radialPositionFrom-radialPositionTo)
                if (separationRadial < radialSeparationMaximum) then
                   cosTheta=min(                                                     &
                        &           +1.0d0,                                          &
                        &       max(                                                 &
                        &           -1.0d0,                                          &
                        &           Dot_Product(                                     &
                        &                       fromPoints(:,              iPoint ), &
                        &                       toPoints  (:,neighborIndex(jPoint))  &
                        &                      )                                     &
                        &           /radialPositionFrom                              &
                        &           /radialPositionTo                                &
                        &          )                                                 &
                        &      )
                   separationProjected=(radialPositionFrom+radialPositionTo)*tan(0.5d0*acos(cosTheta))
                   if (separationProjected > 0.0d0 .and. separationRadial > 0.0d0) then
                      iBin=int((log(separationProjected)-separationLogarithmicMinimum)*separationLogarithmicStepInverse)+1
                      jBin=int((log(separationRadial   )-    radialLogarithmicMinimum)*    radialLogarithmicStepInverse)+1
                      if     (                                        &
                           &   iBin > 0 .and. iBin <= separationCount &
                           &  .and.                                   &
                           &   jBin > 0 .and. jBin <=  radialBinCount &
                           & ) pairCount  (iBin,jBin)=                &
                           &    +pairCount(iBin,jBin)                 &
                           &    +1
                   end if

                end if
             end do
          end select
       end do
       !$omp end do
       !$omp barrier
       deallocate(neighborFinder)
       ! Store pair counts.
       select case (iPass)
       case (1)
          ! Data-data pairs.
          dataDataCount    =pairCount
       case (2)
          ! Data-random pairs.
          dataRandomCount  =pairCount
       case (3)
          ! Random-random pairs.
          randomRandomCount=pairCount
       end select
       !$omp end parallel
       call displayCounterClear()
       call displayUnindent('done')
    end do
    ! Normalize counts.
    dataData    =dble(dataDataCount    )/dble(size(dataPosition,dim=2))                                 **2
    dataRandom  =dble(dataRandomCount  )/dble(size(dataPosition,dim=2))/dble(size(randomPosition,dim=2))
    randomRandom=dble(randomRandomCount)                               /dble(size(randomPosition,dim=2))**2
    ! Compute correlation.
    where (randomRandom > 0.0d0)
       correlationTmp =(dataData-2.0d0*dataRandom+randomRandom)/randomRandom
    elsewhere
       correlationTmp=0.0d0
    end where
    ! Compute final correlation.
    if (projectedActual) then
       ! For projected correlation weight each radial bin by the step in radius, then sum them (to integrate over radial
       ! separation). Summation is multiplied by a factor two since the projected correlation function sums over positive and
       ! negative radial separations, unless a half integral is requested.
       forall(jBin=1:radialBinCount)
          correlationTmp        (:,jBin)=                               &
               & +correlationTmp(:,jBin)                                &
               & *radialMinimum                                         &
               & *(                                                     &
               &  +exp((dble(jBin)-0.5d0)/radialLogarithmicStepInverse) &
               &  -exp((dble(jBin)-1.5d0)/radialLogarithmicStepInverse) &
               & )
       end forall
       correlation=2.0d0*sum(correlationTmp,dim=2)
       if (present(halfIntegral).and.halfIntegral) correlation=correlation/2.0d0
    else
       ! For non-projected correlation simply copy.
       correlation=correlationTmp(:,1)
    end if
    ! Deallocate.
    deallocate(separationSquaredMinimum)
    deallocate(separationSquaredMaximum)
    deallocate(dataDataCount           )
    deallocate(randomRandomCount       )
    deallocate(dataRandomCount         )
    deallocate(pairCount               )
    deallocate(dataData                )
    deallocate(randomRandom            )
    deallocate(dataRandom              )
    deallocate(correlationTmp          )
    return
  end subroutine Statistics_Points_Correlation

end module Statistics_Points_Correlations
