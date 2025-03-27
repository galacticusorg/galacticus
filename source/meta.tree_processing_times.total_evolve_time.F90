!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Implements a merger tree processing time estimator that estimates times based on the total time to evolve in the tree.
!!}

  !![
  <metaTreeProcessingTime name="metaTreeProcessingTimeTotalEvolveTime">
   <description>
    A merger tree processing time class 
   </description>
  </metaTreeProcessingTime>
  !!]
  type, extends(metaTreeProcessingTimeClass) :: metaTreeProcessingTimeTotalEvolveTime
     !!{
     A merger tree processing time estimator that estimates times based on the total time to evolve in the tree.
     !!}
     private
     double precision :: exponentTime, updateInterval
   contains
     procedure :: timeRemaining => totalEvolveTimeTimeRemaining
  end type metaTreeProcessingTimeTotalEvolveTime

  interface metaTreeProcessingTimeTotalEvolveTime
     !!{
     Constructors for the {\normalfont \ttfamily totalEvolveTime} merger tree processing time estimator.
     !!}
     module procedure totalEvolveTimeConstructorParameters
     module procedure totalEvolveTimeConstructorInternal
  end interface metaTreeProcessingTimeTotalEvolveTime

contains

  function totalEvolveTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily totalEvolveTime} merger tree processing time estimator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (metaTreeProcessingTimeTotalEvolveTime)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: exponentTime, updateInterval
    
    !![
    <inputParameter>
      <name>exponentTime</name>
      <description>The exponent of cosmic time used in estimating the work associated with evolving a node.</description>
      <defaultValue>0.5d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>updateInterval</name>
      <description>The minimum interval (in seconds) between updates of the estimated time remaining.</description>
      <defaultValue>10.0d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=metaTreeProcessingTimeTotalEvolveTime(exponentTime,updateInterval)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function totalEvolveTimeConstructorParameters

  function totalEvolveTimeConstructorInternal(exponentTime,updateInterval) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily totalEvolveTime} merger tree processing time estimator class.
    !!}
    implicit none
    type            (metaTreeProcessingTimeTotalEvolveTime)                :: self
    double precision                                       , intent(in   ) :: exponentTime, updateInterval
    !![
    <constructorAssign variables="exponentTime, updateInterval"/>
    !!]
    
    return
  end function totalEvolveTimeConstructorInternal

  double precision function totalEvolveTimeTimeRemaining(self,tree,timeFinal) result(timeRemaining)
    !!{
    Estimate the remaining time to process the tree.
    !!}
    use   , intrinsic :: ISO_C_Binding      , only : c_size_t
    !$ use            :: OMP_Lib            , only : OMP_Get_Num_Threads
    use               :: Galacticus_Nodes   , only : treeNode                , nodeComponentBasic
    use               :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class           (metaTreeProcessingTimeTotalEvolveTime), intent(inout) :: self
    type            (mergerTree                           ), intent(inout) :: tree
    double precision                                       , intent(in   ) :: timeFinal
    type            (treeNode                             ), pointer       :: node            , nodeParent
    class           (nodeComponentBasic                   ), pointer       :: basic           , basicParent
    type            (mergerTreeWalkerAllNodes             )                :: treeWalker
    double precision                                                       :: time            , timeCPU           , &
         &                                                                    workTotalCentral, workTotalSatellite
    logical                                                                :: isFirst         , compute
    real                                                                   :: timeCPU_
    integer         (c_size_t                             )                :: countNodes

    ! Assume that we have no estimate by default.
    timeRemaining=-1.0d0
    ! Determine if this is the first estimate to be made for this tree.
    isFirst=.not.tree%properties%exists('metaTreeProcessingTimeLast')
    ! Find the current CPU time and decide if we want to make an estimate of the remaining time.
    call CPU_Time(timeCPU_)
    timeCPU=dble(timeCPU_)
    !$ timeCPU=timeCPU/dble(OMP_Get_Num_Threads())
    if (isFirst) then
       compute=.true.
    else
       compute=timeCPU-tree%properties%value('metaTreeProcessingTimeLast') > self%updateInterval
    end if
    if (.not.compute) return
    ! Estimate the remaining work. We treat centrals and satellites separately (assuming that the work involved in evolving each
    ! could be quite different).
    countNodes        =0_c_size_t
    workTotalCentral  =0.0d0
    workTotalSatellite=0.0d0
    treeWalker        =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! Consider only nodes at the tip of their branch.
       if (associated(node%firstChild)) cycle
       basic      => node %basic()
       time       =  basic%time ()
       countNodes =  countNodes+1_c_size_t
       if (node%isSatellite()) then
          ! Satellite node - accumulate all work to the satellite work.
          workTotalSatellite=workTotalSatellite+max(timeFinal**self%exponentTime-time**self%exponentTime,0.0d0)
       else
          ! Central node - accumulate work to both central and satellite work - split at the time when this node will become a
          ! satellite.
          nodeParent => node
          do while (nodeParent%isPrimaryProgenitor())
             nodeParent => nodeParent%parent
          end do
          basicParent        => nodeParent%basic()
          workTotalCentral   =  workTotalCentral  +max(basicParent%time     ()**self%exponentTime-            time  **self%exponentTime,0.0d0)
          workTotalSatellite =  workTotalSatellite+max(            timeFinal  **self%exponentTime-basicParent%time()**self%exponentTime,0.0d0)
       end if
    end do
    if (.not.isFirst) then
       ! Accumulate time remaining.
       timeRemaining=0.0d0
       if (workTotalCentral   < tree%properties%value('metaTreeProcessingWorkTotalCentralLast'  ))                   &
           & timeRemaining=+timeRemaining                                                                           &
           &               +workTotalCentral                                                                        &
           &               *(-tree%properties%value('metaTreeProcessingTimeLast'              )+timeCPU           ) &
           &               /(+tree%properties%value('metaTreeProcessingWorkTotalCentralLast'  )-workTotalCentral  )
      if (workTotalSatellite < tree%properties%value('metaTreeProcessingWorkTotalSatelliteLast'))                   &
           & timeRemaining=+timeRemaining                                                                           &
           &               +workTotalSatellite                                                                      &
           &               *(-tree%properties%value('metaTreeProcessingTimeLast'              )+timeCPU           ) &
           &               /(+tree%properties%value('metaTreeProcessingWorkTotalSatelliteLast')-workTotalSatellite)      
   end if
   ! Store the CPU time and work for reference at the next estimate.
   call tree%properties%set('metaTreeProcessingTimeLast'              ,timeCPU           )
   call tree%properties%set('metaTreeProcessingWorkTotalCentralLast'  ,workTotalCentral  )
   call tree%properties%set('metaTreeProcessingWorkTotalSatelliteLast',workTotalSatellite)
   return
  end function totalEvolveTimeTimeRemaining

