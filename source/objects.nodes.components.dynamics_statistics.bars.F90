!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements tracking of dynamics statistics related to bars.

module Node_Component_Dynamics_Statistics_Bars
  !% Implements tracking of dynamics statistics related to bars.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Dynamics_Statistics_Bars_Timestep, Node_Component_Dynamics_Statistics_Bars_Output

  !# <component>
  !#  <class>dynamicsStatistics</class>
  !#  <name>bars</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>time</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>barInstabilityTimescale</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>adiabaticRatio</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="record" function="Node_Component_Dynamics_Statistics_Bars_Record" bindsTo="component" description="Record the current bar-dynamical state of the node." returnType="\void" arguments="\doublezero\ time\argin, \doublezero\ barInstabilityTimescale\argin, \doublezero adiabaticRatio\argin" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.dynamics_statistics.bars.bound_functions.inc</functions>
  !# </component>

  ! Module initialization state.
  logical          :: dynamicsStatisticsBarsInitialized=.false.
  ! Frequency for recording state.
  double precision :: dynamicsStatisticsBarsFrequency

contains

  !# <timeStepsTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Timestep</unitName>
  !# </timeStepsTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Timestep(thisNode,timeStep,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Determines the timestep to go to the next tabulation point for galactic bar dynamics storage.
    use Input_Parameters
    use Evolve_To_Time_Reports
    use ISO_Varying_String
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    implicit none
    type            (treeNode                                      ), intent(inout)             , pointer :: thisNode
    procedure       (Node_Component_Dynamics_Statistics_Bars_Record), intent(inout)             , pointer :: End_Of_Timestep_Task
    double precision                                                , intent(inout)                       :: timeStep
    logical                                                         , intent(in   )                       :: report
    type            (treeNode                                      ), intent(inout), optional   , pointer :: lockNode
    type            (varying_string                                ), intent(inout), optional             :: lockType
    type            (treeNode                                      )                            , pointer :: hostNode
    class           (nodeComponentBasic                            )                            , pointer :: thisBasic
    class           (nodeComponentDynamicsStatistics               )                            , pointer :: thisDynamicsStatistics
    class           (darkMatterHaloScaleClass                      )                            , pointer :: darkMatterHaloScale_
    double precision                                                , allocatable  , dimension(:)         :: timeRecord
    double precision                                                                                      :: ourTimeStep

    ! Return immediately if this class is not active or if this galaxy is not a satellite.
    if (.not.(defaultDynamicsStatisticsComponent%barsIsActive().and.thisNode%isSatellite())) return
    ! Initialize if necessary
    if (.not.dynamicsStatisticsBarsInitialized) then
       !$omp critical (dynamicsStatisticsBarsInitialize)
       if (.not.dynamicsStatisticsBarsInitialized) then
          ! Get module parameters.
          !@ <inputParameter>
          !@   <name>dynamicsStatisticsBarsFrequency</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The frequency (in fractions of the host halo dynamical time) at which to record the bar dynamical status of satellite galaxies.
          !@   </description>
          !@   <type>double</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('dynamicsStatisticsBarsFrequency',dynamicsStatisticsBarsFrequency,defaultValue=0.1d0)
          ! Record that initialization is now complete.
          dynamicsStatisticsBarsInitialized=.true.
       end if
       !$omp end critical (dynamicsStatisticsBarsInitialize)
    end if
    ! Determine the allowed timestep.
    thisDynamicsStatistics => thisNode%dynamicsStatistics()
    select type (thisDynamicsStatistics)
    type is (nodeComponentDynamicsStatistics)
       ! Create the component now, and record state immediately.
       ourTimeStep=0.0d0
    class is (nodeComponentDynamicsStatisticsBars)
       ! Set return value if our timestep is smaller than current one.
       timeRecord           =  thisDynamicsStatistics%time  ()
       hostNode             => thisNode              %parent
       thisBasic            => thisNode              %basic ()
       darkMatterHaloScale_ => darkMatterHaloScale          ()
       ourTimestep=                                                   &
            & max(                                                    &
            &      timeRecord(size(timeRecord))                       &
            &     +dynamicsStatisticsBarsFrequency                    &
            &     *darkMatterHaloScale_%dynamicalTimescale(hostNode)  &
            &     -thisBasic%time()                                 , &
            &     0.0d0                                               &
            &    )
    class default
       ourTimestep=0.0d0
       call Galacticus_Error_Report('Node_Component_Dynamics_Statistics_Bars_Timestep','unknown class')
    end select
    ! Check if our timestep is the limiting factor.
    if (ourTimeStep <= timeStep) then
       if (present(lockNode)) lockNode => thisNode
       if (present(lockType)) lockType =  "galactic dynamics statistics (bars)"
       timeStep=ourTimeStep
       End_Of_Timestep_Task => Node_Component_Dynamics_Statistics_Bars_Record
    end if
    if (report) call Evolve_To_Time_Report("galactic dynamics statistics (bars): ",timeStep)
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Timestep

  subroutine Node_Component_Dynamics_Statistics_Bars_Record(thisTree,thisNode,deadlockStatus)
    !% Record the bar dynamical state of a satellite galaxy.
    use Numerical_Interpolation
    use Numerical_Constants_Math
    use Galactic_Dynamics_Bar_Instabilities
    use Satellite_Orbits
    use Kepler_Orbits
    implicit none
    type            (mergerTree                     ), intent(in   )          :: thisTree
    type            (treeNode                       ), intent(inout), pointer :: thisNode
    integer                                          , intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic             )               , pointer :: thisBasic
    class           (nodeComponentDisk              )               , pointer :: thisDisk
    class           (nodeComponentSatellite         )               , pointer :: thisSatellite
    class           (nodeComponentDynamicsStatistics)               , pointer :: thisDynamicsStatistics
    type            (treeNode                       )               , pointer :: hostNode
    type            (keplerOrbit                    )                         :: thisOrbit
    double precision                                                          :: barInstabilityTimescale, barInstabilityExternalDrivingSpecificTorque, &
         &                                                                       adiabaticRatio         , velocityPericenter                         , &
         &                                                                       radiusPericenter
    !GCC$ attributes unused :: thisTree, deadlockStatus
    
    ! Get components.
    thisBasic              => thisNode%basic             (                 )
    thisDynamicsStatistics => thisNode%dynamicsStatistics(autoCreate=.true.)
    ! Record the state.
    select type (thisDynamicsStatistics)
    class is (nodeComponentDynamicsStatisticsBars)
       thisDisk      => thisNode     %disk       ()
       thisSatellite => thisNode     %satellite  ()
       hostNode      => thisNode     %parent
       thisOrbit     =  thisSatellite%virialOrbit()
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(hostNode,thisOrbit,extremumPericenter,radiusPericenter,velocityPericenter)
       call Bar_Instability_Timescale(thisNode,barInstabilityTimescale,barInstabilityExternalDrivingSpecificTorque)
       if (thisDisk%radius() > 0.0d0) then
          adiabaticRatio=(radiusPericenter/velocityPericenter)/(2.0d0*Pi*thisDisk%radius()/thisDisk%velocity())
       else
          adiabaticRatio=-1.0d0
       end if
       call thisDynamicsStatistics%record(thisBasic%time(),barInstabilityTimescale,adiabaticRatio)
    end select
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Record

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store the dynamical histories of galaxies to \glc\ output file.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_HDF5
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    use String_Handling
    use ISO_Varying_String
    use IO_HDF5
    implicit none
    type            (treeNode                       ), intent(inout), pointer      :: thisNode
    integer         (kind=c_size_t                  ), intent(in   )               :: iOutput
    integer         (kind=kind_int8                 ), intent(in   )               :: treeIndex
    logical                                          , intent(in   )               :: nodePassesFilter
    class           (nodeComponentDynamicsStatistics)               , pointer      :: thisDynamicsStatistics
    double precision                                 , allocatable  , dimension(:) :: thisData
    type            (varying_string                 )                              :: outputGroupName          , treeGroupName       , &
         &                                                                            timeDatasetName          , timescaleDatasetName, &
         &                                                                            adiabaticRatioDatasetName
    type            (hdf5Object                     )                              :: outputs                  , output              , &
         &                                                                            dynamics                 , tree                , &
         &                                                                            thisDataset
    !GCC$ attributes unused :: nodePassesFilter
    
    ! Output the history data if and only if any has been collated.
    if (dynamicsStatisticsBarsInitialized) then
       ! Get the dynamics statistics component.
       thisDynamicsStatistics => thisNode%dynamicsStatistics()
       select type (thisDynamicsStatistics)
       class is (nodeComponentDynamicsStatisticsBars)
          ! Open groups for writing.
          outputGroupName="Output"
          outputGroupName=outputGroupName//iOutput
          treeGroupName  ="tree"
          treeGroupName=treeGroupName//treeIndex
          outputs =galacticusOutputFile%openGroup('Outputs'                                                                                  )
          output  =outputs             %openGroup(char(outputGroupName)                                                                      )
          dynamics=output              %openGroup("dynamicsStatistics" ,"Group containing statistics of satellite galaxy bar dynamical state")
          tree    =dynamics            %openGroup(char(treeGroupName)  ,"Group containing statistics of satellite galaxy bar dynamical state")
          !@ <outputProperty>
          !@   <name>time</name>
          !@   <datatype>double</datatype>
          !@   <cardinality>0..</cardinality>
          !@   <description>The times at which bar dynamical state data are stored.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          !@ <outputProperty>
          !@   <name>timeScale</name>
          !@   <datatype>double</datatype>
          !@   <cardinality>0..</cardinality>
          !@   <description>The bar instability timescale.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          !@ <outputProperty>
          !@   <name>adiabaticRatio</name>
          !@   <datatype>double</datatype>
          !@   <cardinality>0..</cardinality>
          !@   <description>The adiabatic ratio.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          timeDatasetName          ="time"
          timescaleDatasetName     ="timeScale"
          adiabaticRatioDatasetName="adiabaticRatio"
          timeDatasetName          =timeDatasetName          //thisNode%index()
          timescaleDatasetName     =timescaleDatasetName     //thisNode%index()
          adiabaticRatioDatasetName=adiabaticRatioDatasetName//thisNode%index()
          thisData            =thisDynamicsStatistics%time                   ()
          call tree       %writeDataset  (thisData,char(timeDatasetName          ),"Time [Gyr]"      ,datasetReturned=thisDataset)
          call thisDataset%writeAttribute(gigaYear                                ,"unitsInSI"                                   )
          call thisDataset%close         (                                                                                       )
          thisData            =thisDynamicsStatistics%barInstabilityTimescale()
          call tree       %writeDataset  (thisData,char(timescaleDatasetName     ),"Time scale [Gyr]",datasetReturned=thisDataset)
          call thisDataset%writeAttribute(gigaYear                                ,"unitsInSI"                                   )
          call thisDataset%close         (                                                                                       )
          thisData            =thisDynamicsStatistics%adiabaticRatio         ()
          call tree       %writeDataset  (thisData,char(adiabaticRatioDatasetName),"[]"              ,datasetReturned=thisDataset)
          call thisDataset%writeAttribute(gigaYear                                ,"unitsInSI"                                   )
          call thisDataset%close         (                                                                                       )
          ! Close groups. 
          call tree    %close()
          call dynamics%close()
          call output  %close()
          call outputs %close()
       end select
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Output

end module Node_Component_Dynamics_Statistics_Bars
