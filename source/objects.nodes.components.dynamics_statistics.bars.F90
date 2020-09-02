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

!% Contains a module which implements tracking of dynamics statistics related to bars.

module Node_Component_Dynamics_Statistics_Bars
  !% Implements tracking of dynamics statistics related to bars.
  use :: Dark_Matter_Halo_Scales            , only : darkMatterHaloScaleClass
  use :: Galactic_Dynamics_Bar_Instabilities, only : galacticDynamicsBarInstabilityClass
  implicit none
  private
  public :: Node_Component_Dynamics_Statistics_Bars_Rate_Compute       , Node_Component_Dynamics_Statistics_Bars_Output    , &
       &    Node_Component_Dynamics_Statistics_Bars_Thread_Initialize  , Node_Component_Dynamics_Statistics_Bars_Initialize, &
       &    Node_Component_Dynamics_Statistics_Bars_Thread_Uninitialize

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

  ! Objects used by this component.
  class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
  class(galacticDynamicsBarInstabilityClass), pointer :: galacticDynamicsBarInstability_
  !$omp threadprivate(darkMatterHaloScale_,galacticDynamicsBarInstability_)

  ! Module initialization state.
  logical          :: dynamicsStatisticsBarsInitialized=.false.
  ! Frequency for recording state.
  double precision :: dynamicsStatisticsBarsFrequency

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Initialize(parameters_)
    !% Initializes the tree node standard disk methods module.
    use :: Galacticus_Nodes, only : defaultDynamicsStatisticsComponent
    use :: Input_Parameters, only : inputParameter                    , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultDynamicsStatisticsComponent%barsIsActive()) then
       !# <inputParameter>
       !#   <name>dynamicsStatisticsBarsFrequency</name>
       !#   <defaultValue>0.1d0</defaultValue>
       !#   <description>The frequency (in fractions of the host halo dynamical time) at which to record the bar dynamical status of satellite galaxies.</description>
       !#   <source>parameters_</source>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Initialize(parameters_)
    !% Initializes the tree node very simple disk profile module.
    use :: Galacticus_Nodes, only : defaultDynamicsStatisticsComponent
    use :: Input_Parameters, only : inputParameter                    , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultDynamicsStatisticsComponent%barsIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters_"/>
       !# <objectBuilder class="galacticDynamicsBarInstability" name="galacticDynamicsBarInstability_" source="parameters_"/>
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Uninitialize()
    !% Uninitializes the tree node very simple disk profile module.
    use :: Galacticus_Nodes, only : defaultDynamicsStatisticsComponent
    implicit none

    if (defaultDynamicsStatisticsComponent%barsIsActive()) then
       !# <objectDestructor name="darkMatterHaloScale_"           />
       !# <objectDestructor name="galacticDynamicsBarInstability_"/>
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Uninitialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute the standard disk node mass rate of change.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultDynamicsStatisticsComponent , interruptTask, nodeComponentBasic, nodeComponentDynamicsStatistics, &
          &                         nodeComponentDynamicsStatisticsBars, treeNode
    implicit none
    type            (treeNode                       ), intent(inout), pointer      :: node
    logical                                          , intent(in   )               :: odeConverged
    logical                                          , intent(inout)               :: interrupt
    procedure       (interruptTask                  ), intent(inout), pointer      :: interruptProcedure
    integer                                          , intent(in   )               :: propertyType
    type            (treeNode                       )               , pointer      :: hostNode
    class           (nodeComponentBasic             )               , pointer      :: basic
    class           (nodeComponentDynamicsStatistics)               , pointer      :: dynamicsStatistics
    double precision                                 , allocatable  , dimension(:) :: timeRecord
    double precision                                                               :: time
    !$GLC attributes unused :: odeConverged, propertyType

    ! Do not compute rates if this component is not active.
    if (.not.defaultDynamicsStatisticsComponent%barsIsActive().or..not.node%isSatellite()) return
    ! Determine the allowed timestep.
    dynamicsStatistics => node%dynamicsStatistics()
    select type (dynamicsStatistics)
    type is (nodeComponentDynamicsStatistics)
       ! Create the component now, and record state immediately.
       interrupt          =  .true.
       interruptProcedure => Node_Component_Dynamics_Statistics_Bars_Record
       return
    class is (nodeComponentDynamicsStatisticsBars)
       ! Set return value if our timestep is smaller than current one.
       timeRecord =  dynamicsStatistics%time  ()
       hostNode   =>  node%parent
       basic      =>  node%basic ()
       time       =  +timeRecord(size(timeRecord))                      &
            &        +dynamicsStatisticsBarsFrequency                   &
            &        *darkMatterHaloScale_%dynamicalTimescale(hostNode)
       ! Check if our timestep is the limiting factor.
       if (basic%time() >= time) then
          interrupt          =  .true.
          interruptProcedure => Node_Component_Dynamics_Statistics_Bars_Record
       end if
    class default
       call Galacticus_Error_Report('unknown class'//{introspection:location})
    end select
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Rate_Compute

  subroutine Node_Component_Dynamics_Statistics_Bars_Record(node)
    !% Record the bar dynamical state of a satellite galaxy.
    use :: Galacticus_Nodes        , only : nodeComponentBasic                              , nodeComponentDisk , nodeComponentDynamicsStatistics, nodeComponentDynamicsStatisticsBars, &
          &                                 nodeComponentSatellite                          , treeNode
    use :: Kepler_Orbits           , only : keplerOrbit
    use :: Numerical_Constants_Math, only : Pi
    use :: Satellite_Orbits        , only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumPericenter
    implicit none
    type            (treeNode                       ), intent(inout), target  :: node
    class           (nodeComponentBasic             )               , pointer :: basic
    class           (nodeComponentDisk              )               , pointer :: disk
    class           (nodeComponentSatellite         )               , pointer :: satellite
    class           (nodeComponentDynamicsStatistics)               , pointer :: dynamicsStatistics
    type            (treeNode                       )               , pointer :: hostNode
    type            (keplerOrbit                    )                         :: orbit
    double precision                                                          :: barInstabilityTimescale, barInstabilityExternalDrivingSpecificTorque, &
         &                                                                       adiabaticRatio         , velocityPericenter                         , &
         &                                                                       radiusPericenter

    ! Get components.
    basic              => node%basic             (                 )
    dynamicsStatistics => node%dynamicsStatistics(autoCreate=.true.)
    ! Record the state.
    select type (dynamicsStatistics)
    class is (nodeComponentDynamicsStatisticsBars)
       disk      => node     %disk       ()
       satellite => node     %satellite  ()
       hostNode  => node     %parent
       orbit     =  satellite%virialOrbit()
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(hostNode,orbit,extremumPericenter,radiusPericenter,velocityPericenter)
       call galacticDynamicsBarInstability_%timescale(node,barInstabilityTimescale,barInstabilityExternalDrivingSpecificTorque)
       if (disk%radius() > 0.0d0) then
          adiabaticRatio=(radiusPericenter/velocityPericenter)/(2.0d0*Pi*disk%radius()/disk%velocity())
       else
          adiabaticRatio=-1.0d0
       end if
       call dynamicsStatistics%record(basic%time(),barInstabilityTimescale,adiabaticRatio)
    end select
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Record

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Store the dynamical histories of galaxies to \glc\ output file.
    use            :: Galacticus_HDF5                 , only : galacticusOutputFile
    use            :: Galacticus_Nodes                , only : nodeComponentDynamicsStatistics, nodeComponentDynamicsStatisticsBars, treeNode
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : assignment(=)                  , char                               , varying_string
    use            :: Kind_Numbers                    , only : kind_int8
    use            :: Numerical_Constants_Astronomical, only : gigaYear
    use            :: String_Handling                 , only : operator(//)
    implicit none
    type            (treeNode                       ), intent(inout), pointer      :: node
    integer         (kind=c_size_t                  ), intent(in   )               :: iOutput
    integer         (kind=kind_int8                 ), intent(in   )               :: treeIndex
    logical                                          , intent(in   )               :: nodePassesFilter
    class           (nodeComponentDynamicsStatistics)               , pointer      :: dynamicsStatistics
    double precision                                 , allocatable  , dimension(:) :: dataValues
    type            (varying_string                 )                              :: outputGroupName          , treeGroupName       , &
         &                                                                            timeDatasetName          , timescaleDatasetName, &
         &                                                                            adiabaticRatioDatasetName
    type            (hdf5Object                     )                              :: outputs                  , output              , &
         &                                                                            dynamics                 , tree                , &
         &                                                                            dataset
    !$GLC attributes unused :: nodePassesFilter

    ! Output the history data if and only if any has been collated.
    if (dynamicsStatisticsBarsInitialized) then
       ! Get the dynamics statistics component.
       dynamicsStatistics => node%dynamicsStatistics()
       select type (dynamicsStatistics)
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
          timeDatasetName          ="time"
          timescaleDatasetName     ="timeScale"
          adiabaticRatioDatasetName="adiabaticRatio"
          timeDatasetName          =timeDatasetName          //node%index()
          timescaleDatasetName     =timescaleDatasetName     //node%index()
          adiabaticRatioDatasetName=adiabaticRatioDatasetName//node%index()
          dataValues               =dynamicsStatistics%time                   ()
          call tree   %writeDataset  (dataValues,char(timeDatasetName          ),"Time [Gyr]"      ,datasetReturned=dataset)
          call dataset%writeAttribute(gigaYear                                  ,"unitsInSI"                               )
          call dataset%close         (                                                                                     )
          dataValues               =dynamicsStatistics%barInstabilityTimescale()
          call tree   %writeDataset  (dataValues,char(timescaleDatasetName     ),"Time scale [Gyr]",datasetReturned=dataset)
          call dataset%writeAttribute(gigaYear                                  ,"unitsInSI"                               )
          call dataset%close         (                                                                                     )
          dataValues               =dynamicsStatistics%adiabaticRatio         ()
          call tree   %writeDataset  (dataValues,char(adiabaticRatioDatasetName),"[]"              ,datasetReturned=dataset)
          call dataset%writeAttribute(gigaYear                                  ,"unitsInSI"                               )
          call dataset%close         (                                                                                     )
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
