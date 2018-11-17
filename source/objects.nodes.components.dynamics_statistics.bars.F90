!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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
  use Dark_Matter_Halo_Scales
  use Galactic_Dynamics_Bar_Instabilities
  implicit none
  private
  public :: Node_Component_Dynamics_Statistics_Bars_Rate_Compute     , Node_Component_Dynamics_Statistics_Bars_Output    , &
       &    Node_Component_Dynamics_Statistics_Bars_Thread_Initialize, Node_Component_Dynamics_Statistics_Bars_Initialize

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
  subroutine Node_Component_Dynamics_Statistics_Bars_Initialize(parameters)
    !% Initializes the tree node standard disk methods module.
    use Input_Parameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    
    if (defaultDynamicsStatisticsComponent%barsIsActive()) then
       !# <inputParameter>
       !#   <name>dynamicsStatisticsBarsFrequency</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>0.1d0</defaultValue>
       !#   <description>The frequency (in fractions of the host halo dynamical time) at which to record the bar dynamical status of satellite galaxies.</description>
       !#   <group>timeStepping</group>
       !#   <source>parameters</source>
       !#   <type>double</type>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Initialize

  !# <mergerTreeEvolveThreadInitialize>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Thread_Initialize</unitName>
  !# </mergerTreeEvolveThreadInitialize>
  subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Initialize(parameters)
    !% Initializes the tree node very simple disk profile module.
    use Input_Parameters
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (defaultDynamicsStatisticsComponent%barsIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
       !# <objectBuilder class="galacticDynamicsBarInstability" name="galacticDynamicsBarInstability_" source="parameters"/>
    end if
    return
  end subroutine Node_Component_Dynamics_Statistics_Bars_Thread_Initialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dynamics_Statistics_Bars_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dynamics_Statistics_Bars_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute the standard disk node mass rate of change.
    use Galacticus_Error
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
    !GCC$ attributes unused :: odeConverged, propertyType

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
    use Numerical_Interpolation
    use Numerical_Constants_Math
    use Satellite_Orbits
    use Kepler_Orbits
    implicit none
    type            (treeNode                       ), intent(inout), pointer :: node
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
    use, intrinsic :: ISO_C_Binding
    use Galacticus_HDF5
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    use String_Handling
    use ISO_Varying_String
    use IO_HDF5
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
    !GCC$ attributes unused :: nodePassesFilter
    
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
