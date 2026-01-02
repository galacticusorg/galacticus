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
  Implements a node operator class that triggers destruction of satellites based on their bound mass.
  !!}

  !![
  <nodeOperator name="nodeOperatorSatelliteDestructionMassThreshold">
   <description>A node operator class that triggers destruction of satellites based on their bound mass.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteDestructionMassThreshold
     !!{
     A node operator class that triggers destruction of satellites based on their mass.
     !!}
     private
     double precision :: massDestructionAbsolute        , massDestructionMassInfallFraction, &
          &              massDestructionMassTreeFraction
     logical          :: mergeOnDestruction
   contains
     !![
     <methods>
       <method description="Compute the mass at which the satellite will be destroyed." method="massDestroy" />
     </methods>
     !!]
     procedure :: differentialEvolution => satelliteDestructionMassThresholdDifferentialEvolution
     procedure :: massDestroy           => satelliteDestructionMassThresholdMassDestroy
  end type nodeOperatorSatelliteDestructionMassThreshold
  
  interface nodeOperatorSatelliteDestructionMassThreshold
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteDestructionMassThreshold} node operator class.
     !!}
     module procedure satelliteDestructionMassThresholdConstructorParameters
     module procedure satelliteDestructionMassThresholdConstructorInternal
  end interface nodeOperatorSatelliteDestructionMassThreshold

  ! Submodule-scope pointer to self, used in callback functions.
  class(nodeOperatorSatelliteDestructionMassThreshold), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function satelliteDestructionMassThresholdConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteDestructionMassThreshold} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteDestructionMassThreshold)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: massDestructionAbsolute        , massDestructionMassInfallFraction, &
         &                                                                            massDestructionMassTreeFraction
    logical                                                                        :: mergeOnDestruction

    !![
    <inputParameter>
      <name>massDestructionAbsolute</name>
      <defaultValue>0.00d0</defaultValue>
      <description>The absolute mass below which satellites are destroyed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massDestructionMassInfallFraction</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The fraction of the infall mass below which satellites are destroyed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massDestructionMassTreeFraction</name>
      <defaultValue>0.00d0</defaultValue>
      <description>The fraction of the tree mass below which satellites are destroyed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mergeOnDestruction</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, destruction actually triggers the satellite to be merged with its central galaxy instead of being destroyed.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorSatelliteDestructionMassThreshold(massDestructionAbsolute,massDestructionMassInfallFraction,massDestructionMassTreeFraction,mergeOnDestruction)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteDestructionMassThresholdConstructorParameters

  function satelliteDestructionMassThresholdConstructorInternal(massDestructionAbsolute,massDestructionMassInfallFraction,massDestructionMassTreeFraction,mergeOnDestruction) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteDestructionMassThreshold} node operator class.
    !!}
    implicit none
    type            (nodeOperatorSatelliteDestructionMassThreshold)                :: self
    double precision                                               , intent(in   ) :: massDestructionAbsolute        , massDestructionMassInfallFraction, &
         &                                                                            massDestructionMassTreeFraction
    logical                                                        , intent(in   ) :: mergeOnDestruction
    !![
    <constructorAssign variables="massDestructionAbsolute, massDestructionMassInfallFraction, massDestructionMassTreeFraction, mergeOnDestruction"/>
    !!]
    
    return
  end function satelliteDestructionMassThresholdConstructorInternal

  subroutine satelliteDestructionMassThresholdDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Trigger destruction of a satellite halo based on its bound mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class           (nodeOperatorSatelliteDestructionMassThreshold), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    logical                                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                                ), intent(inout), pointer :: functionInterrupt
    integer                                                        , intent(in   )          :: propertyType
    class           (nodeComponentSatellite                       )               , pointer :: satellite
    double precision                                                                        :: massSatellite
    !$GLC attributes unused :: propertyType
    
    if (.not.node%isSatellite()) return
    satellite     => node     %satellite()
    massSatellite =  satellite%boundMass()
    if     (                                        &
         &   massSatellite > 0.0d0                  &
         &  .and.                                   &
         &   massSatellite < self%massDestroy(node) &
         & ) then
       ! Destruction criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => destructionTrigger
       self_             => self
       return
    end if
    return
  end subroutine satelliteDestructionMassThresholdDifferentialEvolution
  
  subroutine destructionTrigger(node,timeEnd)
    !!{
    Trigger destruction of the satellite by setting the time until destruction to zero.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic, treeNode
    use :: Display         , only : displayGreen          , displayReset
    use :: Error           , only : Error_Report
    use :: String_Handling , only : stringXMLFormat
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentSatellite)               , pointer  :: satellite
    class           (nodeComponentBasic    )               , pointer  :: basic
    !$GLC attributes unused :: timeEnd
    
    satellite => node%satellite()
    if (satellite%boundMass() < self_%massDestroy(node)) then
       if (self_%mergeOnDestruction) then
          ! Destruction is really just merging. Set the merging time to now to trigger merging to occur.
          basic => node%basic()
          call satellite%timeOfMergingSet(basic%time())
       else
          ! Destruction really mean destruction. Set the destruction time to zero to trigger destruction to occur.
          if (satellite%destructionTime() >= 0.0d0)                                                                                                                   &
               & call Error_Report(                                                                                                                                   &
               &                   'satellite was previously triggered for destruction - but still exists'                                               //char(10)// &
               &                   displayGreen()//'  HELP:'//displayReset()//' destruction requires the following timestepper to be included:'//char(10)//char(10)// &
               &                   stringXMLFormat('mergerTreeEvolveTimestep value="satelliteDestruction"/>')                                                      // &
               &                   {introspection:location}                                                                                                           &
               &                  )
          call satellite%destructionTimeSet(0.0d0)
       end if
    end if
    return
  end subroutine destructionTrigger

  double precision function satelliteDestructionMassThresholdMassDestroy(self,node)
    !!{
    Compute the destruction mass for a node.
    !!}
    use :: Galacticus_Nodes, only : treeNode, nodeComponentBasic
    implicit none
    class           (nodeOperatorSatelliteDestructionMassThreshold), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    class           (nodeComponentBasic                           ), pointer       :: basic
    double precision                                                               :: massDestructionAbsolute, massDestructionInfall, &
         &                                                                            massDestructionTree
    
    massDestructionAbsolute=self%massDestructionAbsolute
    massDestructionInfall  =0.0d0
    massDestructionTree    =0.0d0
    if (self%massDestructionMassInfallFraction > 0.0d0) then
       basic                 =>  node                   %basic                            ()
       massDestructionInfall =  +basic                  %mass                             () &
            &                   *self                   %massDestructionMassInfallFraction
    end if
    if (self%massDestructionMassTreeFraction > 0.0d0)  then
       basic                 =>  node %hostTree%nodeBase%basic                            ()
       massDestructionTree   =  +basic                  %mass                             () &
            &                   *self                   %massDestructionMassTreeFraction
    end if
    satelliteDestructionMassThresholdMassDestroy=max(                         &
         &                                           massDestructionAbsolute, &
         &                                           massDestructionInfall  , &
         &                                           massDestructionTree      &
         &                                          )
    return
  end function satelliteDestructionMassThresholdMassDestroy
