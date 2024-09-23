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
  Implements a node operator class that triggers destruction of satellites based on their density profiles.
  !!}
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorSatelliteDestructionDensityProfileThreshold">
   <description>A node operator class that triggers destruction of satellites based on their density profiles.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteDestructionDensityProfileThreshold
     !!{
     A node operator class that triggers destruction of satellites based on their density profiles.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_                => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_                 => null()
     double precision                                     :: fractionDensityProfileVirialFraction
   contains
     !![
     <methods>
       <method description="Return true if the node should be destroyed." method="shouldDestroy"/>
     </methods>
     !!]
     final     ::                          satelliteDestructionDensityProfileThresholdDestructor
     procedure :: differentialEvolution => satelliteDestructionDensityProfileThresholdDiffEvoltn
     procedure :: shouldDestroy         => satelliteDestructionDensityProfileThresholdShouldDestroy
  end type nodeOperatorSatelliteDestructionDensityProfileThreshold
  
  interface nodeOperatorSatelliteDestructionDensityProfileThreshold
     !!{
     Constructors for the {\normalfont \ttfamily satelliteDestructionDensityProfileThreshold} node operator class.
     !!}
     module procedure satelliteDestructionDensityProfileThresholdCnstrctrPrmtrs
     module procedure satelliteDestructionDensityProfileThresholdCnstrctrIntrnl
  end interface nodeOperatorSatelliteDestructionDensityProfileThreshold

  ! Submodule-scope pointer to self, used in callback functions.
  class(nodeOperatorSatelliteDestructionDensityProfileThreshold), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function satelliteDestructionDensityProfileThresholdCnstrctrPrmtrs(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteDestructionDensityProfileThreshold} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteDestructionDensityProfileThreshold)                :: self
    type            (inputParameters                                        ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                               ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                              ), pointer       :: darkMatterProfileDMO_
    double precision                                                                         :: fractionDensityProfileVirialFraction

    !![
    <inputParameter>
      <name>fractionDensityProfileVirialFraction</name>
      <description>The absolute mass below which satellites are destroyed.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteDestructionDensityProfileThreshold(fractionDensityProfileVirialFraction,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function satelliteDestructionDensityProfileThresholdCnstrctrPrmtrs

  function satelliteDestructionDensityProfileThresholdCnstrctrIntrnl(fractionDensityProfileVirialFraction,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily satelliteDestructionDensityProfileThreshold} node operator class.
    !!}
    implicit none
    type            (nodeOperatorSatelliteDestructionDensityProfileThreshold)                        :: self
    class           (darkMatterHaloScaleClass                               ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                              ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                                         , intent(in   )         :: fractionDensityProfileVirialFraction
    !![
    <constructorAssign variables="fractionDensityProfileVirialFraction, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]
    
    return
  end function satelliteDestructionDensityProfileThresholdCnstrctrIntrnl

  subroutine satelliteDestructionDensityProfileThresholdDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily satelliteDestructionDensityProfileThreshold} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteDestructionDensityProfileThreshold), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine satelliteDestructionDensityProfileThresholdDestructor

  logical function satelliteDestructionDensityProfileThresholdShouldDestroy(self,node) result(shouldDestroy)
    !!{
    Determine if the satellite should be destroyed.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorSatelliteDestructionDensityProfileThreshold), intent(inout), target  :: self
    type            (treeNode                                               ), intent(inout), target  :: node
    class           (nodeComponentDarkMatterProfile                         )               , pointer :: darkMatterProfile
    double precision                                                                                  :: radiusScale

    darkMatterProfile =>   node                                   %darkMatterProfile                   (                )
    radiusScale       =   +darkMatterProfile                      %scale                               (                )
    shouldDestroy     =   +self             %darkMatterProfileDMO_%density                             (node,radiusScale) &
         &               <                                                                                                &
         &                +self                                   %fractionDensityProfileVirialFraction                   &
         &                *self             %darkMatterHaloScale_ %densityMean                         (node            )
    return
  end function satelliteDestructionDensityProfileThresholdShouldDestroy
  
  subroutine satelliteDestructionDensityProfileThresholdDiffEvoltn(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Trigger destruction of a satellite halo based on its density profile.
    !!}
    implicit none
    class    (nodeOperatorSatelliteDestructionDensityProfileThreshold), intent(inout), target  :: self
    type     (treeNode                                               ), intent(inout), target  :: node
    logical                                                           , intent(inout)          :: interrupt
    procedure(interruptTask                                          ), intent(inout), pointer :: functionInterrupt
    integer                                                           , intent(in   )          :: propertyType
    !$GLC attributes unused :: propertyType
    
    if (.not.node%isSatellite()) return
    if (self%shouldDestroy(node)) then
       ! Destruction criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => destructionTrigger
       self_             => self
       return
    end if
    return
  end subroutine satelliteDestructionDensityProfileThresholdDiffEvoltn
  
  subroutine destructionTrigger(node,timeEnd)
    !!{
    Trigger destruction of the satellite by setting the time until destruction to zero.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    use :: Display         , only : displayBlue           , displayYellow, displayGreen, displayReset
    use :: Error           , only : Error_Report
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentSatellite)               , pointer  :: satellite
    !$GLC attributes unused :: timeEnd
    
    satellite => node%satellite()
    if (self_%shouldDestroy(node)) then
       if (satellite%destructionTime() >= 0.0d0)                                                                                                                                                                                &
            call Error_Report(                                                                                                                                                                                                  &
            &                 'satellite was previously triggered for destruction - but still exists'                                                                                                              //char(10)// &
            &                 displayGreen()//'  HELP:'//displayReset()//' destruction requires the following timestepper to be included:'                                                                         //char(10)// &
            &                 '    <'//displayBlue()//'mergerTreeEvolveTimestep'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"satelliteDestruction"'//displayReset()//'>'          // &
            &                 {introspection:location}                                                                                                                                                                          &
            &                )
       call satellite%destructionTimeSet(0.0d0)
    end if
    return
  end subroutine destructionTrigger
