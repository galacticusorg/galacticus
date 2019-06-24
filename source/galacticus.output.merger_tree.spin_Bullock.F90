!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of node virial data to the \glc\ output file.

module Galacticus_Output_Trees_Spin_Bullock
  !% Handles outputting of node virial data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Spin_Bullock, Galacticus_Output_Tree_Spin_Bullock_Property_Count, Galacticus_Output_Tree_Spin_Bullock_Names

  ! Flag indicating whether or not Bullock-style spin parameters are to be output.
  logical :: outputSpinBullock

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputSpinBullockInitialized=.false.

  ! Record of whether vector spin information is available.
  logical :: vectorSpinAvailable
  
contains

  subroutine Galacticus_Output_Tree_Spin_Bullock_Initialize()
    !% Initializes the module by determining whether or not virial data should be output.
    use Input_Parameters
    use Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (.not.outputSpinBullockInitialized) then
       !$omp critical(Galacticus_Output_Tree_Spin_Bullock_Initialize)
       if (.not.outputSpinBullockInitialized) then
          !# <inputParameter>
          !#   <name>outputSpinBullock</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not Bullock-style spin parameters should be included in the output.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Determine if vector spin information is available.
          vectorSpinAvailable=defaultSpinComponent%spinVectorIsGettable()
          ! Flag that module is now initialized.
          outputSpinBullockInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Spin_Bullock_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Bullock_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Spin_Bullock_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Bullock</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Spin_Bullock_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of virial properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: thisNode
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: thisNode, time, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Spin_Bullock_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputSpinBullock) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spinBullock'
       doublePropertyComments(doubleProperty)='Spin parameter of the halo under the Bullock et al. (2001) definition [].'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
       if (vectorSpinAvailable) then
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spinBullockX'
          doublePropertyComments(doubleProperty)='x-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].'
          doublePropertyUnitsSI (doubleProperty)=0.0d0
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spinBullockY'
          doublePropertyComments(doubleProperty)='y-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].'
          doublePropertyUnitsSI (doubleProperty)=0.0d0
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spinBullockZ'
          doublePropertyComments(doubleProperty)='z-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].'
          doublePropertyUnitsSI (doubleProperty)=0.0d0
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Bullock_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Spin_Bullock_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Bullock</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Spin_Bullock_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of virial properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: thisNode, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Spin_Bullock_Initialize()

    ! Increment property count if we are outputting virial data.
    if (outputSpinBullock) then
       doublePropertyCount=doublePropertyCount+1
       if (vectorSpinAvailable) doublePropertyCount=doublePropertyCount+3
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Bullock_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Spin_Bullock</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Bullock</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Spin_Bullock(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store virial properties in the \glc\ output file buffers.
    use Galacticus_Nodes        , only : treeNode                 , nodeComponentBasic  , nodeComponentSpin
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass, darkMatterProfileDMO
    use Kind_Numbers
    use Multi_Counters
    implicit none
    double precision                           , intent(in   ) :: time
    type            (treeNode                 ), intent(inout) :: thisNode
    integer                                    , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                        integerProperty
    integer         (kind=kind_int8           ), intent(inout) :: integerBuffer    (:,:)
    double precision                           , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter             ), intent(inout) :: instance
    class           (darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass), pointer       :: darkMatterProfileDMO_
    class           (nodeComponentBasic       ), pointer       :: basic
    class           (nodeComponentSpin        ), pointer       :: spin
    double precision                           , dimension(3)  :: spinVectorUnit
    double precision                                           :: spinBullock
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Spin_Bullock_Initialize()

    ! Store property data if we are outputting virial data.
    if (outputSpinBullock) then
       darkMatterHaloScale_  => darkMatterHaloScale       ()
       darkMatterProfileDMO_ => darkMatterProfileDMO      ()
       basic                 => thisNode            %basic()
       spinBullock           =  +Dark_Matter_Halo_Angular_Momentum  (thisNode,darkMatterProfileDMO_) &
            &                   /sqrt(2.0d0)                                                         &
            &                   /basic               %mass          (                              ) &
            &                   /darkMatterHaloScale_%virialRadius  (thisNode                      ) &
            &                   /darkMatterHaloScale_%virialVelocity(thisNode                      )
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=spinBullock
       if (vectorSpinAvailable) then
          spin                                                              =>  thisNode%spin      ()
          spinVectorUnit                                                    =  +spin    %spinVector() &
               &                                                               /spin    %spin      ()
          doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+3) =  +spinBullock           &
               &                                                               *spinVectorUnit
          doubleProperty                                                    =  +doubleProperty        &
               &                                                               +3
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Bullock

end module Galacticus_Output_Trees_Spin_Bullock
