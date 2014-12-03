!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a class implementing accretion of gas from the \gls{igm} onto halos.

module Accretion_Halos
  !% Implements a class implementing accretion of gas from the \gls{igm} onto halos.
  use ISO_Varying_String
  use Galacticus_Nodes
  use Abundances_Structure
  use Chemical_Abundances_Structure
  !# <include directive="accretionHalo" type="functionModules" >
  include 'accretionHalo.functionModules.inc'
  !# </include>
  private
  public :: Accretion_Halos_Hot_Halo_Output, Accretion_Halos_Hot_Halo_Output_Count, Accretion_Halos_Hot_Halo_Output_Names

  ! Options controlling whether hot, cold, or total accretion is required.
  integer, public, parameter :: accretionModeTotal=0
  integer, public, parameter :: accretionModeHot  =1
  integer, public, parameter :: accretionModeCold =2

  ! Options controlling output.
  logical                    :: outputHaloAccretionMode, accretionHalosOutputInitialized

  !# <include directive="accretionHalo" type="function" >
  !#  <descriptiveName>Accretion Onto Halos</descriptiveName>
  !#  <description>Class providing rates of accretion of gas from the \gls{igm} onto halos.</description>
  !#  <default>simple</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>no</stateful>
  !#  <method name="accretionRate" >
  !#   <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of mass from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="accretedMass" >
  !#   <description>Returns the mass (in units of $M_\odot$) of accreted from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}. Used to initialize nodes.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="failedAccretionRate" >
  !#   <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of failed accretion of mass from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="failedAccretedMass" >
  !#   <description>Returns the mass (in units of $M_\odot$) of that failed to accrete from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}. Used to initialize nodes.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="accretionRateMetals" >
  !#   <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of metals from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}.</description>
  !#   <type>type(abundances)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="accretedMassMetals" >
  !#   <description>Returns the mass of metals (in units of $M_\odot$) of accreted from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}. Used to initialize nodes.</description>
  !#   <type>type(abundances)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="accretionRateChemicals" >
  !#   <description>Returns the rate (in units of $M_\odot$ Gyr$^{-1}$) of accretion of chemicals from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}.</description>
  !#   <type>type(chemicalAbundances)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  !#  <method name="accretedMassChemicals" >
  !#   <description>Returns the mass of chemicals (in units of $M_\odot$) of accreted from the \gls{igm} onto {\tt node} in the given {\tt accretionMode}. Used to initialize nodes.</description>
  !#   <type>type(chemicalAbundances)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>integer       , intent(in   )          :: accretionMode</argument>
  !#  </method>
  include 'accretionHalo.type.inc'
  !# </include>

  subroutine Accretion_Halos_Output_Initialize()
    !% Initialize output in the halo accretion module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.accretionHalosOutputInitialized) then
       !$omp critical(Accretion_Halos_Output_Initialization)
       if (.not.accretionHalosOutputInitialized) then
          ! Get options controlling output.
          !@ <inputParameter>
          !@   <name>outputHaloAccretionMode</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Determines whether or not halo accretion rates are output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHaloAccretionMode',outputHaloAccretionMode,defaultValue=.false.)
          accretionHalosOutputInitialized=.true.
       end if
       !$omp end critical(Accretion_Halos_Output_Initialization)
    end if
    return
  end subroutine Accretion_Halos_Output_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Accretion_Halos_Hot_Halo_Output_Names</unitName>
  !#  <sortName>Accretion_Halos_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Accretion_Halos_Hot_Halo_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode            )              , intent(inout), pointer :: thisNode
    double precision                                    , intent(in   )          :: time
    integer                                             , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*               ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                                          integerPropertyComments, integerPropertyNames
    double precision                      , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Accretion_Halos_Output_Initialize()

    if (outputHaloAccretionMode) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>haloAccretionHotModeFraction</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Fraction of halo accretion rate occuring via the hot mode.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>hotHalo</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='haloAccretionHotModeFraction'
       doublePropertyComments(doubleProperty)='Fraction of halo accretion rate occuring via the hot mode.'
       doublePropertyUnitsSI (doubleProperty)=1.0d0
    end if
    return
  end subroutine Accretion_Halos_Hot_Halo_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Accretion_Halos_Hot_Halo_Output_Count</unitName>
  !#  <sortName>Accretion_Halos_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Accretion_Halos_Hot_Halo_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of hot halo cooling properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    double precision                      , intent(in   )          :: time
    integer                               , intent(inout)          :: doublePropertyCount  , integerPropertyCount
    integer                               , parameter              :: propertyCount      =1

    ! Initialize the module.
    call Accretion_Halos_Output_Initialize()

    if (outputHaloAccretionMode) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Accretion_Halos_Hot_Halo_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Accretion_Halos_Hot_Halo_Output</unitName>
  !#  <sortName>Accretion_Halos_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Accretion_Halos_Hot_Halo_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision                    , intent(in   )          :: time
    type            (treeNode          ), intent(inout), pointer :: thisNode
    integer                             , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                          integerProperty
    integer         (kind=kind_int8    ), intent(inout)          :: integerBuffer    (:,:)
    double precision                    , intent(inout)          :: doubleBuffer     (:,:)
    class           (accretionHaloClass)               , pointer :: accretionHalo_
    double precision                                             :: accretionRateHot      ,accretionRateTotal

    ! Initialize the module.
    call Accretion_Halos_Output_Initialize()
    if (outputHaloAccretionMode) then
       doubleProperty=doubleProperty+1
       accretionHalo_ => accretionHalo()
       accretionRateHot  =accretionHalo_%accretionRate(thisNode,accretionModeHot  )
       accretionRateTotal=accretionHalo_%accretionRate(thisNode,accretionModeTotal)
       if (accretionRateTotal /= 0.0d0) then
          doubleBuffer(doubleBufferCount,doubleProperty)=accretionRateHot/accretionRateTotal
       else
          doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
       end if
    end if
    return
  end subroutine Accretion_Halos_Hot_Halo_Output

end module Accretion_Halos
