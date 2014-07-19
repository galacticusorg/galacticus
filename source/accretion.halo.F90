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

!% Contains a module which implements calculations of baryonic accretion into halos.

module Accretion_Halos
  !% Implements calculations of baryonic accretion into halos.
  use ISO_Varying_String
  use Galacticus_Nodes
  use Abundances_Structure
  use Chemical_Abundances_Structure
  implicit none
  private
  public :: Halo_Baryonic_Accretion_Rate, Halo_Baryonic_Accreted_Mass, Halo_Baryonic_Failed_Accretion_Rate,&
       & Halo_Baryonic_Failed_Accreted_Mass, Halo_Baryonic_Accretion_Rate_Abundances, Halo_Baryonic_Accreted_Abundances,&
       & Halo_Baryonic_Accretion_Rate_Chemicals, Halo_Baryonic_Accreted_Chemicals, Accretion_Halos_Hot_Halo_Output,&
       & Accretion_Halos_Hot_Halo_Output_Count, Accretion_Halos_Hot_Halo_Output_Names

  ! Flag to indicate if this module has been initialized.
  logical                                                     :: accretionHalosInitialized                  =.false.

  ! Name of mass movement method used.
  type     (varying_string                       )            :: accretionHalosMethod

  ! Pointers to functions that return baryonic mass accretion rates/masses.
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer  :: Halo_Baryonic_Accretion_Rate_Get           =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer  :: Halo_Baryonic_Accreted_Mass_Get            =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer  :: Halo_Baryonic_Failed_Accretion_Rate_Get    =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer  :: Halo_Baryonic_Failed_Accreted_Mass_Get     =>null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer  :: Halo_Baryonic_Accreted_Abundances_Get      =>null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer  :: Halo_Baryonic_Accreted_Chemicals_Get       =>null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer  :: Halo_Baryonic_Accretion_Rate_Abundances_Get=>null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer  :: Halo_Baryonic_Accretion_Rate_Chemicals_Get =>null()
  abstract interface
     double precision function Halo_Baryonic_Accretion_Get_Template(thisNode,accretionMode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
       integer       , intent(in   )          :: accretionMode
     end function Halo_Baryonic_Accretion_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Abundances_Get_Template(thisNode,accretedAbundances,accretionMode)
       import treeNode
       import abundances
       type(treeNode  ), intent(inout), pointer :: thisNode
       type(abundances), intent(inout)          :: accretedAbundances
       integer         , intent(in   )          :: accretionMode
     end subroutine Halo_Baryonic_Abundances_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Chemicals_Get_Template(thisNode,accretedChemicals,accretionMode)
       import treeNode
       import chemicalAbundances
       type(treeNode          ), intent(inout), pointer :: thisNode
       type(chemicalAbundances), intent(inout)          :: accretedChemicals
       integer                 , intent(in   )          :: accretionMode
     end subroutine Halo_Baryonic_Chemicals_Get_Template
  end interface

  ! Option controlling whether accretion fractions are output.
  logical :: outputHaloAccretionMode, accretionHalosOutputInitialized=.false.

contains

  subroutine Accretion_Halos_Initialize
    !% Initalize the accretion disk module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="accretionHalosMethod" type="moduleUse">
    include 'accretion.halos.modules.inc'
    !# </include>
    implicit none

    if (.not.accretionHalosInitialized) then
       !$omp critical(accretionHalosInitialize)
       if (.not.accretionHalosInitialized) then
          ! Get the halo accretion method parameter.
          !@ <inputParameter>
          !@   <name>accretionHalosMethod</name>
          !@   <defaultValue>simple</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Selects which method should be used for accretion onto halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('accretionHalosMethod',accretionHalosMethod,defaultValue='simple')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="accretionHalosMethod" type="functionCall" functionType="void">
          !#  <functionArgs>accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get,Halo_Baryonic_Accreted_Abundances_Get,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accretion_Rate_Chemicals_Get,Halo_Baryonic_Accreted_Chemicals_Get</functionArgs>
          include 'accretion.halos.inc'
          !# </include>
          if     (.not.(     associated(Halo_Baryonic_Accretion_Rate_Get           ) &
               &        .and.associated(Halo_Baryonic_Accreted_Mass_Get            ) &
               &        .and.associated(Halo_Baryonic_Failed_Accretion_Rate_Get    ) &
               &        .and.associated(Halo_Baryonic_Failed_Accreted_Mass_Get     ) &
               &        .and.associated(Halo_Baryonic_Accretion_Rate_Abundances_Get) &
               &        .and.associated(Halo_Baryonic_Accreted_Abundances_Get      ) &
               &        .and.associated(Halo_Baryonic_Accretion_Rate_Chemicals_Get ) &
               &        .and.associated(Halo_Baryonic_Accreted_Chemicals_Get       ) &
               &       )                                                             &
               & ) call Galacticus_Error_Report('Accretion_Halos_Initialize','method ' //char(accretionHalosMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          accretionHalosInitialized=.true.
       end if
       !$omp end critical(accretionHalosInitialize)
    end if
    return
  end subroutine Accretion_Halos_Initialize

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

  double precision function Halo_Baryonic_Accretion_Rate(thisNode,accretionMode)
    !% Computes the rate of baryonic mass accretion (in $M_\odot$/Gyr) onto {\tt thisNode} from the intergalactic medium.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer       , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Accretion_Rate=Halo_Baryonic_Accretion_Rate_Get(thisNode,accretionMode)

    return
  end function Halo_Baryonic_Accretion_Rate

  double precision function Halo_Baryonic_Accreted_Mass(thisNode,accretionMode)
    !% Computes the mass of baryons accreted (in $M_\odot$) into {\tt thisNode} from the intergalactic medium. Used to initialize
    !% nodes.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer       , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Accreted_Mass=Halo_Baryonic_Accreted_Mass_Get(thisNode,accretionMode)

    return
  end function Halo_Baryonic_Accreted_Mass

  double precision function Halo_Baryonic_Failed_Accretion_Rate(thisNode,accretionMode)
    !% Computes the rate of failed baryonic mass accretion (in $M_\odot$/Gyr) onto {\tt thisNode} from the intergalactic medium.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer       , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accretion_Rate=Halo_Baryonic_Failed_Accretion_Rate_Get(thisNode,accretionMode)

    return
  end function Halo_Baryonic_Failed_Accretion_Rate

  double precision function Halo_Baryonic_Failed_Accreted_Mass(thisNode,accretionMode)
    !% Computes the mass of baryons that failed to accrete (in $M_\odot$) into {\tt thisNode} from the intergalactic medium. Used to initialize
    !% nodes.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer       , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accreted_Mass=Halo_Baryonic_Failed_Accreted_Mass_Get(thisNode,accretionMode)

    return
  end function Halo_Baryonic_Failed_Accreted_Mass

  subroutine Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances,accretionMode)
    !% Compute the rate of mass accretion of abundances (in $M_\odot/$Gyr) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretionRateAbundances
    integer         , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Abundances_Get(thisNode,accretionRateAbundances,accretionMode)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances

  subroutine Halo_Baryonic_Accreted_Abundances(thisNode,accretedAbundances,accretionMode)
    !% Compute the mass of abundances (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretedAbundances
    integer         , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Abundances_Get(thisNode,accretedAbundances,accretionMode)

    return
  end subroutine Halo_Baryonic_Accreted_Abundances

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals(thisNode,accretionRateChemicals,accretionMode)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretionRateChemicals
    integer                 , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Chemicals_Get(thisNode,accretionRateChemicals,accretionMode)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals

  subroutine Halo_Baryonic_Accreted_Chemicals(thisNode,accretedChemicals,accretionMode)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretedChemicals
    integer                 , intent(in   )          :: accretionMode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Chemicals_Get(thisNode,accretedChemicals,accretionMode)

    return
  end subroutine Halo_Baryonic_Accreted_Chemicals

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
    use Accretion_Halos_Options
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)
    double precision                                         :: accretionRateHot      ,accretionRateTotal
    ! Initialize the module.
    call Accretion_Halos_Output_Initialize()

    if (outputHaloAccretionMode) then
       doubleProperty=doubleProperty+1
       accretionRateHot  =Halo_Baryonic_Accretion_Rate(thisNode,accretionModeHot  )
       accretionRateTotal=Halo_Baryonic_Accretion_Rate(thisNode,accretionModeTotal)
       if (accretionRateTotal /= 0.0d0) then
          doubleBuffer(doubleBufferCount,doubleProperty)=accretionRateHot/accretionRateTotal
       else
          doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
       end if
    end if
    return
  end subroutine Accretion_Halos_Hot_Halo_Output

end module Accretion_Halos
