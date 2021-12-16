!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which extends the standard implementation of basic component to track the
Bertschinger mass.
!!}

module Node_Component_Basic_Standard_Extended
  !!{
  Extends the standard implementation of basic component to track the Bertschinger mass.
  !!}
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass
  implicit none
  private
  public :: Node_Component_Basic_Standard_Extended_Initialize  , Node_Component_Basic_Standard_Extended_Node_Merger , &
       &    Node_Component_Basic_Standard_Extended_Scale_Set   , Node_Component_Basic_Extended_Bindings             , &
       &    Node_Component_Basic_Extended_Thread_Initialize    , Node_Component_Basic_Extended_Thread_Uninitialize  , &
       &    Node_Component_Basic_Standard_Extended_Rate_Compute, Node_Component_Basic_Extended_State_Store          , &
       &    Node_Component_Basic_Extended_State_Restore

  !![
  <component>
   <class>basic</class>
   <name>standardExtended</name>
   <extends>
    <class>basic</class>
    <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>massBertschinger</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
      <classDefault>-1.0d0</classDefault>
      <output unitsInSI="massSolar" comment="Bertschinger mass of the node, assuming universal baryon fraction."/>
    </property>
    <property>
      <name>accretionRateBertschinger</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>radiusTurnaround</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
      <classDefault>-1.0d0</classDefault>
    </property>
    <property>
      <name>radiusTurnaroundGrowthRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>-1.0d0</classDefault>
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(cosmologyParametersClass ), pointer :: cosmologyParameters_
  class(cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_
  class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_
  !$omp threadprivate(cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_)

  ! Options controlling spherical collapse model to use.
  integer                                        :: nodeComponentBasicExtendedSphericalCollapseType                 , nodeComponentBasicExtendedSphericalCollapseEnergyFixedAt
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeLambda         =0
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeDE             =1
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998=2

  ! Virial density contrast object.
  logical                                          :: virialDensityContrastInitialized                            =.false.
  class  (virialDensityContrastClass), allocatable :: virialDensityContrast_
  !$omp threadprivate(virialDensityContrast_,virialDensityContrastInitialized)


contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Basic_Extended_Bindings</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Basic_Extended_Bindings(parameters_)
    !!{
    Initializes the ``extended'' implementation of the basic component.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasicStandardExtended
    use :: ISO_Varying_String        , only : var_str                                       , varying_string , char
    use :: Input_Parameters          , only : inputParameter                                , inputParameters
    use :: Spherical_Collapse_Solvers, only : enumerationCllsnlssMttrDarkEnergyFixedAtEncode
    implicit none
    type(inputParameters                   ), intent(inout) :: parameters_
    type(varying_string                    )                :: nodeComponentBasicExtendedSphericalCollapseTypeText, nodeComponentBasicExtendedSphericalCollapseEnergyFixedAtText
    type(nodeComponentBasicStandardExtended)                :: basic

    !![
    <inputParameter>
      <name>nodeComponentBasicExtendedSphericalCollapseType</name>
      <defaultValue>var_str('matterLambda')</defaultValue>
      <description>The type of spherical collapse model to assume in the extended basic node component class.</description>
      <source>parameters_</source>
      <variable>nodeComponentBasicExtendedSphericalCollapseTypeText</variable>
    </inputParameter>
    !!]
    select case (char(nodeComponentBasicExtendedSphericalCollapseTypeText))
    case ('matterLambda'    )
       nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeLambda
    case ('matterDarkEnergy')
       nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeDE
       !![
       <inputParameter>
         <name>nodeComponentBasicExtendedSphericalCollapseEnergyFixedAt</name>
         <defaultValue>var_str('turnaround')</defaultValue>
         <description>Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
           ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
           \citealt{percival_cosmological_2005}; \S8.).</description>
         <source>parameters_</source>
         <variable>nodeComponentBasicExtendedSphericalCollapseEnergyFixedAtText</variable>
       </inputParameter>
       !!]
       nodeComponentBasicExtendedSphericalCollapseEnergyFixedAt=enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(nodeComponentBasicExtendedSphericalCollapseEnergyFixedAtText),includesPrefix=.false.)
    case ('bryanNorman')
       nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998
    end select
    ! Bind deferred functions.
    call basic%massBertschingerFunction(Node_Component_Basic_Extended_Mass_Bertschinger)
    call basic%radiusTurnaroundFunction(Node_Component_Basic_Extended_Radius_Turnaround)
    return
  end subroutine Node_Component_Basic_Extended_Bindings

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Basic_Extended_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Basic_Extended_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent   , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultBasicComponent
    use :: Input_Parameters, only : inputParameter       , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultBasicComponent%standardExtendedIsActive()) then
       !![
       <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters_"/>
       <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters_"/>
       <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters_"/>
       !!]
       call nodePromotionEvent%attach(defaultBasicComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentBasicExtended')
    end if
    return
  end subroutine Node_Component_Basic_Extended_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Basic_Extended_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Basic_Extended_Thread_Uninitialize()
    !!{
    Uninitializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none

    if (defaultBasicComponent%standardExtendedIsActive()) then
       !![
       <objectDestructor name="cosmologyParameters_" />
       <objectDestructor name="cosmologyFunctions_"  />
       <objectDestructor name="darkMatterProfileDMO_"/>
       !!]
       if (nodePromotionEvent%isAttached(defaultBasicComponent,nodePromotion)) call nodePromotionEvent%detach(defaultBasicComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Basic_Extended_Thread_Uninitialize

  subroutine Node_Component_Basic_Extended_Bertschinger_Solver(self)
    !!{
    Compute the Bertschinger mass and turnaround radii
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasicStandardExtended  , treeNode
    use :: Virial_Density_Contrast             , only : virialDensityContrastBryanNorman1998, virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy
    implicit none
    class           (nodeComponentBasicStandardExtended), intent(inout) :: self
    type            (treeNode                          ), pointer       :: selfNode
    double precision                                                    :: radiusVirial

    ! Initialize virial density contrast objects.
    if (.not.virialDensityContrastInitialized) then
       select case (nodeComponentBasicExtendedSphericalCollapseType)
       case (nodeComponentBasicExtendedSphericalCollapseTypeLambda         )
          allocate(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)
             virialDensityContrast_=virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(.true.                                                         ,cosmologyFunctions_)
          end select
       case (nodeComponentBasicExtendedSphericalCollapseTypeDE                    )
          allocate(virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy      :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy     )
             virialDensityContrast_=virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy     (.true.,nodeComponentBasicExtendedSphericalCollapseEnergyFixedAt,cosmologyFunctions_)
          end select
       case (nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998)
          allocate(virialDensityContrastBryanNorman1998                           :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastBryanNorman1998                           )
             virialDensityContrast_=virialDensityContrastBryanNorman1998                          (       cosmologyParameters_                                    ,cosmologyFunctions_)
          end select
       end select
       virialDensityContrastInitialized=.true.
    end if
    ! Compute Bertschinger mass and turnaround radius.
    selfNode => self%hostNode
    call self%massBertschingerSet(                                                                                                                           &
         &                        Dark_Matter_Profile_Mass_Definition(                                                                                       &
         &                                                                                   selfNode                                                      , &
         &                                                                                   virialDensityContrast_%densityContrast(                         &
         &                                                                                                                          self%mass            (), &
         &                                                                                                                          self%timeLastIsolated()  &
         &                                                                                                                         )                       , &
         &                                                            radius                =radiusVirial                                                  , &
         &                                                            cosmologyParameters_  =cosmologyParameters_                                          , &
         &                                                            cosmologyFunctions_   =cosmologyFunctions_                                           , &
         &                                                            darkMatterProfileDMO_ =darkMatterProfileDMO_                                         , &
         &                                                            virialDensityContrast_=virialDensityContrast_                                          &
         &                                                           )                                                                                       &
         &                       )
    call self%radiusTurnaroundSet(virialDensityContrast_%turnAroundOverVirialRadii(mass=self%mass(),time=self%time())*radiusVirial)
    return
  end subroutine Node_Component_Basic_Extended_Bertschinger_Solver

  double precision function Node_Component_Basic_Extended_Mass_Bertschinger(self)
    !!{
    Return the Bertschinger mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasicStandardExtended
    implicit none
    class(nodeComponentBasicStandardExtended), intent(inout) :: self

    if (self%massBertschingerValue() <= 0.0d0) call Node_Component_Basic_Extended_Bertschinger_Solver(self)
    Node_Component_Basic_Extended_Mass_Bertschinger=self%massBertschingerValue()
    return
  end function Node_Component_Basic_Extended_Mass_Bertschinger

  double precision function Node_Component_Basic_Extended_Radius_Turnaround(self)
    !!{
    Return the turnaround radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasicStandardExtended
    implicit none
    class(nodeComponentBasicStandardExtended), intent(inout) :: self

    if (self%radiusTurnaroundValue() <= 0.0d0) call Node_Component_Basic_Extended_Bertschinger_Solver(self)
    Node_Component_Basic_Extended_Radius_Turnaround=self%radiusTurnaroundValue()
    return
  end function Node_Component_Basic_Extended_Radius_Turnaround

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Basic_Standard_Extended_Initialize</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Basic_Standard_Extended_Initialize(node)
    !!{
    Set the mass accretion rate for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandardExtended, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    type            (treeNode          )               , pointer :: nodeChild          , nodeParent
    class           (nodeComponentBasic)               , pointer :: basic              , basicChild, &
         &                                                          basicParent
    double precision                                             :: massUnresolved     , deltaTime , &
         &                                                          progenitorMassTotal

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard extended class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Determine if this node has a descendent.
       if (.not.associated(node%parent)) then
          ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
          ! progenitor, if it has one.
          nodeChild => node%firstChild
          if (associated(nodeChild)) then
             ! Get the basic component of the child node.
             basicChild => nodeChild%basic()
             ! Ensure the child has a mass growth rate computed.
             call Node_Component_Basic_Standard_Extended_Initialize(nodeChild)
             ! Get the growth rate of the child.
             call basic% accretionRateBertschingerSet(basicChild% accretionRateBertschinger())
             call basic%radiusTurnaroundGrowthRateSet(basicChild%radiusTurnaroundGrowthRate())
          else
             ! Parentless node has no child - set a zero growth rate.
             call basic%accretionRateBertschingerSet (0.0d0                                  )
             call basic%radiusTurnaroundGrowthRateSet(0.0d0                                  )
          end if
       else
          ! Get the parent node.
          nodeParent => node%parent
          ! Get the basic component of the parent node.
          basicParent => nodeParent%basic()
          ! Compute the unresolved mass.
          massUnresolved=Node_Component_Basic_Standard_Extended_Unresolved_Mass(nodeParent)
          if (massUnresolved > 0.0d0) then
             ! Positive mass growth - assume this occurs entirely in the main progenitor.
             if (node%isPrimaryProgenitor()) then
                ! Main progenitor - compute required growth rate.
                deltaTime=basicParent%time()-basic%time()
                if (deltaTime > 0.0d0) then
                   call basic% accretionRateBertschingerSet(massUnresolved/deltaTime)
                   call basic%radiusTurnaroundGrowthRateSet((basicParent%radiusTurnaround()-basic%radiusTurnaround())/deltaTime)
                end if
             else
                ! Non-main progenitor - assume zero growth rate.
                call basic%             accretionRateSet(0.0d0)
                call basic%radiusTurnaroundGrowthRateSet(0.0d0)
             end if
          else
             ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
             ! Compute the total mass in progenitors.
             progenitorMassTotal=basicParent%massBertschinger()-massUnresolved
             ! Compute the time available for accretion.
             deltaTime=basicParent%time()-basic%time()
             ! Compute mass growth rate.
             if (deltaTime > 0.0d0) then
                call basic% accretionRateBertschingerSet((massUnresolved/deltaTime)*(basic%massBertschinger()/progenitorMassTotal))
                call basic%radiusTurnaroundGrowthRateSet((basicParent%radiusTurnaround()-basic%radiusTurnaround())/deltaTime)
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Initialize

  double precision function Node_Component_Basic_Standard_Extended_Unresolved_Mass(node)
    !!{
    Return the unresolved mass for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: child
    class(nodeComponentBasic)               , pointer :: basicChild, basic

    ! Get the basic component.
    basic => node%basic()
    ! Initialize the unresolved mass to the mass of the current node's basic component.
    Node_Component_Basic_Standard_Extended_Unresolved_Mass=basic%massBertschinger()
    ! Remove the mass of all child nodes.
    child => node%firstChild
    do while (associated(child))
       basicChild                                             => child%basic()
       Node_Component_Basic_Standard_Extended_Unresolved_Mass =  Node_Component_Basic_Standard_Extended_Unresolved_Mass-basicChild%massBertschinger()
       child                                                  => child%sibling
    end do
    return
  end function Node_Component_Basic_Standard_Extended_Unresolved_Mass

  !![
  <rateComputeTask>
   <unitName>Node_Component_Basic_Standard_Extended_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Basic_Standard_Extended_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute rates of change of properties in the standard implementation of the basic component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandardExtended, propertyTypeInactive, treeNode
    implicit none
    type     (treeNode          ), intent(inout)          :: node
    logical                      , intent(inout)          :: interrupt
    procedure(                  ), intent(inout), pointer :: interruptProcedure
    integer                      , intent(in   )          :: propertyType
    class    (nodeComponentBasic)               , pointer :: basic
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Mass rate of change is set to the accretion rate.
       call basic%massBertschingerRate(basic%accretionRateBertschinger ())
       ! Set the radius growth rate.
       call basic%radiusTurnaroundRate(basic%radiusTurnaroundGrowthRate())
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Basic_Standard_Extended_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Basic_Standard_Extended_Scale_Set(node)
    !!{
    Set scales for properties in the standard implementation of the basic component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandardExtended, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    double precision                    , parameter              :: scaleLengthRelative=1.0d-6
    double precision                    , parameter              :: scaleMassRelative  =1.0d-6
    class           (nodeComponentBasic)               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Set scale for mass.
       call basic%massBertschingerScale(basic%massBertschinger()*scaleMassRelative  )
       ! Set scale for raadius.
       call basic%radiusTurnaroundScale(basic%radiusTurnaround()*scaleLengthRelative)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Scale_Set

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Basic_Standard_Extended_Node_Merger</unitName>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Basic_Standard_Extended_Node_Merger(node)
    !!{
    Switch off accretion of new mass onto this node once it becomes a satellite.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandardExtended, treeNode
    implicit none
    type (treeNode          ), intent(inout) :: node
    class(nodeComponentBasic), pointer       :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandardExtended)
       ! Shut down mass accretion onto the halo and growth now that it is a satellite.
       call basic%accretionRateBertschingerSet (0.0d0)
       call basic%radiusTurnAroundGrowthRateSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Extended_Node_Merger
  
  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply
    update the mass of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentBasicStandardExtended, treeNode
    implicit none
    class(*                 ), intent(inout) :: self
    type (treeNode          ), intent(inout) :: node
    type (treeNode          ), pointer       :: nodeParent
    class(nodeComponentBasic), pointer       :: basicParent, basic
    !$GLC attributes unused :: self
    
    basic       => node      %basic ()
    nodeParent  => node      %parent
    basicParent => nodeParent%basic ()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time())                                                                  &
         & call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust mass, turnaround radius, and their growth rates to those of the parent node.
    call basic%massBertschingerSet          (basicParent%massBertschinger          ())
    call basic%accretionRateBertschingerSet (basicParent%accretionRateBertschinger ())
    call basic%radiusTurnaroundSet          (basicParent%radiusTurnaround          ())
    call basic%radiusTurnaroundGrowthRateSet(basicParent%radiusTurnaroundGrowthRate())
    return
  end subroutine nodePromotion

  !![
  <galacticusStateStoreTask>
   <unitName>Node_Component_Basic_Extended_State_Store</unitName>
  </galacticusStateStoreTask>
  !!]
  subroutine Node_Component_Basic_Extended_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentBasic -> extended',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="cosmologyParameters_ cosmologyFunctions_ darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine Node_Component_Basic_Extended_State_Store

  !![
  <galacticusStateRetrieveTask>
   <unitName>Node_Component_Basic_Extended_State_Restore</unitName>
  </galacticusStateRetrieveTask>
  !!]
  subroutine Node_Component_Basic_Extended_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentBasic -> extended',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="cosmologyParameters_ cosmologyFunctions_ darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine Node_Component_Basic_Extended_State_Restore

end module Node_Component_Basic_Standard_Extended
