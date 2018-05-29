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

!% Contains a module which extends the standard implementation of basic component to track the
!% Bertschinger mass.

module Node_Component_Basic_Standard_Extended
  !% Extends the standard implementation of basic component to track the Bertschinger mass.
  use Virial_Density_Contrast
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Basic_Standard_Extended_Initialize, Node_Component_Basic_Standard_Extended_Node_Merger , &
       &    Node_Component_Basic_Standard_Extended_Promote   , Node_Component_Basic_Standard_Extended_Rate_Compute, &
       &    Node_Component_Basic_Standard_Extended_Scale_Set , Node_Component_Basic_Extended_Bindings

  !# <component>
  !#  <class>basic</class>
  !#  <name>standardExtended</name>
  !#  <extends>
  !#   <class>basic</class>
  !#   <name>standard</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massBertschinger</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" /> 
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="Bertschinger mass of the node, assuming universal baryon fraction."/>
  !#   </property>
  !#   <property>
  !#     <name>accretionRateBertschinger</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>radiusTurnaround</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>radiusTurnaroundGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Options controlling spherical collapse model to use.
  integer                                        :: nodeComponentBasicExtendedSphericalCollapseType
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeLambda         =0
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeDE             =1
  integer                            , parameter :: nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998=2

  ! Initialization state.
  logical                                          :: moduleInitialized                                  =.false.

  ! Virial density contrast object.
  logical                                          :: virialDensityContrastInitialized                   =.false.
  class  (virialDensityContrastClass), allocatable :: virialDensityContrast_
  !$omp threadprivate(virialDensityContrast_,virialDensityContrastInitialized)


contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Basic_Extended_Bindings</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Basic_Extended_Bindings()
    !% Initializes the ``extended'' implementation of the basic component.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string                    ) :: nodeComponentBasicExtendedSphericalCollapseTypeText
    type(nodeComponentBasicStandardExtended) :: basic

    ! Initialize the bindings.
    if (.not.moduleInitialized) then
       !$omp critical (Node_Component_Basic_Extended_Bindings)
       if (.not.moduleInitialized) then
          !# <inputParameter>
          !#   <name>nodeComponentBasicExtendedSphericalCollapseType</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('matterLambda')</defaultValue>
          !#   <description>The type of spherical collapse model to assume in the extended basic node component class.</description>
          !#   <group>cosmology</group>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>nodeComponentBasicExtendedSphericalCollapseTypeText</variable>
          !# </inputParameter>
          select case (char(nodeComponentBasicExtendedSphericalCollapseTypeText))
          case ('matterLambda'    )
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeLambda
          case ('matterDarkEnergy')
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeDE
          case ('bryanNorman')
             nodeComponentBasicExtendedSphericalCollapseType=nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998
          end select
          ! Bind deferred functions.
          call basic%massBertschingerFunction(Node_Component_Basic_Extended_Mass_Bertschinger)
          call basic%radiusTurnaroundFunction(Node_Component_Basic_Extended_Radius_Turnaround)
          ! Record that the module is now initialize.
          moduleInitialized=.true.
       end if
       !$omp end critical (Node_Component_Basic_Extended_Bindings)
    end if
    return
  end subroutine Node_Component_Basic_Extended_Bindings
  
  subroutine Node_Component_Basic_Extended_Bertschinger_Solver(self)
    !% Compute the Bertschinger mass and turnaround radii
    use Dark_Matter_Profile_Mass_Definitions
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profile_Scales
    use Cosmology_Functions
    use Cosmology_Parameters
    implicit none
    class           (nodeComponentBasicStandardExtended), intent(inout) :: self
    type            (treeNode                          ), pointer       :: selfNode
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: radiusVirial

    ! Initialize virial density contrast objects.
    if (.not.virialDensityContrastInitialized) then
       cosmologyParameters_ => cosmologyParameters()
       cosmologyFunctions_  => cosmologyFunctions ()
       select case (nodeComponentBasicExtendedSphericalCollapseType)
       case (nodeComponentBasicExtendedSphericalCollapseTypeLambda         )
          allocate(virialDensityContrastSphericalCollapseMatterLambda :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastSphericalCollapseMatterLambda)
             virialDensityContrast_=virialDensityContrastSphericalCollapseMatterLambda(                     cosmologyFunctions_)
          end select
       case (nodeComponentBasicExtendedSphericalCollapseTypeDE             )
          allocate(virialDensityContrastSphericalCollapseMatterDE     :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastSphericalCollapseMatterDE    )
             virialDensityContrast_=virialDensityContrastSphericalCollapseMatterDE    (                     cosmologyFunctions_)
          end select
       case (nodeComponentBasicExtendedSphericalCollapseTypeBryanNorman1998)
          allocate(virialDensityContrastBryanNorman1998               :: virialDensityContrast_)
          select type (virialDensityContrast_)
          type is (virialDensityContrastBryanNorman1998              )
             virialDensityContrast_=virialDensityContrastBryanNorman1998              (cosmologyParameters_,cosmologyFunctions_)
          end select
       end select
       virialDensityContrastInitialized=.true.
    end if
    ! Compute Bertschinger mass and turnaround radius.
    selfNode => self%hostNode
    call self%massBertschingerSet(                                                                                                           &
         &                        Dark_Matter_Profile_Mass_Definition(                                                                       &
         &                                                                   selfNode                                                      , &
         &                                                                   virialDensityContrast_%densityContrast(                         &
         &                                                                                                          self%mass            (), &
         &                                                                                                          self%timeLastIsolated()  &
         &                                                                                                         )                       , &
         &                                                            radius=radiusVirial                                                    &
         &                                                           )                                                                       &
         &                       )
    call self%radiusTurnaroundSet(virialDensityContrast_%turnAroundOverVirialRadii(time=self%time())*radiusVirial)
    return
  end subroutine Node_Component_Basic_Extended_Bertschinger_Solver
  
  double precision function Node_Component_Basic_Extended_Mass_Bertschinger(self)
    !% Return the Bertschinger mass.
    implicit none
    class(nodeComponentBasicStandardExtended), intent(inout) :: self
    
    if (self%massBertschingerValue() <= 0.0d0) call Node_Component_Basic_Extended_Bertschinger_Solver(self)
    Node_Component_Basic_Extended_Mass_Bertschinger=self%massBertschingerValue()
    return
  end function Node_Component_Basic_Extended_Mass_Bertschinger

  double precision function Node_Component_Basic_Extended_Radius_Turnaround(self)
    !% Return the turnaround radius.
    implicit none
    class(nodeComponentBasicStandardExtended), intent(inout) :: self
    
    if (self%radiusTurnaroundValue() <= 0.0d0) call Node_Component_Basic_Extended_Bertschinger_Solver(self)
    Node_Component_Basic_Extended_Radius_Turnaround=self%radiusTurnaroundValue()
    return
  end function Node_Component_Basic_Extended_Radius_Turnaround

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Standard_Extended_Initialize(node)
    !% Set the mass accretion rate for {\normalfont \ttfamily node}.
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
    !% Return the unresolved mass for {\normalfont \ttfamily node}.
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

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Basic_Standard_Extended_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    implicit none
    type     (treeNode          ), intent(inout), pointer :: node
    logical                      , intent(in   )          :: odeConverged
    logical                      , intent(inout)          :: interrupt
    procedure(                  ), intent(inout), pointer :: interruptProcedure
    integer                      , intent(in   )          :: propertyType
    class    (nodeComponentBasic)               , pointer :: basic
    !GCC$ attributes unused :: interrupt, interruptProcedure, odeConverged
    
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

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Basic_Standard_Extended_Scale_Set(node)
    !% Set scales for properties in the standard implementation of the basic component.
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

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Basic_Standard_Extended_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Basic_Standard_Extended_Node_Merger(node)
    !% Switch off accretion of new mass onto this node once it becomes a satellite.
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    class(nodeComponentBasic)               , pointer :: basic

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

   !# <nodePromotionTask>
   !#  <unitName>Node_Component_Basic_Standard_Extended_Promote</unitName>
   !# </nodePromotionTask>
   subroutine Node_Component_Basic_Standard_Extended_Promote(node)
     !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply
     !% update the mass of {\normalfont \ttfamily node} to be that of its parent.
     use Galacticus_Error
     implicit none
     type (treeNode          ), intent(inout), pointer :: node
     type (treeNode          )               , pointer :: nodeParent
     class(nodeComponentBasic)               , pointer :: basicParent, basic

     ! Get the basic component.
     basic => node%basic()
     ! Ensure that it is of the standard class.
     select type (basic)
     class is (nodeComponentBasicStandardExtended)
        ! Get the parent node and its basic component.
        nodeParent  => node      %parent
        basicParent => nodeParent%basic ()
        ! Ensure the two halos exist at the same time.
        if (basic%time() /= basicParent%time())                                                                  &
             & call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
        ! Adjust the mass to that of the parent node.
        call basic%massBertschingerSet          (basicParent%massBertschinger          ())
        ! Adjust the accretion rate to that of the parent node.
        call basic%accretionRateBertschingerSet (basicParent%accretionRateBertschinger ())
        ! Adjust the turnaround radius to that of the parent node.
        call basic%radiusTurnaroundSet          (basicParent%radiusTurnaround          ())
        ! Adjust the accretion rate to that of the parent node.
        call basic%radiusTurnaroundGrowthRateSet(basicParent%radiusTurnaroundGrowthRate())
     end select
     return
   end subroutine Node_Component_Basic_Standard_Extended_Promote

end module Node_Component_Basic_Standard_Extended
