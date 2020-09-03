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

!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Node_Component_Dark_Matter_Profile_Scale
  !% Implements a dark matter profile method that provides a scale radius.
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
  use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfileScale
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Rate_Compute     , Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize    , &
       &    Node_Component_Dark_Matter_Profile_Scale_Scale_Set        , Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize, Node_Component_Dark_Matter_Profile_Scale_Plausibility       , &
       &    Node_Component_Dark_Matter_Profile_Scale_Initialize

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scale</name>
  !#  <isDefault>true</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>scale</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>scaleGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>scaleIsLimited</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>.true.</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Objects used by this component.
  class           (darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
  class           (darkMatterProfileScaleRadiusClass  ), pointer :: darkMatterProfileScaleRadius_
  !$omp threadprivate(darkMatterHaloScale_,darkMatterProfileScaleRadius_)

  ! Parameters of the method.
  double precision                                               :: darkMatterProfileMaximumConcentration                  , darkMatterProfileMinimumConcentration

  ! Flag indicating whether scale radius data should be output when full merger trees are output.
  logical                                                        :: mergerTreeStructureOutputDarkMatterProfileScale

  ! Queriable dark matter profile object.
  type            (nodeComponentDarkMatterProfileScale)          :: darkMatterProfile

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize(parameters_)
    !% Initializes the ``scale'' implementation of the dark matter halo profile component.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    ! Get parameters.
    !# <inputParameter>
    !#   <name>darkMatterProfileMinimumConcentration</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>4.0d0</defaultValue>
    !#   <description>The minimum concentration allowed for dark matter profiles.</description>
    !#   <source>parameters_</source>
    !#   <type>double</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>darkMatterProfileMaximumConcentration</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>100.0d0</defaultValue>
    !#   <description>The maximum concentration allowed for dark matter profiles.</description>
    !#   <source>parameters_</source>
    !#   <type>double</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mergerTreeStructureOutputDarkMatterProfileScale</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Determines whether or not dark matter halo scale radius is included in outputs of merger trees.</description>
    !#   <source>parameters_</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    ! Bind the scale get function.
    call darkMatterProfile%scaleFunction(Node_Component_Dark_Matter_Profile_Scale_Scale)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters_"/>
       !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters_"/>
       call nodePromotionEvent%attach(defaultDarkMatterProfileComponent,nodePromotion,openMPThreadBindingAtLevel,label="nodeComponentDarkMatterProfileScale")
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       !# <objectDestructor name="darkMatterHaloScale_"         />
       !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
       call nodePromotionEvent%detach(defaultDarkMatterProfileComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize

  double precision function Node_Component_Dark_Matter_Profile_Scale_Scale(self)
    !% Return the scale radius in the dark matter halo profile.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    class           (nodeComponentDarkMatterProfileScale), intent(inout) :: self
    type            (treeNode                           ), pointer       :: selfNode
    double precision                                                     :: scaleLengthMaximum  , scaleLengthMinimum, &
         &                                                                  radiusVirial

    ! Set the scale if it isn't already set.
    selfNode => self%host()
    if (self%scaleValue() < 0.0d0) call self%scaleSet(darkMatterProfileScaleRadius_%radius(selfNode))
    if (self%scaleIsLimited()) then
       radiusVirial      =darkMatterHaloScale_%virialRadius(selfNode)
       scaleLengthMaximum=radiusVirial/darkMatterProfileMinimumConcentration
       scaleLengthMinimum=radiusVirial/darkMatterProfileMaximumConcentration
       Node_Component_Dark_Matter_Profile_Scale_Scale= &
            & min(                                     &
            &     scaleLengthMaximum    ,              &
            &     max(                                 &
            &         scaleLengthMinimum,              &
            &         self%scaleValue()                &
            &        )                                 &
            &    )
    else
       Node_Component_Dark_Matter_Profile_Scale_Scale=self%scaleValue()
    end if
    return
  end function Node_Component_Dark_Matter_Profile_Scale_Scale

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the rate of change of the scale radius.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, propertyTypeInactive, treeNode
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (                              ), intent(inout), pointer :: interruptProcedure
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: concentration       , growthRate
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
     ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScale)
       ! Find the concentration of this halo.
       concentration=darkMatterHaloScale_%virialRadius(node)/darkMatterProfile%scale()
       ! Find the growth rate and limit to ensure minimum and maximum concentrations are not exceeded.
       growthRate=darkMatterProfile%scaleGrowthRate()
       if (concentration <= darkMatterProfileMinimumConcentration) growthRate=min(growthRate,darkMatterHaloScale_%virialRadiusGrowthRate(node)/darkMatterProfileMinimumConcentration)
       if (concentration >= darkMatterProfileMaximumConcentration) growthRate=max(growthRate,darkMatterHaloScale_%virialRadiusGrowthRate(node)/darkMatterProfileMaximumConcentration)
       call darkMatterProfile%scaleRate(growthRate)
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Rate_Compute

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Plausibility</unitName>
  !#  <after>Node_Component_Basic_Standard_Plausibility</after>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility(node)
    !% Determines whether the dark matter profile is physically plausible for radius solving tasks.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type   (treeNode                      ), intent(inout) :: node
    class  (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    ! Return immediately if already non-plausible.
    if (.not.(node%isPhysicallyPlausible.and.node%isSolvable)) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScale)
       if (darkMatterProfile%scale() <= 0.0d0) node%isPhysicallyPlausible=.false.
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize(node)
    !% Initialize the scale radius of {\normalfont \ttfamily node}.
    use :: Galacticus_Error   , only : Galacticus_Error_Report
    use :: Galacticus_Nodes   , only : defaultDarkMatterProfileComponent, nodeComponentBasic, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, &
          &                            treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfileParent, darkMatterProfile, &
         &                                                                      darkMatterProfileWork
    class           (nodeComponentBasic            )               , pointer :: basicParent            , basic
    type            (treeNode                      )               , pointer :: nodeWork
    type            (mergerTreeWalkerAllNodes      )                         :: treeWalker
    double precision                                                         :: deltaTime              , radiusScale

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       ! Get the dark matter profile component - creating if if necessary.
       darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
       ! If the scale has already been initialized, no need to repeat this part.
       select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScale)
          if (darkMatterProfile%scaleValue() < 0.0d0) then
             ! Perform our own depth-first tree walk to set scales in all nodes of the tree. This is necessary as we require access
             ! to the parent scale to set scale growth rates, but must initialize scales in a strictly depth-first manner as some
             ! algorithms rely on knowing the progenitor structure of the tree to compute scale radii.
             treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.false.)
             do while (treeWalker%next(nodeWork))
                ! Get the scale radius - this will initialize the radius if necessary.
                darkMatterProfileWork => nodeWork             %darkMatterProfile(autoCreate=.true.)
                radiusScale           =  darkMatterProfileWork%scale            (                 )
             end do
          end if
       class default
          call Galacticus_Error_Report('unexpected class'//{introspection:location})
       end select
       ! Check if this node is the primary progenitor.
       if (node%isPrimaryProgenitor()) then
          ! It is, so compute the scale radius growth rate.
          ! Now compute the growth rate.
          basic       =>  node              %basic()
          basicParent =>  node       %parent%basic()
          deltaTime   =  +basicParent       %time () &
               &         -basic             %time ()
          if (deltaTime > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile(autoCreate=.true.)
             call darkMatterProfile%scaleGrowthRateSet(                                  &
                  &                                    (                                 &
                  &                                      darkMatterProfileParent%scale() &
                  &                                     -darkMatterProfile      %scale() &
                  &                                    )                                 &
                  &                                    /deltaTime                        &
                  &                                   )
          else
             call darkMatterProfile%scaleGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set scale radius growth rate to zero.
          call    darkMatterProfile%scaleGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Tree_Initialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the growth rate of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    class(*                             ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfileParent, darkMatterProfile
    class(nodeComponentBasic            ), pointer       :: basicParent            , basic
    !$GLC attributes unused :: self    

    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basic                   => node       %basic            ()
    basicParent             => node%parent%basic            ()
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the scale radius and its growth rate to that of the parent node.
    call darkMatterProfile%scaleSet          (darkMatterProfileParent%scale          ())
    call darkMatterProfile%scaleGrowthRateSet(darkMatterProfileParent%scaleGrowthRate())
    return
  end subroutine nodePromotion
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileScale)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(darkMatterProfile%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale
