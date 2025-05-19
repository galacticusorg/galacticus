!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a merger tree build controller class which builds constrained trees.
!!}

  use :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field        , only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Linear_Growth                     , only : linearGrowthClass 
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerConstrained">
   <description>A merger tree build controller class which builds constrained trees.</description>
  </mergerTreeBuildController>
  !!]

  ! Enumeration for different fitting function types.
  !![
  <enumeration>
   <name>constructionOption</name>
   <description>Specifies option for constructing merger tree.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="constrainedBranchOnly"       />
   <entry label="constrainedAndMainBranchOnly"/>
   <entry label="allBranches"                 />
  </enumeration>
  !!]

  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerConstrained
     !!{     
     A merger tree build controller class which builds constrained trees.
     !!}
     private
     class           (cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_                          => null()
     class           (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbabilityUnconstrained_ => null(), mergerTreeBranchingProbabilityConstrained_ => null()
     class           (linearGrowthClass                  ), pointer :: linearGrowth_                                => null()
     class           (criticalOverdensityClass           ), pointer :: criticalOverdensity_                         => null()
     class           (cosmologicalMassVarianceClass      ), pointer :: cosmologicalMassVariance_                    => null()
     class           (mergerTreeMassResolutionClass      ), pointer :: mergerTreeMassResolution_                    => null()
     integer                                                        :: isConstrainedID
     type            (enumerationConstructionOptionType  )          :: constructionOption
     double precision                                               :: criticalOverdensityConstrained                        , varianceConstrained                                 , &
          &                                                            timeConstrained                                       , massConstrained                                     , &
          &                                                            redshiftConstrained
     type            (varying_string                     )          :: label                                                 , labelDescription
     integer                                                        :: labelID
   contains
     final     ::                               constrainedDestructor
     procedure :: control                    => constrainedControl
     procedure :: branchingProbabilityObject => constrainedBranchingProbabilityObject
     procedure :: nodesInserted              => constrainedNodesInserted
     procedure :: timeMaximum                => constrainedTimeMaximum
     procedure :: controlTimeMaximum         => constrainedControlTimeMaximum
  end type mergerTreeBuildControllerConstrained
  
  interface mergerTreeBuildControllerConstrained
     !!{
     Constructors for the {\normalfont \ttfamily constrained} merger tree build controller class.
     !!}
     module procedure constrainedConstructorParameters
     module procedure constrainedConstructorInternal
  end interface mergerTreeBuildControllerConstrained

contains

  function constrainedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily constrained} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters

    implicit none
    type            (mergerTreeBuildControllerConstrained)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    type            (varying_string                      )                :: constructionOption
    class           (mergerTreeBranchingProbabilityClass ), pointer       :: mergerTreeBranchingProbabilityUnconstrained_, mergerTreeBranchingProbabilityConstrained_
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                   ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass            ), pointer       :: criticalOverdensity_ 
    class           (cosmologicalMassVarianceClass       ), pointer       :: cosmologicalMassVariance_
    class           (mergerTreeMassResolutionClass       ), pointer       :: mergerTreeMassResolution_
    double precision                                                      :: criticalOverdensityConstrained              , varianceConstrained                       , &
         &                                                                   timeConstrained                             , massConstrained                           , &
         &                                                                   timePresent                                 , redshiftConstrained                       , &
         &                                                                   expansionFactor
    type            (varying_string                      )                :: label                                       , labelDescription
  
    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <defaultValue>var_str(' ')</defaultValue>
      <description>A label to apply to the constrained node.</description>
    </inputParameter>
    !!]
    if (label == '') label=' '
    if (trim(label) /= '') then
       !![
       <inputParameter>
         <name>labelDescription</name>
         <source>parameters</source>
         <description>A description of the label.</description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbabilityUnconstrained_" parameterName="mergerTreeBranchingProbabilityUnconstrained" source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbabilityConstrained_"   parameterName="mergerTreeBranchingProbabilityConstrained"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"                                                                                      source="parameters"/>
    <objectBuilder class="linearGrowth"                   name="linearGrowth_"                                                                                            source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"                                                                                     source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"                                                                                source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"       name="mergerTreeMassResolution_"                                                                                source="parameters"/>
     <inputParameter>
      <name>constructionOption</name>
      <source>parameters</source>
      <description>Controls which branches of the tree to build.</description>
    </inputParameter>
    !!]
    timePresent=cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    if      (parameters%isPresent('criticalOverdensityConstrained')) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('varianceConstrained')                                                                                                                       &
            & ) call Error_Report('both "criticalOverdensityConstrained" and "varianceConstrained" must be provided'                                      //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('redshiftConstrained'           )                                                                                                            &
            &  .or.                                                                                                                                                                   &
            &       parameters%isPresent('massConstrained'               )                                                                                                            &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>criticalOverdensityConstrained</name>
         <source>parameters</source>
         <description>The critical overdensity at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>varianceConstrained</name>
         <source>parameters</source>
         <description>The variance at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       massConstrained=self%cosmologicalMassVariance_%mass          (time               =timePresent                   ,rootVariance=sqrt(varianceConstrained))
       timeConstrained=self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=criticalOverdensityConstrained,mass        =     massConstrained     )
    else if (parameters%isPresent('redshiftConstrained           ')) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('massConstrained'    )                                                                                                                       &
            & ) call Error_Report('both "redshiftConstrained" and "massConstrained" must be provided'                                                     //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('criticalOverdensityConstrained')                                                                                                            &
            &  .or.                                                                                                                                                                   &
            &       parameters%isPresent('varianceConstrained'           )                                                                                                            &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>redshiftConstrained</name>
         <source>parameters</source>
         <description>The redshift at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>massConstrained</name>
         <source>parameters</source>
         <description>The halo mass at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       expansionFactor               =+cosmologyFunctions_      %expansionFactorFromRedshift(redshift       =redshiftConstrained                 )
       timeConstrained               =+cosmologyFunctions_      %cosmicTime                 (expansionFactor=expansionFactor                     )
       criticalOverdensityConstrained=+criticalOverdensity_     %value                      (time           =timeConstrained,mass=massConstrained)    &
            &                         /linearGrowth_            %value                      (time           =timeConstrained                     )
       varianceConstrained           =+cosmologicalMassVariance_%rootVariance               (time           =timePresent    ,mass=massConstrained)**2
    else
       criticalOverdensityConstrained=0.0d0
       varianceConstrained           =0.0d0
       timeConstrained               =0.0d0
       massConstrained               =0.0d0
       call Error_Report('must provide either [criticalOverdensityConstrained] and [varianceConstrained], or [timeConstrained] and [massConstrained]')
    end if
    self=mergerTreeBuildControllerConstrained(criticalOverdensityConstrained,varianceConstrained,enumerationConstructionOptionEncode(char(constructionOption),includesPrefix=.false.),label,labelDescription,mergerTreeBranchingProbabilityUnconstrained_,mergerTreeBranchingProbabilityConstrained_,cosmologyFunctions_,linearGrowth_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeMassResolution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbabilityUnconstrained_"/>
    <objectDestructor name="mergerTreeBranchingProbabilityConstrained_"  />
    <objectDestructor name="cosmologyFunctions_"                         />
    <objectDestructor name="linearGrowth_"                               />
    <objectDestructor name="criticalOverdensity_"                        />
    <objectDestructor name="cosmologicalMassVariance_"                   />
    <objectDestructor name="mergerTreeMassResolution_"                   />
    !!]
    return
  end function constrainedConstructorParameters

  function constrainedConstructorInternal(criticalOverdensityConstrained,varianceConstrained,constructionOption,label,labelDescription,mergerTreeBranchingProbabilityUnconstrained_,mergerTreeBranchingProbabilityConstrained_,cosmologyFunctions_,linearGrowth_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeMassResolution_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily constrained} merger tree build controller class.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type            (mergerTreeBuildControllerConstrained)                        :: self
    type            (enumerationConstructionOptionType   ), intent(in   )         :: constructionOption
    type            (varying_string                      ), intent(in   )         :: label                                       , labelDescription
    class           (mergerTreeBranchingProbabilityClass ), intent(in   ), target :: mergerTreeBranchingProbabilityUnconstrained_, mergerTreeBranchingProbabilityConstrained_
    class           (cosmologyFunctionsClass             ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                   ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass            ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass       ), intent(in   ), target :: cosmologicalMassVariance_
    class           (mergerTreeMassResolutionClass       ), intent(in   ), target :: mergerTreeMassResolution_
    double precision                                      , intent(in   )         :: criticalOverdensityConstrained              , varianceConstrained
    double precision                                                              :: timePresent                                 , expansionFactor
    !![
    <constructorAssign variables="criticalOverdensityConstrained, varianceConstrained, constructionOption, label, labelDescription, *mergerTreeBranchingProbabilityUnconstrained_, *mergerTreeBranchingProbabilityConstrained_, *cosmologyFunctions_, *linearGrowth_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeMassResolution_"/>
    !!]

    ! Find mass and time corresponding to the constraint point.
    timePresent             =self%cosmologyFunctions_      %cosmicTime                 (expansionFactor    =1.0d0                                                                          )
    self%massConstrained    =self%cosmologicalMassVariance_%mass                       (time               =timePresent                        ,rootVariance=sqrt(self%varianceConstrained))
    self%timeConstrained    =self%criticalOverdensity_     %timeOfCollapse             (criticalOverdensity=self%criticalOverdensityConstrained,mass        =     self%massConstrained     )
    expansionFactor         =self%cosmologyFunctions_      %expansionFactor            (time               =self%timeConstrained                                                           )
    self%redshiftConstrained=self%cosmologyFunctions_      %redshiftFromExpansionFactor(expansionFactor    =expansionFactor                                                                )
    ! Here we add a "meta-property" to the "basic" component of each node to store the status of whether a node is on the
    ! constrained branch or not.    
    !![
    <addMetaProperty component="basic" name="isConstrained" type="integer" id="self%isConstrainedID" isCreator="yes"/>
    !!]
    ! Register a label if required.
    if (trim(label) /= '') then
       self%labelID=nodeLabelRegister(char(label),char(labelDescription))
    else
       self%labelID=-1
    end if
    return
  end function constrainedConstructorInternal

  subroutine constrainedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily constrained} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerConstrained), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbabilityUnconstrained_"/>
    <objectDestructor name="self%mergerTreeBranchingProbabilityConstrained_"  />
    <objectDestructor name="self%cosmologyFunctions_"                         />
    <objectDestructor name="self%linearGrowth_"                               />
    <objectDestructor name="self%criticalOverdensity_"                        />
    <objectDestructor name="self%cosmologicalMassVariance_"                   />
    <objectDestructor name="self%mergerTreeMassResolution_"                   />
    !!]
    return
  end subroutine constrainedDestructor

  logical function constrainedControl(self,node,treeWalker_)
    !!{
    Apply control to merger tree building.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(mergerTreeBuildControllerConstrained), intent(inout)           :: self
    type (treeNode                            ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass               ), intent(inout), optional :: treeWalker_
    class(nodeComponentBasic                  )               , pointer  :: basic

    ! Always return true as we never want to halt tree building.
    constrainedControl=.true.
    ! Mark the root node as being on the constrained branch.
    if (.not.associated(node%parent)) then
       basic => node%basic()
       call basic%integerRank0MetaPropertySet(self%isConstrainedID,1)
    end if
    ! Enforce that the mass on the constrained branch can not be below the constrained mass at times after the constrained time.
    basic => node%basic()
    if     (                                                                                                &
         &   basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1                                   &
         &  .and.                                                                                           &
         &   basic%time                       (                    ) <  self%criticalOverdensityConstrained &
         &  .and.                                                                                           &
         &   basic%mass                       (                    ) <  self%massConstrained                &
         & ) call basic%massSet(self%massConstrained)
    ! Continue with this node if it meets the criteria.
    select case (self%constructionOption%ID)
    case (constructionOptionConstrainedBranchOnly       %ID)
       basic => node%basic()
       do while (constrainedControl.and.associated(node%parent).and. basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 0                                 )
          if (present(treeWalker_)) then
             constrainedControl=treeWalker_%next(node)
          else
             constrainedControl=.false.
          end if
          if (constrainedControl) basic => node%basic()
       end do
    case (constructionOptionConstrainedAndMainBranchOnly%ID)
       basic => node%basic()
       do while (constrainedControl.and.associated(node%parent).and.(basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 0 .and. .not.node%isOnMainBranch()))
          if (present(treeWalker_)) then
             constrainedControl=treeWalker_%next(node)
          else
             constrainedControl=.false.
          end if
          if (constrainedControl) basic => node%basic()
       end do
    case (constructionOptionAllBranches                 %ID)
       constrainedControl=.true.
    end select
    return
  end function constrainedControl

  double precision function constrainedTimeMaximum(self,node,massBranch,criticalOverdensityBranch,timeReference,insertNode)
    !!{
    Return the maximum allowed time for this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (mergerTreeBuildControllerConstrained), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: massBranch   , criticalOverdensityBranch, &
         &                                                                   timeReference
    logical                                               , intent(  out) :: insertNode
    class           (nodeComponentBasic                  ), pointer       :: basic
    logical                                                               :: isConstrained
    !$GLC attributes unused :: timeReference

    insertNode    =  .false.
    basic         => node %basic                      (                   )
    isConstrained =  basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    if (isConstrained .and. criticalOverdensityBranch < self%criticalOverdensityConstrained) then
       constrainedTimeMaximum=self%criticalOverdensityConstrained
    else
       constrainedTimeMaximum=huge(0.0d0)
    end if
    return
  end function constrainedTimeMaximum
  
  logical function constrainedControlTimeMaximum(self,node,massBranch,criticalOverdensityBranch,nodeIndex)
    !!{
    Control when the maximum time is reached.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: Kind_Numbers    , only : kind_int8
    use :: Nodes_Labels    , only : nodeLabelSet
    implicit none
    class           (mergerTreeBuildControllerConstrained), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: massBranch                      , criticalOverdensityBranch
    integer         (kind=kind_int8                      ), intent(inout)         :: nodeIndex
    type            (treeNode                            ), pointer               :: nodeParent                      , nodeConstrained              , &
         &                                                                           nodeUnconstrained
    class           (nodeComponentBasic                  ), pointer               :: basic                           , basicConstrained             , &
         &                                                                           basicUnconstrained
    double precision                                      , parameter             :: toleranceRelative        =1.0d-3
    double precision                                                              :: criticalOverdensityParent       , criticalOverdensityProgenitor, &
         &                                                                           massParent                      , massSecondary

    ! The constraint time has been reached. We need to insert the merger to the constrained halo mass. First, check if the halo
    ! has any progenitors - this is unlikely but could happen if another merger happened precisely at this time.
    if (associated(node%firstChild)) then
       ! Node has children - find the one on the constrained branch.
       nodeParent => node      %firstChild
       basic      => nodeParent%basic     ()
       do while (basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 0)
          nodeParent => nodeParent%sibling
          if (.not.associated(nodeParent)) call Error_Report('unable to locate constrained branch'//{introspection:location})
          basic      => nodeParent%basic ()
       end do
       criticalOverdensityParent =  basic%time()
       massParent                =  basic%mass()
    else
       ! Node has no children, so it is the parent for our merger.
       nodeParent                => node
       criticalOverdensityParent =  criticalOverdensityBranch
       massParent                =  massBranch
    end if
    ! Determine the time at which to insert the progenitors.
    criticalOverdensityProgenitor=max(self%criticalOverdensityConstrained,criticalOverdensityBranch*(1.0d0+toleranceRelative))
    ! Create progenitors for the parent.
    nodeIndex               =  nodeIndex+1
    nodeConstrained         => treeNode(nodeIndex,nodeParent%hostTree)
    basicConstrained        => nodeConstrained%basic(autoCreate=.true.)
    nodeConstrained %parent => nodeParent
    call basicConstrained%massSet                    (           self%massConstrained                )
    call basicConstrained%timeSet                    (                criticalOverdensityProgenitor  )
    call basicConstrained%integerRank0MetaPropertySet(self%isConstrainedID                         ,1)
    if (self%labelID > 0) call nodeLabelSet(self%labelID,nodeConstrained)
    massSecondary=+     massParent      &
         &        -self%massConstrained
    if (massSecondary > self%mergerTreeMassResolution_%resolution(node%hostTree)) then
       nodeIndex                 =  nodeIndex+1
       nodeUnconstrained         => treeNode(nodeIndex,nodeParent%hostTree)
       basicUnconstrained        => nodeUnconstrained%basic(autoCreate=.true.)
       nodeUnconstrained %parent => nodeParent
       call basicUnconstrained%massSet(massParent-self%massConstrained              )
       call basicUnconstrained%timeSet(                criticalOverdensityProgenitor)
       if (self%massConstrained > massSecondary) then
          ! Constrained node is the main progenitor.
          nodeParent       %firstChild => nodeConstrained
          nodeConstrained  %sibling    => nodeUnconstrained
       else
          ! Constrained node is the secondary progenitor.
          nodeParent       %firstChild => nodeUnconstrained
          nodeUnconstrained%sibling    => nodeConstrained
       end if
    else
       nodeParent%firstChild => nodeConstrained
    end if
    ! Return false indicating that the current node is finished, so building should continue from its progenitor nodes.
    constrainedControlTimeMaximum=.false.
    return
  end function constrainedControlTimeMaximum
  
  function constrainedBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (mergerTreeBranchingProbabilityClass ), pointer       :: mergerTreeBranchingProbability_
    class  (mergerTreeBuildControllerConstrained), intent(inout) :: self
    type   (treeNode                            ), intent(inout) :: node
    class  (nodeComponentBasic                  ), pointer       :: basic
    logical                                                      :: isConstrained

    basic         => node %basic                      (                   )
    isConstrained =  basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    if (isConstrained) then
       mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbabilityConstrained_
    else
       mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbabilityUnconstrained_
    end if
    return
  end function constrainedBranchingProbabilityObject

  subroutine constrainedNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (mergerTreeBuildControllerConstrained), intent(inout)           :: self
    type   (treeNode                            ), intent(inout)           :: nodeCurrent     , nodeProgenitor1
    type   (treeNode                            ), intent(inout), optional :: nodeProgenitor2
    logical                                      , intent(in   ), optional :: didBranch
    class  (nodeComponentBasic                  ), pointer                 :: basicCurrent    , basicProgenitor1, &
         &                                                                    basicProgenitor2
    logical                                                                :: isConstrained
    !$GLC attributes unused :: didBranch

    basicCurrent     => nodeCurrent    %basic                      (                    )
    basicProgenitor1 => nodeProgenitor1%basic                      (                    )
    isConstrained    =  basicCurrent   %integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    if (basicCurrent%time() >= self%criticalOverdensityConstrained) then
       if (isConstrained) then
          if (present(nodeProgenitor2)) then
             basicProgenitor2 => nodeProgenitor2%basic()
             if (basicProgenitor1%mass() > basicProgenitor2%mass()) then
                call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,1)
                call basicProgenitor2%integerRank0MetaPropertySet(self%isConstrainedID,0)
             else
                call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,0)
                call basicProgenitor2%integerRank0MetaPropertySet(self%isConstrainedID,1)
             end if
          else
             call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,1)
          end if
       else
          call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,0)
          if (present(nodeProgenitor2)) then
             basicProgenitor2 => nodeProgenitor2%basic()
             ! Need to mark this secondary progenitor as not on the constrained branch, e.g.:
             call basicProgenitor2%integerRank0MetaPropertySet(self%isConstrainedID,0)
          end if
       end if
    else
       if (isConstrained) then
          ! Parent is on the constrained branch, so this progenitor also is - mark it as such using the "integerRank0MetaPropertySet" function.
          call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,1)
       else
          ! Parent is not on the constrained branch, so this progenitor also is not - mark it as such using the "integerRank0MetaPropertySet" function.
          call basicProgenitor1%integerRank0MetaPropertySet(self%isConstrainedID,0)
       end if
       ! If the second progenitor is present, mark it as not on the constrained branch.
       if (present(nodeProgenitor2)) then
          basicProgenitor2 => nodeProgenitor2%basic()
          ! Need to mark this secondary progenitor as not on the constrained branch, e.g.:
          call basicProgenitor2%integerRank0MetaPropertySet(self%isConstrainedID,0)
       end if
    end if
    return
  end subroutine constrainedNodesInserted
