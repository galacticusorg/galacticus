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
  Implements a node operator class that applies heuristics to look for unphysical behavior in merger trees.
  !!}
  
  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <nodeOperator name="nodeOperatorAssemblyHistoryHeuristics">
  <description>
 
      A node operator class that applies heuristics to look for unphysical behavior in merger trees. Specifically, rapid changes
      in the halo mass from a node to its parent\footnote{This is only considered for primary progenitors (the jump in mass from a
      non-primary progenitor to its parent can, of course, be very large), including cases of subhalo promotions, and for the
      bound mass history of subhalos.} are flagged as unphysical.

      Halos for which the change in mass is less than {\normalfont \ttfamily [sigmaThreshold]} times the uncertainty in the
      difference of the halo mass between child and parent are treated as being physical. For other halos, the unphysicalness is
      judged by the heuristic:
      \begin{equation}
       \left| \log \frac{M_2}{M_1} \right| > \alpha \left| \log \frac{t_2}{t_1} \right|,
      \end{equation}
      where $\alpha=${\normalfont \ttfamily [exponentGrowth]}, $M_1$ is the mass of a primary progenitor halo, $M_2$ is the mass
      of its parent \emph{minus any mass resolved in non-primary progenitors}. Essentially, this says that we expect masses to
      grow (or shrink) no faster that a power-law in cosmic time, with exponent $\alpha$.

      </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorAssemblyHistoryHeuristics
     !!{
     A node operator class that applies heuristics to look for unphysical behavior in merger trees.
     !!}
     private
     class           (nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
     double precision                                   :: exponentGrowth               , sigmaThreshold
     integer                                            :: labelID
   contains
     final     ::                       assemblyHistoryHeuristicsDestructor
     procedure :: nodeTreeInitialize => assemblyHistoryHeuristicsNodeTreeInitialize
     procedure :: nodePromote        => assemblyHistoryHeuristicsNodePromote
     procedure :: galaxiesMerge      => assemblyHistoryHeuristicsGalaxiesMerge
  end type nodeOperatorAssemblyHistoryHeuristics
  
  interface nodeOperatorAssemblyHistoryHeuristics
     !!{
     Constructors for the \refClass{nodeOperatorAssemblyHistoryHeuristics} node operator class.
     !!}
     module procedure assemblyHistoryHeuristicsConstructorParameters
     module procedure assemblyHistoryHeuristicsConstructorInternal
  end interface nodeOperatorAssemblyHistoryHeuristics
  
contains

  function assemblyHistoryHeuristicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorAssemblyHistoryHeuristics} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorAssemblyHistoryHeuristics)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (nbodyHaloMassErrorClass              ), pointer       :: nbodyHaloMassError_
    double precision                                                       :: exponentGrowth     , sigmaThreshold

    !![
    <inputParameter>
      <name>exponentGrowth</name>
      <defaultValue>100.0d0</defaultValue>
      <description>The maximum plausible growth exponent.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigmaThreshold</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The number of $\sigma$ in halo mass uncertainty below which we ignore changes in halo mass.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    !!]
    self=nodeOperatorAssemblyHistoryHeuristics(exponentGrowth,sigmaThreshold,nbodyHaloMassError_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyHaloMassError_"/>
    !!]
    return
  end function assemblyHistoryHeuristicsConstructorParameters

  function assemblyHistoryHeuristicsConstructorInternal(exponentGrowth,sigmaThreshold,nbodyHaloMassError_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorAssemblyHistoryHeuristics} node operator class.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type            (nodeOperatorAssemblyHistoryHeuristics)                        :: self
    class           (nbodyHaloMassErrorClass              ), target, intent(in   ) :: nbodyHaloMassError_
    double precision                                               , intent(in   ) :: exponentGrowth     , sigmaThreshold
    !![
    <constructorAssign variables="exponentGrowth, sigmaThreshold, *nbodyHaloMassError_"/>
    !!]
    
    self%labelID=nodeLabelRegister('nodeAssemblyIsGood')
    return
  end function assemblyHistoryHeuristicsConstructorInternal

  subroutine assemblyHistoryHeuristicsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorAssemblyHistoryHeuristics} node operator class.
    !!}
    implicit none
    type(nodeOperatorAssemblyHistoryHeuristics), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_"/>
    !!]
    return
  end subroutine assemblyHistoryHeuristicsDestructor

  subroutine assemblyHistoryHeuristicsNodeTreeInitialize(self,node)
    !!{
    Initialize the maximum host mass of this node.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentSatellite, nodeEvent, nodeEventSubhaloPromotion
    use :: Nodes_Labels        , only : nodeLabelSet      , nodeLabelUnset
    use :: Numerical_Comparison, only : Values_Agree
    use :: Histories           , only : history
    implicit none
    class           (nodeOperatorAssemblyHistoryHeuristics), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    type            (treeNode                             )               , pointer :: nodeSibling, nodeParent        , &
         &                                                                             nodeWork
    class           (nodeComponentBasic                   )               , pointer :: basic      , basicParent       , &
         &                                                                             basicSibling
    class           (nodeComponentSatellite               )               , pointer :: satellite
    class           (nodeEvent                            )               , pointer :: event
    type            (history                              )                         :: massBound
    double precision                                                                :: massRatio  , timeRatio         , &
         &                                                                             massPrimary, massUncertainty
    logical                                                                         :: isGood     , isPrimary         , &
         &                                                                             foundParent, isSubhaloPromotion
    integer                                                                         :: i

    ! Initialize state.
    isGood             =  .true.
    isSubhaloPromotion =  .false.
    foundParent        =  .false.
    nodeParent         => null()
    ! Look for events attached to this node.
    event => node%event
    do while (associated(event))
       select type (event)
       type is (nodeEventSubhaloPromotion)
          ! Subhalo promotion event.
          nodeParent  => event     %node
          basic       => node      %basic()
          basicParent => nodeParent%basic()
          ! Only consider the child->parent event, not the inverse event.
          if (basicParent%time() > basic%time()) then
             ! Check if we are the primary progenitor.
             if (associated(nodeParent%parent)) then
                basic       => nodeParent       %basic()
                basicParent => nodeParent%parent%basic()
                isPrimary   =  Values_Agree(basic%time(),basicParent%time(),relTol=2.0d-6)
             else
                isPrimary=.false.
             end if
             ! If we are the primary progenitor, set appropriate pointers to our parent and siblings.
             if (isPrimary) then
                nodeSibling        => nodeParent%sibling
                foundParent        =  .true.
                isSubhaloPromotion =  .true.
             else
                nodeParent         => null()
                nodeSibling        => null()
             end if
          end if
          exit          
       end select
       event => event%next
    end do
    ! If no parent has yet been found, check for a direct parent, of which our node is the primary progenitor.
    if (.not.foundParent) then
       if (associated(node%parent).and.associated(node%parent%firstChild,node)) then
          nodeParent  => node%parent
          nodeSibling => node%sibling
          foundParent =  .true.
       end if
    end if
    ! If a parent was found, apply heuristics.    
    if (foundParent) then
       basic       => node      %basic()
       basicParent => nodeParent%basic()
       ! Find the uncertainty in the mass difference.
       massUncertainty=+sqrt(                                                                                   &
            &                +(self%nbodyHaloMassError_%errorFractional(node           )*basic      %mass())**2 &
            &                +(self%nbodyHaloMassError_%errorFractional(     nodeParent)*basicParent%mass())**2 &
            &                - self%nbodyHaloMassError_%errorFractional(node           )*basic      %mass()     &
            &                * self%nbodyHaloMassError_%errorFractional(     nodeParent)*basicParent%mass()     &
            &                * self%nbodyHaloMassError_%correlation    (node,nodeParent)                        &
            &               )
       ! If the difference in masses is below the mass uncertainty threshold, accept this halo as valid.
       if (abs(basicParent%mass()-basic%mass()) > self%sigmaThreshold*massUncertainty) then
          ! Otherwise, apply heuristics.
          if (basicParent%mass() > basic%mass()) then
             ! The parent mass is greater than the child mass. Find the mass of the parent excluding resolved mergers.
             massPrimary =  basicParent%mass()
             nodeWork    => nodeSibling
             do while (associated(nodeWork))
                basicSibling =>  nodeWork    %basic      ()
                massPrimary  =  +             massPrimary   &
                     &          -basicSibling%mass       ()
                nodeWork     =>  nodeWork    %sibling
             end do
          else
             ! The parent mass is less than the mass of the child. In this case, we compare the two masses directly.
             massPrimary     =   basicParent %mass       ()
          end if
          ! Check the rate of mass growth.
          massRatio=+      massPrimary   &
               &    /basic%mass       ()
          if (massRatio > 0.0d0) then
             timeRatio=+basicParent%time() &
                  &    /basic      %time()
             if (abs(log(massRatio)) > self%exponentGrowth*abs(log(timeRatio))) then
                isGood=.false.
                call heuristicWarn(node,nodeParent,nodeSibling,'rapid mass change',isSubhaloPromotion)
             end if
          else
             isGood=.false.
             call heuristicWarn(node,nodeParent,nodeSibling,'significant negative unresolved mass',isSubhaloPromotion)
          end if
       end if
    end if
    ! Check any subhalo bound mass history.
    satellite => node%satellite()
    if (satellite%boundMassHistoryIsGettable()) then       
       massBound=satellite%boundMassHistory()
       if (massBound%exists()) then
          do i=1,size(massBound%time)-1
             massRatio=massBound%data(i+1,1)/massBound%data(i,1)
             timeRatio=massBound%time(i+1  )/massBound%time(i  )
             if (abs(log(massRatio)) > self%exponentGrowth*abs(log(timeRatio))) then
                isGood      = .false.
                nodeParent  => null()
                nodeSibling => null()
                call heuristicWarn(node,nodeParent,nodeSibling,'rapid mass change in subhalo bound mass',isSubhaloPromotion,massBound,i)
             end if
          end do
       end if
    end if
    ! Set the label appropriately.
    if (isGood) then
       call nodeLabelSet  (self%labelID,node)
    else
       call nodeLabelUnset(self%labelID,node)
    end if
    return
  end subroutine assemblyHistoryHeuristicsNodeTreeInitialize

  subroutine assemblyHistoryHeuristicsNodePromote(self,node)
    !!{
    Propagate bad status into parent node on promotion.
    !!}
    use :: Nodes_Labels, only : nodeLabelIsPresent, nodeLabelUnset    
    implicit none
    class(nodeOperatorAssemblyHistoryHeuristics), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    if (.not.nodeLabelIsPresent(self%labelID,node%parent)) call nodeLabelUnset(self%labelID,node)
    return
  end subroutine assemblyHistoryHeuristicsNodePromote

  subroutine assemblyHistoryHeuristicsGalaxiesMerge(self,node)
    !!{
    Propagate bad status into merge target when two galaxies merge.
    !!}
    use :: Nodes_Labels, only : nodeLabelIsPresent, nodeLabelUnset
    implicit none
    class(nodeOperatorAssemblyHistoryHeuristics), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    type (treeNode                             ), pointer       :: nodeHost

    nodeHost => node%mergesWith()
    if (.not.nodeLabelIsPresent(self%labelID,node)) call nodeLabelUnset(self%labelID,nodeHost)
    return
  end subroutine assemblyHistoryHeuristicsGalaxiesMerge

  subroutine heuristicWarn(node,nodeParent,nodeSibling,description,isSubhaloPromotion,massBound,indexBound)
    !!{
    Warn about violations of assembly history heuristics.
    !!}
    use :: Display           , only : displayMessage    , displayIndent     , displayUnindent, displayMagenta, &
         &                            displayReset      , verbosityLevelWarn
    use :: Galacticus_Nodes  , only : nodeComponentBasic, treeNode
    use :: Histories         , only : history
    use :: ISO_Varying_String, only : varying_string    , assignment(=)     , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (treeNode          ), intent(inout), target   :: node
    type     (treeNode          ), intent(inout), pointer  :: nodeParent       , nodeSibling
    character(len=*             ), intent(in   )           :: description
    logical                      , intent(in   )           :: isSubhaloPromotion
    type     (history           ), intent(in   ), optional :: massBound
    integer                      , intent(in   ), optional :: indexBound
    type     (treeNode          )               , pointer  :: nodeWork
    class    (nodeComponentBasic)               , pointer  :: basic
    type     (varying_string    )                          :: message
    character(len=30            )                          :: label

    message =  description
    if (isSubhaloPromotion) message=message//" (subhalo promotion)"
    basic   => node%basic()
    write (label,'(a2,e12.6,a2,e12.6)') ", ",basic%time(),", ",basic%mass()
    message   =message//char(10)//"   node: [index, time, mass] = "//node      %index()//trim(label)
    if (associated(nodeParent)) then
       basic => nodeParent%basic()
       write (label,'(a2,e12.6,a2,e12.6)') ", ",basic%time(),", ",basic%mass()
       message=message//char(10)//" parent: [index, time, mass] = "//nodeParent%index()//trim(label)
    end if
    nodeWork => nodeSibling
    do while (associated(nodeWork))
       basic => nodeWork%basic()
       write (label,'(a2,e12.6,a2,e12.6)') ", ",basic%time(),", ",basic%mass()
       message=message//char(10)//"sibling: [index, time, mass] = "//nodeWork  %index()//trim(label)
       nodeWork => nodeWork%sibling
    end do    
    if (present(massBound).and.present(indexBound)) then
       write (label,'(a2,e12.6,a2,e12.6)') ", ",massBound%time(indexBound  ),", ",massBound%data(indexBound  ,1)
       message=message//char(10)//"subhalo: [index, time, mass] = "//indexBound//trim(label)
       write (label,'(a2,e12.6,a2,e12.6)') ", ",massBound%time(indexBound+1),", ",massBound%data(indexBound+1,1)
       message=message//char(10)//"subhalo: [index, time, mass] = "//indexBound+1//trim(label)
    end if
    call displayIndent  (displayMagenta()//"WARN: "//displayReset()//"node fails assembly history heuristics",verbosityLevelWarn)
    call displayMessage (message                                                                             ,verbosityLevelWarn)
    call displayUnindent(""                                                                                  ,verbosityLevelWarn)
    return
  end subroutine heuristicWarn
