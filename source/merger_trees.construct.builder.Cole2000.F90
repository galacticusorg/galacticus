!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of a merger tree builder using the algorithm of \cite{cole_hierarchical_2000}.
  use Merger_Trees_Build_Mass_Resolution

  !# <mergerTreeBuilder name="mergerTreeBuilderCole2000">
  !#  <description>Merger trees are built using the algorithm of \cite{cole_hierarchical_2000}.</description>
  !# </mergerTreeBuilder>
  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderCole2000
     !% A merger tree builder class using the algorithm of \cite{cole_hierarchical_2000}.
     private
     class           (mergerTreeMassResolutionClass), pointer :: mergerTreeMassResolution_
     ! Variables controlling merger tree accuracy.
     double precision                                         :: accretionLimit      , timeEarliest    , &
          &                                                      mergeProbability
     ! Option controlling random number sequences.
     logical                                                  :: randomSeedsFixed
     
     ! Random number sequence variables
     type            (fgsl_rng                              ) :: clonedPseudoSequence, pseudoSequence
     logical                                                  :: reset               , ompThreadOffset , &
          &                                                      resetSnapshot
     integer                                                  :: incrementSeed    
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeBuilderCole2000</object>
     !@   <objectMethod>
     !@     <method>shouldAbort</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textless type(mergerTree)\textgreater\ tree\argin</arguments>
     !@     <description>Return true if construction of the merger tree should be aborted.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>shouldFollowBranch</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textless type(mergerTree)\textgreater\ tree\argin, \textless type(treeNode)\textgreater\ node\argin</arguments>
     !@     <description>Return true if the branch should be followed.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: build              => cole2000Build
     procedure :: shouldAbort        => cole2000ShouldAbort
     procedure :: shouldFollowBranch => cole2000ShouldFollowBranch
     procedure :: timeEarliestSet    => cole2000TimeEarliestSet
     procedure :: stateStore         => cole2000StateStore 
     procedure :: stateRestore       => cole2000StateRestore
     procedure :: stateSnapshot      => cole2000StateSnapshot
  end type mergerTreeBuilderCole2000

  interface mergerTreeBuilderCole2000
     !% Constructors for the {\normalfont \ttfamily cole2000} merger tree builder class.
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerTreeBuilderCole2000
  
contains

  function cole2000ConstructorParameters(parameters)
    !% Constructor for the \cite{cole_hierarchical_2000} merger tree building class which reads parameters from a provided parameter list.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    implicit none
    type            (mergerTreeBuilderCole2000)                :: cole2000ConstructorParameters    
    type            (inputParameters          ), intent(in   ) :: parameters
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    double precision                                           :: redshiftMaximum
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>mergeProbability</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%mergeProbability</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The largest probability of branching allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>accretionLimit</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%accretionLimit</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The largest fractional mass change due to subresolution accretion allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d5</defaultValue>
    !#   <description>The highest redshift to which merger trees will be built in the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomSeedsFixed</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%randomSeedsFixed</variable>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description></description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="mergerTreeMassResolution" name="cole2000ConstructorParameters%mergerTreeMassResolution_" source="parameters"/>
    ! Convert maximum redshift to earliest time.
    cosmologyFunctions_ => cosmologyFunctions()
    cole2000ConstructorParameters%timeEarliest=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    ! Initialize state.
    cole2000ConstructorParameters%reset           =.true.
    cole2000ConstructorParameters%ompThreadOffset =.true.
    cole2000ConstructorParameters%incrementSeed   =0
    if (FGSL_Well_Defined(cole2000ConstructorParameters%      pseudoSequence))                &
         & call FGSL_Obj_C_Ptr(cole2000ConstructorParameters%      pseudoSequence,C_Null_Ptr)
    if (FGSL_Well_Defined(cole2000ConstructorParameters%clonedPseudoSequence))                &
         & call FGSL_Obj_C_Ptr(cole2000ConstructorParameters%clonedPseudoSequence,C_Null_Ptr)
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(mergeProbability,accretionLimit,timeEarliest,randomSeedsFixed,mergerTreeMassResolution_)
    !% Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    implicit none
    type            (mergerTreeBuilderCole2000    )                        :: cole2000ConstructorInternal
    double precision                               , intent(in   )         :: mergeProbability           , accretionLimit, &
         &                                                                    timeEarliest
    logical                                        , intent(in   )         :: randomSeedsFixed
    class           (mergerTreeMassResolutionClass), intent(in   ), target :: mergerTreeMassResolution_
    
    ! Store options.
    cole2000ConstructorInternal%mergeProbability          =  mergeProbability
    cole2000ConstructorInternal%accretionLimit            =  accretionLimit
    cole2000ConstructorInternal%timeEarliest              =  timeEarliest
    cole2000ConstructorInternal%randomSeedsFixed          =  randomSeedsFixed
    cole2000ConstructorInternal%mergerTreeMassResolution_ => mergerTreeMassResolution_
    ! Initialize state.
    cole2000ConstructorInternal%reset           =.true.
    cole2000ConstructorInternal%ompThreadOffset =.true.
    cole2000ConstructorInternal%incrementSeed   =0
    if (FGSL_Well_Defined(cole2000ConstructorInternal%      pseudoSequence))                &
         & call FGSL_Obj_C_Ptr(cole2000ConstructorInternal%      pseudoSequence,C_Null_Ptr)
    if (FGSL_Well_Defined(cole2000ConstructorInternal%clonedPseudoSequence))                &
         & call FGSL_Obj_C_Ptr(cole2000ConstructorInternal%clonedPseudoSequence,C_Null_Ptr)
    return
  end function cole2000ConstructorInternal

  subroutine cole2000Build(self,tree)
    !% Build a merger tree.
    use Galacticus_Nodes
    use Galacticus_Error
    use Critical_Overdensity
    use Merger_Tree_Branching
    use Pseudo_Random
    use Kind_Numbers
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout)         :: self
    type            (mergerTree               ), intent(inout), target :: tree
    type            (treeNode                 ), pointer               :: newNode1         , newNode2     , thisNode
    class           (nodeComponentBasic       ), pointer               :: newBasic1        , newBasic2    , thisBasic           , &
         &                                                                parentBasic
    integer         (kind=kind_int8           )                        :: nodeIndex
    double precision                                                   :: accretionFraction, baseNodeTime , branchingProbability, &
         &                                                                collapseTime     , deltaCritical, deltaCritical1      , &
         &                                                                deltaCritical2   , deltaW       , nodeMass1           , &
         &                                                                nodeMass2        , time         , uniformRandom       , &
         &                                                                massResolution
    logical                                                            :: doBranch

    nodeIndex =  1                   ! Initialize the node index counter to unity.
    thisNode  => tree    %baseNode   ! Point to the base node.
    thisBasic => thisNode%basic   () ! Get the basic component of the node.
    ! Restart the random number sequence.
    if (self%randomSeedsFixed) then
       if (.not.self%reset) call Pseudo_Random_Free(self%pseudoSequence)
       self%reset          =.true.
       self%incrementSeed  =int(tree%index)
       self%ompThreadOffset=.false.
    end if
    ! Get the mass resolution for this tree.
    massResolution=self%mergerTreeMassResolution_%resolution(tree)
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this module).
    baseNodeTime =thisBasic%time()
    deltaCritical=Critical_Overdensity_for_Collapse(time=thisBasic%time(),mass=thisBasic%mass())
    call thisBasic%timeSet(deltaCritical)
    ! Begin tree build loop.
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasic => thisNode%basic()
       ! If halo is above the resolution limit, then evolve it.
       if     (                                                                                                   &
            &                                                              thisBasic%mass()  > massResolution     &
            &  .and.                                                                                              &
            &   Time_of_Collapse(criticalOverdensity=thisBasic%time(),mass=thisBasic%mass()) > self%timeEarliest  &
            &  .and.                                                                                              &
            &   self%shouldFollowBranch(tree,thisNode)                                                            &
            & ) then
          ! Find branching probability.
          branchingProbability=Tree_Branching_Probability (thisBasic%mass(),thisBasic%time(),massResolution)
          ! Find accretion rate.
          accretionFraction   =Tree_Subresolution_Fraction(thisBasic%mass(),thisBasic%time(),massResolution)
          ! A negative accretion fraction indicates that the node is so close to the resolution limit that
          ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
          ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
          if (accretionFraction < 0.0d0) then
             thisNode => thisNode%walkTreeUnderConstruction()
             cycle
          end if
          ! Finding maximum allowed step in w.
          deltaW=min(self%accretionLimit/accretionFraction,Tree_Maximum_Step(thisBasic%mass(),thisBasic%time(),massResolution))
          deltaW=Tree_Maximum_Step(thisBasic%mass(),thisBasic%time(),massResolution)
          if (accretionFraction    > 0.0d0) deltaW=min(deltaW,self%accretionLimit  /accretionFraction   )
          if (branchingProbability > 0.0d0) deltaW=min(deltaW,self%mergeProbability/branchingProbability)
          ! Scale values to the determined timestep.
          branchingProbability=branchingProbability*deltaW
          accretionFraction   =accretionFraction   *deltaW
          ! Accretion fraction must be less than unity. Reduce timestep (and branching probability and accretion fraction) by
          ! factors of two until this condition is satisfied.
          do while (accretionFraction >= 1.0d0)
             deltaW              =deltaW              *0.5d0
             branchingProbability=branchingProbability*0.5d0
             accretionFraction   =accretionFraction   *0.5d0
          end do
          ! Decide if a branching occurs.
          if (branchingProbability > 0.0d0) then
             uniformRandom=Pseudo_Random_Get(self%pseudoSequence,reset=self%reset,ompThreadOffset=self%ompThreadOffset,incrementSeed=self%incrementSeed)
             doBranch=(uniformRandom <= branchingProbability)
          else
             doBranch=.false.
          end if
          ! Determine the critical overdensity for collapse for the new halo(s).
          deltaCritical=thisBasic%time()+deltaW
          ! Create new nodes.
          select case (doBranch)
          case (.true.)
             ! Branching occurs - create two progenitors.
             nodeIndex=nodeIndex+1
             newNode1  => treeNode(nodeIndex,tree)
             newBasic1 => newNode1%basic(autoCreate=.true.)
             ! Compute mass of the new nodes. First convert the realized probability back to a rate.
             branchingProbability=uniformRandom/deltaW
             nodeMass1=Tree_Branch_Mass(thisBasic%mass(),thisBasic%time(),massResolution,branchingProbability)
             nodeMass2=thisBasic%mass()-nodeMass1
             nodeMass1=nodeMass1*(1.0d0-accretionFraction)
             nodeMass2=nodeMass2*(1.0d0-accretionFraction)            
             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=thisBasic%mass())
             ! Set properties of first new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call newBasic1%massSet(nodeMass1     )
             call newBasic1%timeSet(deltaCritical1)
             ! Create second progenitor.
             nodeIndex=nodeIndex+1
             newNode2  => treeNode(nodeIndex,tree)
             newBasic2 => newNode2%basic(autoCreate=.true.)
             ! Set properties of second new node.
             deltaCritical2=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass2)
             call newBasic2%massSet(nodeMass2     )
             call newBasic2%timeSet(deltaCritical2)
             ! Create links from old to new nodes and vice-versa. (Ensure that child node is the more massive progenitor.)
             if (nodeMass2 > nodeMass1) then
                thisNode%firstChild => newNode2
                newNode2%sibling    => newNode1
             else
                thisNode%firstChild => newNode1
                newNode1%sibling    => newNode2
             end if
             newNode1%parent        => thisNode
             newNode2%parent        => thisNode
          case (.false.)
             ! No branching occurs - create one progenitor.
             nodeIndex=nodeIndex+1
             newNode1  => treeNode(nodeIndex,tree)
             newBasic1 => newNode1%basic(autoCreate=.true.)
             ! Compute new mass accounting for sub-resolution accretion.
             nodeMass1=thisBasic%mass()*(1.0d0-accretionFraction)
             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=thisBasic%mass())
             ! Set properties of the new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call newBasic1%massSet(nodeMass1     )
             call newBasic1%timeSet(deltaCritical1)
             ! Create links from old to new node and vice-versa.
             thisNode%firstChild => newNode1
             newNode1%parent     => thisNode
          end select
       end if
       ! Check if tree should be aborted.
       if (self%shouldAbort(tree)) then
          thisNode => null()
       else
          ! Walk to the next node.
          thisNode => thisNode%walkTreeUnderConstruction()
       end if
    end do
    ! Walk the tree and convert w to time.
    thisNode => tree%baseNode
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasic => thisNode%basic()
       ! Compute the collapse time.
       collapseTime=Time_of_Collapse(criticalOverdensity=thisBasic%time(),mass=thisBasic%mass())
       call thisBasic%timeSet(collapseTime)
       thisNode => thisNode%walkTree()
    end do
    thisBasic => tree%baseNode%basic()
    call thisBasic%timeSet(baseNodeTime)
    ! Check for well-ordering in time.
    thisNode => tree%baseNode
    do while (associated(thisNode))       
       if (associated(thisnode%parent)) then
          thisBasic   => thisNode       %basic()
          parentBasic => thisNode%parent%basic()
          if (parentBasic%time() <= thisBasic%time()) call Galacticus_Error_Report('cole2000Build','branch is not well-ordered in time')
       end if
       thisNode => thisNode%walkTree()
    end do       
    return
  end subroutine cole2000Build

  logical function cole2000ShouldAbort(self,tree)
    !% Return {\normalfont \ttfamily true} if tree construction should be aborted. In the {\normalfont \ttfamily cole2000} tree
    !% builder we never abort.
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout) :: self
    type (mergerTree               ), intent(in   ) :: tree

    cole2000ShouldAbort=.false.
    return
  end function cole2000ShouldAbort
  
  logical function cole2000ShouldFollowBranch(self,tree,node)
    !% Return {\normalfont \ttfamily true} if tree construction should continue to follow the current branch. In the {\normalfont
    !% \ttfamily cole2000} tree builder we always continue.
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout)          :: self
    type (mergerTree               ), intent(in   )          :: tree
    type (treeNode                 ), intent(inout), pointer :: node

    cole2000ShouldFollowBranch=.true.
    return
  end function cole2000ShouldFollowBranch
  
  subroutine cole2000TimeEarliestSet(self,timeEarliest)
    !% Set the earliest time for the tree builder.
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout) :: self
    double precision                           , intent(in   ) :: timeEarliest

    self%timeEarliest=timeEarliest
    return
  end subroutine cole2000TimeEarliestSet
  
  subroutine cole2000StateSnapshot(self)
    !% Store a snapshot of the random number generator internal state.
    use Pseudo_Random
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout) :: self
    
    if (.not.self%reset) then
       if (FGSL_Well_Defined(self%clonedPseudoSequence)) call Pseudo_Random_Free(self%clonedPseudoSequence)
       self%clonedPseudoSequence=FGSL_Rng_Clone(self%pseudoSequence)
    end if
    self%resetSnapshot=self%reset
    return
  end subroutine cole2000StateSnapshot

  subroutine cole2000StateStore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (mergerTreeBuilderCole2000), intent(inout) :: self
    integer                           , intent(in   ) :: stateFile
    type   (fgsl_file                ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetSnapshot
    if (.not.self%resetSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequence,fgslStateFile)
    return
  end subroutine cole2000StateStore

  subroutine cole2000StateRestore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (mergerTreeBuilderCole2000), intent(inout) :: self
    integer                           , intent(in   ) :: stateFile
    type   (fgsl_file                ), intent(in   ) :: fgslStateFile

    read (stateFile) self%reset
    if (.not.self%reset) call Pseudo_Random_Retrieve(self%pseudoSequence,fgslStateFile)
    return
  end subroutine cole2000StateRestore
