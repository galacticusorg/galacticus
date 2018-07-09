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
  
  !% Implements a merger tree constructor class which builds merger trees after drawing masses at random from a mass distribution.
  
  use Cosmology_Functions
  use Cosmology_Parameters
  use Halo_Mass_Functions
  use Merger_Trees_Build_Masses
  use Merger_Trees_Builders

  !# <mergerTreeConstructor name="mergerTreeConstructorBuild">
  !#  <description>Merger tree constructor class which builds merger trees.</description>
  !# </mergerTreeConstructor>
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorBuild
     !% A class implementing merger tree construction by building trees.
     private
     class           (cosmologyParametersClass  ), pointer                   :: cosmologyParameters_
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_
     class           (mergerTreeBuildMassesClass), pointer                   :: mergerTreeBuildMasses_
     class           (mergerTreeBuilderClass    ), pointer                   :: mergerTreeBuilder_
     class           (haloMassFunctionClass     ), pointer                   :: haloMassFunction_
     ! Variables giving the mass range and sampling frequency for mass function sampling.
     double precision                                                        :: timeBase
     integer                                                                 :: treeBeginAt
     ! Direction in which to process trees. 
     logical                                                                 :: processDescending
     ! Array of halo masses to use. 
     integer         (c_size_t                  )                            :: treeCount              , treeNumberOffset
     double precision                            , allocatable, dimension(:) :: treeMass               , treeWeight      , &
          &                                                                     treeMassMinimum        , treeMassMaximum
     integer         (c_size_t                  ), allocatable, dimension(:) :: treeMassCount          , rankMass
     logical                                                                 :: computeTreeWeights
   contains
     final     ::              buildDestructor
     procedure :: construct => buildConstruct
  end type mergerTreeConstructorBuild

  interface mergerTreeConstructorBuild
     !% Constructors for the {\normalfont \ttfamily build} merger tree constructor class.
     module procedure buildConstructorParameters
     module procedure buildConstructorInternal
  end interface mergerTreeConstructorBuild

contains
  
  function buildConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily augment} merger tree operator class which takes a parameter set as input.
    use Input_Parameters
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (mergerTreeConstructorBuild)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class           (mergerTreeBuilderClass    ), pointer       :: mergerTreeBuilder_
    class           (haloMassFunctionClass     ), pointer       :: haloMassFunction_
    class           (mergerTreeBuildMassesClass), pointer       :: mergerTreeBuildMasses_
    double precision                                            :: redshiftBase
    integer                                                     :: treeBeginAt
    logical                                                     :: processDescending

    !# <inputParameter>
    !#   <name>redshiftBase</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift at which to plant the base node when building merger trees.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>treeBeginAt</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0</defaultValue>
    !#   <description>The index (in order of increasing base halo mass) of the tree at which to begin when building merger trees. A value of ``0'' means to begin with tree number 1 (if processing trees in ascending order), or equal to the number of trees (otherwise).</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>processDescending</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true, causes merger trees to be processed in order of decreasing mass.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="mergerTreeBuilder"     name="mergerTreeBuilder_"     source="parameters"/>
    !# <objectBuilder class="haloMassFunction"      name="haloMassFunction_"      source="parameters"/>
    !# <objectBuilder class="mergerTreeBuildMasses" name="mergerTreeBuildMasses_" source="parameters"/>
    self=mergerTreeConstructorBuild(                                                                                                        &
         &                          cosmologyFunctions_    %cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftBase     )), &
         &                                                                                                             treeBeginAt        , &
         &                                                                                                             processDescending  , &
         &                          cosmologyParameters_                                                                                  , &
         &                          cosmologyFunctions_                                                                                   , &
         &                          mergerTreeBuildMasses_                                                                                , &
         &                          mergerTreeBuilder_                                                                                    , &
         &                          haloMassFunction_                                                                                       &
         &                         )
    !# <inputParametersValidate source="parameters"/>
    return
  end function buildConstructorParameters

  function buildConstructorInternal(timeBase,treeBeginAt,processDescending,cosmologyParameters_,cosmologyFunctions_,mergerTreeBuildMasses_,mergerTreeBuilder_,haloMassFunction_) result(self)
    !% Initializes the merger tree building module.
    use Memory_Management
    use Sort
    use Galacticus_Error
    implicit none
    type            (mergerTreeConstructorBuild)                        :: self
    double precision                            , intent(in   )         :: timeBase
    integer                                     , intent(in   )         :: treeBeginAt
    logical                                     , intent(in   )         :: processDescending
    class           (cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    class           (mergerTreeBuilderClass    ), intent(in   ), target :: mergerTreeBuilder_
    class           (haloMassFunctionClass     ), intent(in   ), target :: haloMassFunction_
    class           (mergerTreeBuildMassesClass), intent(in   ), target :: mergerTreeBuildMasses_
    integer         (c_size_t                  )                        :: iTreeFirst            , iTreeLast
    !# <constructorAssign variables="timeBase, treeBeginAt, processDescending, *cosmologyParameters_, *cosmologyFunctions_, *mergerTreeBuildMasses_, *mergerTreeBuilder_, *haloMassFunction_"/>

    ! Set offset for tree numbers.
    if (self%treeBeginAt == 0) then
       self%treeNumberOffset=0_c_size_t
    else
       self%treeNumberOffset=self%treeBeginAt-1_c_size_t
    end if
    ! Generate set of merger tree masses
    call self%mergerTreeBuildMasses_%construct(self%timeBase,self%treeMass,self%treeMassMinimum,self%treeMassMaximum,self%treeWeight)
    self%treeCount=size(self%treeMass,kind=c_size_t)
    ! Sort halos by mass.
    allocate(self%rankMass(self%treeCount))
    self%rankMass=Sort_Index_Do(self%treeMass)
    ! Compute the weight (number of trees per unit volume) for each tree if weights were not supplied.
    self%computeTreeWeights=.not.allocated(self%treeWeight)
    if (.not.self%computeTreeWeights) then
       if     (                                 &
            &   allocated(self%treeMassMinimum) &
            &  .or.                             &
            &   allocated(self%treeMassMaximum) &
            & ) call Galacticus_Error_Report('mass interval should not be set if tree weights are provided'//{introspection:location})
    else
       call allocateArray(self%treeWeight   ,[self%treeCount])
       call allocateArray(self%treeMassCount,[self%treeCount])
       iTreeFirst=0
       do while (iTreeFirst < self%treeCount)
          iTreeFirst=iTreeFirst+1
          ! Find the last tree with the same mass.
          iTreeLast=iTreeFirst
          if (iTreeLast < self%treeCount) then
             do while (self%treeMass(self%rankMass(iTreeLast+1)) == self%treeMass(self%rankMass(iTreeFirst)))
                iTreeLast=iTreeLast+1
                if (iTreeLast == self%treeCount) exit
             end do
          end if
          ! Store the number of trees at this mass.
          self%treeMassCount(iTreeFirst:iTreeLast)=iTreeLast-iTreeFirst+1
          ! Update to the last tree processed.
          iTreeFirst=iTreeLast
       end do
    end if
    return
  end function buildConstructorInternal

  subroutine buildDestructor(self)
    !% Destructor for the {\normalfont \ttfamily build} merger tree constructor class.
    implicit none
    type(mergerTreeConstructorBuild), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%mergerTreeBuilder_"    />
    !# <objectDestructor name="self%haloMassFunction_"     />
    !# <objectDestructor name="self%mergerTreeBuildMasses_"/>
    return
  end subroutine buildDestructor
  
  function buildConstruct(self,treeNumber) result(tree)
    !% Build a merger tree.
    use    Galacticus_State
    use    Kind_Numbers
    use    String_Handling
    use    Merger_Tree_State_Store
    use    Pseudo_Random
    !$ use OMP_Lib
    implicit none
    type            (mergerTree                ), pointer       :: tree
    class           (mergerTreeConstructorBuild), intent(inout) :: self
    integer         (c_size_t                  ), intent(in   ) :: treeNumber
    class           (nodeComponentBasic        ), pointer       :: basicBase
    integer         (kind_int8                 ), parameter     :: baseNodeIndex=1
    integer         (kind_int8                 )                :: treeIndex
    type            (varying_string            )                :: message
    double precision                                            :: uniformRandom

    ! Prepare to store/restore internal state.
    treeStateStoreSequence=-1_c_size_t
    ! Retrieve stored internal state if possible.
    call Galacticus_State_Retrieve()
    if (treeStateStoreSequence > 0_c_size_t) then
       if (self%processDescending) then
          self%treeNumberOffset=self%treeCount-treeStateStoreSequence
       else
          self%treeNumberOffset=              +treeStateStoreSequence-1_c_size_t
       end if
    end if    
    ! Determine the index of the tree to process.
    if (self%processDescending) then
       ! Processing trees in descending order, so begin from the final index and work back.
       treeIndex=self%treeCount+1-(treeNumber+self%treeNumberOffset)
    else
       ! Processing trees in ascending order, to just use treeNumber as the index of the tree to process.
       treeIndex=                +(treeNumber+self%treeNumberOffset)
    end if
    if     (                             &
         &   treeIndex >  0_kind_int8    &
         &  .and.                        &
         &   treeIndex <= self%treeCount &
         & ) then
       ! Allocate the tree.
       allocate(tree)
       ! Give the tree an index.
       tree%index=treeIndex
       ! Restart the random number sequence.
       tree%randomNumberGenerator=pseudoRandom()
       uniformRandom=tree%randomNumberGenerator%uniformSample(ompThreadOffset=.false.,mpiRankOFfset=.false.,incrementSeed=int(tree%index))
       ! Store the internal state.
       if (treeStateStoreSequence == -1_c_size_t) treeStateStoreSequence=treeNumber
       message=var_str('Storing state for tree #')//treeNumber
       call Galacticus_State_Store(message)
        ! Initialize.
       tree%event            => null()
       tree%initializedUntil =  0.0d0
       call tree%properties%initialize()
       ! Create the base node.
       tree%baseNode => treeNode(baseNodeIndex,tree)
       ! Get the basic component of the base node.
       basicBase     => tree%baseNode%basic(autoCreate=.true.)
       ! Assign a mass to it.
       call basicBase%massSet(self%treeMass(self%rankMass(treeIndex)))
       ! Assign a time.
       call basicBase%timeSet(self%timeBase           )
       ! Assign a weight to the tree, computing it if necessary.
       if (self%computeTreeWeights) then
          ! The weight is computed by finding the total mass in halos per unit volume within the mass range represented by this
          ! tree, and dividing that by the tree mass. This ensures that we account for all mass that should be in halos. (The
          ! alternative, just computing the number density of halos within the mass range by integrating the mass function, gives
          ! a slightly different answer as the contribution from each mass in the range is then not mass-weighted.) We finally
          ! divide through by the number of trees with this mass.
          self%treeWeight(self%rankMass(treeIndex))=+self%haloMassFunction_   %massFraction   (                                                &
               &                                                                               self%timeBase                                 , &
               &                                                                               self%treeMassMinimum(self%rankMass(treeIndex)), &
               &                                                                               self%treeMassMaximum(self%rankMass(treeIndex)), &
               &                                                                               tree%baseNode                                   &
               &                                                                              )                                                &
               &                                    *self%cosmologyParameters_%densityCritical(                                                &
               &                                                                              )                                                &
               &                                    *self%cosmologyParameters_%OmegaMatter    (                                                &
               &                                                                              )                                                &
               &                                    /                                          self%treeMass       (self%rankMass(treeIndex))  &
               &                                    /dble                                     (                                                &
               &                                                                               self%treeMassCount  (              treeIndex)   &
               &                                                                              )
       end if
       tree%volumeWeight=self%treeWeight(self%rankMass(treeIndex))
       ! Build the tree.
       call self%mergerTreeBuilder_%build(tree)
    else
       nullify(tree)
    end if
    return
  end function buildConstruct
