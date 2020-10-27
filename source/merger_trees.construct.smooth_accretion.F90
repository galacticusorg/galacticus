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

  !% Implements a merger tree constructor class which builds merger trees assuming smooth accretion.

  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
  use :: Numerical_Random_Numbers                 , only : randomNumberGeneratorClass

  !# <mergerTreeConstructor name="mergerTreeConstructorSmoothAccretion">
  !#  <description>
  !#   A merger tree constructor class which builds a branchless merger tree with a smooth accretion history using the selected
  !#   \refPhysics{darkMatterHaloMassAccretionHistory} class. The tree has a final mass of {\normalfont \ttfamily massHalo} (in
  !#   units of $M_\odot$) at redshift {\normalfont \ttfamily redshiftBase} and is continued back in time by decreasing the halo
  !#   mass by a factor {\normalfont \ttfamily massHaloDeclineFactor} at each new \gls{node} until a specified {\normalfont
  !#   \ttfamily massHaloResolution} (in units of $M_\odot$) is reached.
  !#  </description>
  !# </mergerTreeConstructor>
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorSmoothAccretion
     !% A class implementing merger tree construction by building trees assuming smooth accretion.
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class           (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     class           (randomNumberGeneratorClass             ), pointer :: randomNumberGenerator_              => null()
     double precision                                                   :: redshiftBase                                 , massHalo          , &
          &                                                                massHaloDeclineFactor                        , massHaloResolution
   contains
     final     ::              smoothAccretionDestructor
     procedure :: construct => smoothAccretionConstruct
  end type mergerTreeConstructorSmoothAccretion

  interface mergerTreeConstructorSmoothAccretion
     !% Constructors for the {\normalfont \ttfamily smoothAccretion} merger tree constructor class.
     module procedure smoothAccretionConstructorParameters
     module procedure smoothAccretionConstructorInternal
  end interface mergerTreeConstructorSmoothAccretion

contains

  function smoothAccretionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily augment} merger tree operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeConstructorSmoothAccretion   )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_
    class           (randomNumberGeneratorClass             ), pointer       :: randomNumberGenerator_
    double precision                                                         :: redshiftBase                       , massHalo          , &
         &                                                                      massHaloDeclineFactor              , massHaloResolution

    !# <inputParameter>
    !#   <name>massHalo</name>
    !#   <defaultValue>1.0d12</defaultValue>
    !#   <description>The final mass of the merger tree base halo to consider when building a smoothly accreting merger tree, in units of $M_\odot$.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massHaloResolution</name>
    !#   <defaultValue>1.0d9</defaultValue>
    !#   <description>The final mass of the merger tree base halo to consider when building a smoothly accreting merger tree, in units of $M_\odot$.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massHaloDeclineFactor</name>
    !#   <defaultValue>0.9d0</defaultValue>
    !#   <description>The factor by which halo mass should decrease in each step back in time building a smoothly accreting merger tree, in units of $M_\odot$.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftBase</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift at which to plant the base node when building the smoothly accreting merger tree.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    !# <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    !# <objectBuilder class="randomNumberGenerator"              name="randomNumberGenerator_"              source="parameters"/>
    self=mergerTreeConstructorSmoothAccretion(redshiftBase,massHalo,massHaloDeclineFactor,massHaloResolution,cosmologyFunctions_,darkMatterHaloMassAccretionHistory_,randomNumberGenerator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"                />
    !# <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !# <objectDestructor name="randomNumberGenerator_"             />
    return
  end function smoothAccretionConstructorParameters

  function smoothAccretionConstructorInternal(redshiftBase,massHalo,massHaloDeclineFactor,massHaloResolution,cosmologyFunctions_,darkMatterHaloMassAccretionHistory_,randomNumberGenerator_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily augment} merger tree operator class.
    implicit none
    type            (mergerTreeConstructorSmoothAccretion   )                        :: self
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    class           (randomNumberGeneratorClass             ), intent(in   ), target :: randomNumberGenerator_
    double precision                                         , intent(in   )         :: redshiftBase                       , massHalo          , &
         &                                                                              massHaloDeclineFactor              , massHaloResolution
    !# <constructorAssign variables="redshiftBase, massHalo, massHaloDeclineFactor, massHaloResolution, *cosmologyFunctions_, *darkMatterHaloMassAccretionHistory_, *randomNumberGenerator_"/>

    return
  end function smoothAccretionConstructorInternal

  subroutine smoothAccretionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily smoothAccretion} merger tree constructor class.
    implicit none
    type(mergerTreeConstructorSmoothAccretion), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"                />
    !# <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    !# <objectDestructor name="self%randomNumberGenerator_"             />
    return
  end subroutine smoothAccretionDestructor

  function smoothAccretionConstruct(self,treeNumber) result(tree)
    !% Build a merger tree with a smooth mass accretion history using the fitting function of \cite{wechsler_concentrations_2002}.
    use            :: Galacticus_Nodes, only : mergerTree  , nodeComponentBasic, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    class           (mergerTreeConstructorSmoothAccretion   ), intent(inout) :: self
    type            (mergerTree                             ), pointer       :: tree
    integer         (c_size_t                               ), intent(in   ) :: treeNumber
    type            (treeNode                               ), pointer       :: nodeCurrent        , nodeNew
    class           (nodeComponentBasic                     ), pointer       :: basicBase          , basicNew
    integer         (kind=kind_int8                         )                :: indexNode
    double precision                                                         :: expansionFactorBase, timeBase, &
         &                                                                      massNode           , timeNode

    ! Build the merger tree.
    if (treeNumber == 1_c_size_t) then
       ! Create the tree.
       allocate(tree)
       ! Give the tree an index.
       tree%index=1
       ! Create the base node.
       indexNode      =  1
       tree%firstTree => tree
       tree%baseNode  => treeNode               (indexNode,tree   )
       basicBase      => tree    %baseNode%basic(autoCreate=.true.)
       ! Assign an arbitrary weight to the tree.
       tree%volumeWeight     =  1.0
       tree%event            => null()
       tree%initializedUntil =  0.0d0
       call tree%properties%initialize()
       ! Restart the random number sequence.
       allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
       !$omp critical(mergerTreeConstructSmoothAccretionDeepCopyReset)
       !# <deepCopyReset variables="self%randomNumberGenerator_"/>
       !# <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
       !$omp end critical(mergerTreeConstructSmoothAccretionDeepCopyReset)
       call tree%randomNumberGenerator_%seedSet(seed=tree%index,offset=.true.)
       ! Assign a mass to the base node.
       call basicBase%massSet(self%massHalo)
       ! Find the cosmic time at which the tree is based.
       expansionFactorBase=self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftBase       )
       timeBase           =self%cosmologyFunctions_%cosmicTime                 (     expansionFactorBase)
       ! Assign a time to the base node.
       call basicBase%            timeSet(timeBase)
       call basicBase%timeLastIsolatedSet(timeBase)
       ! Get a pointer to the current node (i.e. the base node).
       nodeCurrent => tree%baseNode
       ! Initialize current node mass.
       massNode    =  self%massHalo
       ! Step backwards, creating nodes until a sufficiently low mass has been reached.
       do while (massNode > self%massHaloResolution)
          ! Increment node index.
          indexNode=indexNode+1
          ! Create a node.
          nodeNew  => treeNode     (indexNode,tree   )
          basicNew => nodeNew%basic(autoCreate=.true.)
          ! Adjust the mass by the specified factor.
          massNode=massNode*self%massHaloDeclineFactor
          ! Set the mass of the node.
          call basicNew%massSet(massNode)
          ! Find the time corresponding to this expansion factor.
          timeNode=self%darkMatterHaloMassAccretionHistory_%time(tree%baseNode,massNode)
          ! Set the time for the new node.
          call basicNew%            timeSet(timeNode)
          call basicNew%timeLastIsolatedSet(timeNode)
          ! Create parent and child links.
          nodeCurrent%firstChild => nodeNew
          nodeNew    %parent     => nodeCurrent
          ! Move the current node to the new node.
          nodeCurrent            => nodeNew
       end do
    else
       nullify(tree)
    end if
    return
  end function smoothAccretionConstruct
