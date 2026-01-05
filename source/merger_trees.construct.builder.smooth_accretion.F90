!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a merger tree builder class which builds merger trees assuming smooth accretion.
  !!}

  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
  use :: Merger_Trees_Build_Mass_Resolution       , only : mergerTreeMassResolutionClass

  !![
  <mergerTreeBuilder name="mergerTreeBuilderSmoothAccretion">
   <description>
    A merger tree builder class which builds a branchless merger tree with a smooth accretion history using the selected
    \refPhysics{darkMatterHaloMassAccretionHistory} class. The tree has a final mass of {\normalfont \ttfamily massHalo} (in
    units of $M_\odot$) at redshift {\normalfont \ttfamily redshiftBase} and is continued back in time by decreasing the halo
    mass by a factor {\normalfont \ttfamily massHaloDeclineFactor} at each new \gls{node} until a specified {\normalfont
    \ttfamily massHaloResolution} (in units of $M_\odot$) is reached.
   </description>
  </mergerTreeBuilder>
  !!]
  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderSmoothAccretion
     !!{
     A class implementing merger tree building assuming smooth accretion.
     !!}
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class           (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     class           (mergerTreeMassResolutionClass          ), pointer :: mergerTreeMassResolution_           => null()
     double precision                                                   :: massHaloDeclineFactor                        , timeEarliest, &
          &                                                                redshiftEarliest
   contains
     final     ::          smoothAccretionDestructor
     procedure :: build => smoothAccretionBuild
  end type mergerTreeBuilderSmoothAccretion

  interface mergerTreeBuilderSmoothAccretion
     !!{
     Constructors for the \refClass{mergerTreeBuilderSmoothAccretion} merger tree constructor class.
     !!}
     module procedure smoothAccretionConstructorParameters
     module procedure smoothAccretionConstructorInternal
  end interface mergerTreeBuilderSmoothAccretion

contains

  function smoothAccretionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuilderSmoothAccretion} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuilderSmoothAccretion       )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_
    class           (mergerTreeMassResolutionClass          ), pointer       :: mergerTreeMassResolution_
    double precision                                                         :: massHaloDeclineFactor              , redshiftEarliest, &
         &                                                                      timeEarliest

    !![
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"           name="mergerTreeMassResolution_"           source="parameters"/>
    <inputParameter>
      <name>massHaloDeclineFactor</name>
      <defaultValue>0.9d0</defaultValue>
      <description>The factor by which halo mass should decrease in each step back in time building a smoothly accreting merger tree.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('redshiftEarliest')) then
       !![
       <inputParameter>
         <name>redshiftEarliest</name>
         <description>The earliest redshift to which to build a smoothly accreting merger tree.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       timeEarliest=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftEarliest))
    else
       timeEarliest=0.0d0
    end if
    self=mergerTreeBuilderSmoothAccretion(massHaloDeclineFactor,timeEarliest,cosmologyFunctions_,darkMatterHaloMassAccretionHistory_,mergerTreeMassResolution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="mergerTreeMassResolution_"          />
    !!]
    return
  end function smoothAccretionConstructorParameters

  function smoothAccretionConstructorInternal(massHaloDeclineFactor,timeEarliest,cosmologyFunctions_,darkMatterHaloMassAccretionHistory_,mergerTreeMassResolution_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuilderSmoothAccretion} merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeBuilderSmoothAccretion       )                        :: self
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    class           (mergerTreeMassResolutionClass          ), intent(in   ), target :: mergerTreeMassResolution_
    double precision                                         , intent(in   )         :: massHaloDeclineFactor              , timeEarliest
    !![
    <constructorAssign variables="massHaloDeclineFactor, timeEarliest, *cosmologyFunctions_, *darkMatterHaloMassAccretionHistory_, *mergerTreeMassResolution_"/>
    !!]

    if (timeEarliest > 0.0d0) then
       self%redshiftEarliest=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeEarliest))
    else
       self%redshiftEarliest=huge(0.0d0)
    end if
    return
  end function smoothAccretionConstructorInternal

  subroutine smoothAccretionDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuilderSmoothAccretion} merger tree constructor class.
    !!}
    implicit none
    type(mergerTreeBuilderSmoothAccretion), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%mergerTreeMassResolution_"          />
    !!]
    return
  end subroutine smoothAccretionDestructor

  subroutine smoothAccretionBuild(self,tree)
    !!{
    Build a merger tree with a smooth mass accretion history.
    !!}
    use            :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    class           (mergerTreeBuilderSmoothAccretion), intent(inout), target :: self
    type            (mergerTree                      ), intent(inout), target :: tree
    type            (treeNode                        ), pointer               :: nodeCurrent   , nodeNew
    class           (nodeComponentBasic              ), pointer               :: basic
    integer         (kind=kind_int8                  )                        :: indexNode
    double precision                                                          :: massNode      , timeNode, &
         &                                                                       massResolution

    ! Get the mass resolution for this tree.
    massResolution =  self%       mergerTreeMassResolution_%resolution(tree)
    ! Extract the base node properties.
    nodeCurrent    => tree                                 %nodeBase
    basic          => nodeCurrent                          %basic     (    )
    indexNode      =  nodeCurrent                          %index     (    )
    ! Begin building. Step backward in time, creating nodes until a sufficiently low mass has been reached.
    nodeCurrent    => tree                                 %nodeBase
    massNode       =  basic                                %mass      (    )
    timeNode       =  basic                                %time      (    )
    do while (massNode > massResolution .and. timeNode > self%timeEarliest)
       indexNode =  +indexNode   &
            &       +1_kind_int8
       nodeNew   =>  treeNode     (index     =indexNode,hostTree=tree)
       basic     =>  nodeNew%basic(autoCreate=.true.                 )
       ! Adjust the mass by the specified factor.
       massNode  =  +      massNode             &
            &       *self%massHaloDeclineFactor
       ! Find the time corresponding to this mass.
       timeNode  =  +self%darkMatterHaloMassAccretionHistory_%time(tree%nodeBase,massNode)
       ! Set the properties of the new node.
       call basic%            massSet(massNode)
       call basic%            timeSet(timeNode)
       call basic%timeLastIsolatedSet(timeNode)
       ! Create parent and child links.
       nodeCurrent%firstChild => nodeNew
       nodeNew    %parent     => nodeCurrent
       ! Move the current node to the new node.
       nodeCurrent            => nodeNew
    end do
    return
  end subroutine smoothAccretionBuild
