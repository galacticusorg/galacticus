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
Implements a node property extractor for the most massive progenitor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMostMassiveProgenitor">
   <description>
    A node property extractor class which extracts a value of $1$ for the most massive progenitor node in a tree at each output
    time and $0$ for all other nodes.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorMostMassiveProgenitor
     !!{
     A node property extractor for most massive progenitor status.
     !!}
     private
     double precision            :: timePrevious
     integer         (kind_int8) :: uniqueIdMatched, uniqueIdPrevious
   contains
     procedure :: extract     => mostMassiveProgenitorExtract
     procedure :: name        => mostMassiveProgenitorName
     procedure :: description => mostMassiveProgenitorDescription
  end type nodePropertyExtractorMostMassiveProgenitor

  interface nodePropertyExtractorMostMassiveProgenitor
     !!{
     Constructors for the {\normalfont \ttfamily mostMassiveProgenitor} output analysis class.
     !!}
     module procedure mostMassiveProgenitorConstructorParameters
     module procedure mostMassiveProgenitorConstructorInternal
  end interface nodePropertyExtractorMostMassiveProgenitor

contains

  function mostMassiveProgenitorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily mostMassiveProgenitor} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMostMassiveProgenitor)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorMostMassiveProgenitor()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mostMassiveProgenitorConstructorParameters

  function mostMassiveProgenitorConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily mostMassiveProgenitor} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMostMassiveProgenitor) :: self

    self%timePrevious    =-1.0d0
    self%uniqueIDPrevious=-1_kind_int8
    self%uniqueIDMatched =-1_kind_int8
    return
  end function mostMassiveProgenitorConstructorInternal

  function mostMassiveProgenitorExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily mostMassiveProgenitor} node property extractor.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentBasic           , treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    integer         (kind_int8                                 )                          :: mostMassiveProgenitorExtract
    class           (nodePropertyExtractorMostMassiveProgenitor), intent(inout)           :: self
    type            (treeNode                                  ), intent(inout), target   :: node
    double precision                                            , intent(in   )           :: time
    type            (multiCounter                              ), intent(inout), optional :: instance
    type            (treeNode                                  ), pointer                 :: nodeCurrent
    class           (nodeComponentBasic                        ), pointer                 :: basicCurrent
    type            (mergerTreeWalkerIsolatedNodes             )                          :: treeWalker
    double precision                                                                      :: massMostMassiveProgenitor
    !$GLC attributes unused :: self, instance

    ! Find the root node in the tree.
    nodeCurrent => node
    do while (associated(nodeCurrent%parent))
       nodeCurrent => nodeCurrent%parent
    end do
    ! Check if this is the same tree, at the same time as on the previous call.
    if (time /= self%timePrevious .or. nodeCurrent%uniqueId() /= self%uniqueIdPrevious) then
       ! It is not, so record the new tree root unique ID and the new time.
       self%timePrevious    =time
       self%uniqueIdPrevious=nodeCurrent%uniqueId()
       ! Find the most massive progenitor in the tree at this time.
       massMostMassiveProgenitor=0.0d0
       treeWalker=mergerTreeWalkerIsolatedNodes(nodeCurrent%hostTree)
       do while (treeWalker%next(nodeCurrent))
          basicCurrent => nodeCurrent%basic()
          if (basicCurrent%time() == time .and. basicCurrent%mass() > massMostMassiveProgenitor) then
             self%uniqueIdMatched     =nodeCurrent %uniqueId()
             massMostMassiveProgenitor=basicCurrent%mass    ()
          end if
       end do
    end if
    if (node%uniqueId() == self%uniqueIdMatched) then
       mostMassiveProgenitorExtract=1
    else
       mostMassiveProgenitorExtract=0
    end if
    return
  end function mostMassiveProgenitorExtract


  function mostMassiveProgenitorName(self)
    !!{
    Return the name of the mostMassiveProgenitor property.
    !!}
    implicit none
    type (varying_string                            )                :: mostMassiveProgenitorName
    class(nodePropertyExtractorMostMassiveProgenitor), intent(inout) :: self
    !$GLC attributes unused :: self

    mostMassiveProgenitorName=var_str('isMostMassiveProgenitor')
    return
  end function mostMassiveProgenitorName

  function mostMassiveProgenitorDescription(self)
    !!{
    Return a description of the mostMassiveProgenitor property.
    !!}
    implicit none
    type (varying_string                            )                :: mostMassiveProgenitorDescription
    class(nodePropertyExtractorMostMassiveProgenitor), intent(inout) :: self
    !$GLC attributes unused :: self

    mostMassiveProgenitorDescription=var_str('Flag indicating if this node is the most massive progenitor in its tree at this time.')
    return
  end function mostMassiveProgenitorDescription
