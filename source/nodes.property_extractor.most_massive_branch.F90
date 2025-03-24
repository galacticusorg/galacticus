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
  Contains a node property extractor which reports if a node is on the most massive branch of its merger tree.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorBranchMostMassive">
   <description>
    A node property extractor class which indicates if a node is on the most massive branch of its tree. The status will be
    extracted as {\normalfont \ttfamily nodeIsOnMostMassiveBranch}, with a value of 1 indicating that the node is on the most
    massive branch and a value of 0 indicating that it is not.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorBranchMostMassive
     !!{
     A node property extractor class which indicates if a node is on the most massive branch of its tree. The status will be
     extracted as {\normalfont \ttfamily nodeIsOnMostMassiveBranch}, with a value of 1 indicating that the node is on the most
     massive branch and a value of 0 indicating that it is not.
     !!}
     private
     integer :: isMostMassiveBranchID
   contains
     procedure :: extract     => branchMostMassiveExtract
     procedure :: name        => branchMostMassiveName
     procedure :: description => branchMostMassiveDescription
  end type nodePropertyExtractorBranchMostMassive

  interface nodePropertyExtractorBranchMostMassive
     !!{
     Constructors for the {\normalfont \ttfamily branchMostMassive} output analysis class.
     !!}
     module procedure branchMostMassiveConstructorParameters
     module procedure branchMostMassiveConstructorInternal
  end interface nodePropertyExtractorBranchMostMassive

contains

  function branchMostMassiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily branchMostMassive} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBranchMostMassive)               :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=nodePropertyExtractorBranchMostMassive()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function branchMostMassiveConstructorParameters

  function branchMostMassiveConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily branchMostMassive} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBranchMostMassive) :: self

    !![
    <addMetaProperty component="basic" name="isMostMassiveBranch" type="integer" id="self%isMostMassiveBranchID" isCreator="no"/>
    !!]
    return
  end function branchMostMassiveConstructorInternal

  function branchMostMassiveExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily branchMostMassive} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                             )                          :: branchMostMassiveExtract
    class           (nodePropertyExtractorBranchMostMassive), intent(inout)           :: self
    type            (treeNode                              ), intent(inout), target   :: node
    double precision                                        , intent(in   )           :: time
    type            (multiCounter                          ), intent(inout), optional :: instance
    class           (nodeComponentBasic                    )               , pointer  :: basic
    logical                                                                           :: isMostMassiveBranch
    !$GLC attributes unused :: instance, time

    basic               => node %basic                      (                          )
    isMostMassiveBranch =  basic%integerRank0MetaPropertyGet(self%isMostMassiveBranchID) == 1
    if (isMostMassiveBranch) then
       branchMostMassiveExtract=1
    else
       branchMostMassiveExtract=0
    end if
    return
  end function branchMostMassiveExtract

  function branchMostMassiveName(self)
    !!{
    Return the name of the branchMostMassive property.
    !!}
    implicit none
    type (varying_string                        )                :: branchMostMassiveName
    class(nodePropertyExtractorBranchMostMassive), intent(inout) :: self
    !$GLC attributes unused :: self
    
    branchMostMassiveName=var_str('nodeIsOnMostMassiveBranch')
    return
  end function branchMostMassiveName
  
  function branchMostMassiveDescription(self)
    !!{
    Return a description of the branchMostMassive property.
    !!}
    implicit none
    type (varying_string                        )                :: branchMostMassiveDescription
    class(nodePropertyExtractorBranchMostMassive), intent(inout) :: self
    !$GLC attributes unused :: self

    branchMostMassiveDescription=var_str('Indicates if the node is on the most massive branch of its tree (0|1).')
    return
  end function branchMostMassiveDescription
