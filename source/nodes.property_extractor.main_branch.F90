!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Contains a module which implements a node property extractor which reports if a node is on the main branch of its merger
  tree.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMainBranchStatus">
   <description>
    A node property extractor class which extracts the status of each node with respect to the main branch of its merger
    tree. The status will be extracted as {\normalfont \ttfamily nodeIsOnMainBranch}, with a value of 1 indicating that the
    node is a primary progenitor of the final halo (i.e. is on the main branch of the tree) and a value of 0 indicating that it
    is not.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorMainBranchStatus
     !!{
     A stelalr mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => mainBranchStatusExtract
     procedure :: type        => mainBranchStatusType
     procedure :: name        => mainBranchStatusName
     procedure :: description => mainBranchStatusDescription
  end type nodePropertyExtractorMainBranchStatus

  interface nodePropertyExtractorMainBranchStatus
     !!{
     Constructors for the ``mainBranchStatus'' output analysis class.
     !!}
     module procedure mainBranchStatusConstructorParameters
  end interface nodePropertyExtractorMainBranchStatus

contains

  function mainBranchStatusConstructorParameters(parameters)
    !!{
    Constructor for the {\normalfont \ttfamily mainBranchStatus} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMainBranchStatus)                :: mainBranchStatusConstructorParameters
    type(inputParameters                      ), intent(inout) :: parameters

    mainBranchStatusConstructorParameters=nodePropertyExtractorMainBranchStatus()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mainBranchStatusConstructorParameters

  function mainBranchStatusExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily mainBranchStatus} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                            )                          :: mainBranchStatusExtract
    class           (nodePropertyExtractorMainBranchStatus), intent(inout)           :: self
    type            (treeNode                             ), intent(inout), target   :: node
    double precision                                       , intent(in   )           :: time
    type            (multiCounter                         ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance, time

    if (node%isOnMainBranch()) then
       mainBranchStatusExtract=1
    else
       mainBranchStatusExtract=0
    end if
    return
  end function mainBranchStatusExtract

  integer function mainBranchStatusType(self)
    !!{
    Return the type of the stellar mass property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorMainBranchStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    mainBranchStatusType=outputAnalysisPropertyTypeLinear
    return
  end function mainBranchStatusType

  function mainBranchStatusName(self)
    !!{
    Return the name of the mainBranchStatus property.
    !!}
    implicit none
    type (varying_string                       )                :: mainBranchStatusName
    class(nodePropertyExtractorMainBranchStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    mainBranchStatusName=var_str('nodeIsOnMainBranch')
    return
  end function mainBranchStatusName

  function mainBranchStatusDescription(self)
    !!{
    Return a description of the mainBranchStatus property.
    !!}
    implicit none
    type (varying_string                       )                :: mainBranchStatusDescription
    class(nodePropertyExtractorMainBranchStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    mainBranchStatusDescription=var_str('Indicates if the node is on the main branch of the merger tree (0|1).')
    return
  end function mainBranchStatusDescription
