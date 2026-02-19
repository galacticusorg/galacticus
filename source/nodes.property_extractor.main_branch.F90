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
  Implements a node property extractor which reports if a node is on the main branch of its merger
  tree.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMainBranchStatus">
   <description>
    A node property extractor class which extracts the status of each node with respect to the main branch of its merger
    tree. The status will be extracted as {\normalfont \ttfamily nodeIsOnMainBranch}, with a value of 1 indicating that the
    node is a primary progenitor of the final halo (i.e. is on the main branch of the tree) and a value of 0 indicating that it
    is not.

    If {\normalfont \ttfamily [includeSubhalos]} is set to true then subhalos of the main branch halo are also assigned a value of
    1 (with subhalos of non-main branch halos assigned a value of 0). Otherwise, all subhalos are assigned a value of 0.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorMainBranchStatus
     !!{
     A node property extractor class which extracts the status of each node with respect to the main branch of its merger
     tree. The status will be extracted as {\normalfont \ttfamily nodeIsOnMainBranch}, with a value of 1 indicating that the
     node is a primary progenitor of the final halo (i.e. is on the main branch of the tree) and a value of 0 indicating that it
     is not.

     If {\normalfont \ttfamily [includeSubhalos]} is set to true then subhalos of the main branch halo are also assigned a value of
     1 (with subhalos of non-main branch halos assigned a value of 0). Otherwise, all subhalos are assigned a value of 0.
      !!}
     private
     logical :: includeSubhalos
   contains
     procedure :: extract     => mainBranchStatusExtract
     procedure :: name        => mainBranchStatusName
     procedure :: description => mainBranchStatusDescription
  end type nodePropertyExtractorMainBranchStatus

  interface nodePropertyExtractorMainBranchStatus
     !!{
     Constructors for the \refClass{nodePropertyExtractorMainBranchStatus} output analysis class.
     !!}
     module procedure mainBranchStatusConstructorParameters
     module procedure mainBranchStatusConstructorInternal
  end interface nodePropertyExtractorMainBranchStatus

contains

  function mainBranchStatusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMainBranchStatus} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodePropertyExtractorMainBranchStatus)                :: self
    type   (inputParameters                      ), intent(inout) :: parameters
    logical                                                       :: includeSubhalos

    !![
    <inputParameter>
      <name>includeSubhalos</name>
      <description>
	If set to true then subhalos of the main branch halo are also assigned a value of 1 (with subhalos of non-main branch
	halos assigned a value of 0). Otherwise, all subhalos are assigned a value of 0.
      </description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    !!]
    self=nodePropertyExtractorMainBranchStatus(includeSubhalos)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mainBranchStatusConstructorParameters

  function mainBranchStatusConstructorInternal(includeSubhalos) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMainBranchStatus} node property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorMainBranchStatus)                :: self
    logical                                       , intent(in   ) :: includeSubhalos
    !![
    <constructorAssign variables="includeSubhalos"/>
    !!]
   
    return
  end function mainBranchStatusConstructorInternal

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
    type            (treeNode                             )               , pointer  :: nodeHost
    !$GLC attributes unused :: self, instance, time

    if (node%isSatellite()) then
       if (self%includeSubhalos) then
          nodeHost => node
          do while (nodeHost%isSatellite())
             nodeHost => nodeHost%parent
          end do
       else
          mainBranchStatusExtract=0
          return  
       end if
    else
       nodeHost => node
    end if
    if (nodeHost%isOnMainBranch()) then
       mainBranchStatusExtract=1
    else
       mainBranchStatusExtract=0
    end if
    return
  end function mainBranchStatusExtract

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
