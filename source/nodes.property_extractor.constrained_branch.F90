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
  Implements a node property extractor which reports if a node is drawn from the constrained branching rate solution.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorConstrainedStatus">
   <description>
    A node property extractor class which extracts the constrained excursion set solution status of each node. The status will be
    extracted as {\normalfont \ttfamily nodeIsConstrained}, with a value of 1 indicating that the node follows the constrained
    branching rate solution and a value of 0 indicating that it does not.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorConstrainedStatus
     !!{
     A node property extractor class which extracts the constrained excursion set solution status of each node. The status will be
     extracted as {\normalfont \ttfamily nodeIsConstrained}, with a value of 1 indicating that the node follows the constrained
     branching rate solution and a value of 0 indicating that it does not.
     !!}
     private
     integer :: isConstrainedID
   contains
     procedure :: extract     => constrainedStatusExtract
     procedure :: name        => constrainedStatusName
     procedure :: description => constrainedStatusDescription
  end type nodePropertyExtractorConstrainedStatus

  interface nodePropertyExtractorConstrainedStatus
     !!{
     Constructors for the \refClass{nodePropertyExtractorConstrainedStatus} output analysis class.
     !!}
     module procedure constrainedStatusConstructorParameters
     module procedure constrainedStatusConstructorInternal
  end interface nodePropertyExtractorConstrainedStatus

contains

  function constrainedStatusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorConstrainedStatus} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorConstrainedStatus)               :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=nodePropertyExtractorConstrainedStatus()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function constrainedStatusConstructorParameters

  function constrainedStatusConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorConstrainedStatus} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorConstrainedStatus) :: self

    ! Use a matched "meta-property" as we used in the constrained tree build controller. This allows us to recover the stored
    ! "isConstrained" state of each node that was written by that object.
    !![
    <addMetaProperty component="basic" name="isConstrained" type="integer" id="self%isConstrainedID" isCreator="no"/>
    !!]
    return
  end function constrainedStatusConstructorInternal

  function constrainedStatusExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily constrainedStatus} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                             )                          :: constrainedStatusExtract
    class           (nodePropertyExtractorConstrainedStatus), intent(inout)           :: self
    type            (treeNode                              ), intent(inout), target   :: node
    double precision                                        , intent(in   )           :: time
    type            (multiCounter                          ), intent(inout), optional :: instance
    class           (nodeComponentBasic                    )               , pointer  :: basic
    logical                                                                           :: isConstrained
    !$GLC attributes unused :: instance, time


    basic         => node %basic                      (                    )
    isConstrained =  basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    if (isConstrained) then
       constrainedStatusExtract=1
    else
       constrainedStatusExtract=0
    end if
    return
  end function constrainedStatusExtract

  function constrainedStatusName(self)
    !!{
    Return the name of the constrainedStatus property.
    !!}
    implicit none
    type (varying_string                        )                :: constrainedStatusName
    class(nodePropertyExtractorConstrainedStatus), intent(inout) :: self
    !$GLC attributes unused :: self
    
    constrainedStatusName=var_str('nodeIsConstrained')
    return
  end function constrainedStatusName
  
  function constrainedStatusDescription(self)
    !!{
    Return a description of the constrainedStatus property.
    !!}
    implicit none
    type (varying_string                        )                :: constrainedStatusDescription
    class(nodePropertyExtractorConstrainedStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    constrainedStatusDescription=var_str('Indicates if the node is on the constrained branch (0|1).')
    return
  end function constrainedStatusDescription
