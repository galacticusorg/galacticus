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
  Contains a module which implements a node property extractor which reports the formation channel for the central black hole of a given node.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleFormationChannel">
   <description>
    A node property extractor class which extracts the formation channel for black hole seeds, with a value of 1 indicating the \cite{Vergara_2023} black hole seeding scenario
    and a value of 0 indicating the fixed black hole seeding scenario.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorBlackHoleFormationChannel
     !!{
     A node property extractor class which extracts the formation channel for black hole seeds, with a value of 1 indicating the \cite{Vergara_2023} black hole seeding scenario
    and a value of 0 indicating the fixed black hole seeding scenario.
     !!}
     private
     integer :: blackHoleSeedsFormationChannelID
   contains
     procedure :: extract     => blackHoleFormationChannelExtract
     procedure :: name        => blackHoleFormationChannelName
     procedure :: description => blackHoleFormationChannelDescription
  end type nodePropertyExtractorBlackHoleFormationChannel

  interface nodePropertyExtractorBlackHoleFormationChannel
     !!{
     Constructors for the ``blackHoleFormationChannel'' node property extractor class.
     !!}
     module procedure blackHoleFormationChannelConstructorParameters
     module procedure blackHoleFormationChannelConstructorInternal
  end interface nodePropertyExtractorBlackHoleFormationChannel

contains

  function blackHoleFormationChannelConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleFormationChannel} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleFormationChannel)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleFormationChannel()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleFormationChannelConstructorParameters

  function blackHoleFormationChannelConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationChannel} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleFormationChannel) :: self
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedsFormationChannel" type="integer" id="self%blackHoleSeedsFormationChannelID" isCreator="no"/>
    !!]
    return
  end function blackHoleFormationChannelConstructorInternal

  function blackHoleFormationChannelExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily blackHoleFormationChannel} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentBlackHoleStandard
    implicit none
    integer         (kind_int8                                     )                          :: blackHoleFormationChannelExtract
    class           (nodePropertyExtractorBlackHoleFormationChannel), intent(inout)           :: self
    type            (treeNode                                      ), intent(inout), target   :: node
    double precision                                                , intent(in   )           :: time
    type            (multiCounter                                  ), intent(inout), optional :: instance
    class           (nodeComponentBlackHole                        )               , pointer  :: blackHole
    !$GLC attributes unused :: instance, time

    blackHole => node%blackHole()
    select type (blackHole)
      class is (nodeComponentBlackHoleStandard)
        ! Standard class, extract the integer.
        blackHoleFormationChannelExtract= blackHole%integerRank0MetaPropertyGet(self%blackHoleSeedsFormationChannelID)
    end select
    return
  end function blackHoleFormationChannelExtract

  function blackHoleFormationChannelName(self)
    !!{
    Return the name of the blackHoleFormationChannel property.
    !!}
    implicit none
    type (varying_string                                )                :: blackHoleFormationChannelName
    class(nodePropertyExtractorBlackHoleFormationChannel), intent(inout) :: self
    !$GLC attributes unused :: self
  
    blackHoleFormationChannelName=var_str('blackHoleSeedsFormationChannel')
    return
  end function blackHoleFormationChannelName
  
  function blackHoleFormationChannelDescription(self)
    !!{
    Return a description of the blackHoleFormationChannel property.
    !!}
    implicit none
    type (varying_string                                )                :: blackHoleFormationChannelDescription
    class(nodePropertyExtractorBlackHoleFormationChannel), intent(inout) :: self
    !$GLC attributes unused :: self

    blackHoleFormationChannelDescription=var_str('Indicates the formation channel of the black hole seed 0:fixed, 1: stellar collisions in nuclear star clusters.')
    return
  end function blackHoleFormationChannelDescription
