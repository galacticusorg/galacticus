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
  Implements a node operator class that creates black hole seeds in nodes.
  !!}

  use :: Black_Hole_Seeds, only : blackHoleSeedsClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesSeed">
   <description>A node operator class that create black hole seeds in nodes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesSeed
     !!{
     A node operator class that creates black hole seeds in nodes. 
     !!}
     private
     class(blackHoleSeedsClass), pointer :: blackHoleSeeds_                  => null()
     integer                             :: blackHoleSeedsFormationChannelID
   contains
     final     ::                          blackHolesSeedDestructor
     procedure :: nodeInitialize        => blackHolesSeedNodeInitialize
     procedure :: differentialEvolution => blackHoleSeedDifferentialEvolution
  end type nodeOperatorBlackHolesSeed
  
  interface nodeOperatorBlackHolesSeed
     !!{
     Constructors for the {\normalfont \ttfamily blackHolesSeed} node operator class.
     !!}
     module procedure blackHolesSeedConstructorParameters
     module procedure blackHolesSeedConstructorInternal
  end interface nodeOperatorBlackHolesSeed

  ! Submodule-scope variable used in create function.
  class           (nodeOperatorBlackHolesSeed), pointer :: self_     => null()
  double precision                                      :: massSeed_          , spinSeed_
  !$omp threadprivate(self_,massSeed_,spinSeed_)

contains

  function blackHolesSeedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesSeed)                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(blackHoleSeedsClass       ), pointer       :: blackHoleSeeds_
    
    !![
    <objectBuilder class="blackHoleSeeds" name="blackHoleSeeds_" source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesSeed(blackHoleSeeds_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleSeeds_"/>
    !!]
    return
  end function blackHolesSeedConstructorParameters

  function blackHolesSeedConstructorInternal(blackHoleSeeds_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHolesSeed} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesSeed)                        :: self
    class(blackHoleSeedsClass       ), intent(in   ), target :: blackHoleSeeds_
    !![
    <constructorAssign variables="*blackHoleSeeds_"/>
    !!]

    !![
    <addMetaProperty component="blackHole" type="integer" name="blackHoleSeedsFormationChannel" id="self%blackHoleSeedsFormationChannelID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function blackHolesSeedConstructorInternal

  subroutine blackHolesSeedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily blackHolesSeed} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesSeed), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleSeeds_"/>
    !!]
    return
  end subroutine blackHolesSeedDestructor

  subroutine blackHolesSeedNodeInitialize(self,node)
    !!{
    Create any initial black hole seeds.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    use :: Black_Hole_Seeds, only : enumerationBlackHoleFormationChannelType
    implicit none
    class           (nodeOperatorBlackHolesSeed              ), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBlackHole                  )               , pointer :: blackHole
    double precision                                                                   :: massSeed
    type            (enumerationBlackHoleFormationChannelType)                         :: formationChannel

    massSeed=self%blackHoleSeeds_%mass(node)
    ! Create a black hole component only if the seed mass is non-zero.
    if (massSeed > 0.0d0) then
       blackHole        => node                %blackHole       (autoCreate=.true.)
       formationChannel =  self%blackHoleSeeds_%formationChannel(           node  )
       call        blackHole%massSet                    (                                                                       massSeed  )
       call        blackHole%integerRank0MetaPropertySet(self%blackHoleSeedsFormationChannelID,formationChannel                %ID        )
       if (blackHole%spinIsSettable()) &
            & call blackHole%spinSet                    (                                      self            %blackHoleSeeds_%spin(node))
    end if
    return
  end subroutine blackHolesSeedNodeInitialize
  
  subroutine blackHoleSeedDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
      Create any initial black hole seeds via interrupt.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeOperatorBlackHolesSeed), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    logical                                     , intent(inout),         :: interrupt
    procedure       (interruptTask             ), intent(inout), pointer :: functionInterrupt
    integer                                     , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole    )               , pointer :: blackHole
    double precision                                                     :: massSeed         , spinSeed
    
    blackHole => node%blackHole()

    select type (blackHole)
      type is (nodeComponentBlackHole)
        massSeed=self%blackHoleSeeds_%mass(node)
        spinSeed=self%blackHoleSeeds_%spin(node)
        ! Create a black hole component only if the seed mass is non-zero and the type is non-standard.
        if (massSeed > 0.0d0) then
          self_     => self
          massSeed_ =  massSeed
          spinSeed_ =  spinSeed
          interrupt =  .true.
          functionInterrupt => blackHoleCreate
        end if
      class default
        ! A black hole already exists - nothing to do.
    end select
    return
  end subroutine blackHoleSeedDifferentialEvolution 

  subroutine blackHoleCreate(node,timeEnd)
  !!{
      Creates the black hole via interrupt.
  !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    type            (treeNode              ), intent(inout), target  :: node
    double precision                        , intent(in   ), optional:: timeEnd
    class           (nodeComponentBlackHole),                pointer :: blackHole 
    !$GLC attributes unused :: timeEnd

    blackHole=> node%blackHole(autoCreate=.true.)
    call        blackHole%massSet(massSeed_)
    call        blackHole%integerRank0MetaPropertySet(self_%blackHoleSeedsFormationChannelID,self_%blackHoleSeeds_%formationChannel(node))
    if (blackHole%spinIsSettable()) &
         & call blackHole%spinSet(spinSeed_)
    return 
  end subroutine blackHoleCreate
