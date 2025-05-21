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
  Implements a node operator class that promotes sub-sub-halos.
  !!}

  !![
  <nodeOperator name="nodeOperatorSubsubhaloPromotion">
   <description>A node operator class that promotes sub-sub-halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSubsubhaloPromotion
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
   contains
     procedure :: differentialEvolution => subsubhaloPromotionDifferentialEvolution
  end type nodeOperatorSubsubhaloPromotion
  
  interface nodeOperatorSubsubhaloPromotion
     !!{
     Constructors for the \refClass{nodeOperatorSubsubhaloPromotion} node operator class.
     !!}
     module procedure subsubhaloPromotionConstructorParameters
     module procedure subsubhaloPromotionConstructorInternal
  end interface nodeOperatorSubsubhaloPromotion
  
contains
  
  function subsubhaloPromotionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSubsubhaloPromotion} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSubsubhaloPromotion)                :: self
    type (inputParameters                ), intent(inout) :: parameters

    self=nodeOperatorSubsubhaloPromotion()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function subsubhaloPromotionConstructorParameters

  function subsubhaloPromotionConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSubsubhaloPromotion} node operator class.
    !!}
    use:: Error           , only : Component_List           , Error_Report
    use:: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none
    type(nodeOperatorSubsubhaloPromotion) :: self

    if (.not.defaultSatelliteComponent%positionIsGettable())                                                             &
            & call Error_Report                                                                                          &
            &      (                                                                                                     &
            &       'sub-subhalo promotion requires that the position property of the satellite component be gettable'// &
            &       Component_List(                                                                                      &
            &                      'satellite'                                                                        ,  &
            &                       defaultSatelliteComponent%positionAttributeMatch(requireGettable=.true.)             &
            &                     )                                                                                   // &
            &       {introspection:location}                                                                             &
            &      )
    return
  end function subsubhaloPromotionConstructorInternal

  subroutine subsubhaloPromotionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Determine if sub-sub-halos should be promoted.
    !!}
    use :: Galacticus_Nodes          , only : propertyInactive           , nodeComponentSatellite
    use :: Galactic_Structure_Options, only : componentTypeDarkMatterOnly, massTypeDark
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nodeOperatorSubsubhaloPromotion), intent(inout), target  :: self
    type            (treeNode                       ), intent(inout), target  :: node
    logical                                          , intent(inout)          :: interrupt
    procedure       (interruptTask                  ), intent(inout), pointer :: functionInterrupt
    integer                                          , intent(in   )          :: propertyType
    class           (nodeComponentSatellite         )               , pointer :: satellite        , satelliteHost
    class           (massDistributionClass          )               , pointer :: massDistribution_
    double precision                                 , dimension(3)           :: positionSatellite
    double precision                                                          :: radiusSatellite  , massBoundHost, &
         &                                                                       massEnclosedHost
    !$GLC attributes unused :: self
    
    ! Return if inactive variables are requested.
    if (propertyInactive(propertyType)            ) return
    ! Return if this is not a sub-sub-halo.
    if (.not.associated(node%parent              )) return
    if (.not.           node%parent%isSatellite() ) return
    ! Find orbital radius of this node.
    satellite         => node     %satellite()
    positionSatellite =  satellite%position ()
    radiusSatellite   =  sqrt(sum(positionSatellite**2))
    ! Determine the mass enclosed by this orbit and the bound mass of the host.
    massDistribution_ => node             %parent%massDistribution    (componentTypeDarkMatterOnly,massTypeDark)
    satelliteHost     => node             %parent%satellite           (                                        )
    massBoundHost     =  satelliteHost           %boundMass           (                                        )
    massEnclosedHost  =  massDistribution_       %massEnclosedBySphere(radiusSatellite                         )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! If the satellite is within the radius enclosing the total bound mass of the host it will not be promoted.
    if (massEnclosedHost <= massBoundHost) return
    ! The satellite is outside the current bound radius of the host, so should be promoted. Trigger an interrupt.
    interrupt         =  .true.
    functionInterrupt => subsubhaloPromotionPromote
    return
  end subroutine subsubhaloPromotionDifferentialEvolution
  
  subroutine subsubhaloPromotionPromote(node,timeEnd)
    !!{
    Promote a sub-sub-halo into its host's host.
    !!}
    use :: Satellite_Promotion, only : Satellite_Move_To_New_Host
    implicit none
    type            (treeNode), intent(inout), target   :: node
    double precision          , intent(in   ), optional :: timeEnd
    !$GLC attributes unused :: timeEnd
    
    ! Move the sub-sub-halo into its host's host.
    call Satellite_Move_To_New_Host(node,node%parent%parent)
    return
  end subroutine subsubhaloPromotionPromote
