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

  !+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

  !!{
  Implements a node operator class that applies tidal mass loss to orbiting satellite halos.
  !!}

  use :: Satellite_Tidal_Heating, only : satelliteTidalHeatingRateClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorSatelliteTidalHeating">
   <description>A node operator class that applies tidal mass loss to orbiting satellite halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteTidalHeating
     !!{
     A node operator class that applies tidal mass loss to orbiting satellite halos.
     !!}
     private
     class           (satelliteTidalHeatingRateClass), pointer :: satelliteTidalHeatingRate_ => null()
     class           (darkMatterHaloScaleClass      ), pointer :: darkMatterHaloScale_       => null()
     double precision                                          :: efficiencyDecay
     logical                                                   :: applyPreInfall
   contains
     final     ::                          satelliteTidalHeatingRateDestructor
     procedure :: differentialEvolution => satelliteTidalHeatingRateDifferentialEvolution
  end type nodeOperatorSatelliteTidalHeating
  
  interface nodeOperatorSatelliteTidalHeating
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteTidalHeating} node operator class.
     !!}
     module procedure satelliteTidalHeatingRateConstructorParameters
     module procedure satelliteTidalHeatingRateConstructorInternal
  end interface nodeOperatorSatelliteTidalHeating
  
contains

  function satelliteTidalHeatingRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteTidalHeating} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Error           , only : Error_Report
    implicit none
    type            (nodeOperatorSatelliteTidalHeating)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (satelliteTidalHeatingRateClass   ), pointer       :: satelliteTidalHeatingRate_
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    double precision                                                   :: efficiencyDecay
    logical                                                            :: applyPreInfall

    !![
    <inputParameter>
      <name>efficiencyDecay</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Efficiency of the decay of the tidal tensor integral.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>applyPreInfall</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, tidal heating is applied pre-infall.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteTidalHeatingRate" name="satelliteTidalHeatingRate_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    !!]
    self=nodeOperatorSatelliteTidalHeating(efficiencyDecay,applyPreInfall,satelliteTidalHeatingRate_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalHeatingRate_"/>
    <objectDestructor name="darkMatterHaloScale_"      />
    !!]
    return
  end function satelliteTidalHeatingRateConstructorParameters

  function satelliteTidalHeatingRateConstructorInternal(efficiencyDecay,applyPreInfall,satelliteTidalHeatingRate_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteTidalHeating} node operator class.
    !!}
    implicit none
    type            (nodeOperatorSatelliteTidalHeating)                        :: self
    class           (satelliteTidalHeatingRateClass   ), intent(in   ), target :: satelliteTidalHeatingRate_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                   , intent(in   )         :: efficiencyDecay
    logical                                            , intent(in   )         :: applyPreInfall
    !![
    <constructorAssign variables="efficiencyDecay, applyPreInfall, *satelliteTidalHeatingRate_, *darkMatterHaloScale_"/>
    !!]

    return
  end function satelliteTidalHeatingRateConstructorInternal

  subroutine satelliteTidalHeatingRateDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteTidalHeating} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteTidalHeating), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalHeatingRate_"/>
    <objectDestructor name="self%darkMatterHaloScale_"      />
    !!]
    return
  end subroutine satelliteTidalHeatingRateDestructor
  
  subroutine satelliteTidalHeatingRateDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass loss from a satellite due to tidal stripping.
    !!}
    use :: Coordinates                     , only : coordinateCartesian   , assignment(=)
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, nodeComponentBasic
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear              , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Tensors                         , only : assignment(=)         , operator(*)   , tensorRank2Dimension3Symmetric
    use :: Vectors                         , only : Vector_Magnitude      , Vector_Product
    implicit none
    class           (nodeOperatorSatelliteTidalHeating), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentBasic               )               , pointer :: basic                         , basicHost
    class           (nodeComponentSatellite           )               , pointer :: satellite
    type            (treeNode                         )               , pointer :: nodeHost
    class           (massDistributionClass            )               , pointer :: massDistribution_
    double precision                                   , dimension(3)           :: position                      , velocity
    double precision                                   , parameter              :: frequencyFractionalTiny=1.0d-6
    double precision                                                            :: radius                        , orbitalPeriod            , &
         &                                                                         radialFrequency               , angularFrequency         , &
         &                                                                         timescaleDynamical            , frequencyOrbital
    type            (tensorRank2Dimension3Symmetric)                            :: tidalTensor                   , tidalTensorPathIntegrated
    type            (coordinateCartesian           )                            :: coordinates
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if     (                                   &
         &          node%isOnMainBranch     () &
         &  .or.                               &
         &   (                                 &
         &     .not.self%applyPreInfall        &
         &    .and.                            &
         &     .not.node%isSatellite        () &
         &   )                                 &
         & ) return
    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%mergesWith()
    else
       ! Walk up the branch to find the node with which this node will merge.
       nodeHost => node
       do while (nodeHost%isPrimaryProgenitor())
          nodeHost => nodeHost%parent
       end do
       nodeHost  => nodeHost%parent%firstChild
       ! Follow that node back through descendants until we find the node at the corresponding time.
       basic     => node    %basic()
       basicHost => nodeHost%basic()
       do while (basicHost%time() > basic%time())
          if (.not.associated(nodeHost%firstChild)) exit
          nodeHost  => nodeHost%firstChild
          basicHost => nodeHost%basic     ()
       end do
    end if
    ! Get required quantities from satellite and host nodes.
    satellite                 => node     %satellite                (        )
    position                  =  satellite%position                 (        )
    velocity                  =  satellite%velocity                 (        )
    tidalTensorPathIntegrated =  satellite%tidalTensorPathIntegrated(        )
    radius                    =  Vector_Magnitude                   (position)
    if (radius <= 0.0d0) return ! Do not compute rates at zero radius.
    ! Calculate tidal tensor and rate of change of integrated tidal tensor.
    coordinates       =  position
    massDistribution_ => nodeHost         %massDistribution(           )
    tidalTensor       =  massDistribution_%tidalTensor     (coordinates)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the orbital frequency. We use the larger of the angular and radial frequencies to avoid numerical problems for
    ! purely radial or purely circular orbits.
    angularFrequency=+Vector_Magnitude(Vector_Product(position,velocity)) &
         &           /radius**2                                           &
         &           *kilo                                                &
         &           *gigaYear                                            &
         &           /megaParsec
    radialFrequency =+abs             (   Dot_Product(position,velocity)) &
         &           /radius**2                                           &
         &           *kilo                                                &
         &           *gigaYear                                            &
         &           /megaParsec
    frequencyOrbital=max(                  &
         &               angularFrequency, &
         &               radialFrequency   &
         &              )
    ! Find the orbital period.
    timescaleDynamical=self%darkMatterHaloScale_%timescaleDynamical(nodeHost)
    if (frequencyOrbital > frequencyFractionalTiny/timescaleDynamical) then
       orbitalPeriod=+2.0d0                 &
            &        *Pi                    &
            &        /max(                  &
            &             angularFrequency, &
            &             radialFrequency   &
            &            )
    else
       orbitalPeriod=+timescaleDynamical
    end if
    ! Calculate integrated tidal tensor, and heating rates.
    call satellite%tidalTensorPathIntegratedRate(                                                   &
         &                                       +     tidalTensor                                  &
         &                                       -self%efficiencyDecay                              &
         &                                       *     tidalTensorPathIntegrated                    &
         &                                       /     orbitalPeriod                                &
         &                                      )
    call satellite%tidalHeatingNormalizedRate   (                                                   &
         &                                       +self%satelliteTidalHeatingRate_%heatingRate(node) &
         &                                      )
    return
  end subroutine satelliteTidalHeatingRateDifferentialEvolution
  
