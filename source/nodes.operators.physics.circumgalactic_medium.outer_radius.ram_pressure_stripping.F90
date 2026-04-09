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
  Implements a node operator class that evolves the \gls{cgm} outer radius in response to ram pressure stripping.
  !!}

  use :: Cosmology_Parameters                      , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales                   , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Outflows_Reincorporations        , only : hotHaloOutflowReincorporationClass
  use :: Hot_Halo_Ram_Pressure_Stripping           , only : hotHaloRamPressureStrippingClass
  use :: Hot_Halo_Ram_Pressure_Stripping_Timescales, only : hotHaloRamPressureTimescaleClass
  use :: Hot_Halo_Outflows_Stripping               , only : hotHaloOutflowStrippingClass

  !![
  <nodeOperator name="nodeOperatorCGMOuterRadiusRamPressureStripping">
   <description>
    A node operator class that evolves the \gls{cgm} outer radius in response to ram pressure stripping.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMOuterRadiusRamPressureStripping
     !!{
     A node operator class that evolves the \gls{cgm} outer radius in response to ram pressure stripping.
     !!}
     private
     class(cosmologyParametersClass           ), pointer :: cosmologyParameters_           => null()
     class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_           => null()
     class(hotHaloRamPressureStrippingClass   ), pointer :: hotHaloRamPressureStripping_   => null()
     class(hotHaloRamPressureTimescaleClass   ), pointer :: hotHaloRamPressureTimescale_   => null()
     class(hotHaloOutflowReincorporationClass ), pointer :: hotHaloOutflowReincorporation_ => null()
     class(hotHaloOutflowStrippingClass       ), pointer :: hotHaloOutflowStripping_       => null()
   contains
     final     ::                          cgmOuterRadiusRamPressureStrippingDestructor
     procedure :: differentialEvolution => cgmOuterRadiusRamPressureStrippingDifferentialEvolution
  end type nodeOperatorCGMOuterRadiusRamPressureStripping
  
  interface nodeOperatorCGMOuterRadiusRamPressureStripping
     !!{
     Constructors for the \refClass{nodeOperatorCGMOuterRadiusRamPressureStripping} node operator class.
     !!}
     module procedure cgmOuterRadiusRamPressureStrippingConstructorParameters
     module procedure cgmOuterRadiusRamPressureStrippingConstructorInternal
  end interface nodeOperatorCGMOuterRadiusRamPressureStripping
  
contains
  
  function cgmOuterRadiusRamPressureStrippingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMOuterRadiusRamPressureStripping} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorCGMOuterRadiusRamPressureStripping)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    class(darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_
    class(hotHaloRamPressureStrippingClass              ), pointer       :: hotHaloRamPressureStripping_
    class(hotHaloRamPressureTimescaleClass              ), pointer       :: hotHaloRamPressureTimescale_
    class(hotHaloOutflowReincorporationClass            ), pointer       :: hotHaloOutflowReincorporation_
    class(hotHaloOutflowStrippingClass                  ), pointer       :: hotHaloOutflowStripping_

    !![
    <objectBuilder class="cosmologyParameters"           name="cosmologyParameters_"           source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"           name="darkMatterHaloScale_"           source="parameters"/>
    <objectBuilder class="hotHaloRamPressureStripping"   name="hotHaloRamPressureStripping_"   source="parameters"/>
    <objectBuilder class="hotHaloRamPressureTimescale"   name="hotHaloRamPressureTimescale_"   source="parameters"/>
    <objectBuilder class="hotHaloOutflowReincorporation" name="hotHaloOutflowReincorporation_" source="parameters"/>
    <objectBuilder class="hotHaloOutflowStripping"       name="hotHaloOutflowStripping_"       source="parameters"/>
    !!]
    self=nodeOperatorCGMOuterRadiusRamPressureStripping(cosmologyParameters_,darkMatterHaloScale_,hotHaloRamPressureStripping_,hotHaloRamPressureTimescale_,hotHaloOutflowReincorporation_,hotHaloOutflowStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"          />
    <objectDestructor name="darkMatterHaloScale_"          />
    <objectDestructor name="hotHaloRamPressureStripping_"  />
    <objectDestructor name="hotHaloRamPressureTimescale_"  />
    <objectDestructor name="hotHaloOutflowReincorporation_"/>
    <objectDestructor name="hotHaloOutflowStripping_"      />
    !!]
    return
  end function cgmOuterRadiusRamPressureStrippingConstructorParameters

  function cgmOuterRadiusRamPressureStrippingConstructorInternal(cosmologyParameters_,darkMatterHaloScale_,hotHaloRamPressureStripping_,hotHaloRamPressureTimescale_,hotHaloOutflowReincorporation_,hotHaloOutflowStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMOuterRadiusRamPressureStripping} node operator class.
    !!}
    implicit none
    type (nodeOperatorCGMOuterRadiusRamPressureStripping)                        :: self
    class(cosmologyParametersClass                      ), intent(in   ), target :: cosmologyParameters_
    class(darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    class(hotHaloRamPressureStrippingClass              ), intent(in   ), target :: hotHaloRamPressureStripping_
    class(hotHaloRamPressureTimescaleClass              ), intent(in   ), target :: hotHaloRamPressureTimescale_
    class(hotHaloOutflowReincorporationClass            ), intent(in   ), target :: hotHaloOutflowReincorporation_
    class(hotHaloOutflowStrippingClass                  ), intent(in   ), target :: hotHaloOutflowStripping_
    !![
    <constructorAssign variables="*cosmologyParameters_, *darkMatterHaloScale_, *hotHaloRamPressureStripping_, *hotHaloRamPressureTimescale_, *hotHaloOutflowReincorporation_, *hotHaloOutflowStripping_"/>
    !!]

    return
  end function cgmOuterRadiusRamPressureStrippingConstructorInternal

  subroutine cgmOuterRadiusRamPressureStrippingDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMOuterRadiusRamPressureStripping} node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMOuterRadiusRamPressureStripping), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"          />
    <objectDestructor name="self%darkMatterHaloScale_"          />
    <objectDestructor name="self%hotHaloRamPressureStripping_"  />
    <objectDestructor name="self%hotHaloRamPressureTimescale_"  />
    <objectDestructor name="self%hotHaloOutflowReincorporation_"/>
    <objectDestructor name="self%hotHaloOutflowStripping_"      />
    !!]
    return
  end subroutine cgmOuterRadiusRamPressureStrippingDestructor
  
  subroutine cgmOuterRadiusRamPressureStrippingDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Compute the rate of change of the outer radius of the \gls{cgm} in response to ram pressure
    stripping (and return of gas).
    !!}
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentHotHalo
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , componentTypeColdHalo, massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    class           (nodeOperatorCGMOuterRadiusRamPressureStripping), intent(inout), target  :: self
    type            (treeNode                                      ), intent(inout), target  :: node
    logical                                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                                 ), intent(inout), pointer :: functionInterrupt
    integer                                                         , intent(in   )          :: propertyType
    class           (nodeComponentBasic                            )               , pointer :: basic
    class           (nodeComponentHotHalo                          )               , pointer :: hotHalo
    class           (massDistributionClass                         )               , pointer :: massDistribution_
    double precision                                                , parameter              :: radiusOuterOverRadiusVirialMinimum=1.0d-3
    type            (coordinateSpherical                           )                         :: coordinates
    double precision                                                                         :: massReturnRate                           , radiusOuter         , &
         &                                                                                      radiusVirial                             , densityAtOuterRadius, &
         &                                                                                      densityMinimum                           , massGas             , &
         &                                                                                      radiusRamPressure                        , massLossRate        , &
         &                                                                                      radiusOuterGrowthRate
    
    ! Get the hot halo component to work with.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! First handle return of mass to the circumgalactic medium, which makes the outer radius expand.    
       ! Get the rate at which mass is being returned to the circumgalactic medium.
       massReturnRate=self%hotHaloOutflowReincorporation_%rate(node)
       ! The outer radius must be increased as the halo fills up with gas.
       radiusOuter =hotHalo                     %outerRadius (    )
       radiusVirial=self   %darkMatterHaloScale_%radiusVirial(node)
       if (radiusOuter < radiusVirial) then
          coordinates          =  [radiusOuter,0.0d0,0.0d0]
          basic                => node             %basic           (                                    )
          massDistribution_    => node             %massDistribution(componentTypeHotHalo,massTypeGaseous)
          densityAtOuterRadius =  massDistribution_%density         (coordinates                         )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (radiusOuter > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             densityMinimum=+self %cosmologyParameters_%omegaBaryon ()    &
                  &         /self %cosmologyParameters_%omegaMatter ()    &
                  &         *basic                     %mass        ()    &
                  &         /                           radiusVirial  **3 &
                  &         /4.0d0                                        &
                  &         /Pi
             call hotHalo%outerRadiusRate(                           &
                  &                       +massReturnRate            &
                  &                       /4.0d0                     &
                  &                       /Pi                        &
                  &                       /radiusOuter**2            &
                  &                       /max(                      &
                  &                            densityAtOuterRadius, &
                  &                            densityMinimum        &
                  &                           )                      &
                  &                      )
             ! Otherwise, if we have a positive rate of mass return, simply grow the radius at the virial velocity.
          else if (massReturnRate > 0.0d0) then
             ! Force some growth here so the radius is not trapped at zero.
             call hotHalo%outerRadiusRate(                       &
                  &                       +massReturnRate        &
                  &                       /basic         %mass() &
                  &                       *radiusVirial          &
                  &                      )
          end if
       end if
       ! Now handle ram pressure stripping.
       if (node%isSatellite()) then
          ! For satellites, first compute the outer radius growth rate due to ram pressure stripping.
          radiusRamPressure=self%hotHaloRamPressureStripping_%radiusStripped(node)
          ! Test whether the ram pressure radius is smaller than the current outer radius of the hot gas profile.
          if     (                                           &
               &           radiusRamPressure   < radiusOuter &
               &  .and.                                      &
               &   hotHalo%angularMomentum  () > 0.0d0       &
               & ) then
             ! The ram pressure stripping radius is within the outer radius. Cause the outer radius to shrink to the ram pressure
             ! stripping radius on the ram pressure stripping timescale.
             radiusOuterGrowthRate=+(                                                 &
                  &                  +radiusRamPressure                               &
                  &                  -radiusOuter                                     &
                  &                 )                                                 &
                  &                /self%hotHaloRamPressureTimescale_%timescale(node)
          else
             radiusOuterGrowthRate=+0.0d0
          end if
          ! Now apply this growth rate and the associated rates of change of mass.
          if     (                                                                                                           &
               &   radiusOuterGrowthRate   /= 0.0d0                                                                          &
               &  .and.                                                                                                      &
               &   hotHalo%mass         () >  0.0d0                                                                          &
               &  .and.                                                                                                      &
               &   radiusOuter             <=                                   self%darkMatterHaloScale_%radiusVirial(node) &
               &  .and.                                                                                                      &
               &   radiusOuter             > radiusOuterOverRadiusVirialMinimum*self%darkMatterHaloScale_%radiusVirial(node) &
               & ) then
             coordinates          =  [radiusOuter,0.0d0,0.0d0]
             massDistribution_    => node             %massDistribution(componentTypeHotHalo,massTypeGaseous)
             densityAtOuterRadius =  massDistribution_%density         (coordinates                         )
             !![
             <objectDestructor name="massDistribution_"/>
             !!]
             massGas=hotHalo%mass()
             if (massGas > 0.0d0) then
                massLossRate=+4.0d0                    &
                     &       *Pi                       &
                     &       *densityAtOuterRadius     &
                     &       *radiusOuter          **2 &
                     &       *radiusOuterGrowthRate
                call    hotHalo%       outerRadiusRate(                           +radiusOuterGrowthRate,interrupt,functionInterrupt)
                call    hotHalo%              massRate(                           +massLossRate         ,interrupt,functionInterrupt)
                call    hotHalo%   angularMomentumRate(hotHalo%angularMomentum()*(+massLossRate/massGas),interrupt,functionInterrupt)
                call    hotHalo%        abundancesRate(hotHalo%abundances     ()*(+massLossRate/massGas),interrupt,functionInterrupt)
                call    hotHalo%         chemicalsRate(hotHalo%chemicals      ()*(+massLossRate/massGas),interrupt,functionInterrupt)
                if (.not.self%hotHaloOutflowStripping_%neverStripped(node)) then
                   call hotHalo%      strippedMassRate(                           -massLossRate         ,interrupt,functionInterrupt)
                   call hotHalo%strippedAbundancesRate(hotHalo%abundances     ()*(-massLossRate/massGas),interrupt,functionInterrupt)
                   call hotHalo% strippedChemicalsRate(hotHalo%chemicals      ()*(-massLossRate/massGas),interrupt,functionInterrupt)
                end if
             end if
             ! Handle cold mode gas.
             if (hotHalo%massColdIsSettable()) then
                ! Compute the density at the outer radius for the cold mode component.
                coordinates          =  [radiusOuter,0.0d0,0.0d0]
                massDistribution_    => node             %massDistribution(componentType=componentTypeColdHalo,massType=massTypeGaseous)
                densityAtOuterRadius =  massDistribution_%density         (              coordinates                                   )
                !![
		<objectDestructor name="massDistribution_"/>
	        !!]
                massGas=hotHalo%massCold()
                if (massGas > 0.0d0) then
                   ! Compute the mass loss rate.
                   massLossRate=+4.0d0                    &
                        &       *Pi                       &
                        &       *densityAtOuterRadius     &
                        &       *radiusOuter          **2 &
                        &       *radiusOuterGrowthRate
                   ! Adjust the rates.
                   call hotHalo%           massColdRate(+                              massLossRate        ,interrupt,functionInterrupt)
                   call hotHalo%angularMomentumColdRate(+hotHalo%angularMomentumCold()*massLossRate/massGas,interrupt,functionInterrupt)
                   call hotHalo%     abundancesColdRate(+hotHalo%abundancesCold     ()*massLossRate/massGas,interrupt,functionInterrupt)
                   call hotHalo%       strippedMassRate(-                              massLossRate        ,interrupt,functionInterrupt)
                   call hotHalo% strippedAbundancesRate(-hotHalo%abundancesCold     ()*massLossRate/massGas,interrupt,functionInterrupt)
                end if
             end if
          end if
       else
          ! For isolated halos, the outer radius should grow with the virial radius.
          call hotHalo%outerRadiusRate(self%darkMatterHaloScale_%radiusVirialGrowthRate(node),interrupt,functionInterrupt)
       end if
    end select
    return
  end subroutine cgmOuterRadiusRamPressureStrippingDifferentialEvolution
