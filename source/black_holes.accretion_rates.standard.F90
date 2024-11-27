!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements the standard black hole accretion rate calculation.
  !!}

  use :: Black_Hole_Binary_Separations, only : blackHoleBinarySeparationGrowthRateClass
  use :: Accretion_Disks              , only : accretionDisksClass
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass
  use :: Cooling_Radii                , only : coolingRadiusClass
  use :: Dark_Matter_Halo_Scales      , only : darkMatterHaloScaleClass

  !![
  <blackHoleAccretionRate name="blackHoleAccretionRateStandard">
   <description>
    The standard black hole accretion rate calculation.
   </description>
  </blackHoleAccretionRate>
  !!]
  type, extends(blackHoleAccretionRateClass) :: blackHoleAccretionRateStandard
     !!{
     The standard black hole accretion rate calculation.      
     !!}
     private
     class           (blackHoleBinarySeparationGrowthRateClass), pointer :: blackHoleBinarySeparationGrowthRate_   => null()
     class           (accretionDisksClass                     ), pointer :: accretionDisks_                        => null()
     class           (hotHaloTemperatureProfileClass          ), pointer :: hotHaloTemperatureProfile_             => null()
     class           (coolingRadiusClass                      ), pointer :: coolingRadius_                         => null()
     class           (darkMatterHaloScaleClass                ), pointer :: darkMatterHaloScale_                   => null()
     ! Enhancement factors for the accretion rate.
     double precision                                                    :: bondiHoyleAccretionEnhancementHotHalo           , bondiHoyleAccretionEnhancementSpheroid, &
      &                                                                     bondiHoyleAccretionEnhancementNSC
     ! Temperature of accreting gas.
     double precision                                                    :: bondiHoyleAccretionTemperatureSpheroid          , bondiHoyleAccretionTemperatureNSC
     ! Control for hot mode only accretion.
     logical                                                             :: bondiHoyleAccretionHotModeOnly
     ! Record of whether cold mode is explicitly tracked.
     logical                                                             :: coldModeTracked
     ! Control for Accretion
     logical                                                             :: stopAccretion  

   contains
     !![
     <methods>
       <method method="hotModeFraction" description="Compute the fraction of a halo which is in the hot accretion mode."/>
     </methods>
     !!]
     final     ::                    standardDestructor
     procedure :: rateAccretion   => standardRateAccretion
     procedure :: hotModeFraction => standardHotModeFraction
  end type blackHoleAccretionRateStandard

  interface blackHoleAccretionRateStandard
     !!{
     Constructors for the {\normalfont \ttfamily standard} black hole accretion rate class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface blackHoleAccretionRateStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily standard} black hole accretion rate class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleAccretionRateStandard          )                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (blackHoleBinarySeparationGrowthRateClass), pointer       :: blackHoleBinarySeparationGrowthRate_
    class           (accretionDisksClass                     ), pointer       :: accretionDisks_
    class           (hotHaloTemperatureProfileClass          ), pointer       :: hotHaloTemperatureProfile_
    class           (coolingRadiusClass                      ), pointer       :: coolingRadius_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    double precision                                                          :: bondiHoyleAccretionEnhancementHotHalo , bondiHoyleAccretionEnhancementSpheroid, &
         &                                                                       bondiHoyleAccretionEnhancementNSC     , bondiHoyleAccretionTemperatureSpheroid, &
         &                                                                       bondiHoyleAccretionTemperatureNSC
    logical                                                                   :: bondiHoyleAccretionHotModeOnly        , stopAccretion

    !![
    <inputParameter>
      <name>stopAccretion</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines if the SMBH initialized in each galaxy accretes or not.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionEnhancementSpheroid</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The factor by which the Bondi-Hoyle accretion rate of spheroid gas onto black holes in enhanced.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionEnhancementHotHalo</name>
      <defaultValue>6.0d0</defaultValue>
      <description>The factor by which the Bondi-Hoyle accretion rate of hot halo gas onto black holes in enhanced.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionEnhancementNSC</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The factor by which the Bondi-Hoyle accretion rate of NSC gas onto black holes in enhanced.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionHotModeOnly</name>
      <defaultValue>.true.</defaultValue>
      <description>Determines whether accretion from the hot halo should only occur if the halo is in the hot accretion mode.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionTemperatureSpheroid</name>
      <defaultValue>1.0d2</defaultValue>
      <description>The assumed temperature (in Kelvin) of gas in the spheroid when computing Bondi-Hoyle accretion rates onto black holes.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionTemperatureNSC</name>
      <defaultValue>1.0d2</defaultValue>
      <description>The assumed temperature (in Kelvin) of gas in the NSC when computing Bondi-Hoyle accretion rates onto black holes.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="accretionDisks"                      name="accretionDisks_"                      source="parameters"/>
    <objectBuilder class="blackHoleBinarySeparationGrowthRate" name="blackHoleBinarySeparationGrowthRate_" source="parameters"/>
    <objectBuilder class="hotHaloTemperatureProfile"           name="hotHaloTemperatureProfile_"           source="parameters"/>
    <objectBuilder class="coolingRadius"                       name="coolingRadius_"                       source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"                 name="darkMatterHaloScale_"                 source="parameters"/>
    !!]
    self=blackHoleAccretionRateStandard(stopAccretion,bondiHoyleAccretionEnhancementHotHalo,bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementNSC,bondiHoyleAccretionTemperatureSpheroid,bondiHoyleAccretionTemperatureNSC,bondiHoyleAccretionHotModeOnly,blackHoleBinarySeparationGrowthRate_,hotHaloTemperatureProfile_,accretionDisks_,coolingRadius_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionDisks_"                     />
    <objectDestructor name="blackHoleBinarySeparationGrowthRate_"/>
    <objectDestructor name="hotHaloTemperatureProfile_"          />
    <objectDestructor name="coolingRadius_"                      />
    <objectDestructor name="darkMatterHaloScale_"                />
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(stopAccretion,bondiHoyleAccretionEnhancementHotHalo,bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementNSC,bondiHoyleAccretionTemperatureSpheroid,bondiHoyleAccretionTemperatureNSC,bondiHoyleAccretionHotModeOnly,blackHoleBinarySeparationGrowthRate_,hotHaloTemperatureProfile_,accretionDisks_,coolingRadius_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily standard} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none
    type            (blackHoleAccretionRateStandard          )                        :: self
    class           (blackHoleBinarySeparationGrowthRateClass), target, intent(in   ) :: blackHoleBinarySeparationGrowthRate_
    class           (accretionDisksClass                     ), target, intent(in   ) :: accretionDisks_
    class           (hotHaloTemperatureProfileClass          ), target, intent(in   ) :: hotHaloTemperatureProfile_
    class           (coolingRadiusClass                      ), target, intent(in   ) :: coolingRadius_
    class           (darkMatterHaloScaleClass                ), target, intent(in   ) :: darkMatterHaloScale_
    double precision                                                  , intent(in   ) :: bondiHoyleAccretionEnhancementHotHalo , bondiHoyleAccretionEnhancementSpheroid, &
         &                                                                               bondiHoyleAccretionEnhancementNSC     , bondiHoyleAccretionTemperatureNSC     , &
         &                                                                               bondiHoyleAccretionTemperatureSpheroid
    logical                                                           , intent(in   ) :: bondiHoyleAccretionHotModeOnly        , stopAccretion
    !![
    <constructorAssign variables="stopAccretion,bondiHoyleAccretionEnhancementHotHalo,bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementNSC,bondiHoyleAccretionTemperatureSpheroid,bondiHoyleAccretionTemperatureNSC,bondiHoyleAccretionHotModeOnly, *blackHoleBinarySeparationGrowthRate_, *hotHaloTemperatureProfile_, *accretionDisks_, *coolingRadius_, *darkMatterHaloScale_"/>
    !!]

    ! Check if cold mode is explicitly tracked.
    self%coldModeTracked=defaultHotHaloComponent%massColdIsGettable()
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the critical overdensity standard set barrier class.
    !!}
    implicit none
    type(blackHoleAccretionRateStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionDisks_"                     />
    <objectDestructor name="self%blackHoleBinarySeparationGrowthRate_"/>
    <objectDestructor name="self%hotHaloTemperatureProfile_"          />
    <objectDestructor name="self%coolingRadius_"                      />
    <objectDestructor name="self%darkMatterHaloScale_"                />
    !!]                                                                                                                                                                                                               
    return
  end subroutine standardDestructor

  subroutine standardRateAccretion(self,blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNSC)
    !!{
    Compute the accretion rate onto a black hole.
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Eddington_Accretion_Rate
    use :: Bondi_Hoyle_Lyttleton_Accretion , only : Bondi_Hoyle_Lyttleton_Accretion_Radius, Bondi_Hoyle_Lyttleton_Accretion_Rate
    use :: Galactic_Structure_Options      , only : componentTypeColdHalo                 , componentTypeHotHalo                , componentTypeSpheroid, componentTypeNSC, &
          &                                         coordinateSystemCylindrical           , massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole                , nodeComponentHotHalo                , nodeComponentSpheroid, nodeComponentNSC, &
          &                                         treeNode
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Jeans_Length                , Ideal_Gas_Sound_Speed
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr               , gigaYear                            , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Mass_Distributions              , only : massDistributionClass                 , kinematicsDistributionClass
    use :: Coordinates                     , only : coordinateSpherical                   , assignment(=)
    implicit none
    class           (blackHoleAccretionRateStandard), intent(inout) :: self
    class           (nodeComponentBlackHole        ), intent(inout) :: blackHole
    double precision                                , intent(  out) :: rateMassAccretionSpheroid      ,rateMassAccretionHotHalo, rateMassAccretionNSC
    type            (treeNode                      ), pointer       :: node
    class           (nodeComponentSpheroid         ), pointer       :: spheroid
    class           (nodeComponentHotHalo          ), pointer       :: hotHalo
    class           (nodeComponentNSC              ), pointer       :: NSC
    class           (massDistributionClass         ), pointer       :: massDistributionSpheroid_      , massDistributionHotHalo_, &
         &                                                             massDistributionColdHalo_      , massDistributionNSC_
    class           (kinematicsDistributionClass   ), pointer       :: kinematicsDistribution_
    type            (coordinateSpherical           )                :: coordinates
    ! Lowest gas density to consider when computing accretion rates onto black hole (in units of M☉/Mpc³).
    double precision                                , parameter     :: densityGasMinimum        =1.0d0
    double precision                                                :: radiusAccretion                , rateAccretionMaximum    , &
         &                                                             massBlackHole                  , densityGas              , &
         &                                                             temperatureHotHalo             , fractionHotMode         , &
         &                                                             lengthJeans                    , fractionColdMode        , &
         &                                                             efficiencyRadiative            , velocityRelative
         
    ! Get the host node.
    node          => blackHole%host()
    ! Get black hole mass.
    massBlackHole =  blackHole%mass()

    ! First, we check the formation channel of the BH.
    select case (blackHole%NSCChannel())
    ! If the NSC was created via NSC, the accretion is allowed.
    case (.true.)
      self%stopAccretion = .false.
    ! if the BH wasn't form via NSC collapse we need to check if the StopAccretion is set to true of false in the parameter file.
    case (.false.)
      select case (self%stopAccretion)
    ! In case that StopAccretion=true in the parameter file, this stop the accretion for all the BHs, so the rest accretion rate is 0.0.
      case (.true.)
        self%stopAccretion = .true.

    ! In case that StopAccretion=false in the parameter file, the accretion rate is computed as normal.
      case (.false.)
        self%stopAccretion = .false.
      end select
    end select

    select case (self%stopAccretion)
    case(.true.)
        rateMassAccretionSpheroid=0.0d0
        rateMassAccretionHotHalo =0.0d0
        rateMAssAccretionNSC     =0.0d0
    case(.false.)
    ! Check black hole mass is positive.
    if (massBlackHole > 0.0d0) then
       ! Compute the relative velocity of black hole and gas. We assume that relative motion arises only from the radial
       ! migration of the black hole.
       velocityRelative=+self%blackHoleBinarySeparationGrowthRate_%growthRate(blackHole) &
            &           *Mpc_per_km_per_s_To_Gyr
       ! Contribution from spheroid:
       ! Get the accretion radius. We take this to be the larger of the Bondi-Hoyle radius and the current radius position of
       ! the black hole.
       radiusAccretion=max(                                                                                                   &
            &               Bondi_Hoyle_Lyttleton_Accretion_Radius(massBlackHole,self%bondiHoyleAccretionTemperatureSpheroid) &
            &              ,blackHole%radialPosition()                                                                        &
            &             )
       ! Set the position.
       coordinates               =  [radiusAccretion,0.0d0,0.0d0]
       ! Get density of gas at the galactic center.
       massDistributionSpheroid_ => node                     %massDistribution(componentTypeSpheroid,massTypeGaseous)
       densityGas                =  massDistributionSpheroid_%density         (coordinates                          )
       ! Check if we have a non-negligible gas density.
       if (densityGas > densityGasMinimum) then
          ! Get the spheroid component.
          spheroid => node%spheroid()
          ! Get the Jeans length scale.
          lengthJeans=Ideal_Gas_Jeans_Length(self%bondiHoyleAccretionTemperatureSpheroid,densityGas)
          ! Limit the smoothing scale to the scale of the spheroid.
          lengthJeans=min(lengthJeans,spheroid%radius())
          ! If the Jeans length exceeds the Bondi-Hoyle-Lyttleton accretion radius, then recompute gas density for a larger
          ! radius, as the gas should be smoothly distributed on scales below the Jeans length.
          if (lengthJeans > radiusAccretion) then
             ! Set the position.
             coordinates=[lengthJeans,0.0d0,0.0d0]
             ! Get density of gas at the galactic center.
             densityGas =massDistributionSpheroid_%density(coordinates)
          end if
          ! Compute the accretion rate.
          rateMassAccretionSpheroid=max(                                                                                                                              &
               &                        +self%bondiHoyleAccretionEnhancementSpheroid                                                                                  &
               &                        *Bondi_Hoyle_Lyttleton_Accretion_Rate(massBlackHole,densityGas,velocityRelative,self%bondiHoyleAccretionTemperatureSpheroid), &
               &                        +0.0d0                                                                                                                        &
               &                       )
          ! Get the radiative efficiency of the accretion.
          efficiencyRadiative=self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionSpheroid)
          ! Limit the accretion rate to the Eddington limit.
          if (efficiencyRadiative > 0.0d0) rateMassAccretionSpheroid=min(rateMassAccretionSpheroid,Black_Hole_Eddington_Accretion_Rate(blackHole) /efficiencyRadiative)
       else
          ! Gas density is negative - set zero accretion rate.
          rateMassAccretionSpheroid=0.0d0
       end if
        ! Get the nuclear star cluster component
         NSC => node%NSC()

         radiusAccretion=max(                                                                                               &
               &              Bondi_Hoyle_Lyttleton_Accretion_Radius(massBlackHole,self%bondiHoyleAccretionTemperatureNSC)  &
               &             ,blackHole%radialPosition()                                                                    &
               &            )
         ! Set the position.
         coordinates         = [radiusAccretion,0.0d0,0.0d0]
         massDistributionNSC_=> node              %massDistribution(componentTypeNSC,massTypeGaseous)
         ! Get density of gas at the galactic center.
         densityGas          =massDistributionNSC_%density         (coordinates                     )
         ! Check if we have a non-negligible gas density.
         if (densityGas > densityGasMinimum) then
            lengthJeans=Ideal_Gas_Jeans_Length(self%bondiHoyleAccretionTemperatureNSC,densityGas)
            ! Limit the smoothing scale to the scale of the nuclear star cluster.
            lengthJeans=min(lengthJeans,NSC%radius())
            ! If the Jeans length exceeds the Bondi-Hoyle-Lyttleton accretion radius, then recompute gas density for a larger
            ! radius, as the gas should be smoothly distributed on scales below the Jeans length.
            if (lengthJeans > radiusAccretion) then
               ! Set the position.
               coordinates=[lengthJeans,0.0d0,0.0d0]
               ! Get density of gas at the galactic center.
               densityGas =massDistributionNSC_%density(coordinates)
            end if
            ! Compute the accretion rate.
            rateMassAccretionNSC=max(self%bondiHoyleAccretionEnhancementNSC*Bondi_Hoyle_Lyttleton_Accretion_Rate(massBlackHole&
                 &,densityGas ,velocityRelative,self%bondiHoyleAccretionTemperatureNSC),0.0d0)
            ! Get the radiative efficiency of the accretion.
            efficiencyRadiative=self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionNSC)
            ! Limit the accretion rate to the Eddington limit.
            if (efficiencyRadiative > 0.0d0) rateMassAccretionNSC=min(rateMassAccretionNSC    &
                 &,Black_Hole_Eddington_Accretion_Rate(blackHole) /efficiencyRadiative)
         else
            ! Gas density is negative - set zero accretion rate.
            rateMassAccretionNSC=0.0d0
         end if
       ! Contribution from hot halo:
       ! Get the hot halo component.
       hotHalo => node%hotHalo()
       ! Get halo gas temperature.
       massDistributionHotHalo_ => node                    %massDistribution      (componentType=componentTypeHotHalo,massType=massTypeGaseous)
       kinematicsDistribution_  => massDistributionHotHalo_%kinematicsDistribution(                                                           )
       coordinates              =  [0.0d0,0.0d0,0.0d0]
       temperatureHotHalo       =  kinematicsDistribution_%temperature(coordinates)
       ! Get the accretion radius.
       radiusAccretion=Bondi_Hoyle_Lyttleton_Accretion_Radius(massBlackHole,temperatureHotHalo)
       radiusAccretion=min(radiusAccretion,hotHalo%outerRadius())
       ! Set the position.
       coordinates=[radiusAccretion,0.0d0,0.0d0]
       ! Find the fraction of gas in the halo which is in the hot mode. Set this to unity if hot/cold mode is not to be
       ! considered.
       select case (self%bondiHoyleAccretionHotModeOnly)
       case (.true.)
          if (self%coldModeTracked) then
             fractionHotMode=1.0d0
          else
             fractionHotMode=self%hotModeFraction(node)
          end if
          fractionColdMode=0.0d0
       case (.false.)
          fractionHotMode=1.0d0
          if (self%coldModeTracked) then
             fractionColdMode=1.0d0
          else
             fractionColdMode=0.0d0
          end if
       end select
       ! Get density of gas at the galactic center - scaled by the fraction in the hot accretion mode.
       densityGas                   =  +fractionHotMode                                                                   &
            &                          *massDistributionHotHalo_%density          (coordinates                          )
       if (self%coldModeTracked.and.fractionColdMode > 0.0d0) then
          massDistributionColdHalo_ =>  node                     %massDistribution(componentTypeColdHalo,massTypeGaseous)
          densityGas                =  +densityGas                                                                        &
            &                          +fractionColdMode                                                                  &
            &                          *massDistributionColdHalo_%density         (coordinates                          )
          !![
          <objectDestructor name="massDistributionColdHalo_"/>
          !!]
       end if
       ! Check if we have a non-zero gas density.
       if (densityGas > densityGasMinimum) then
          ! Compute the accretion rate.
          rateMassAccretionHotHalo=max(                                                                                                                     &
               &                       +self%bondiHoyleAccretionEnhancementHotHalo                                                                          &
               &                       *Bondi_Hoyle_Lyttleton_Accretion_Rate(massBlackHole,densityGas,velocityRelative,temperatureHotHalo,radiusAccretion), &
               &                       +0.0d0                                                                                                               &
               &                      )
          ! Limit the accretion rate to the total mass of the hot halo, divided by the sound crossing time.
          rateAccretionMaximum=max(                                              &
               &                   +hotHalo%mass()                               &
               &                   /(                                            &
               &                     +hotHalo%outerRadius()                      &
               &                     /Ideal_Gas_Sound_Speed(temperatureHotHalo)  &
               &                     *Mpc_per_km_per_s_To_Gyr                    &
               &                    )                                          , &
               &                   +0.0d0                                        &
               &                  )
          rateMassAccretionHotHalo=min(rateMassAccretionHotHalo,rateAccretionMaximum)
          ! Get the radiative efficiency of the accretion.
          efficiencyRadiative=self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionHotHalo)
          ! Limit the accretion rate to the Eddington limit.
          if (efficiencyRadiative > 0.0d0)                                                     &
               & rateMassAccretionHotHalo=min(                                                 &
               &                              +rateMassAccretionHotHalo                      , &
               &                              +Black_Hole_Eddington_Accretion_Rate(blackHole)  &
               &                              /efficiencyRadiative                             &
               &                             )
       else
          ! No gas density, so zero accretion rate.
          rateMassAccretionHotHalo=0.0d0
       end if
       !![
       <objectDestructor name="massDistributionSpheroid_"/>
       <objectDestructor name="massDistributionHotHalo_" />
       <objectDestructor name="massDistributionNSC_"     />
       <objectDestructor name="kinematicsDistribution_"  />
       !!]
    else
       rateMassAccretionSpheroid=0.0d0
       rateMassAccretionHotHalo =0.0d0
       rateMassAccretionNSC     =0.0d0
    end if
    end select
    return
  end subroutine standardRateAccretion

  double precision function standardHotModeFraction(self,node)
    !!{
    A simple interpolating function which is used as a measure of the fraction of a halo which is in the hot accretion mode.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class           (blackHoleAccretionRateStandard), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , parameter     :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision                                , parameter     :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                                                :: coolingRadiusFractional                       , x

    coolingRadiusFractional=+self%coolingRadius_      %      radius(node) &
         &                  /self%darkMatterHaloScale_%radiusVirial(node)
    if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
       standardHotModeFraction=1.0d0
    else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
       standardHotModeFraction=0.0d0
    else
       x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
            & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
       standardHotModeFraction=x**2*(2.0d0*x-3.0d0)+1.0d0
    end if
    return
  end function standardHotModeFraction
