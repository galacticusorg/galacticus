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

  !!{RST
  Implementation of a cooling radius class for :math:`\beta`-profile halos, assuming collisional ionization equilibrium such that cooling time scales as inverse density.
  !!}

  use :: Cooling_Times          , only : coolingTimeClass
  use :: Cooling_Times_Available, only : coolingTimeAvailableClass
  use :: Cosmology_Functions    , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <coolingRadius name="coolingRadiusBetaProfile" docformat="rst">
   <description>
   A cooling radius class for :math:`\beta`-profile halos. Computes the cooling radius by assuming that the hot gas density profile is a :math:`\beta`-profile (:math:`\rho(r) \propto [r^2+r_\mathrm{c}^2]^{-1}`), and that the cooling rate scales as density squared, :math:`\dot{E}\propto \rho^2`, such that the cooling time scales as inverse density, :math:`t_\mathrm{cool} \propto \rho^{-1}`. Consequently, the cooling radius is given by

   .. math::

      r_\mathrm{cool} = r_\mathrm{virial} \left( \left[ {t_\mathrm{avail} \over t_0} - 1 \right] \left[ {t_\mathrm{virial} \over t_0} - 1 \right]^{-1} \right)^{1/2},

   where :math:`t_0`, and :math:`t_\mathrm{virial}` are the cooling times at zero radius and the virial radius respectively.
   </description>
  </coolingRadius>
  !!]
  type, extends(coolingRadiusCoolingTime) :: coolingRadiusBetaProfile
     !!{RST
     Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available for cooling equals the cooling time.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                     betaProfileDestructor
     procedure :: radius           => betaProfileRadius
     procedure :: radiusGrowthRate => betaProfileRadiusGrowthRate
  end type coolingRadiusBetaProfile

  interface coolingRadiusBetaProfile
     !!{RST
     Constructors for the betaProfile cooling radius class.
     !!}
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface coolingRadiusBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :math:`\beta`-profile cooling radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRadiusBetaProfile      )                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(coolingTimeAvailableClass     ), pointer       :: coolingTimeAvailable_
    class(coolingTimeClass              ), pointer       :: coolingTime_
    class(darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_" source="parameters"/>
    <objectBuilder class="coolingTime"          name="coolingTime_"          source="parameters"/>
    !!]
    self=coolingRadiusBetaProfile(cosmologyFunctions_,darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="coolingTimeAvailable_"/>
    <objectDestructor name="coolingTime_"         />
    !!]
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_) result(self)
    !!{RST
    Internal constructor for the :math:`\beta`-profile cooling radius class.
    !!}
    implicit none
    type (coolingRadiusBetaProfile )                        :: self
    class(cosmologyFunctionsClass  ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    class(coolingTimeAvailableClass), intent(in   ), target :: coolingTimeAvailable_
    class(coolingTimeClass         ), intent(in   ), target :: coolingTime_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_, *coolingTimeAvailable_, *coolingTime_"/>
    !!]

    ! Initialize the state shared with other cooling time-based classes.
    call self%initialize()
    return
  end function betaProfileConstructorInternal

  subroutine betaProfileDestructor(self)
    !!{RST
    Destructor for the :math:`\beta`-profile cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(coolingRadiusBetaProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%coolingTimeAvailable_"/>
    <objectDestructor name="self%coolingTime_"         />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%radiation"            />
    !!]
    if (calculationResetEvent%isAttached(self,coolingTimeCalculationReset)) call calculationResetEvent%detach(self,coolingTimeCalculationReset)
    return
  end subroutine betaProfileDestructor

  double precision function betaProfileRadiusGrowthRate(self,node)
    !!{RST
    Returns the cooling radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    !!}
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                   , nodeComponentHotHalo       , treeNode
    use :: Mass_Distributions               , only : massDistributionClass                , kinematicsDistributionClass
    use :: Coordinates                      , only : coordinateSpherical                  , assignment(=)
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo                 , massTypeGaseous
    implicit none
    class           (coolingRadiusBetaProfile   ), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    class           (nodeComponentBasic         ), pointer       :: basic
    class           (nodeComponentHotHalo       ), pointer       :: hotHalo
    class           (massDistributionClass      ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass), pointer       :: kinematicsDistribution_
    type            (coordinateSpherical        )                :: coordinates
    double precision                                             :: coolingTimeZero        , timeAvailable          , &
         &                                                          densityZero            , massToDensityConversion, &
         &                                                          temperature            , outerRadius            , &
         &                                                          densityOuter           , coolingTimeOuter
    type            (abundances                 )                :: hotAbundances
    type            (chemicalAbundances         )                :: chemicalFractions      , chemicalMasses

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())

    ! Check if cooling radius growth rate is already computed.
    if (.not.self%radiusGrowthRateComputed) then
       ! Get the time available for cooling in node.
       timeAvailable              =  self%coolingTimeAvailable_%timeAvailable(node)
       ! Get the abundances for this node.
       hotHalo                    => node   %hotHalo   ()
       hotAbundances              =  hotHalo%abundances()
       call hotAbundances%massToMassFraction(hotHalo%mass())
       ! Get the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalFractions)
          if (hotHalo%mass() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
          else
             massToDensityConversion=0.0d0
          end if          
          ! Convert to number density per unit total mass density.
          chemicalFractions=chemicalFractions*massToDensityConversion
       end if
       ! Set epoch for radiation field.
       basic => node%basic()
       call self%radiation%timeSet(basic%time())
       ! Get the outer radius.
       outerRadius=hotHalo%outerRadius()
       ! Get the mass distribution.
       massDistribution_       => node             %massDistribution      (componentTypeHotHalo,massTypeGaseous)
       kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                    )
       ! Get the temperature.
       coordinates             =  [outerRadius,0.0d0,0.0d0]
       temperature             =  kinematicsDistribution_       %temperature(coordinates                                                                              )
       ! Compute density and cooling time at outer radius and zero radius.
       densityOuter            =  massDistribution_             %density    (coordinates                                                                              )
       coordinates             =  [0.0d0      ,0.0d0,0.0d0]
       densityZero             =  massDistribution_             %density    (coordinates                                                                              )
       coolingTimeZero         =  self             %coolingTime_%time       (node,temperature,densityZero ,hotAbundances,chemicalFractions*densityZero ,self%radiation)
       coolingTimeOuter        =  self             %coolingTime_%time       (node,temperature,densityOuter,hotAbundances,chemicalFractions*densityOuter,self%radiation)
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]          
       if (coolingTimeOuter < timeAvailable .or. coolingTimeZero > timeAvailable) then
          ! Cooling radius is static.
          self%radiusGrowthRateStored=0.0d0
       else
          ! Cooling radius is between zero and virial radii.
          self%radiusGrowthRateStored=+0.5d0                                                      &
               &                      *outerRadius                                                &
               &                      /coolingTimeZero                                            &
               &                      /sqrt(                                                      &
               &                            +(timeAvailable   /coolingTimeZero-1.0d0)             &
               &                            *(coolingTimeOuter/coolingTimeZero-1.0d0)             &
               &                           )                                                      &
               &                      *self%coolingTimeAvailable_%timeAvailableIncreaseRate(node)
       end if
       ! Flag that cooling radius is now computed.
       self%radiusGrowthRateComputed=.true.
    end if
    ! Return the stored value.
    betaProfileRadiusGrowthRate=self%radiusGrowthRateStored
    return
  end function betaProfileRadiusGrowthRate

  double precision function betaProfileRadius(self,node)
    !!{RST
    Return the cooling radius in the :math:`\beta`-profile model.
    !!}
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                   , nodeComponentHotHalo       , treeNode
    use :: Mass_Distributions               , only : massDistributionClass                , kinematicsDistributionClass
    use :: Coordinates                      , only : coordinateSpherical                  , assignment(=)
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo                 , massTypeGaseous
    implicit none
    class           (coolingRadiusBetaProfile   ), intent(inout), target  :: self
    type            (treeNode                   ), intent(inout), target  :: node
    class           (nodeComponentBasic         )               , pointer :: basic
    class           (nodeComponentHotHalo       )               , pointer :: hotHalo
    class           (massDistributionClass      ), pointer                :: massDistribution_
    class           (kinematicsDistributionClass), pointer                :: kinematicsDistribution_
    type            (coordinateSpherical        )                         :: coordinates
    double precision                                                      :: coolingTimeZero        , timeAvailable          , &
         &                                                                   densityZero            , massToDensityConversion, &
         &                                                                   temperature            , outerRadius            , &
         &                                                                   densityOuter           , coolingTimeOuter
    type            (abundances                 )                         :: hotAbundances
    type            (chemicalAbundances         )                         :: chemicalFractions      , chemicalMasses

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if cooling radius is already computed.
    if (.not.self%radiusComputed) then
       ! Get the time available for cooling in node.
       timeAvailable              =  self%coolingTimeAvailable_%timeAvailable(node)
       ! Get the abundances for this node.
       hotHalo                    => node   %hotHalo   ()
       hotAbundances              =  hotHalo%abundances()
       call hotAbundances%massToMassFraction(hotHalo%mass())
       ! Get the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalFractions)
          if (hotHalo%mass() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
          else
             massToDensityConversion=0.0d0
          end if          
          ! Convert to number density per unit total mass density.
          chemicalFractions=chemicalFractions*massToDensityConversion
       end if
       ! Set epoch for radiation field.
       basic => node%basic()
       call self%radiation%timeSet(basic%time())
       ! Get the outer radius.
       outerRadius=hotHalo%outerRadius()
       ! Get the mass distribution.
       massDistribution_       => node             %massDistribution      (componentTypeHotHalo,massTypeGaseous)
       kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                    )
       ! Get the temperature.
       coordinates             =  [outerRadius,0.0d0,0.0d0]
       temperature             =  kinematicsDistribution_       %temperature(coordinates                                                                              )
       ! Compute density and cooling time at outer radius and zero radius.
       densityOuter            =  massDistribution_             %density    (coordinates                                                                              )
       coordinates             =  [0.0d0      ,0.0d0,0.0d0]
       densityZero             =  massDistribution_             %density    (coordinates                                                                              )
       coolingTimeZero         =  self             %coolingTime_%time       (node,temperature,densityZero ,hotAbundances,chemicalFractions*densityZero ,self%radiation)
       coolingTimeOuter        =  self             %coolingTime_%time       (node,temperature,densityOuter,hotAbundances,chemicalFractions*densityOuter,self%radiation)
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]          
       if (coolingTimeOuter < timeAvailable) then
          ! Cooling time available exceeds cooling time at virial radius, return virial radius.
          self%radiusStored=outerRadius
       else if (coolingTimeZero > timeAvailable) then
          ! Gas at zero radius can not cool.
          self%radiusStored=0.0d0
       else
          ! Cooling radius is between zero and virial radii.
          self%radiusStored=+outerRadius                                     &
               &            *sqrt(                                          &
               &                  +(timeAvailable   /coolingTimeZero-1.0d0) &
               &                  /(coolingTimeOuter/coolingTimeZero-1.0d0) &
               &                 )
       end if
       ! Flag that cooling radius is now computed.
       self%radiusComputed=.true.
    end if
    betaProfileRadius=self%radiusStored
    return
  end function betaProfileRadius
