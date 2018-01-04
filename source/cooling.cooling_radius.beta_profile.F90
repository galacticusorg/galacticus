!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a cooling radius class for $\beta$-profile halos, assuming collisional ionization equilibrium such that cooling
  !% time scales as inverse density.
  
  use Kind_Numbers
  use Dark_Matter_Halo_Scales
  use Cooling_Times_Available
  use Cooling_Times

  !# <coolingRadius name="coolingRadiusBetaProfile" defaultThreadPrivate="yes">
  !#  <description>
  !#   A cooling radius class for $\beta$-profile halos. Computes the cooling radius by assuming that the hot gas density profile is a
  !#   $\beta$-profile ($\rho(r) \propto [r^2+r_{\rm c}^2]^{-1}$), and that the cooling rate scales as density squared, $\dot{E}\propto
  !#   \rho^2$, such that the cooling time scales as inverse density, $t_{\mathrm cool} \propto \rho^{-1}$. Consequently, the cooling radius is given by
  !#   \begin{equation}
  !#    r_{\rm cool} = r_{\rm virial} \left( \left[ {t_{\rm avail} \over t_0} - 1 \right] \left[ {t_{\rm virial} \over t_0} - 1 \right]^{-1} \right)^{1/2},
  !#   \end{equation}
  !#   where $t_0$, and $t_{\rm virial}$ are the cooling times at zero radius and the virial radius respectively.
  !#  </description>
  !# </coolingRadius>
  type, extends(coolingRadiusClass) :: coolingRadiusBetaProfile
     !% Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available
     !% for cooling equals the cooling time.
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_
     class           (coolingTimeAvailableClass), pointer :: coolingTimeAvailable_
     class           (coolingTimeClass         ), pointer :: coolingTime_
     integer         (kind=kind_int8           )          :: lastUniqueID                  =-1
     integer                                              :: abundancesCount                  , chemicalsCount
     ! Stored values of cooling radius.
     logical                                              :: radiusComputed                   , radiusGrowthRateComputed
     double precision                                     :: radiusGrowthRateStored           , radiusStored
   contains
     final     ::                     betaProfileDestructor
     procedure :: radius           => betaProfileRadius
     procedure :: radiusGrowthRate => betaProfileRadiusGrowthRate
     procedure :: calculationReset => betaProfileCalculationReset
  end type coolingRadiusBetaProfile

  interface coolingRadiusBetaProfile
     !% Constructors for the betaProfile cooling radius class.
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface coolingRadiusBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !% Constructor for the $\beta$-profile cooling radius class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type (coolingRadiusBetaProfile )                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(coolingTimeAvailableClass), pointer       :: coolingTimeAvailable_
    class(coolingTimeClass         ), pointer       :: coolingTime_
    class(darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_

    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !# <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_" source="parameters"/>
    !# <objectBuilder class="coolingTime"          name="coolingTime_"          source="parameters"/>
    self=coolingRadiusBetaProfile(darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_) result(self)
    !% Internal constructor for the $\beta$-profile cooling radius class.
    use ISO_Varying_String
    use Galacticus_Error
    use Array_Utilities
    use String_Handling
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Hot_Halo_Mass_Distributions
    use Hot_Halo_Temperature_Profiles    
    implicit none
    type (coolingRadiusBetaProfile      )                         :: self
    class(darkMatterHaloScaleClass      ), intent(in   ), target  :: darkMatterHaloScale_
    class(coolingTimeAvailableClass     ), intent(in   ), target  :: coolingTimeAvailable_
    class(coolingTimeClass              ), intent(in   ), target  :: coolingTime_
    class(hotHaloMassDistributionClass  )               , pointer :: hotHaloMassDistribution_
    class(hotHaloTemperatureProfileClass)               , pointer :: hotHaloTemperatureProfile_
    !# <constructorAssign variables="*darkMatterHaloScale_, *coolingTimeAvailable_, *coolingTime_"/>

    ! Initial state of stored solutions.
    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Check that required components are gettable.
    if     (                                                                                                                        &
         &  .not.(                                                                                                                  &
         &         defaultHotHaloComponent%       massIsGettable() .and.                                                            &
         &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                            &
         &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                            &
         &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  self%chemicalsCount == 0)                                  &
         &       )                                                                                                                  &
         & ) call Galacticus_Error_Report                                                                                           &
         & (                                                                                                                        &
         &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                                 &
         &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                                 &
         &  Galacticus_Component_List(                                                                                              &
         &                            'hotHalo'                                                                                  ,  &
         &                             defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=self%chemicalsCount > 0)   &
         &                           )                                                                                           // &
         &  {introspection:location}                                                                                                &
         & )
    ! Check that assumptions are valid. Note that currently we can not check all assumptions. We check that:
    !  * Hot halo temperature profile is isothermal;
    !  * Hot halo density profile is a β-profile.
    ! We do not check that:
    !  * β=2/3;
    !  * Cooling function is always proportional to ρ².
    hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
    hotHaloMassDistribution_   => hotHaloMassDistribution  ()
    select type (hotHaloTemperatureProfile_)
    class is (hotHaloTemperatureProfileVirial)
       ! An isothermal profile - this is acceptable.
    class default
       call Galacticus_Error_Report('assumption of isothermal hot halo temperature profile is not met'//{introspection:location})
    end select
    select type (hotHaloMassDistribution_)
    class is (hotHaloMassDistributionBetaProfile)
       ! An beta-model profile - this is acceptable.
    class default
       call Galacticus_Error_Report('assumption of β-model hot halo mass distribution is not met'     //{introspection:location})
    end select
    return
  end function betaProfileConstructorInternal
  
  subroutine betaProfileDestructor(self)
    !% Destructor for the $\beta$-profile cooling radius class.
    implicit none
    type(coolingRadiusBetaProfile), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%coolingTimeAvailable_"/> 
    !# <objectDestructor name="self%coolingTime_"         />
   return
  end subroutine betaProfileDestructor

  subroutine betaProfileCalculationReset(self,node)
    !% Reset the cooling radius calculation.
    implicit none
    class(coolingRadiusBetaProfile), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    self%lastUniqueID            =node%uniqueID()
    return
  end subroutine betaProfileCalculationReset

  double precision function betaProfileRadiusGrowthRate(self,node)
    !% Returns the cooling radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    use Radiation_Structure
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Hot_Halo_Mass_Distributions
    use Hot_Halo_Temperature_Profiles    
    implicit none
    class           (coolingRadiusBetaProfile      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout)          :: node
    class           (nodeComponentHotHalo          )               , pointer :: hotHalo
    class           (hotHaloMassDistributionClass  )               , pointer :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass)               , pointer :: hotHaloTemperatureProfile_
    double precision                                                         :: coolingTimeZero           , timeAvailable          , &
         &                                                                      densityZero               , massToDensityConversion, &
         &                                                                      temperature               , outerRadius            , &
         &                                                                      densityOuter              , coolingTimeOuter
    type            (abundances                    )                         :: hotAbundances
    type            (chemicalAbundances            )                         :: chemicalDensities         , chemicalMasses
    type            (radiationStructure            )                         :: radiation

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)

    ! Check if cooling radius growth rate is already computed.
    if (.not.self%radiusGrowthRateComputed) then
       hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
       hotHaloMassDistribution_   => hotHaloMassDistribution  ()
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
          call chemicalMasses%massToNumber(chemicalDensities)
          if (hotHalo%outerRadius() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
          else
             massToDensityConversion=0.0d0
          end if
          ! Convert to number density.
          chemicalDensities=chemicalDensities*massToDensityConversion
       end if
       ! Set the radiation field.
       call radiation%define([radiationTypeCMB])
       call radiation%set(node)
       ! Get the outer radius.
       outerRadius=hotHalo%outerRadius()
       ! Get the temperature.
       temperature=hotHaloTemperatureProfile_%temperature(node,outerRadius)
       ! Compute density and cooling time at outer radius and zero radius.
       densityZero     =     hotHaloMassDistribution_  %density(node,0.0d0      )
       densityOuter    =     hotHaloMassDistribution_  %density(node,outerRadius)
       coolingTimeZero =self%coolingTime_              %time   (temperature,densityZero ,hotAbundances,chemicalDensities,radiation)
       coolingTimeOuter=self%coolingTime_              %time   (temperature,densityOuter,hotAbundances,chemicalDensities,radiation)
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
    !% Return the cooling radius in the $\beta$-profile model.
    use Radiation_Structure
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Hot_Halo_Mass_Distributions
    use Hot_Halo_Temperature_Profiles    
    implicit none
    class           (coolingRadiusBetaProfile      ), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    class           (nodeComponentHotHalo          )               , pointer :: hotHalo
    class           (hotHaloMassDistributionClass  )               , pointer :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass)               , pointer :: hotHaloTemperatureProfile_
    double precision                                                         :: coolingTimeZero           , timeAvailable          , &
         &                                                                      densityZero               , massToDensityConversion, &
         &                                                                      temperature               , outerRadius            , &
         &                                                                      densityOuter              , coolingTimeOuter
    type            (abundances                    )                         :: hotAbundances
    type            (chemicalAbundances            )                         :: chemicalDensities         , chemicalMasses
    type            (radiationStructure            )                         :: radiation

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if cooling radius is already computed.
    if (.not.self%radiusComputed) then
       ! Get required objects.
       hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
       hotHaloMassDistribution_   => hotHaloMassDistribution  ()
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
          call chemicalMasses%massToNumber(chemicalDensities)
          if (hotHalo%outerRadius() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
          else
             massToDensityConversion=0.0d0
          end if
          ! Convert to number density.
          chemicalDensities=chemicalDensities*massToDensityConversion
       end if
       ! Set the radiation field.
       call radiation%define([radiationTypeCMB])
       call radiation%set(node)
       ! Get the outer radius.
       outerRadius=hotHalo%outerRadius()
       ! Get the temperature.
       temperature=hotHaloTemperatureProfile_%temperature(node,outerRadius)
       ! Compute density and cooling time at outer radius and zero radius.
       densityZero     =     hotHaloMassDistribution_  %density(node       ,0.0d0                                                 )
       densityOuter    =     hotHaloMassDistribution_  %density(node       ,outerRadius                                           )
       coolingTimeZero =self%coolingTime_              %time   (temperature,densityZero ,hotAbundances,chemicalDensities,radiation)
       coolingTimeOuter=self%coolingTime_              %time   (temperature,densityOuter,hotAbundances,chemicalDensities,radiation)
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
