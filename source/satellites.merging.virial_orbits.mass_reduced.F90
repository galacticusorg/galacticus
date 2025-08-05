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
  An implementation of virial orbits that modifies another \refClass{virialOrbitClass} by accounting for any reduction in mass of
  the primary halo below its ``\gls{dmou}'' value.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <virialOrbit name="virialOrbitMassReduced">
   <description>
    A virial orbit class that modifies another \refClass{virialOrbitClass} by accounting for any reduction in mass of the primary
    halo below its ``\gls{dmou}'' value.
   </description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitMassReduced
     !!{
     A virial orbit class that modifies another \refClass{virialOrbitClass} by accounting for any reduction in mass of the primary
     halo below its ``\gls{dmou}'' value.
     !!}
     private
     class(virialOrbitClass          ), pointer :: virialOrbit_           => null()
     class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
     class(cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class(cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class(darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
   contains
     final     ::                                 massReducedDestructor
     procedure :: orbit                        => massReducedOrbit
     procedure :: densityContrastDefinition    => massReducedDensityContrastDefinition
     procedure :: velocityTotalRootMeanSquared => massReducedVelocityTotalRootMeanSquared
     procedure :: energyMean                   => massReducedEnergyMean
  end type virialOrbitMassReduced

  interface virialOrbitMassReduced
     !!{
     Constructors for the \refClass{virialOrbitMassReduced} virial orbit class.
     !!}
     module procedure massReducedConstructorParameters
     module procedure massReducedConstructorInternal
  end interface virialOrbitMassReduced

contains

  function massReducedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitMassReduced} satellite virial orbit class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialOrbitMassReduced    )                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(virialOrbitClass          ), pointer       :: virialOrbit_
    class(cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class(virialDensityContrastClass), pointer       :: virialDensityContrast_
    class(darkMatterProfileDMOClass ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="virialOrbit"           name="virialOrbit_"           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=virialOrbitMassReduced(virialOrbit_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"          />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function massReducedConstructorParameters

  function massReducedConstructorInternal(virialOrbit_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitMassReduced} virial orbits class.
    !!}
    implicit none
    type (virialOrbitMassReduced    )                        :: self
    class(virialOrbitClass          ), intent(in   ), target :: virialOrbit_
    class(darkMatterProfileDMOClass ), intent(in   ), target :: darkMatterProfileDMO_
    class(cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    class(virialDensityContrastClass), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*virialOrbit_, *cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function massReducedConstructorInternal

  subroutine massReducedDestructor(self)
    !!{
    Destructor for the \refClass{virialOrbitMassReduced} virial orbits class.
    !!}
    implicit none
    type(virialOrbitMassReduced), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"          />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%virialDensityContrast_"/>
    <objectDestructor name="self%cosmologyParameters_"  />
    !!]
    return
  end subroutine massReducedDestructor

  function massReducedOrbit(self,node,host,acceptUnboundOrbits) result(orbit)
    !!{
    Return massReduced orbital parameters for a satellite.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Kepler_Orbits                       , only : keplerOrbitPhi                     , keplerOrbitRadius, keplerOrbitTheta, keplerOrbitVelocityTangential
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    use :: Galactic_Structure_Options          , only : componentTypeAll                   , massTypeAll
    use :: Virial_Density_Contrast             , only : virialDensityContrastClass
    implicit none
    type            (keplerOrbit               )                        :: orbit
    class           (virialOrbitMassReduced    ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                           , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    integer                                     , parameter             :: iterationMaximum          =1000
    class           (nodeComponentBasic        ), pointer               :: basicHost
    class           (virialDensityContrastClass), pointer               :: densityContrastDefinition_
    class           (massDistributionClass     ), pointer               :: massDistribution_
    double precision                                                    :: velocityHost                   , massHost             , &
         &                                                                 radiusHost                     , massHostDMO          , &
         &                                                                 massSatellite                  , velocityRadialSquared
    integer                                                             :: iteration
    logical                                                             :: acceptOrbit
    
    ! Find dark matter-only mass, radius, and velocity in the host corresponding to the virial density contrast definition.
    !![
    <referenceAcquire target="densityContrastDefinition_" source="self%virialOrbit_%densityContrastDefinition()"/>
    !!]
    basicHost   => host%basic()
    massHostDMO =  Dark_Matter_Profile_Mass_Definition(                                                                                                                  &
         &                                                                    host                                                                                     , &
         &                                                                    densityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                    radiusHost                                                                               , &
         &                                                                    velocityHost                                                                             , &
         &                                             cosmologyParameters_  =self%cosmologyParameters_                                                                , &
         &                                             cosmologyFunctions_   =self%cosmologyFunctions_                                                                 , &
         &                                             virialDensityContrast_=self%virialDensityContrast_                                                              , &
         &                                             darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                 &
         &                                            )
    !![
    <objectDestructor name="densityContrastDefinition_"/>
    !!]
    ! Find the mass of the host including baryons.
    massDistribution_ => host             %massDistribution    (massType=massTypeAll,componentType=componentTypeAll)
    massHost          =  massDistribution_%massEnclosedBySphere(                                   radiusHost      )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Iterate until an acceptable orbit is found.
    acceptOrbit=.false.
    iteration  =0
    do while (.not.acceptOrbit .and. iteration < iterationMaximum)
       acceptOrbit=.true.
       iteration  =iteration+1
       ! Get the dark matter-only orbit.
       orbit=self%virialOrbit_%orbit(node,host,acceptUnboundOrbits)
       ! Keeping the angular momentum unchanged, adjust the energy of the orbit.
       velocityRadialSquared=+orbit%velocityRadial()**2      &
            &                +2.0d0                          &
            &                *gravitationalConstant_internal &
            &                *(                              &
            &                  +massHost                     &
            &                  -massHostDMO                  &
            &                 )                              &
            &                /radiusHost
       if (velocityRadialSquared < 0.0d0) then
          acceptOrbit=.false.
          cycle
       end if
       massSatellite=(1.0d0/orbit%specificReducedMass()-1.0d0)*massHost
       call orbit%reset            (                                     &
            &                       keep=[                               &
            &                             keplerOrbitRadius            , &
            &                             keplerOrbitVelocityTangential, &
            &                             keplerOrbitTheta             , &
            &                             keplerOrbitPhi                 &
            &                            ]                               &
            &                      )
       call orbit%massesSet        (                                     &
            &                             massSatellite                , &
            &                             massHost                       &
            &                      )
       call orbit%velocityRadialSet(                                     &
            &                             -sqrt(+velocityRadialSquared)  &
            &                      )
       if (orbit%energy() >= 0.0d0 .and. .not.acceptUnboundOrbits) &
            & acceptOrbit=.false.
    end do
    if (.not.acceptOrbit) call Error_Report('no acceptable orbit found'//{introspection:location})
    return
  end function massReducedOrbit

  function massReducedDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of massReduced virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: massReducedDensityContrastDefinition
    class(virialOrbitMassReduced    ), intent(inout) :: self

    massReducedDensityContrastDefinition => self%virialOrbit_%densityContrastDefinition()
    return
  end function massReducedDensityContrastDefinition

  double precision function massReducedVelocityTotalRootMeanSquared(self,node,host) result (velocityRootMeanSquared)
    !!{
    Return the root mean squared of the total velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    use :: Galactic_Structure_Options          , only : componentTypeAll                   , massTypeAll
    use :: Virial_Density_Contrast             , only : virialDensityContrastClass
    implicit none
    class           (virialOrbitMassReduced    ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                      , host
    class           (nodeComponentBasic        ), pointer       :: basicHost
    class           (virialDensityContrastClass), pointer       :: densityContrastDefinition_
    class           (massDistributionClass     ), pointer       :: massDistribution_
    double precision                                            :: massHost                  , radiusHost , &
         &                                                         velocityHost              , massHostDMO
    !$GLC attributes unused :: node

    !![
    <referenceAcquire target="densityContrastDefinition_" source="self%virialOrbit_%densityContrastDefinition()"/>
    !!]
    basicHost               =>  host%basic()
    massHostDMO             =  +Dark_Matter_Profile_Mass_Definition(                                                                                                                  &
         &                                                                                 host                                                                                     , &
         &                                                                                 densityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                                 radiusHost                                                                               , &
         &                                                                                 velocityHost                                                                             , &
         &                                                          cosmologyParameters_  =self%cosmologyParameters_                                                                , &
         &                                                          cosmologyFunctions_   =self%cosmologyFunctions_                                                                 , &
         &                                                          virialDensityContrast_=self%virialDensityContrast_                                                              , &
         &                                                          darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                 &
         &                                                        )
    !![
    <objectDestructor name="densityContrastDefinition_"/>
    !!]
    massDistribution_ => host             %massDistribution    (massType=massTypeAll,componentType=componentTypeAll)
    massHost          =  massDistribution_%massEnclosedBySphere(                                   radiusHost      )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    velocityRootMeanSquared =  +sqrt(                                                              &
         &                           +self%virialOrbit_%velocityTotalRootMeanSquared(node,host)**2 &
         &                           +2.0d0                                                        &
         &                           *gravitationalConstant_internal                               &
         &                           *(                                                            &
         &                             +massHost                                                   &
         &                             -massHostDMO                                                &
         &                            )                                                            &
         &                           /radiusHost                                                   &
         &                          )
    return
  end function massReducedVelocityTotalRootMeanSquared

  double precision function massReducedEnergyMean(self,node,host) result(energyMean)
    !!{
    Return the mean energy of the orbits.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    use :: Galactic_Structure_Options          , only : componentTypeAll                   , massTypeAll
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Virial_Density_Contrast             , only : virialDensityContrastClass
    implicit none
    class           (virialOrbitMassReduced    ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                      , host
    class           (nodeComponentBasic        ), pointer       :: basicHost
    class           (virialDensityContrastClass), pointer       :: densityContrastDefinition_
    class           (massDistributionClass     ), pointer               :: massDistribution_
    double precision                                            :: massHost                  , radiusHost , &
         &                                                         velocityHost              , massHostDMO

    !![
    <referenceAcquire target="densityContrastDefinition_" source="self%virialOrbit_%densityContrastDefinition()"/>
    !!]
    basicHost   =>  host%basic()
    massHostDMO =  +Dark_Matter_Profile_Mass_Definition(                                                                                                                  &
         &                                                                     host                                                                                     , &
         &                                                                     densityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                     radiusHost                                                                               , &
         &                                                                     velocityHost                                                                             , &
         &                                              cosmologyParameters_  =self%cosmologyParameters_                                                                , &
         &                                              cosmologyFunctions_   =self%cosmologyFunctions_                                                                 , &
         &                                              virialDensityContrast_=self%virialDensityContrast_                                                              , &
         &                                              darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                 &
         &                                             )
    !![
    <objectDestructor name="densityContrastDefinition_"/>
    !!]
    massDistribution_ => host             %massDistribution    (massType=massTypeAll,componentType=componentTypeAll)
    massHost          =  massDistribution_%massEnclosedBySphere(                                   radiusHost      )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    energyMean=+self%virialOrbit_%energyMean(node,host) &
         &     +gravitationalConstant_internal          &
         &     *(                                       &
         &       +massHost                              &
         &       -massHostDMO                           &
         &      )                                       &
         &     /radiusHost
    return
  end function massReducedEnergyMean
