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
  An implementation of virial orbits using the \cite{benson_orbital_2005} orbital parameter distribution.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass, virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt

  !![
  <virialOrbit name="virialOrbitBenson2005">
   <description>
    A virial orbits class which selects orbital parameters randomly from the distribution given by
    \cite{benson_orbital_2005}. If the virial density contrast definition differs from that used by \cite{benson_orbital_2005}
    then the orbit is assigned based on \cite{benson_orbital_2005}'s definition and then propagated to the virial radius
    relevant to the current definition of density contrast.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </stateStorable>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitBenson2005
     !!{
     A virial orbit class using the \cite{benson_orbital_2005} orbital parameter distribution.
     !!}
     private
     class(darkMatterHaloScaleClass                                      ), pointer :: darkMatterHaloScale_             => null()
     class(virialDensityContrastClass                                    ), pointer :: virialDensityContrast_           => null()
     class(cosmologyParametersClass                                      ), pointer :: cosmologyParameters_             => null()
     class(cosmologyFunctionsClass                                       ), pointer :: cosmologyFunctions_              => null()
     class(darkMatterProfileDMOClass                                     ), pointer :: darkMatterProfileDMO_            => null()
     type (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: virialDensityContrastDefinition_ => null()
   contains
     final     ::                                    benson2005Destructor
     procedure :: orbit                           => benson2005Orbit
     procedure :: densityContrastDefinition       => benson2005DensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => benson2005VelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => benson2005VelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => benson2005AngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => benson2005AngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => benson2005VelocityTotalRootMeanSquared
     procedure :: energyMean                      => benson2005EnergyMean
  end type virialOrbitBenson2005

  interface virialOrbitBenson2005
     !!{
     Constructors for the \refClass{virialOrbitBenson2005} virial orbit class.
     !!}
     module procedure benson2005ConstructorParameters
     module procedure benson2005ConstructorInternal
  end interface virialOrbitBenson2005

contains

  function benson2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitBenson2005} virial orbits class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialOrbitBenson2005     )                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class(darkMatterProfileDMOClass ), pointer       :: darkMatterProfileDMO_
    class(virialDensityContrastClass), pointer       :: virialDensityContrast_
    
    !![
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=virialOrbitBenson2005(darkMatterHaloScale_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function benson2005ConstructorParameters

  function benson2005ConstructorInternal(darkMatterHaloScale_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitBenson2005} virial orbits class.
    !!}
    use :: Virial_Density_Contrast, only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
    implicit none
    type (virialOrbitBenson2005     )                        :: self
    class(darkMatterHaloScaleClass  ), intent(in   ), target :: darkMatterHaloScale_
    class(cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    class(virialDensityContrastClass), intent(in   ), target :: virialDensityContrast_
    class(darkMatterProfileDMOClass ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, *darkMatterProfileDMO_"/>
    !!]

    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="virialDensityContrastDefinition_" constructor="virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(.true.,cosmologyFunctions_)"/>
    !!]
    return
  end function benson2005ConstructorInternal

  subroutine benson2005Destructor(self)
    !!{
    Destructor for the \refClass{virialOrbitBenson2005} virial orbits class.
    !!}
    implicit none
    type(virialOrbitBenson2005), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine benson2005Destructor

  function benson2005Orbit(self,node,host,acceptUnboundOrbits)
    !!{
    Return benson2005 orbital parameters for a satellite.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    type            (keplerOrbit          )                        :: benson2005Orbit
    class           (virialOrbitBenson2005), intent(inout), target :: self
    type            (treeNode             ), intent(inout)         :: host                   , node
    logical                                , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic   ), pointer               :: basicHost              , basic
    double precision                       , parameter             :: pMax                   =1.96797d0                              , &
         &                                                            velocityMax            =3.00000d0
    double precision                       , parameter             :: a                   (9)=[                                        &
         &                                                                                     0.390052d+01,0.247973d+01,0.102373d+02, &
         &                                                                                     0.683922d+00,0.353953d+00,0.107716d+01, &
         &                                                                                     0.509837d+00,0.206204d+00,0.314641d+00  &
         &                                                                                    ]
    double precision                       , parameter             :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                                               :: b1                             , b2                        , &
         &                                                            distributionFunction           , energyInternal            , &
         &                                                            uniformRandom                  , velocityRadialInternal    , &
         &                                                            velocityHost                   , velocityTangentialInternal, &
         &                                                            massHost                       , radiusHost                , &
         &                                                            massSatellite                  , radiusHostSelf
    logical                                                        :: foundOrbit

    ! Get basic components.
    basic     => node%basic()
    basicHost => host%basic()
    ! Find mass, radius, and velocity in the host corresponding to the Benson (2005) virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHostSelf                                                                                      , &
         &                                                                   velocityHost                                                                                        , &
         &                                            cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                            cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                            virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                            darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                           )
    massSatellite=Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   node                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(    basic%mass(),    basic%timeLastIsolated()), &
         &                                            cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                            cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                            virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                            darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                           )
    ! Handle zero masses gracefully - simply set to the basic mass as no physically-meaningful result should depend on what we do here.
    if (massHost       <= 0.0d0) massHost      =basicHost                     %mass        (    )
    if (massSatellite  <= 0.0d0) massSatellite =basic                         %mass        (    )
    if (radiusHostSelf <= 0.0d0) radiusHostSelf=self     %darkMatterHaloScale_%radiusVirial(host)
    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Reset the orbit.
       call benson2005Orbit%reset()
       ! Set basic properties of the orbit.
       call benson2005Orbit%massesSet(massSatellite,massHost      )
       call benson2005Orbit%radiusSet(              radiusHostSelf)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =node%hostTree%randomNumberGenerator_%uniformSample()*velocityMax
       velocityTangentialInternal=node%hostTree%randomNumberGenerator_%uniformSample()*velocityMax
       ! Evaluate distribution function for these parameters.
       b1                  =+a(3)                                                &
            &               *exp(-a (4)*( velocityTangentialInternal-a (5))**2)
       b2                  =+a(6)                                                &
            &               *exp(-a (7)*( velocityTangentialInternal-a (8))**2)
       distributionFunction=+a(1)                                                &
            &               *velocityTangentialInternal                          &
            &               *exp(-a (2)*((velocityTangentialInternal-a (9))**2)) &
            &               *exp(-b1   *( velocityRadialInternal    -b2   )**2)
       if (distributionFunction > pMax) call Error_Report('distribution function exceeds expected peak value'//{introspection:location})
       uniformRandom=pMax*node%hostTree%randomNumberGenerator_%uniformSample()
       if (uniformRandom <= distributionFunction) then
          foundOrbit=.true.
          ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
          ! bound that later rounding errors will not make it appear unbound.
          if (.not.acceptUnboundOrbits) then
             energyInternal=-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2)*benson2005Orbit%specificReducedMass()
             foundOrbit=(energyInternal < -boundTolerance)
          end if
       end if
       if (.not.foundOrbit) cycle
       call benson2005Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call benson2005Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%radiusVirial(host)
       if (benson2005Orbit%radiusApocenter() >= radiusHost .and. benson2005Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call benson2005Orbit%propagate(radiusHost  ,infalling=.true.)
          call benson2005Orbit%massesSet(basic%mass(),basicHost%mass())
       end if
    end do
    return
  end function benson2005Orbit

  function benson2005DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of \cite{benson_orbital_2005} virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: benson2005DensityContrastDefinition
    class(virialOrbitBenson2005     ), intent(inout) :: self

    benson2005DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function benson2005DensityContrastDefinition

  double precision function benson2005VelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                             , host
    class           (nodeComponentBasic   ), pointer       :: basicHost
    ! The mean magnitude of tangential velocity. This was by numerical integration over the velocity distribution fitting function.
    double precision                       , parameter     :: velocityTangentialMean=0.748205d0
    double precision                                       :: massHost                         , radiusHost, &
         &                                                    velocityHost
    !$GLC attributes unused :: node

    basicHost                                 =>  host%basic()
    massHost                                  =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    benson2005VelocityTangentialMagnitudeMean =  +velocityTangentialMean &
         &                                       *velocityHost
    return
  end function benson2005VelocityTangentialMagnitudeMean

  function benson2005VelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                       , dimension(3)  :: benson2005VelocityTangentialVectorMean
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                                  , host
    !$GLC attributes unused :: self, node, host

    benson2005VelocityTangentialVectorMean=0.0d0
    call Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function benson2005VelocityTangentialVectorMean

  double precision function benson2005AngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , basicHost
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                                  =>  node%basic()
    basicHost                              =>  host%basic()
    massHost                               =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    benson2005AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                    *radiusHost                                      &
         &                                    /(                                               & ! Account for reduced mass.
         &                                      +1.0d0                                         &
         &                                      +basic    %mass()                              &
         &                                      /basicHost%mass()                              &
         &                                     )
    return
  end function benson2005AngularMomentumMagnitudeMean

  function benson2005AngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector angular momentum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                       , dimension(3)  :: benson2005AngularMomentumVectorMean
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                               , host
    !$GLC attributes unused :: self, node, host

    benson2005AngularMomentumVectorMean=0.0d0
    call Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function benson2005AngularMomentumVectorMean

  double precision function benson2005VelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                                  , host
    class           (nodeComponentBasic   ), pointer       :: basicHost
    ! The root mean squared total velocity. This was by numerical integration over the velocity distribution fitting function.
    double precision                       , parameter     :: velocityTotalRootMeanSquared=1.25534d0
    double precision                                       :: massHost                              , radiusHost, &
         &                                                    velocityHost
    !$GLC attributes unused :: node

    basicHost                              =>  host%basic()
    massHost                               =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    benson2005VelocityTotalRootMeanSquared =  +velocityTotalRootMeanSquared &
         &                                    *velocityHost
    return
  end function benson2005VelocityTotalRootMeanSquared

  double precision function benson2005EnergyMean(self,node,host)
    !!{
    Return the mean energy of the orbits.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    implicit none
    class           (virialOrbitBenson2005), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , basicHost
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                =>  node%basic()
    basicHost            =>  host%basic()
    massHost             =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    benson2005EnergyMean =  +0.5d0                                           &
         &                  *self%velocityTotalRootMeanSquared(node,host)**2 &
         &                  /(                                               & ! Account for reduced mass.
         &                    +1.0d0                                         &
         &                    +basic    %mass()                              &
         &                    /basicHost%mass()                              &
         &                   )                                               &
         &                  -gravitationalConstant_internal                  &
         &                  *massHost                                        &
         &                  /radiusHost
    return
  end function benson2005EnergyMean
