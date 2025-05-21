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
  An implementation of virial orbits using the \cite{li_orbital_2020} orbital parameter distribution.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass       , criticalOverdensityClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastBryanNorman1998, virialDensityContrastClass

  !![
  <virialOrbit name="virialOrbitLi2020">
   <description>Virial orbits using the \cite{li_orbital_2020} orbital parameter distribution.</description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </stateStorable>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitLi2020
     !!{
     A virial orbit class using the \cite{li_orbital_2020} orbital parameter distribution.
     !!}
     private
     class           (darkMatterHaloScaleClass            ), pointer :: darkMatterHaloScale_             => null()
     class           (cosmologyParametersClass            ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass             ), pointer :: cosmologyFunctions_              => null()
     class           (criticalOverdensityClass            ), pointer :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass       ), pointer :: cosmologicalMassVariance_        => null()
     class           (darkMatterProfileDMOClass           ), pointer :: darkMatterProfileDMO_            => null()
     class           (virialDensityContrastClass          ), pointer :: virialDensityContrast_           => null()
     type            (virialDensityContrastBryanNorman1998), pointer :: virialDensityContrastDefinition_ => null()
     double precision                                                :: mu1                                       , mu2   , &
          &                                                             a0                                        , a1    , &
          &                                                             a2                                        , a3    , &
          &                                                             b1                                        , b2    , &
          &                                                             c                                         , sigma1
     logical                                                         :: propagateOrbits
   contains
     !![
     <methods>
       <method description="Evaluate the $\eta$ parameter of the \cite{li_orbital_2020} virial orbit distribution function." method="eta"/>
     </methods>
     !!]
     final     ::                                    li2020Destructor
     procedure :: orbit                           => li2020Orbit
     procedure :: velocityDistributionFunction    => li2020VelocityDistributionFunction
     procedure :: densityContrastDefinition       => li2020DensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => li2020VelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => li2020VelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => li2020AngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => li2020AngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => li2020VelocityTotalRootMeanSquared
     procedure :: energyMean                      => li2020EnergyMean
     procedure :: eta                             => li2020Eta
  end type virialOrbitLi2020

  interface virialOrbitLi2020
     !!{
     Constructors for the \refClass{virialOrbitLi2020} virial orbit class.
     !!}
     module procedure li2020ConstructorParameters
     module procedure li2020ConstructorInternal
  end interface virialOrbitLi2020

contains

  function li2020ConstructorParameters(parameters) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitLi2020} virial orbits class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitLi2020            )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class           (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    class           (virialDensityContrastClass   ), pointer       :: virialDensityContrast_
    double precision                                               :: mu1                      , mu2   , &
         &                                                            a0                       , a1    , &
         &                                                            a2                       , a3    , &
         &                                                            b1                       , b2    , &
         &                                                            c                        , sigma1
    logical                                                        :: propagateOrbits
    
    !![
    <inputParameter>
      <name>mu1</name>
      <defaultValue>1.20d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $\mu_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>mu2</name>
      <defaultValue>1.04d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $\mu_2$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>sigma1</name>
      <defaultValue>0.20d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $\sigma_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>a0</name>
      <defaultValue>0.89d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $a_0$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>a1</name>
      <defaultValue>0.30d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $a_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>a2</name>
      <defaultValue>-3.33d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $a_2$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>a3</name>
      <defaultValue>0.56d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $a_3$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>b1</name>
      <defaultValue>-1.44d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $b_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>b2</name>
      <defaultValue>9.60d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $b_2$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>c</name>
      <defaultValue>0.43d0</defaultValue>
      <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
      <source>parameters</source>
      <description>Values of the $c$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>propagateOrbits</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If true, orbits will be propagated to the virial radius of the host halo. Otherwise, no propagation is performed and orbital velocities will correspond to precisely those from the \cite{li_orbital_2020} orbital velocity distribution.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
     !!]
    self=virialOrbitLi2020(mu1,mu2,sigma1,a0,a1,a2,a3,b1,b2,c,propagateOrbits,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function li2020ConstructorParameters

  function li2020ConstructorInternal(mu1,mu2,sigma1,a0,a1,a2,a3,b1,b2,c,propagateOrbits,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitLi2020} virial orbits class.
    !!}
    implicit none
    type            (virialOrbitLi2020            )                             :: self
    class           (darkMatterHaloScaleClass     ), intent(in   ), target      :: darkMatterHaloScale_
    class           (cosmologyParametersClass     ), intent(in   ), target      :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target      :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), intent(in   ), target      :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), intent(in   ), target      :: cosmologicalMassVariance_
    class           (virialDensityContrastClass   ), intent(in   ), target      :: virialDensityContrast_
    class           (darkMatterProfileDMOClass    ), intent(in   ), target      :: darkMatterProfileDMO_
    double precision                               , intent(in   )              :: mu1                      , mu2   , &
         &                                                                         a0                       , a1    , &
         &                                                                         a2                       , a3    , &
         &                                                                         b1                       , b2    , &
         &                                                                         c                        , sigma1
    logical                                        , intent(in   )              :: propagateOrbits
    !![
    <constructorAssign variables="mu1, mu2, sigma1, a0, a1, a2, a3, b1, b2, c, propagateOrbits, *darkMatterHaloScale_, *cosmologyParameters_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *virialDensityContrast_, *darkMatterProfileDMO_"/>
    !!]
    
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="virialDensityContrastDefinition_">
      <constructor>
        virialDensityContrastBryanNorman1998(                                                     &amp;
         &amp;                               allowUnsupportedCosmology=     .true.              , &amp;
         &amp;                               cosmologyParameters_     =self%cosmologyParameters_, &amp;
         &amp;                               cosmologyFunctions_      =self%cosmologyFunctions_   &amp;
         &amp;                              )
      </constructor>
    </referenceConstruct>
    !!]
    return
  end function li2020ConstructorInternal

  subroutine li2020Destructor(self)
    !!{
    Destructor for the \refClass{virialOrbitLi2020} virial orbits class.
    !!}
    implicit none
    type(virialOrbitLi2020), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"            />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine li2020Destructor

  function li2020Orbit(self,node,host,acceptUnboundOrbits)
    !!{
    Return li2020 orbital parameters for a satellite.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    type            (keplerOrbit               )                        :: li2020Orbit
    class           (virialOrbitLi2020         ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                                   , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: basicHost                              , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrastDefinition_
    integer                                     , parameter             :: attemptsMaximum                 =10000
    double precision                            , parameter             :: boundTolerance                  =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                                                    :: velocityHost                           , radiusHost                , &
         &                                                                 massHost                               , massSatellite             , &
         &                                                                 energyInternal                         , radiusHostSelf            , &
         &                                                                 velocityRadialInternal                 , velocityTangentialInternal, &
         &                                                                 velocityTotalInternal                  , eta
    logical                                                             :: foundOrbit
    integer                                                             :: attempts

    ! Get basic components.
    basic     => node%basic()
    basicHost => host%basic()
    ! Find virial density contrast under Li et al. (2020) definition.
    !![
    <referenceAcquire target="virialDensityContrastDefinition_" source="self%densityContrastDefinition()"/>    
    !!]
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Li et al. (2020) virial density contrast
    ! definition.
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
    !![
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    ! Handle cases of zero mass.
    if (massHost <= 0.0d0 .or. massSatellite <= 0.0d0) then
       massHost      =basicHost                     %mass          (    )
       massSatellite =basic                         %mass          (    )
       radiusHostSelf=self     %darkMatterHaloScale_%radiusVirial  (host)
       velocityHost  =self     %darkMatterHaloScale_%velocityVirial(host)
    end if
    ! Select an orbit.
    foundOrbit=.false.
    attempts  =0
    do while (.not.foundOrbit .and. attempts < attemptsMaximum)
       ! Increment number of attempts.
       attempts=attempts+1
       ! Reset the orbit.
       call li2020Orbit%reset()
       ! Set basic properties of the orbit.
       call li2020Orbit%massesSet(min(massSatellite,massHost),massHost      )
       call li2020Orbit%radiusSet(                            radiusHostSelf)
       ! Select total orbital velocity from a log-normal distribution.
       velocityTotalInternal=exp(self%sigma1*node%hostTree%randomNumberGenerator_%standardNormalSample())*self%mu1
       ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
       ! bound that later rounding errors will not make it appear unbound.
       foundOrbit=.true.
       if (.not.acceptUnboundOrbits) then
          energyInternal=-1.0d0+0.5d0*velocityTotalInternal**2*li2020Orbit%specificReducedMass()
          foundOrbit=(energyInternal < -boundTolerance)
       end if
       if (.not.foundOrbit) cycle
       ! Select the radial velocity.
       eta=self%eta(host,massSatellite,massHost,velocityTotalInternal)
       if (eta > 0.0d0) then
          ! Anisotropic case - exponential distribution.
          velocityRadialInternal=-sqrt(                                                                                      &
               &                       +log((exp(eta)-1.0d0)*node%hostTree%randomNumberGenerator_%uniformSample()+1.0d0)/eta &
               &                      )                                                                                      &
               &                 *velocityTotalInternal
       else
          ! Isotropic case - distribution is uniform in vᵣ².
          velocityRadialInternal=-sqrt(                                                                                      &
               &                                             node%hostTree%randomNumberGenerator_%uniformSample()            &
               &                      )                                                                                      &
               &                 *velocityTotalInternal
       end if
       ! Compute tangential velocity.
       velocityTangentialInternal=sqrt(max(0.0d0,velocityTotalInternal**2-velocityRadialInternal**2))
       call li2020Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call li2020Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       if (self%propagateOrbits) then
          radiusHost=self%darkMatterHaloScale_%radiusVirial(host)
          foundOrbit=.false.
          if (li2020Orbit%radiusApocenter() >= radiusHost .and. li2020Orbit%radiusPericenter() <= radiusHost) then
             foundOrbit=.true.
             call li2020Orbit%propagate(radiusHost  ,infalling=.true.)
             call li2020Orbit%massesSet(basic%mass(),basicHost%mass())
          end if
       end if
    end do
    ! If too many iterations were required to find an orbit, abort.
    if (attempts >= attemptsMaximum) then
       call node%serializeASCII()
       call host%serializeASCII()
       call Error_Report('maximum number of attempts exceeded'//{introspection:location})
    end if
    return
  end function li2020Orbit

  function li2020DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of \cite{li_orbital_2020} virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: li2020DensityContrastDefinition
    class(virialOrbitLi2020         ), intent(inout) :: self

    li2020DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function li2020DensityContrastDefinition

  function li2020VelocityDistributionFunction(self,node,host,velocityRadial,velocityTangential) result(distributionFunction)
    !!{
    Return the orbital velocity distribution function.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Numerical_Constants_Math            , only : Pi
    implicit none
    double precision                                             :: distributionFunction
    class           (virialOrbitLi2020         ), intent(inout)  :: self
    type            (treeNode                  ), intent(inout)  :: host                  , node
    double precision                            , intent(in   )  :: velocityRadial        , velocityTangential
    class           (nodeComponentBasic        ), pointer        :: basic                 , basicHost
    class           (virialDensityContrastClass), pointer        :: virialDensityContrastDefinition_
    double precision                                             :: velocityRadialInternal, velocityTangentialInternal, &
         &                                                          velocityTotalInternal , cosSquaredTheta           , &
         &                                                          massSatellite         , massHost                  , &
         &                                                          radiusHost            , velocityHost              , &
         &                                                          eta

    ! Get basic components.
    basic     => node%basic()
    basicHost => host%basic()
    ! Find virial density contrast under Li et al. (2020) definition.
    !![
    <referenceAcquire target="virialDensityContrastDefinition_" source="self%densityContrastDefinition()"/>    
    !!]
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Li et al. (2020) virial density contrast
    ! definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHost                                                                                          , &
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
    !![
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    ! Compute the total velocity and cos²θ.
    velocityRadialInternal    =velocityRadial    /velocityHost
    velocityTangentialInternal=velocityTangential/velocityHost
    velocityTotalInternal     =sqrt(                               &
         &                          +velocityRadialInternal    **2 &
         &                          +velocityTangentialInternal**2 &
         &                         )
    cosSquaredTheta           =+   (                        &
         &                          +velocityRadialInternal &
         &                          /velocityTotalInternal  &
         &                         )**2
    ! Compute the η parameter.
    eta=self%eta(host,massSatellite,massHost,velocityTotalInternal)
    ! Evaluate the distribution function in terms of scale-free (v_total,cos²θ).
    distributionFunction=+exp (                              &
         &                     -0.5d0                        &
         &                     *(                            &
         &                       +log(                       &
         &                            +velocityTotalInternal &
         &                            /self%mu1              &
         &                           )                       &
         &                       /self%sigma1                &
         &                      )**2                         &
         &                    )                              &
         &               /sqrt(                              &
         &                     +2.0d0                        &
         &                     *Pi                           &
         &                    )                              &
         &               /self%sigma1                        &
         &               /velocityTotalInternal
    if (eta > 0.0d0) then
       ! Anisotropic case - exponential distribution.
       distributionFunction=+distributionFunction             &
            &               *     eta                         &
            &               /(exp(eta                )-1.0d0) &
            &               * exp(eta*cosSquaredTheta)
    else
       ! Isotropic case - distribution is uniform in cos²θ.
       distributionFunction=+distributionFunction
    end if
    ! Transform distribution to (v_r,v_t) coordinates, and make dimensionful.
    distributionFunction=+distributionFunction          &
         &               *2.0d0                         &
         &               *velocityRadialInternal        &
         &               *velocityTangentialInternal    &
         &               /velocityTotalInternal     **3 &
         &               /velocityHost              **2
    return
  end function li2020VelocityDistributionFunction

  double precision function li2020VelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Numerical_Integration               , only : integrator
    implicit none
    class           (virialOrbitLi2020         ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                         , host
    class           (nodeComponentBasic        ), pointer       :: basicHost                    , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrastDefinition_
    double precision                            , parameter     :: extentVelocity        =10.0d0
    double precision                                            :: massHost                     , radiusHost          , &
         &                                                         velocityHost                 , massSatellite       , &
         &                                                         eta                          , velocityTotalMaximum
    type            (integrator                )                :: integratorVelocityTotal

    !![
    <referenceAcquire target="virialDensityContrastDefinition_" source="self%densityContrastDefinition()"/>
    !!]
    basic         => node%basic()
    basicHost     => host%basic()
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHost                                                                                          , &
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
    !![
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    velocityTotalMaximum     =+     self%mu1       &
         &                    *exp(                &
         &                         +extentVelocity &
         &                         *self%sigma1    &
         &                        )
    integratorVelocityTotal              =integrator                       (integrandVelocityTotal,toleranceRelative=1.0d-3)
    li2020VelocityTangentialMagnitudeMean=integratorVelocityTotal%integrate(0.0d0,velocityTotalMaximum)
    return

  contains

    double precision function integrandVelocityTotal(velocityTotal)
      !!{
      Integrand for the total velocity distribution function.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: velocityTotal
      double precision                :: integralCosTheta

      ! Evaluate the integral ∫₀¹ p(cos²θ) u √(1-cos²θ) dcos²θ analytically.
      eta=self%eta(host,massSatellite,massHost,velocityTotal)
      if (eta <= 0.0d0) then
         integralCosTheta=+2.0d0                                                    &
              &           /3.0d0                                                    &
              &           *velocityTotal
      else
         integralCosTheta=+(1.0d0-0.5d0*exp(eta)*sqrt(Pi)*erf(sqrt(eta))/sqrt(eta)) &
              &           /(1.0d0      -exp(eta)                                  ) &
              &           *velocityTotal
      end if
      integrandVelocityTotal=+integralCosTheta                &
           &                 *exp(                            &
           &                      -0.5d0                      &
           &                      *log(                       &
           &                           +     velocityTotal    &
           &                           /self%mu1              &
           &                          )                   **2 &
           &                      /     self%sigma1       **2 &
           &                     )                            &
           &                 /sqrt(2.0d0*Pi)                  &
           &                 /          self%sigma1           &
           &                 /               velocityTotal      
      return
    end function integrandVelocityTotal

  end function li2020VelocityTangentialMagnitudeMean

  function li2020VelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                   , dimension(3)  :: li2020VelocityTangentialVectorMean
    class           (virialOrbitLi2020), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node                              , host
    !$GLC attributes unused :: self, node, host

    li2020VelocityTangentialVectorMean=0.0d0
    call Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function li2020VelocityTangentialVectorMean

  double precision function li2020AngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: node        , host
    class           (nodeComponentBasic), pointer       :: basic       , basicHost
    double precision                                    :: massHost    , radiusHost, &
         &                                                 velocityHost

    basic                                 =>  node%basic()
    basicHost                             =>  host%basic()
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHost                                                                                          , &
         &                                                                   velocityHost                                                                                        , &
         &                                            cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                            cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                            virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                            darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                           )
    if (massHost > 0.0d0) then
       li2020AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
            &                                *radiusHost                                      &
            &                                /(                                               & ! Account for reduced mass.
            &                                  +1.0d0                                         &
            &                                  +basic    %mass()                              &
            &                                  /basicHost%mass()                              &
            &                                 )
    else
       li2020AngularMomentumMagnitudeMean =  +0.0d0
    end if
    return
  end function li2020AngularMomentumMagnitudeMean

  function li2020AngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector angular momentum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                   , dimension(3)  :: li2020AngularMomentumVectorMean
    class           (virialOrbitLi2020), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node                           , host
    !$GLC attributes unused :: self, node, host

    li2020AngularMomentumVectorMean=0.0d0
    call Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function li2020AngularMomentumVectorMean

  double precision function li2020VelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the root mean squared total velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitLi2020         ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                  , host
    class           (nodeComponentBasic        ), pointer       :: basicHost             , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrastDefinition_
    double precision                                            :: massHost              , radiusHost   , &
         &                                                         velocityHost          , massSatellite

    !![
    <referenceAcquire target="virialDensityContrastDefinition_" source="self%densityContrastDefinition()"/>
    !!]
    basic         => node%basic()
    basicHost     => host%basic()
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHost                                                                                          , &
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
    !![
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    li2020VelocityTotalRootMeanSquared=+exp(self%sigma1      **2) &
         &                             *    self%mu1              &
         &                             *         velocityHost
    return
  end function li2020VelocityTotalRootMeanSquared

  double precision function li2020EnergyMean(self,node,host)
    !!{
    Return the mean energy of the orbits.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: node        , host
    class           (nodeComponentBasic), pointer       :: basic       , basicHost
    double precision                                    :: massHost    , radiusHost, &
         &                                                 velocityHost
    
    basic            =>  node%basic()
    basicHost        =>  host%basic()
    massHost     =Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                   host                                                                                                , &
         &                                                                   self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                   radiusHost                                                                                          , &
         &                                                                   velocityHost                                                                                        , &
         &                                            cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                            cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                            virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                            darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                           )
    li2020EnergyMean =  +0.5d0                                           &
         &              *self%velocityTotalRootMeanSquared(node,host)**2 &
         &              /(                                               & ! Account for reduced mass.
         &                +1.0d0                                         &
         &                +basic    %mass()                              &
         &                /basicHost%mass()                              &
         &               )                                               &
         &              -gravitationalConstant_internal                  &
         &              *massHost                                        &
         &              /radiusHost
    return
  end function li2020EnergyMean

  double precision function li2020Eta(self,nodeHost,massSatellite,massHost,velocityTotalInternal)
    !!{
    Evaluate the $\eta$ parameter used in the definition of \cite{li_orbital_2020} virial orbits.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: nodeHost
    double precision                    , intent(in   ) :: massSatellite        , massHost , &
         &                                                 velocityTotalInternal
    class           (nodeComponentBasic), pointer       :: basic
    double precision                                    :: A                    , B        , &
         &                                                 peakHeight           , massRatio, &
         &                                                 time
    
    basic      =>  nodeHost%basic()
    time       =   basic   %time ()
    peakHeight =  +self%criticalOverdensity_     %value       (time=time,mass=massHost) &
         &        /self%cosmologicalMassVariance_%rootVariance(time=time,mass=massHost)
    massRatio  =  min(                &
         &            +massSatellite  &
         &            /massHost     , &
         &            +1.0d0          &
         &           )
    A          =+self%a1*peakHeight                   &
         &      +self%a2           *massRatio**self%c &
         &      +self%a3*peakHeight*massRatio**self%c
    B          =+self%b1                              &
         &      +self%b2           *massRatio**self%c
    li2020Eta  =  +max(                                   &
         &             +0.0d0                           , &
         &             +self%a0                           &
         &             *exp(                              &
         &                  -log(                         &
         &                       +velocityTotalInternal   &
         &                       /self%mu2                &
         &                      )       **2               &
         &                  /2.0d0                        &
         &                  /self%sigma1**2               &
         &                 )                              &
         &              +A*(velocityTotalInternal+1.0d0)  &
         &              +B                                &
         &            )
    return
  end function li2020Eta
