!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An implementation of virial orbits using the \cite{li_orbital_2020} orbital parameter distribution.

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass       , criticalOverdensityClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastBryanNorman1998
  use :: Numerical_Interpolation   , only : interpolator
  
  !# <virialOrbit name="virialOrbitLi2020">
  !#  <description>Virial orbits using the \cite{li_orbital_2020} orbital parameter distribution.</description>
  !#  <deepCopy>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </stateStorable>
  !# </virialOrbit>
  type, extends(virialOrbitClass) :: virialOrbitLi2020
     !% A virial orbit class using the \cite{li_orbital_2020} orbital parameter distribution.
     private
     class           (darkMatterHaloScaleClass            ), pointer :: darkMatterHaloScale_      => null()
     class           (cosmologyParametersClass            ), pointer :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass             ), pointer :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass            ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass       ), pointer :: cosmologicalMassVariance_ => null()
     type            (virialDensityContrastBryanNorman1998), pointer :: virialDensityContrast_    => null()
     double precision                                                :: mu                                 , sigma1  , &
          &                                                             sigma2                             , a0      , &
          &                                                             a1                                 , a2      , &
          &                                                             a3                                 , aMinimum, &
          &                                                             aMaximum
     type            (interpolator                        )          :: interpolatorA
     logical                                                         :: propagateOrbits
   contains
     !# <methods>
     !#   <method description="Evaluate the $\eta$ parameter of the \cite{li_orbital_2020} virial orbit distribution function." method="eta" />
     !# </methods>
     final     ::                                    li2020Destructor
     procedure :: orbit                           => li2020Orbit
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
     !% Constructors for the {\normalfont \ttfamily li2020} virial orbit class.
     module procedure li2020ConstructorParameters
     module procedure li2020ConstructorInternal
  end interface virialOrbitLi2020

contains

  function li2020ConstructorParameters(parameters) result(self)
    !% Generic constructor for the {\normalfont \ttfamily li2020} virial orbits class.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitLi2020            )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                                               :: mu                       , sigma1, &
         &                                                            sigma2                   , a0    , &
         &                                                            a1                       , a2    , &
         &                                                            a3 
    logical                                                        :: propagateOrbits
    
    !# <inputParameter>
    !#   <name>mu</name>
    !#   <defaultValue>1.21d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigma1</name>
    !#   <defaultValue>0.22d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigma2</name>
    !#   <defaultValue>0.40d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma_2$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>a0</name>
    !#   <defaultValue>-0.97d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $a_0$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>a1</name>
    !#   <defaultValue>0.74d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $a_1$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>a2</name>
    !#   <defaultValue>4.80d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $a_2$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>a3</name>
    !#   <defaultValue>0.40d0</defaultValue>
    !#   <defaultSource>\citep[][Table~2]{li_orbital_2020}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>Values of the $a_3$ parameter of the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propagateOrbits</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>If true, orbits will be propagated to the virial radius of the host halo. Otherwise, no propagation is performed and orbital velocities will correspond to precisely those from the \cite{li_orbital_2020} orbital velocity distribution.</description>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    self=virialOrbitLi2020(mu,sigma1,sigma2,a0,a1,a2,a3,propagateOrbits,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"     />
    !# <objectDestructor name="cosmologyParameters_"     />
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="criticalOverdensity_"     />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    return
  end function li2020ConstructorParameters

  function li2020ConstructorInternal(mu,sigma1,sigma2,a0,a1,a2,a3,propagateOrbits,darkMatterHaloScale_,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily li2020} virial orbits class.
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLinear
    implicit none
    type            (virialOrbitLi2020            )                             :: self
    class           (darkMatterHaloScaleClass     ), intent(in   ), target      :: darkMatterHaloScale_
    class           (cosmologyParametersClass     ), intent(in   ), target      :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target      :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), intent(in   ), target      :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), intent(in   ), target      :: cosmologicalMassVariance_
    double precision                               , intent(in   )              :: mu                              , sigma1                           , &
         &                                                                         sigma2                          , a0                               , &
         &                                                                         a1                              , a2                               , &
         &                                                                         a3
    logical                                        , intent(in   )              :: propagateOrbits
    integer                                        , parameter                  :: countTable               =1000
    double precision                               , parameter                  :: peakHeightMaximum        =10.0d0, massRatioMaximum           =2.0d0, &
         &                                                                         extentVelocity           =10.0d0
    double precision                               , dimension(:) , allocatable :: aTable                          , velocityTangentialMeanTable
    type            (integrator                   )                             :: integratorVelocityTotal         , integratorCosSquaredTheta
    double precision                                                               velocityTotal_                  , eta                              , &
         &                                                                         a                               , velocityTotalMaximum
    integer                                                                     :: i
    
    !# <constructorAssign variables="mu, sigma1, sigma2, a0, a1, a2, a3, propagateOrbits, *darkMatterHaloScale_, *cosmologyParameters_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrast_)
    !# <referenceConstruct isResult="yes" owner="self" object="virialDensityContrast_" constructor="virialDensityContrastBryanNorman1998(self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    
    ! Tabulate the mean tangential velocity as a function of a=a₀+a₁ν+a₂ξ^a₃.
    allocate(aTable                     (countTable))
    allocate(velocityTangentialMeanTable(countTable))
    integratorVelocityTotal  = integrator(integrandVelocityTotal  ,toleranceRelative=1.0d-3)
    integratorCosSquaredTheta= integrator(integrandCosSquaredTheta,toleranceRelative=1.0d-3)
    velocityTotalMaximum     =+     self%mu        &
         &                    *exp(                &
         &                         +extentVelocity &
         &                         *self%sigma1    &
         &                        )
    self%aMinimum            =+0.0d0
    self%aMaximum            =+self%a0                            &
         &                    +self%a1*peakHeightMaximum          &
         &                    +self%a2*massRatioMaximum **self%a3
    aTable                   =Make_Range(self%aMinimum,self%aMaximum,countTable,rangeTypeLinear)
    do i=1,countTable
       a                             =aTable(i)
       velocityTangentialMeanTable(i)=integratorVelocityTotal%integrate(0.0d0,velocityTotalMaximum)
    end do
    self%interpolatorA=interpolator(aTable,velocityTangentialMeanTable)
    return

  contains

    double precision function integrandVelocityTotal(velocityTotal)
      !% Integrand for the total velocity distribution function.
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: velocityTotal

      velocityTotal_=velocityTotal
      eta=a*exp(-log(min(velocityTotal,sqrt(2.0d0)))**2/2.0d0/self%sigma2**2)
      integrandVelocityTotal=integratorCosSquaredTheta%integrate(0.0d0,1.0d0)*exp(-0.5d0*log(velocityTotal/self%mu)**2/self%sigma1**2)/sqrt(2.0d0*Pi)/self%sigma1/velocityTotal      
      return
    end function integrandVelocityTotal
    
    double precision function integrandCosSquaredTheta(cosSquaredTheta)
      !% Integrand for the $\cos^2\theta$ distribution, weighted by the tangential velocity.
      implicit none
      double precision, intent(in   ) :: cosSquaredTheta
      double precision, parameter     :: etaSmall       =1.0d-6

      if (eta < etaSmall) then
         ! Series solution for small values of η.
         integrandCosSquaredTheta=+exp(+eta   *cosSquaredTheta) &
              &                   *   (                         &
              &                        +        1.0d0           &
              &                        -eta   / 2.0d0           &
              &                        +eta**2/12.0d0           &
              &                   )
      else
         ! Full solution for larger values of η.
         integrandCosSquaredTheta=+     eta                     &
              &                   *exp(+eta   *cosSquaredTheta) &
              &                   /   (                         &
              &                        +exp(eta)                &
              &                        -1.0d0                   &
              &                       )
      end if
      ! Multiply by the tangential velocity.
      integrandCosSquaredTheta=+integrandCosSquaredTheta*velocityTotal_ &
           &                   *sqrt(                                   &
           &                         +1.0d0                             &
           &                         -cosSquaredTheta                   &
           &                        )
      return
    end function integrandCosSquaredTheta
    
  end function li2020ConstructorInternal

  subroutine li2020Destructor(self)
    !% Destructor for the {\normalfont \ttfamily li2020} virial orbits class.
    implicit none
    type(virialOrbitLi2020), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"     />
    !# <objectDestructor name="self%cosmologyParameters_"     />
    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%criticalOverdensity_"     />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%virialDensityContrast_"   />
    return
  end subroutine li2020Destructor

  function li2020Orbit(self,node,host,acceptUnboundOrbits)
    !% Return li2020 orbital parameters for a satellite.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    type            (keplerOrbit               )                        :: li2020Orbit
    class           (virialOrbitLi2020         ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                          , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: hostBasic                     , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrast_
    integer                                     , parameter             :: attemptsMaximum        =10000
    double precision                            , parameter             :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                                                    :: velocityHost                  , radiusHost                , &
         &                                                                 massHost                      , massSatellite             , &
         &                                                                 energyInternal                , radiusHostSelf            , &
         &                                                                 velocityRadialInternal        , velocityTangentialInternal, &
         &                                                                 velocityTotalInternal         , eta
    logical                                                             :: foundOrbit
    integer                                                             :: attempts

    ! Get basic components.
    basic     => node%basic()
    hostBasic => host%basic()
    ! Find virial density contrast under Li et al. (2020) definition.
    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>    
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Li et al. (2020) virial density contrast
    ! definition.
    massHost     =             Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHostSelf,velocityHost)
    massSatellite=min(massHost,Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                            ))
    !# <objectDestructor name="virialDensityContrast_"/>
    ! Select an orbit.
    foundOrbit=.false.
    attempts  =0
    do while (.not.foundOrbit .and. attempts < attemptsMaximum)
       ! Increment number of attempts.
       attempts=attempts+1
       ! Reset the orbit.
       call li2020Orbit%reset()
       ! Set basic properties of the orbit.
       call li2020Orbit%massesSet(massSatellite,massHost      )
       call li2020Orbit%radiusSet(              radiusHostSelf)
       ! Select total orbital velocity from a log-normal distribution.
       velocityTotalInternal=exp(self%sigma1*node%hostTree%randomNumberGenerator_%standardNormalSample())*self%mu
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
          radiusHost=self%darkMatterHaloScale_%virialRadius(host)
          foundOrbit=.false.
          if (li2020Orbit%radiusApocenter() >= radiusHost .and. li2020Orbit%radiusPericenter() <= radiusHost) then
             foundOrbit=.true.
             call li2020Orbit%propagate(radiusHost  ,infalling=.true.)
             call li2020Orbit%massesSet(basic%mass(),hostBasic%mass())
          end if
       end if
    end do
    ! If too many iterations were required to find an orbit, abort.
    if (attempts >= attemptsMaximum) call Galacticus_Error_Report('maximum number of attempts exceeded'//{introspection:location})
    return
  end function li2020Orbit

  function li2020DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{li_orbital_2020} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: li2020DensityContrastDefinition
    class(virialOrbitLi2020         ), intent(inout) :: self

    li2020DensityContrastDefinition => self%virialDensityContrast_
    return
  end function li2020DensityContrastDefinition

  double precision function li2020VelocityTangentialMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the tangential velocity.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    class           (virialOrbitLi2020         ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                  , host
    class           (nodeComponentBasic        ), pointer       :: hostBasic             , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: massHost              , radiusHost   , &
         &                                                         velocityHost          , massSatellite, &
         &                                                         a

    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>
    basic         => node%basic()
    hostBasic     => host%basic()
    massHost      =  Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    massSatellite =  Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                        )
    !# <objectDestructor name="virialDensityContrast_"/>

    ! Evaluate a=....... - this is just eta(u=1). 
    a=min(max(self%eta(host,massSatellite,massHost,velocityTotalInternal=1.0d0),self%aMinimum),self%aMaximum)
    li2020VelocityTangentialMagnitudeMean=self%interpolatorA%interpolate(a)
    return
  end function li2020VelocityTangentialMagnitudeMean

  function li2020VelocityTangentialVectorMean(self,node,host)
    !% Return the mean of the vector tangential velocity.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                   , dimension(3)  :: li2020VelocityTangentialVectorMean
    class           (virialOrbitLi2020), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node                              , host
    !$GLC attributes unused :: self, node, host

    li2020VelocityTangentialVectorMean=0.0d0
    call Galacticus_Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function li2020VelocityTangentialVectorMean

  double precision function li2020AngularMomentumMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the angular momentum.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: node        , host
    class           (nodeComponentBasic), pointer       :: basic       , hostBasic
    double precision                                    :: massHost    , radiusHost, &
         &                                                 velocityHost

    basic                               =>  node%basic()
    hostBasic                           =>  host%basic()
    massHost                            =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
     li2020AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                 *radiusHost                                      &
         &                                 /(                                               & ! Account for reduced mass.
         &                                   +1.0d0                                         &
         &                                   +basic    %mass()                              &
         &                                   /hostBasic%mass()                              &
         &                                  )
    return
  end function li2020AngularMomentumMagnitudeMean

  function li2020AngularMomentumVectorMean(self,node,host)
    !% Return the mean of the vector angular momentum.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                   , dimension(3)  :: li2020AngularMomentumVectorMean
    class           (virialOrbitLi2020), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node                           , host
    !$GLC attributes unused :: self, node, host

    li2020AngularMomentumVectorMean=0.0d0
    call Galacticus_Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function li2020AngularMomentumVectorMean

  double precision function li2020VelocityTotalRootMeanSquared(self,node,host)
    !% Return the root mean squared total velocity.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitLi2020         ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node                  , host
    class           (nodeComponentBasic        ), pointer       :: hostBasic             , basic
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: massHost              , radiusHost   , &
         &                                                         velocityHost          , massSatellite

    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>
    basic         => node%basic()
    hostBasic     => host%basic()
    massHost      =  Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    massSatellite =  Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                        )
    !# <objectDestructor name="virialDensityContrast_"/>
    li2020VelocityTotalRootMeanSquared=+exp(self%sigma1      **2) &
         &                             *    self%mu               &
         &                             *         velocityHost
    return
  end function li2020VelocityTotalRootMeanSquared

  double precision function li2020EnergyMean(self,node,host)
    !% Return the mean energy of the orbits.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstantGalacticus
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: node        , host
    class           (nodeComponentBasic), pointer       :: basic       , hostBasic
    double precision                                    :: massHost    , radiusHost, &
         &                                                 velocityHost
    
    basic            =>  node%basic()
    hostBasic        =>  host%basic()
    massHost         =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    li2020EnergyMean =  +0.5d0                                           &
         &              *self%velocityTotalRootMeanSquared(node,host)**2 &
         &              /(                                               & ! Account for reduced mass.
         &                +1.0d0                                         &
         &                +basic    %mass()                              &
         &                /hostBasic%mass()                              &
         &               )                                               &
         &              -gravitationalConstantGalacticus                 &
         &              *massHost                                        &
         &              /radiusHost
    return
  end function li2020EnergyMean

  double precision function li2020Eta(self,nodeHost,massSatellite,massHost,velocityTotalInternal)
    !% Evaluate the $\eta$ parameter used in the definition of \cite{li_orbital_2020} virial orbits.
    use :: Galacticus_Nodes, only : nodeComponentBasic 
    implicit none
    class           (virialOrbitLi2020 ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: nodeHost
    double precision                    , intent(in   ) :: massSatellite        , massHost , &
         &                                                 velocityTotalInternal
    class           (nodeComponentBasic), pointer       :: basic
    double precision                                    :: peakHeight           , massRatio, &
         &                                                 time
    
    basic      =>  nodeHost%basic()
    time       =   basic   %time ()
    peakHeight =  +self%criticalOverdensity_     %value       (time=time,mass=massHost) &
         &        /self%cosmologicalMassVariance_%rootVariance(time=time,mass=massHost)
    massRatio  =  +massSatellite &
         &        /massHost
    li2020Eta  =  +max(                                       &
         &             +0.0d0                               , &
         &             +(                                     &
         &               +self%a0                             &
         &               +self%a1*peakHeight                  &
         &               +self%a2*massRatio **self%a3         &
         &              )                                     &
         &             *exp(                                  &
         &                  -log(                             &
         &                       +min(                        &
         &                            +sqrt(2.0d0)          , &
         &                            +velocityTotalInternal  &
         &                           )                        &
         &                      )       **2                   &
         &                  /2.0d0                            &
         &                  /self%sigma2**2                   &
         &                 )                                  &
         &            )
    return
  end function li2020Eta
