!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% An implementation of virial orbits using the \cite{wetzel_orbits_2010} orbital parameter distribution.

  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Root_Finder               , only : rootFinder
  use :: Tables                    , only : table1D                              , table1DLogarithmicLinear
  use :: Virial_Density_Contrast   , only : virialDensityContrastFriendsOfFriends

  !# <virialOrbit name="virialOrbitWetzel2010">
  !#  <description>
  !#   A virial orbits class which selects orbital parameters randomly from the distribution given by \cite{wetzel_orbits_2010},
  !#   including the redshift and mass dependence of the distributions. Note that the parameter $R_1$ can become negative (which
  !#   is unphysical) for certain regimes of mass and redshift according to the fitting function for $R_1$ given by
  !#   \cite{wetzel_orbits_2010}. Therefore, we enforce $R_1>0.05$. Similarly, the parameter $C_1$ can become very large in some
  !#   regimes which is probably an artifact of the fitting function used rather than physically meaningful (and which causes
  !#   numerical difficulties in evaluating the distribution). We therefore prevent $C_1$ from exceeding $9.999999$\footnote{We
  !#   use this value rather than $10$ since the GSL $_2F_1$ hypergeomtric function fails in some cases when $C_1\ge 10$.} If the
  !#   virial density contrast definition differs from that used by \cite{wetzel_orbits_2010} then the orbit is assigned based on
  !#   \cite{wetzel_orbits_2010}'s definition and then propagated to the virial radius relevant to the current definition of
  !#   density contrast.
  !#  </description>
  !#  <deepCopy>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="virialDensityContrast_"/>
  !#  </stateStorable>
  !# </virialOrbit>

  type, extends(virialOrbitClass) :: virialOrbitWetzel2010
     !% A virial orbit class using the \cite{wetzel_orbits_2010} orbital parameter distribution.
     private
     integer                                                     :: pericentricRadiusCount
     type   (table1DLogarithmicLinear             )              :: pericentricRadiusTable
     class  (table1D                              ), allocatable :: pericentricRadiusTableInverse
     type   (rootFinder                           )              :: finder
     class  (darkMatterHaloScaleClass             ), pointer     :: darkMatterHaloScale_          => null()
     class  (cosmologyFunctionsClass              ), pointer     :: cosmologyFunctions_           => null()
     class  (criticalOverdensityClass             ), pointer     :: criticalOverdensity_          => null()
     type   (virialDensityContrastFriendsOfFriends), pointer     :: virialDensityContrast_        => null()
   contains
     final     ::                                    wetzel2010Destructor
     procedure :: orbit                           => wetzel2010Orbit
     procedure :: densityContrastDefinition       => wetzel2010DensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => wetzel2010VelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => wetzel2010VelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => wetzel2010AngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => wetzel2010AngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => wetzel2010VelocityTotalRootMeanSquared
     procedure :: energyMean                      => wetzel2010EnergyMean
  end type virialOrbitWetzel2010

  interface virialOrbitWetzel2010
     !% Constructors for the {\normalfont \ttfamily wetzel2010} virial orbit class.
     module procedure wetzel2010ConstructorParameters
     module procedure wetzel2010ConstructorInternal
  end interface virialOrbitWetzel2010

  ! Table of the cumulative distribution for the pericentric radius.
  integer         , parameter :: wetzel2010PericentricRadiusPointsPerDecade=10
  double precision, parameter :: wetzel2010PericentricRadiusMaximum        =1.0d2     , wetzel2010PericentricRadiusMinimum=1.0d-6
  ! Parameters of the fitting functions.
  double precision, parameter :: wetzel2010CircularityAlpha1               =0.242d0   , wetzel2010CircularityBeta1        =2.360d0 , &
       &                         wetzel2010CircularityGamma1               =0.108d0   , wetzel2010CircularityGamma2       =1.05d0  , &
       &                         wetzel2010CircularityP1                   =0.0d0
  double precision, parameter :: wetzel2010PericenterAlpha1                =0.450d0   , wetzel2010PericenterBeta1         =-0.395d0, &
       &                         wetzel2010PericenterGamma1                =0.109d0   , wetzel2010PericenterGamma2        =0.85d0  , &
       &                         wetzel2010PericenterP1                    =-4.0d0
  double precision, parameter :: wetzel2010C1Maximum                       =9.999999d0, wetzel2010R1Minimum               =0.05d0

  ! Module-scope variables used in root-finding.
  double precision            :: wetzel2010C0                                         , wetzel2010C1                               , &
       &                         wetzel2010UniformDeviate
  !$omp threadprivate(wetzel2010C0,wetzel2010C1,wetzel2010UniformDeviate)

contains

  function wetzel2010ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily wetzel2010} virial orbits class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialOrbitWetzel2010   )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class(criticalOverdensityClass), pointer       :: criticalOverdensity_

    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    self=virialOrbitWetzel2010(darkMatterHaloScale_,cosmologyFunctions_,criticalOverdensity_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="criticalOverdensity_"/>
    return
  end function wetzel2010ConstructorParameters

  function wetzel2010ConstructorInternal(darkMatterHaloScale_,cosmologyFunctions_,criticalOverdensity_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily wetzel2010} virial orbits class.
    use :: Hypergeometric_Functions, only : Hypergeometric_1F1
    implicit none
    type            (virialOrbitWetzel2010   )                        :: self
    class           (darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass), intent(in   ), target :: criticalOverdensity_
    double precision                          , parameter             :: toleranceAbsolute   =0.0d0, toleranceRelative                 =1.0d-2
    integer                                                           :: iRadius
    double precision                                                  :: x                         , xGamma2                                  , &
         &                                                               probabilityCumulative     , probabilityCumulativeNormalization
    !# <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_, *criticalOverdensity_"/>

    ! Initialize root finder.
    self%finder=rootFinder(                                             &
         &                 rootFunction     =wetzel2010CircularityRoot, &
         &                 toleranceAbsolute=toleranceAbsolute        , &
         &                 toleranceRelative=toleranceRelative          &
         &                )
    ! Construct a look-up table for the pericentric radius distribution.
    ! Determine number of points to use in the tabulation.
    self%pericentricRadiusCount=int(log10(wetzel2010PericentricRadiusMaximum/wetzel2010PericentricRadiusMinimum)*dble(wetzel2010PericentricRadiusPointsPerDecade))+1
    ! Construct a range of radii.
    call self%pericentricRadiusTable%destroy()
    call self%pericentricRadiusTable%create(wetzel2010PericentricRadiusMinimum,wetzel2010PericentricRadiusMaximum,self%pericentricRadiusCount)
    ! For each radius, compute the cumulative probability.
    probabilityCumulativeNormalization=0.0d0
    do iRadius=self%pericentricRadiusCount,1,-1
       x      =self%pericentricRadiusTable%x(iRadius)
       xGamma2=x**wetzel2010PericenterGamma2
       probabilityCumulative=                                                                                                                                  &
            &  exp(-xGamma2)                                                                                                                                   &
            &  *x                                                                                                                                              &
            &  *(                                                                                                                                              &
            &     wetzel2010PericenterGamma2                                                                                                                   &
            &    *(                                                                                                                                            &
            &       1.0d0                                                                                                                                      &
            &      +wetzel2010PericenterGamma2                                                                                                                 &
            &      *(1.0d0+xGamma2)                                                                                                                            &
            &     )                                                                                                                                            &
            &    *Hypergeometric_1F1([2.0d0],[(1.0d0+3.0d0*wetzel2010PericenterGamma2)/wetzel2010PericenterGamma2],xGamma2)/(1.0d0+wetzel2010PericenterGamma2) &
            &    +Hypergeometric_1F1([1.0d0],[(1.0d0+3.0d0*wetzel2010PericenterGamma2)/wetzel2010PericenterGamma2],xGamma2)*(1.0d0+wetzel2010PericenterGamma2) &
            &   )
       if (iRadius == self%pericentricRadiusCount) probabilityCumulativeNormalization=probabilityCumulative
       call self%pericentricRadiusTable%populate(probabilityCumulative/probabilityCumulativeNormalization,iRadius)
    end do
    call self%pericentricRadiusTable%reverse(self%pericentricRadiusTableInverse)
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrast_)
    !# <referenceConstruct isResult="yes" owner="self" object="virialDensityContrast_" constructor="virialDensityContrastFriendsOfFriends(linkingLength=0.168d0,densityRatio=4.688d0)"/>
    return
  end function wetzel2010ConstructorInternal

  subroutine wetzel2010Destructor(self)
    !% Destructor for the {\normalfont \ttfamily wetzel2010} virial orbits class.
    implicit none
    type(virialOrbitWetzel2010), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%criticalOverdensity_"  />
    !# <objectDestructor name="self%virialDensityContrast_"/>
    return
  end subroutine wetzel2010Destructor

  function wetzel2010Orbit(self,node,host,acceptUnboundOrbits)
    !% Return wetzel2010 orbital parameters for a satellite.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    type            (keplerOrbit               )                        :: wetzel2010Orbit
    class           (virialOrbitWetzel2010     ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                        , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: hostBasic                   , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrast_
    double precision                            , parameter             :: circularityMaximum    =1.0d0, circularityMinimum    =0.0d0
    double precision                            , parameter             :: redshiftMaximum       =5.0d0, expansionFactorMinimum=1.0d0/(1.0d0+redshiftMaximum)
    double precision                                                    :: R1                          , apocentricRadius                                    , &
         &                                                                 circularity                 , eccentricityInternal                                , &
         &                                                                 expansionFactor             , g1                                                  , &
         &                                                                 massCharacteristic          , pericentricRadius                                   , &
         &                                                                 probabilityTotal            , massSatellite                                       , &
         &                                                                 timeNode                    , velocityHost                                        , &
         &                                                                 radiusHost                  , massHost                                            , &
         &                                                                 radiusHostSelf
    logical                                                     ::         foundOrbit
    !$GLC attributes unused :: acceptUnboundOrbits

    ! Set masses and radius of the orbit.
    basic                => node%basic         ()
    hostBasic            => host%basic         ()
    ! Find virial density contrast under Wetzel (2010) definition.
    !# <referenceAcquire target="virialDensityContrast_" source="self%densityContrastDefinition()"/>
    ! Find mass, radius, and velocity in the host corresponding to the Wetzel (2010) virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHostSelf,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%timeLastIsolated())                            )
    !# <objectDestructor name="virialDensityContrast_"/>
    ! Get the time at which this node exists.
    timeNode=basic%time()
    ! Get the expansion factor.
    expansionFactor=self%cosmologyFunctions_%expansionFactor(timeNode)
    ! Limit the expansion factor to the smallest value considered by Wetzel.
    if (expansionFactor < expansionFactorMinimum) then
       expansionFactor=                               expansionFactorMinimum
       timeNode       =self%cosmologyFunctions_%cosmicTime(expansionFactorMinimum)
    end if
    ! Get the characteristic mass, M*.
    massCharacteristic=self%criticalOverdensity_%collapsingMass(timeNode)
    ! Compute parameter of the circularity fitting function. We limit C1 to a given maximum - the fit is not explored in this
    ! regime and without the truncation we get problems evaluating hypergeometric functions.
    g1              =(1.0d0/expansionFactor)**wetzel2010CircularityP1
    wetzel2010C1    =min(wetzel2010CircularityAlpha1*(1.0d0+wetzel2010CircularityBeta1*(g1*massHost/massCharacteristic)**wetzel2010CircularityGamma1),wetzel2010C1Maximum)
    wetzel2010C0    =1.0d0
    probabilityTotal=wetzel2010CircularityCumulativeProbability(circularityMaximum)
    wetzel2010C0    =1.0d0/probabilityTotal
    ! Compute parameter of the pericentric distance fitting function. Since the fit for R1 can lead to negative pericentric
    ! distances in some cases we force R1 to always be above a specified minimum.
    g1=(1.0d0/expansionFactor)**wetzel2010PericenterP1
    R1=max(wetzel2010PericenterAlpha1*(1.0d0+wetzel2010PericenterBeta1*(g1*massHost/massCharacteristic)**wetzel2010PericenterGamma1),wetzel2010R1Minimum)
    ! Search for an orbit.
    foundOrbit=.false.
    do while (.not.foundOrbit)
       ! Reset the orbit.
       call wetzel2010Orbit%reset()
       ! Set basic properties of the orbit.
       call wetzel2010Orbit%massesSet(massSatellite,massHost      )
       call wetzel2010Orbit%radiusSet(              radiusHostSelf)
       ! Compute pericentric radius by inversion in table.
       wetzel2010UniformDeviate=node%hostTree%randomNumberGenerator_%uniformSample()
       pericentricRadius=R1*self%pericentricRadiusTableInverse%interpolate(wetzel2010UniformDeviate)
       ! Compute circularity by root finding in the cumulative probability distribution.
       wetzel2010UniformDeviate=node%hostTree%randomNumberGenerator_%uniformSample()
       circularity=self%finder%find(rootRange=[circularityMinimum,circularityMaximum])
       ! Check that this is an orbit which actually reaches the virial radius.
       eccentricityInternal=sqrt(1.0d0-circularity**2)
       apocentricRadius    =pericentricRadius*(1.0d0+eccentricityInternal)/(1.0d0-eccentricityInternal)
       foundOrbit          =apocentricRadius >= 1.0d0 .and. pericentricRadius <= 1.0d0
       if (.not.foundOrbit) cycle
       ! Set eccentricity and periapsis.
       call wetzel2010Orbit%eccentricitySet    (sqrt(1.0d0-circularity**2)      )
       call wetzel2010Orbit%radiusPericenterSet(pericentricRadius*radiusHostSelf)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%virialRadius(host)
       if (wetzel2010Orbit%radiusApocenter() >= radiusHost .and. wetzel2010Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call wetzel2010Orbit%propagate(radiusHost  ,infalling=.true.)
          call wetzel2010Orbit%massesSet(basic%mass(),hostBasic%mass())
       end if
    end do
    return
  end function wetzel2010Orbit

  function wetzel2010DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{wetzel_orbits_2010} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: wetzel2010DensityContrastDefinition
    class(virialOrbitWetzel2010     ), intent(inout) :: self

    wetzel2010DensityContrastDefinition => self%virialDensityContrast_
    return
  end function wetzel2010DensityContrastDefinition

  double precision function wetzel2010CircularityRoot(circularity)
    !% Function used in finding the circularity corresponding to a given cumulative probability.
    double precision, intent(in   ) :: circularity
    double precision                :: cumulativeProbability

    cumulativeProbability=wetzel2010CircularityCumulativeProbability(circularity)
    wetzel2010CircularityRoot=cumulativeProbability-wetzel2010UniformDeviate
    return
  end function wetzel2010CircularityRoot

  double precision function wetzel2010CircularityCumulativeProbability(circularity)
    !% The cumulative probability distribution for orbital circularity.
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    implicit none
    double precision, intent(in   ) :: circularity

    wetzel2010CircularityCumulativeProbability=wetzel2010C0*(circularity**(wetzel2010CircularityGamma2+1.0d0))*Hypergeometric_2F1([-wetzel2010C1,1.0d0&
         &+wetzel2010CircularityGamma2],[2.0d0+wetzel2010CircularityGamma2],circularity)/(wetzel2010CircularityGamma2+1.0d0)
    return
  end function wetzel2010CircularityCumulativeProbability

  double precision function wetzel2010VelocityTangentialMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the tangential velocity.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(virialOrbitWetzel2010), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node, host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTangentialMagnitudeMean=0.0d0
    call Galacticus_Error_Report('mean tangential velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTangentialMagnitudeMean

  function wetzel2010VelocityTangentialVectorMean(self,node,host)
    !% Return the mean of the vector tangential velocity.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                       , dimension(3)  :: wetzel2010VelocityTangentialVectorMean
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                                  , host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTangentialVectorMean=0.0d0
    call Galacticus_Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTangentialVectorMean

  double precision function wetzel2010AngularMomentumMagnitudeMean(self,node,host)
    !% Return the mean magnitude of the angular momentum.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , hostBasic
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                                  =>  node%basic()
    hostBasic                              =>  host%basic()
    massHost                               =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    wetzel2010AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                    *radiusHost                                      &
         &                                    /(                                               & ! Account for reduced mass.
         &                                      +1.0d0                                         &
         &                                      +basic    %mass()                              &
         &                                      /hostBasic%mass()                              &
         &                                     )
    return
  end function wetzel2010AngularMomentumMagnitudeMean

  function wetzel2010AngularMomentumVectorMean(self,node,host)
    !% Return the mean of the vector angular momentum.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                       , dimension(3)  :: wetzel2010AngularMomentumVectorMean
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                               , host
    !$GLC attributes unused :: self, node, host

    wetzel2010AngularMomentumVectorMean=0.0d0
    call Galacticus_Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function wetzel2010AngularMomentumVectorMean

  double precision function wetzel2010VelocityTotalRootMeanSquared(self,node,host)
    !% Return the root mean squared total velocity.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(virialOrbitWetzel2010), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node, host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTotalRootMeanSquared=0.0d0
    call Galacticus_Error_Report('root mean squared total velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTotalRootMeanSquared

  double precision function wetzel2010EnergyMean(self,node,host)
    !% Return the mean magnitude of the tangential velocity.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical        , only : gravitationalConstantGalacticus
    implicit none
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , hostBasic
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                =>  node%basic()
    hostBasic            =>  host%basic()
    massHost             =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%timeLastIsolated()),radiusHost,velocityHost)
    wetzel2010EnergyMean =  +0.5d0                                           &
         &                  *self%velocityTotalRootMeanSquared(node,host)**2 &
         &                  /(                                               & ! Account for reduced mass.
         &                    +1.0d0                                         &
         &                    +basic    %mass()                              &
         &                    /hostBasic%mass()                              &
         &                   )                                               &
         &                  -gravitationalConstantGalacticus                 &
         &                  *massHost                                        &
         &                  /radiusHost
    return
  end function wetzel2010EnergyMean
