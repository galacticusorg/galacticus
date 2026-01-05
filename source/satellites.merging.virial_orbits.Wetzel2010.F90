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
  An implementation of virial orbits using the \cite{wetzel_orbits_2010} orbital parameter distribution.
  !!}

  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Root_Finder               , only : rootFinder
  use :: Tables                    , only : table1D                   , table1DLogarithmicLinear
  use :: Virial_Density_Contrast   , only : virialDensityContrastClass, virialDensityContrastFriendsOfFriends

  !![
  <virialOrbit name="virialOrbitWetzel2010">
   <description>
    A virial orbits class which selects orbital parameters randomly from the distribution given by \cite{wetzel_orbits_2010},
    including the redshift and mass dependence of the distributions. Note that the parameter $R_1$ can become negative (which
    is unphysical) for certain regimes of mass and redshift according to the fitting function for $R_1$ given by
    \cite{wetzel_orbits_2010}. Therefore, we enforce $R_1>0.05$. Similarly, the parameter $C_1$ can become very large in some
    regimes which is probably an artifact of the fitting function used rather than physically meaningful (and which causes
    numerical difficulties in evaluating the distribution). We therefore prevent $C_1$ from exceeding $9.999999$\footnote{We
    use this value rather than $10$ since the GSL $_2F_1$ hypergeometric function fails in some cases when $C_1\ge 10$.} If the
    virial density contrast definition differs from that used by \cite{wetzel_orbits_2010} then the orbit is assigned based on
    \cite{wetzel_orbits_2010}'s definition and then propagated to the virial radius relevant to the current definition of
    density contrast.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </stateStorable>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitWetzel2010
     !!{
     A virial orbit class using the \cite{wetzel_orbits_2010} orbital parameter distribution.
     !!}
     private
     integer                                                     :: pericentricRadiusCount
     type   (table1DLogarithmicLinear             )              :: pericentricRadiusTable
     class  (table1D                              ), allocatable :: pericentricRadiusTableInverse
     type   (rootFinder                           )              :: finder
     class  (darkMatterHaloScaleClass             ), pointer     :: darkMatterHaloScale_             => null()
     class  (darkMatterProfileDMOClass            ), pointer     :: darkMatterProfileDMO_            => null()
     class  (cosmologyParametersClass             ), pointer     :: cosmologyParameters_             => null()
     class  (cosmologyFunctionsClass              ), pointer     :: cosmologyFunctions_              => null()
     class  (criticalOverdensityClass             ), pointer     :: criticalOverdensity_             => null()
     class  (virialDensityContrastClass           ), pointer     :: virialDensityContrast_           => null()
     type   (virialDensityContrastFriendsOfFriends), pointer     :: virialDensityContrastDefinition_ => null()
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
     !!{
     Constructors for the \refClass{virialOrbitWetzel2010} virial orbit class.
     !!}
     module procedure wetzel2010ConstructorParameters
     module procedure wetzel2010ConstructorInternal
  end interface virialOrbitWetzel2010

  ! Table of the cumulative distribution for the pericentric radius.
  integer         , parameter :: pericentricRadiusPointsPerDecade=10
  double precision, parameter :: pericentricRadiusMaximum        =1.0d2     , pericentricRadiusMinimum=1.0d-6
  ! Parameters of the fitting functions.
  double precision, parameter :: circularityAlpha1               =0.242d0   , circularityBeta1        =2.360d0 , &
       &                         circularityGamma1               =0.108d0   , circularityGamma2       =1.05d0  , &
       &                         circularityP1                   =0.0d0
  double precision, parameter :: pericenterAlpha1                =0.450d0   , pericenterBeta1         =-0.395d0, &
       &                         pericenterGamma1                =0.109d0   , pericenterGamma2        =0.85d0  , &
       &                         pericenterP1                    =-4.0d0
  double precision, parameter :: c1Maximum                       =9.999999d0, r1Minimum               =0.05d0

  ! Module-scope variables used in root-finding.
  double precision            :: c0                                         , c1                               , &
       &                         uniformDeviate
  !$omp threadprivate(c0,c1,uniformDeviate)

contains

  function wetzel2010ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitWetzel2010} virial orbits class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialOrbitWetzel2010     )                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    class(cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class(criticalOverdensityClass  ), pointer       :: criticalOverdensity_
    class(virialDensityContrastClass), pointer       :: virialDensityContrast_
    class(darkMatterProfileDMOClass ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="criticalOverdensity"   name="criticalOverdensity_"   source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=virialOrbitWetzel2010(darkMatterHaloScale_,cosmologyFunctions_,criticalOverdensity_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="criticalOverdensity_"  />
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function wetzel2010ConstructorParameters

  function wetzel2010ConstructorInternal(darkMatterHaloScale_,cosmologyFunctions_,criticalOverdensity_,cosmologyParameters_,virialDensityContrast_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitWetzel2010} virial orbits class.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_1F1
    implicit none
    type            (virialOrbitWetzel2010     )                        :: self
    class           (darkMatterHaloScaleClass  ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class           (virialDensityContrastClass), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass ), intent(in   ), target :: darkMatterProfileDMO_
    class           (criticalOverdensityClass  ), intent(in   ), target :: criticalOverdensity_
    double precision                            , parameter             :: toleranceAbsolute     =0.0d0, toleranceRelative                 =1.0d-2
    integer                                                             :: iRadius
    double precision                                                    :: x                           , xGamma2                                  , &
         &                                                                 probabilityCumulative       , probabilityCumulativeNormalization
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologyParameters_, *virialDensityContrast_, *darkMatterProfileDMO_"/>
    !!]

    ! Initialize root finder.
    self%finder=rootFinder(                                             &
         &                 rootFunction     =circularityRoot, &
         &                 toleranceAbsolute=toleranceAbsolute        , &
         &                 toleranceRelative=toleranceRelative          &
         &                )
    ! Construct a look-up table for the pericentric radius distribution.
    ! Determine number of points to use in the tabulation.
    self%pericentricRadiusCount=int(log10(pericentricRadiusMaximum/pericentricRadiusMinimum)*dble(pericentricRadiusPointsPerDecade))+1
    ! Construct a range of radii.
    call self%pericentricRadiusTable%destroy()
    call self%pericentricRadiusTable%create(pericentricRadiusMinimum,pericentricRadiusMaximum,self%pericentricRadiusCount)
    ! For each radius, compute the cumulative probability.
    probabilityCumulativeNormalization=0.0d0
    do iRadius=self%pericentricRadiusCount,1,-1
       x      =self%pericentricRadiusTable%x(iRadius)
       xGamma2=x**pericenterGamma2
       probabilityCumulative=                                                                                                    &
            &  exp(-xGamma2)                                                                                                     &
            &  *x                                                                                                                &
            &  *(                                                                                                                &
            &     pericenterGamma2                                                                                               &
            &    *(                                                                                                              &
            &       1.0d0                                                                                                        &
            &      +pericenterGamma2                                                                                             &
            &      *(1.0d0+xGamma2)                                                                                              &
            &     )                                                                                                              &
            &    *Hypergeometric_1F1([2.0d0],[(1.0d0+3.0d0*pericenterGamma2)/pericenterGamma2],xGamma2)/(1.0d0+pericenterGamma2) &
            &    +Hypergeometric_1F1([1.0d0],[(1.0d0+3.0d0*pericenterGamma2)/pericenterGamma2],xGamma2)*(1.0d0+pericenterGamma2) &
            &   )
       if (iRadius == self%pericentricRadiusCount) probabilityCumulativeNormalization=probabilityCumulative
       call self%pericentricRadiusTable%populate(probabilityCumulative/probabilityCumulativeNormalization,iRadius)
    end do
    call self%pericentricRadiusTable%reverse(self%pericentricRadiusTableInverse)
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="virialDensityContrastDefinition_" constructor="virialDensityContrastFriendsOfFriends(linkingLength=0.168d0,densityRatio=4.688d0)"/>
    !!]
    return
  end function wetzel2010ConstructorInternal

  subroutine wetzel2010Destructor(self)
    !!{
    Destructor for the \refClass{virialOrbitWetzel2010} virial orbits class.
    !!}
    implicit none
    type(virialOrbitWetzel2010), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine wetzel2010Destructor

  function wetzel2010Orbit(self,node,host,acceptUnboundOrbits)
    !!{
    Return wetzel2010 orbital parameters for a satellite.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    type            (keplerOrbit               )                        :: wetzel2010Orbit
    class           (virialOrbitWetzel2010     ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                                  , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: basicHost                             , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrastDefinition_
    double precision                            , parameter             :: circularityMaximum              =1.0d0, circularityMinimum    =0.0d0
    double precision                            , parameter             :: redshiftMaximum                 =5.0d0, expansionFactorMinimum=1.0d0/(1.0d0+redshiftMaximum)
    double precision                                                    :: R1                                    , apocentricRadius                                    , &
         &                                                                 circularity                           , eccentricityInternal                                , &
         &                                                                 expansionFactor                       , g1                                                  , &
         &                                                                 massCharacteristic                    , pericentricRadius                                   , &
         &                                                                 probabilityTotal                      , massSatellite                                       , &
         &                                                                 timeNode                              , velocityHost                                        , &
         &                                                                 radiusHost                            , massHost                                            , &
         &                                                                 radiusHostSelf
    logical                                                     ::         foundOrbit
    !$GLC attributes unused :: acceptUnboundOrbits

    ! Set masses and radius of the orbit.
    basic                => node%basic         ()
    basicHost            => host%basic         ()
    ! Find virial density contrast under Wetzel (2010) definition.
    !![
    <referenceAcquire target="virialDensityContrastDefinition_" source="self%densityContrastDefinition()"/>
    !!]
    ! Find mass, radius, and velocity in the host corresponding to the Wetzel (2010) virial density contrast definition.
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
    g1              =(1.0d0/expansionFactor)**circularityP1
    c1    =min(circularityAlpha1*(1.0d0+circularityBeta1*(g1*massHost/massCharacteristic)**circularityGamma1),c1Maximum)
    c0    =1.0d0
    probabilityTotal=circularityCumulativeProbability(circularityMaximum)
    c0    =1.0d0/probabilityTotal
    ! Compute parameter of the pericentric distance fitting function. Since the fit for R1 can lead to negative pericentric
    ! distances in some cases we force R1 to always be above a specified minimum.
    g1=(1.0d0/expansionFactor)**pericenterP1
    R1=max(pericenterAlpha1*(1.0d0+pericenterBeta1*(g1*massHost/massCharacteristic)**pericenterGamma1),r1Minimum)
    ! Search for an orbit.
    foundOrbit=.false.
    do while (.not.foundOrbit)
       ! Reset the orbit.
       call wetzel2010Orbit%reset()
       ! Set basic properties of the orbit.
       call wetzel2010Orbit%massesSet(massSatellite,massHost      )
       call wetzel2010Orbit%radiusSet(              radiusHostSelf)
       ! Compute pericentric radius by inversion in table.
       uniformDeviate=node%hostTree%randomNumberGenerator_%uniformSample()
       pericentricRadius=R1*self%pericentricRadiusTableInverse%interpolate(uniformDeviate)
       ! Compute circularity by root finding in the cumulative probability distribution.
       uniformDeviate=node%hostTree%randomNumberGenerator_%uniformSample()
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
       radiusHost=self%darkMatterHaloScale_%radiusVirial(host)
       if (wetzel2010Orbit%radiusApocenter() >= radiusHost .and. wetzel2010Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call wetzel2010Orbit%propagate(radiusHost  ,infalling=.true.)
          call wetzel2010Orbit%massesSet(basic%mass(),basicHost%mass())
       end if
    end do
    return
  end function wetzel2010Orbit

  function wetzel2010DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of \cite{wetzel_orbits_2010} virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: wetzel2010DensityContrastDefinition
    class(virialOrbitWetzel2010     ), intent(inout) :: self

    wetzel2010DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function wetzel2010DensityContrastDefinition

  double precision function circularityRoot(circularity)
    !!{
    Function used in finding the circularity corresponding to a given cumulative probability.
    !!}
    double precision, intent(in   ) :: circularity
    double precision                :: cumulativeProbability

    cumulativeProbability=circularityCumulativeProbability(circularity)
    circularityRoot=cumulativeProbability-uniformDeviate
    return
  end function circularityRoot

  double precision function circularityCumulativeProbability(circularity)
    !!{
    The cumulative probability distribution for orbital circularity.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    implicit none
    double precision, intent(in   ) :: circularity

    circularityCumulativeProbability=c0*(circularity**(circularityGamma2+1.0d0))*Hypergeometric_2F1([-c1,1.0d0&
         &+circularityGamma2],[2.0d0+circularityGamma2],circularity)/(circularityGamma2+1.0d0)
    return
  end function circularityCumulativeProbability

  double precision function wetzel2010VelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(virialOrbitWetzel2010), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node, host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTangentialMagnitudeMean=0.0d0
    call Error_Report('mean tangential velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTangentialMagnitudeMean

  function wetzel2010VelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                       , dimension(3)  :: wetzel2010VelocityTangentialVectorMean
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                                  , host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTangentialVectorMean=0.0d0
    call Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTangentialVectorMean

  double precision function wetzel2010AngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , basicHost
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                                  =>  node%basic()
    basicHost                              =>  host%basic()
    massHost                               =   Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                                                host                                                                                                , &
         &                                                                                                self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                                                radiusHost                                                                                          , &
         &                                                                                                velocityHost                                                                                        , &
         &                                                                         cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                                         cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                                         virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                                                         darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                                                        )
    wetzel2010AngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                    *radiusHost                                      &
         &                                    /(                                               & ! Account for reduced mass.
         &                                      +1.0d0                                         &
         &                                      +basic    %mass()                              &
         &                                      /basicHost%mass()                              &
         &                                     )
    return
  end function wetzel2010AngularMomentumMagnitudeMean

  function wetzel2010AngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector angular momentum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                       , dimension(3)  :: wetzel2010AngularMomentumVectorMean
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node                               , host
    !$GLC attributes unused :: self, node, host

    wetzel2010AngularMomentumVectorMean=0.0d0
    call Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function wetzel2010AngularMomentumVectorMean

  double precision function wetzel2010VelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the root mean squared total velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(virialOrbitWetzel2010), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node, host
    !$GLC attributes unused :: self, node, host

    wetzel2010VelocityTotalRootMeanSquared=0.0d0
    call Error_Report('root mean squared total velocity is not defined for this class'//{introspection:location})
    return
  end function wetzel2010VelocityTotalRootMeanSquared

  double precision function wetzel2010EnergyMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    implicit none
    class           (virialOrbitWetzel2010), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node        , host
    class           (nodeComponentBasic   ), pointer       :: basic       , basicHost
    double precision                                       :: massHost    , radiusHost, &
         &                                                    velocityHost

    basic                =>  node%basic()
    basicHost            =>  host%basic()
    massHost             =   Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                              host                                                                                                , &
         &                                                                              self%virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()), &
         &                                                                              radiusHost                                                                                          , &
         &                                                                              velocityHost                                                                                        , &
         &                                                       cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                       cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                       virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                                       darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                                      )
    wetzel2010EnergyMean =  +0.5d0                                           &
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
  end function wetzel2010EnergyMean
