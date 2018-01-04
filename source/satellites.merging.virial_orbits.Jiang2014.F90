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

  !% An implementation of virial orbits using the \cite{jiang_orbital_2014} orbital parameter distribution.

  use Statistics_Distributions
  use Tables
  use Cosmology_Functions
  use Dark_Matter_Halo_Scales

  !# <virialOrbit name="virialOrbitJiang2014">
  !#  <description>Virial orbits using the \cite{jiang_orbital_2014} orbital parameter distribution.</description>
  !# </virialOrbit>
  type, extends(virialOrbitClass) :: virialOrbitJiang2014
     !% A virial orbit class using the \cite{jiang_orbital_2014} orbital parameter distribution.
     private
     class           (darkMatterHaloScaleClass), pointer        :: darkMatterHaloScale_ => null()
     class           (cosmologyFunctionsClass ), pointer        :: cosmologyFunctions_  => null()
     double precision                          , dimension(3,3) :: B                             , gamma, &
          &                                                        sigma                         , mu
     type            (table1DLinearLinear     ), dimension(3,3) :: voightDistributions
   contains
     final     ::                              jiang2014Destructor
     procedure :: orbit                     => jiang2014Orbit
     procedure :: densityContrastDefinition => jiang2014DensityContrastDefinition
  end type virialOrbitJiang2014

  interface virialOrbitJiang2014
     !% Constructors for the {\normalfont \ttfamily jiang2014} virial orbit class.
     module procedure jiang2014ConstructorParameters
     module procedure jiang2014ConstructorInternal
  end interface virialOrbitJiang2014

  ! Module-scope variables used in root finding.
  class           (virialOrbitJiang2014), pointer :: jiang2014Self
  double precision                                :: jiang2014XTotal                       , jiang2014XRadial              , &
       &                                             jiang204ProbabilityRadialNormalization, jiang2014VelocityTotalInternal
  integer                                         :: jiang2014I                            , jiang2014J
  !$omp threadprivate(jiang2014Self,jiang2014XTotal,jiang2014XRadial,jiang204ProbabilityRadialNormalization,jiang2014VelocityTotalInternal,jiang2014I,jiang2014J)

contains

  function jiang2014ConstructorParameters(parameters) result(self)
    !% Generic constructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    use Input_Parameters
    implicit none
    type            (virialOrbitJiang2014    )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                          , dimension(3)  :: bRatioLow           , bRatioIntermediate    , bRatioHigh    , &
         &                                                       gammaRatioLow       , gammaRatioIntermediate, gammaRatioHigh, &
         &                                                       sigmaRatioLow       , sigmaRatioIntermediate, sigmaRatioHigh, &
         &                                                       muRatioLow          , muRatioIntermediate   , muRatioHigh
    
    !# <inputParameter>
    !#   <name>bRatioLow</name>
    !#   <defaultValue>[+0.049d0,+0.548d0,+1.229d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bRatioIntermediate</name>
    !#   <defaultValue>[+1.044d0,+1.535d0,+3.396d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bRatioHigh</name>
    !#   <defaultValue>[+2.878d0,+3.946d0,+2.982d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $B$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioLow</name>
    !#   <defaultValue>[+0.109d0,+0.114d0,+0.110d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioIntermediate</name>
    !#   <defaultValue>[+0.098d0,+0.087d0,+0.050d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gammaRatioHigh</name>
    !#   <defaultValue>[+0.071d0,+0.030d0,-0.012d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\gamma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioLow</name>
    !#   <defaultValue>[+0.077d0,+0.094d0,+0.072d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioIntermediate</name>
    !#   <defaultValue>[+0.073d0,+0.083d0,+0.118d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaRatioHigh</name>
    !#   <defaultValue>[+0.091d0,+0.139d0,+0.187d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\sigma$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioLow</name>
    !#   <defaultValue>[+1.220d0,+1.231d0,+1.254d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.0001$--$0.005$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioIntermediate</name>
    !#   <defaultValue>[+1.181d0,+1.201d0,+1.236d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.005$--$0.05$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>muRatioHigh</name>
    !#   <defaultValue>[+1.100d0,+1.100d0,+1.084d0]</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Values of the $\mu$ parameter of the \cite{jiang_orbital_2014} orbital velocity distribution for the three host halo mass ranges, and the $0.05$--$0.5$ mass ratio.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=virialOrbitJiang2014(bRatioLow,bRatioIntermediate,bRatioHigh,gammaRatioLow,gammaRatioIntermediate,gammaRatioHigh,sigmaRatioLow,sigmaRatioIntermediate,sigmaRatioHigh,muRatioLow,muRatioIntermediate,muRatioHigh,darkMatterHaloScale_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function jiang2014ConstructorParameters

  function jiang2014ConstructorInternal(bRatioLow,bRatioIntermediate,bRatioHigh,gammaRatioLow,gammaRatioIntermediate,gammaRatioHigh,sigmaRatioLow,sigmaRatioIntermediate,sigmaRatioHigh,muRatioLow,muRatioIntermediate,muRatioHigh,darkMatterHaloScale_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    implicit none
    type            (virialOrbitJiang2014    )                        :: self
    double precision                          , dimension(3)          :: bRatioLow                , bRatioIntermediate    , bRatioHigh    , &
         &                                                               gammaRatioLow            , gammaRatioIntermediate, gammaRatioHigh, &
         &                                                               sigmaRatioLow            , sigmaRatioIntermediate, sigmaRatioHigh, &
         &                                                               muRatioLow               , muRatioIntermediate   , muRatioHigh
    class           (darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    integer                                   , parameter             :: tableCount          =1000
    integer                                                           :: i                        , j                     , k
    type            (distributionVoight      )                        :: voightDistribution
    !# <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_"/>

    ! Assign parameters of the distribution.
    self%B    (:,1)=    bRatioLow
    self%B    (:,2)=    bRatioIntermediate
    self%B    (:,3)=    bRatioHigh
    self%gamma(:,1)=gammaRatioLow
    self%gamma(:,2)=gammaRatioIntermediate
    self%gamma(:,3)=gammaRatioHigh
    self%sigma(:,1)=sigmaRatioLow
    self%sigma(:,2)=sigmaRatioIntermediate
    self%sigma(:,3)=sigmaRatioHigh
    self%mu   (:,1)=   muRatioLow
    self%mu   (:,2)=   muRatioIntermediate
    self%mu   (:,3)=   muRatioHigh    
    ! Tabulate Voight distribution functions for speed.
    ! Build the distribution function for total velocity.
    do i=1,3
       do j=1,3
          ! Build the distribution.
          voightDistribution=distributionVoight(                       &
               &                                self%gamma(i,j)      , &
               &                                self%mu   (i,j)      , &
               &                                self%sigma(i,j)      , &
               &                                limitLower     =0.0d0, &
               &                                limitUpper     =2.0d0  &
               &                               )
          ! Tabulate the cumulative distribution.
          call self%voightDistributions(i,j)%create(0.0d0,2.0d0,tableCount)
          !$omp parallel do
          do k=1,tableCount
             call self%voightDistributions(i,j)%populate(voightDistribution%cumulative(self%voightDistributions(i,j)%x(k)),k)
          end do
          !$omp end parallel do
       end do
    end do
    return
  end function jiang2014ConstructorInternal

  subroutine jiang2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily jiang2014} virial orbits class.
    implicit none
    type(virialOrbitJiang2014), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%cosmologyFunctions_"  />
    return
  end subroutine jiang2014Destructor

  function jiang2014Orbit(self,node,host,acceptUnboundOrbits)
    !% Return jiang2014 orbital parameters for a satellite.
    use Dark_Matter_Profile_Mass_Definitions
    use Galacticus_Error
    use Root_Finder
    implicit none
    type            (keplerOrbit               )                        :: jiang2014Orbit
    class           (virialOrbitJiang2014      ), intent(inout), target :: self
    type            (treeNode                  ), intent(inout)         :: host                   , node
    logical                                     , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic        ), pointer               :: hostBasic              , basic
    class           (virialDensityContrastClass), pointer               :: virialDensityContrast_
    integer                                     , parameter             :: attemptsMaximum        =10000
    double precision                            , parameter             :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                            , parameter             :: toleranceAbsolute      =0.0d+0
    double precision                            , parameter             :: toleranceRelative      =1.0d-3
    type            (rootFinder                )                        :: totalFinder                   , radialFinder
    double precision                                                    :: velocityHost                  , radiusHost                , &
         &                                                                 massHost                      , massSatellite             , &
         &                                                                 energyInternal                , radiusHostSelf            , &
         &                                                                 velocityRadialInternal        , velocityTangentialInternal
    logical                                                             :: foundOrbit
    integer                                                             :: attempts

    ! Get basic components.
    basic     => node%basic()
    hostBasic => host%basic()
    ! Find virial density contrast under Jiang et al. (2014) definition.
    virialDensityContrast_ => self %densityContrastDefinition()
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Jiang et al. (2014) virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%time()),radiusHostSelf,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%time())                            )
    deallocate(virialDensityContrast_)
    ! Select parameters appropriate for this host-satellite pair.
    if      (massHost < 10.0d0**12.5d0) then
       jiang2014I=1
    else if (massHost < 10.0d0**13.5d0) then
       jiang2014I=2
    else
       jiang2014I=3
    end if
    if      (massSatellite/massHost < 0.005d0) then
       jiang2014J=1
    else if (massSatellite/massHost < 0.050d0) then
       jiang2014J=2
    else
       jiang2014J=3
    end if
    ! Configure finder objects.
    if (.not.totalFinder %isInitialized()) then
       call totalFinder %rootFunction(jiang2014TotalVelocityCDF          )
       call totalFinder %tolerance   (toleranceAbsolute,toleranceRelative)
       call totalFinder %rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    if (.not.radialFinder%isInitialized()) then
       call radialFinder%rootFunction(jiang2014RadialVelocityCDF         )
       call radialFinder%tolerance   (toleranceAbsolute,toleranceRelative)
       call radialFinder%rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    ! Select an orbit.
    foundOrbit    =  .false.
    attempts      =  0
    jiang2014Self => self
    do while (.not.foundOrbit .and. attempts < attemptsMaximum)
       ! Increment number of attempts.
       attempts=attempts+1
       ! Reset the orbit.
       call jiang2014Orbit%reset()
       ! Set basic properties of the orbit.
       call jiang2014Orbit%massesSet(massSatellite,massHost      )
       call jiang2014Orbit%radiusSet(              radiusHostSelf)
       ! Solve for the total velocity.
       jiang2014XTotal               =node%hostTree%randomNumberGenerator%sample()
       jiang2014VelocityTotalInternal=totalFinder %find(rootGuess=1.0d0)
       ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
       ! bound that later rounding errors will not make it appear unbound.
       foundOrbit=.true.
       if (.not.acceptUnboundOrbits) then
          energyInternal=-1.0d0+0.5d0*jiang2014VelocityTotalInternal**2*jiang2014Orbit%specificReducedMass()
          foundOrbit=(energyInternal < -boundTolerance)
       end if
       if (.not.foundOrbit) cycle
       ! Solve for the radial velocity.
       jiang2014XRadial                      =+0.0d0
       jiang204ProbabilityRadialNormalization=+1.0d0
       jiang204ProbabilityRadialNormalization=+1.0d0                                                        &
            &                                 /(                                                            &
            &                                   +jiang2014RadialVelocityCDF(jiang2014VelocityTotalInternal) &
            &                                   -jiang2014RadialVelocityCDF(0.0d0                         ) &
            &                                  )
       jiang2014XRadial                      =node%hostTree%randomNumberGenerator%sample()
       velocityRadialInternal                =radialFinder%find(rootGuess=sqrt(2.0d0)*jiang2014VelocityTotalInternal)
       ! Compute tangential velocity.       
       velocityTangentialInternal=sqrt(max(0.0d0,jiang2014VelocityTotalInternal**2-velocityRadialInternal**2))
       call jiang2014Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call jiang2014Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%virialRadius(host)
       foundOrbit=.false.
       if (jiang2014Orbit%radiusApocenter() >= radiusHost .and. jiang2014Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call jiang2014Orbit%propagate(radiusHost  ,infalling=.true.)
          call jiang2014Orbit%massesSet(basic%mass(),hostBasic%mass())
       end if
    end do
    ! If too many iterations were required to find an orbit, abort.
    if (attempts >= attemptsMaximum) call Galacticus_Error_Report('maximum number of attempts exceeded'//{introspection:location})
    return
  end function jiang2014Orbit
  
  double precision function jiang2014TotalVelocityCDF(velocityTotal)
    !% Cumulative distribution function for the total velocity.
    implicit none
    double precision, intent(in   ) :: velocityTotal

    jiang2014TotalVelocityCDF=jiang2014Self%voightDistributions(jiang2014I,jiang2014J)%interpolate(velocityTotal)-jiang2014XTotal
    return
  end function jiang2014TotalVelocityCDF
  
  double precision function jiang2014RadialVelocityCDF(velocityRadial)
    !% Cumulative distribution function for the radial velocity.
    implicit none
    double precision, intent(in   ) :: velocityRadial
    
    jiang2014RadialVelocityCDF=+jiang204ProbabilityRadialNormalization                         &
         &                     *(                                                              &
         &                       +       jiang2014VelocityTotalInternal                        &
         &                       /       jiang2014Self%B               (jiang2014I,jiang2014J) &
         &                       *(                                                            &
         &                         +exp(                                                       &
         &                              +jiang2014Self%B               (jiang2014I,jiang2014J) &
         &                              *         velocityRadial                               &
         &                              /jiang2014VelocityTotalInternal                        &
         &                             )                                                       &
         &                         -1.0d0                                                      &
         &                        )                                                            &
         &                       -velocityRadial                                               &
         &                      )                                                              &
         &                     -jiang2014XRadial
    return
  end function jiang2014RadialVelocityCDF

  function jiang2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{jiang_orbital_2014} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: jiang2014DensityContrastDefinition
    class(virialOrbitJiang2014      ), intent(inout) :: self

    allocate(virialDensityContrastFixed :: jiang2014DensityContrastDefinition)
    select type (jiang2014DensityContrastDefinition)
    type is (virialDensityContrastFixed)
      jiang2014DensityContrastDefinition=virialDensityContrastFixed(200.0d0,fixedDensityTypeCritical,self%cosmologyFunctions_)
    end select
    return
  end function jiang2014DensityContrastDefinition
  
