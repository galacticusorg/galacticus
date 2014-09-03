!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  use FGSL
  use Tables
  use Root_Finder

  !# <virialOrbit name="virialOrbitWetzel2010">
  !#  <description>Virial orbits using the \cite{wetzel_orbits_2010} orbital parameter distribution.</description>
  !# </virialOrbit>

  type, extends(virialOrbitClass) :: virialOrbitWetzel2010
     !% A virial orbit class using the \cite{wetzel_orbits_2010} orbital parameter distribution.
     private
     type   (fgsl_rng                )              :: clonedPseudoSequenceObject   , pseudoSequenceObject
     logical                                        :: resetSequence                , resetSequenceSnapshot
     integer                                        :: pericentricRadiusCount
     type   (table1DLogarithmicLinear)              :: pericentricRadiusTable
     class  (table1D                 ), allocatable :: pericentricRadiusTableInverse
     type   (rootFinder              )              :: finder
  contains
     procedure :: stateSnapshot             => wetzel2010StateSnapshot
     procedure :: stateStore                => wetzel2010StateStore
     procedure :: stateRestore              => wetzel2010StateRestore
     procedure :: orbit                     => wetzel2010Orbit
     procedure :: densityContrastDefinition => wetzel2010DensityContrastDefinition
  end type virialOrbitWetzel2010

  interface virialOrbitWetzel2010
     !% Constructors for the {\tt wetzel2010} virial orbit class.
     module procedure wetzel2010Constructor
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

  ! Global variables used in root-finding.
  double precision            :: wetzel2010C0                                         , wetzel2010C1                               , &
       &                         wetzel2010UniformDeviate
  !$omp threadprivate(wetzel2010C0,wetzel2010C1,wetzel2010UniformDeviate)

contains

  function wetzel2010Constructor()
    !% Generic constructor for the {\tt wetzel2010} virial orbits class.
    use Input_Parameters
    use Hypergeometric_Functions
    implicit none
    type            (virialOrbitWetzel2010), target    :: wetzel2010Constructor
    double precision                       , parameter :: toleranceAbsolute    =0.0d0, toleranceRelative                 =1.0d-2
    integer                                            :: iRadius
    double precision                                   :: x                          , xGamma2                                  , &
         &                                                probabilityCumulative      , probabilityCumulativeNormalization

    ! Initialize random sequence.
    wetzel2010Constructor%resetSequence=.true.
    ! Initialize root finder.
    call wetzel2010Constructor%finder%rootFunction(wetzel2010CircularityRoot          )
    call wetzel2010Constructor%finder%tolerance   (toleranceAbsolute,toleranceRelative)
    ! Construct a look-up table for the pericentric radius distribution.
    ! Determine number of points to use in the tabulation.
    wetzel2010Constructor%pericentricRadiusCount=int(log10(wetzel2010PericentricRadiusMaximum/wetzel2010PericentricRadiusMinimum)*dble(wetzel2010PericentricRadiusPointsPerDecade))+1
    ! Construct a range of radii.
    call wetzel2010Constructor%pericentricRadiusTable%destroy()
    call wetzel2010Constructor%pericentricRadiusTable%create(wetzel2010PericentricRadiusMinimum,wetzel2010PericentricRadiusMaximum,wetzel2010Constructor%pericentricRadiusCount)
    ! For each radius, compute the cumulative probability.
    do iRadius=wetzel2010Constructor%pericentricRadiusCount,1,-1
       x      =wetzel2010Constructor%pericentricRadiusTable%x(iRadius)
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
       if (iRadius == wetzel2010Constructor%pericentricRadiusCount) probabilityCumulativeNormalization=probabilityCumulative
       call wetzel2010Constructor%pericentricRadiusTable%populate(probabilityCumulative/probabilityCumulativeNormalization,iRadius)
    end do
    call wetzel2010Constructor%pericentricRadiusTable%reverse(wetzel2010Constructor%pericentricRadiusTableInverse)
    return
  end function wetzel2010Constructor

  function wetzel2010Orbit(self,node,host,acceptUnboundOrbits)
    !% Return wetzel2010 orbital parameters for a satellite.
    use Dark_Matter_Profile_Mass_Definitions
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    use Pseudo_Random
    use Critical_Overdensity
    use Cosmology_Functions
    implicit none
    type            (keplerOrbit               )                         :: wetzel2010Orbit
    class           (virialOrbitWetzel2010     ), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: host                   , node
    logical                                     , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentBasic        )               , pointer :: hostBasic             , basic
    class           (cosmologyFunctionsClass   )               , pointer :: cosmologyFunctions_
    class           (virialDensityContrastClass), pointer                :: virialDensityContrast_
    class           (darkMatterHaloScaleClass  )               , pointer :: darkMatterHaloScale_
    double precision                            , parameter              :: circularityMaximum       =1.0d0, circularityMinimum    =0.0d0
    double precision                            , parameter              :: redshiftMaximum          =5.0d0, expansionFactorMinimum=1.0d0/(1.0d0+redshiftMaximum)
    double precision                                                     :: R1                             , apocentricRadius                                    , &
         &                                                                  circularity                    , eccentricityInternal                                , &
         &                                                                  expansionFactor                , g1                                                  , &
         &                                                                  massCharacteristic             , pericentricRadius                                   , &
         &                                                                  probabilityTotal               , radialScale                                         , &
         &                                                                  timeNode                       , velocityHost                                        , &
         &                                                                  radiusHost                     , massHost
    logical                                                              :: foundOrbit

    ! Reset the orbit.
    call wetzel2010Orbit%reset()
    ! Get required objects.
    cosmologyFunctions_  => cosmologyFunctions ()
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Set masses and radius of the orbit.
    basic                => node%basic         ()
    hostBasic            => host%basic         ()
    ! Find virial density contrast under Wetzel (2010) definition.
    virialDensityContrast_          => self                  %densityContrastDefinition(                )
    ! Find mass, radius, and velocity in the host corresponding to the Wetzel (2010) virial density contrast definition.
    massHost=Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%time()),radiusHost,velocityHost)
    deallocate(virialDensityContrast_)
    ! Set basic properties of the orbit.
    call wetzel2010Orbit%massesSet(basic%mass(),hostBasic%mass())
    call wetzel2010Orbit%radiusSet(radiusHost                   )
    ! Get the time at which this node exists.
    timeNode=basic%time()
    ! Get the expansion factor.
    expansionFactor=cosmologyFunctions_%expansionFactor(timeNode)
    ! Limit the expansion factor to the smallest value considered by Wetzel.
    if (expansionFactor < expansionFactorMinimum) then
       expansionFactor=                               expansionFactorMinimum
       timeNode       =cosmologyFunctions_%cosmicTime(expansionFactorMinimum)
    end if
    ! Get the characteristic mass, M*.
    massCharacteristic=Critical_Overdensity_Collapsing_Mass(timeNode)
    ! Compute parameter of the circularity fitting function. We limit C1 to a given maximum - the fit is not explored in this
    ! regime and without the truncation we get problems evaluating hypergeometric functions.
    g1              =(1.0d0/expansionFactor)**wetzel2010CircularityP1
    wetzel2010C1    =min(wetzel2010CircularityAlpha1*(1.0d0+wetzel2010CircularityBeta1*(g1*wetzel2010Orbit%hostMass()/massCharacteristic)**wetzel2010CircularityGamma1),wetzel2010C1Maximum)
    wetzel2010C0    =1.0d0
    probabilityTotal=wetzel2010CircularityCumulativeProbability(circularityMaximum)
    wetzel2010C0    =1.0d0/probabilityTotal
    ! Compute parameter of the pericentric distance fitting function. Since the fit for R1 can lead to negative pericentric
    ! distances in some cases we force R1 to always be above a specified minimum.
    g1=(1.0d0/expansionFactor)**wetzel2010PericenterP1
    R1=max(wetzel2010PericenterAlpha1*(1.0d0+wetzel2010PericenterBeta1*(g1*wetzel2010Orbit%hostMass()/massCharacteristic)**wetzel2010PericenterGamma1),wetzel2010R1Minimum)
    ! Search for an orbit.
    foundOrbit=.false.
    do while (.not.foundOrbit)
       ! Compute pericentric radius by inversion in table.
       wetzel2010UniformDeviate=Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)
       pericentricRadius=R1*self%pericentricRadiusTableInverse%interpolate(wetzel2010UniformDeviate)
       ! Compute circularity by root finding in the cumulative probability distribution.
       wetzel2010UniformDeviate=Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)
       circularity=self%finder%find(rootRange=[circularityMinimum,circularityMaximum])
       ! Check that this is an orbit which actually reaches the virial radius.
       eccentricityInternal=sqrt(1.0d0-circularity**2)
       apocentricRadius    =pericentricRadius*(1.0d0+eccentricityInternal)/(1.0d0-eccentricityInternal)
       foundOrbit          =apocentricRadius >= 1.0d0 .and. pericentricRadius <= 1.0d0
       if (.not.foundOrbit) cycle
       ! Set eccentricity and periapsis.
       call wetzel2010Orbit%eccentricitySet    (sqrt(1.0d0-circularity**2)  )
       call wetzel2010Orbit%radiusPericenterSet(pericentricRadius*radiusHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=darkMatterHaloScale_%virialRadius(host)
       if (wetzel2010Orbit%radiusApocenter() >= radiusHost) then
          foundOrbit=.true.
          call wetzel2010Orbit%propagate(radiusHost,infalling=.true.)
       end if
    end do
    return
  end function wetzel2010Orbit

  function wetzel2010DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{wetzel_orbits_2010} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: wetzel2010DensityContrastDefinition
    class(virialOrbitWetzel2010     ), intent(inout) :: self
    
    allocate(virialDensityContrastFriendsOfFriends :: wetzel2010DensityContrastDefinition)
    select type (wetzel2010DensityContrastDefinition)
    type is (virialDensityContrastFriendsOfFriends)
      wetzel2010DensityContrastDefinition=virialDensityContrastFriendsOfFriends(linkingLength=0.168d0,densityRatio=4.688d0)
    end select
    return
  end function wetzel2010DensityContrastDefinition
  
  subroutine wetzel2010StateSnapshot(self)
    !% Write the tablulation state to file.
    implicit none
    class(virialOrbitWetzel2010), intent(inout) :: self

    if (.not.self%resetSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%pseudoSequenceObject)
    self%resetSequenceSnapshot=self%resetSequence
    return
  end subroutine wetzel2010StateSnapshot

  subroutine wetzel2010StateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitWetzel2010), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetSequenceSnapshot
    if (.not.self%resetSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine wetzel2010StateStore

  subroutine wetzel2010StateRestore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitWetzel2010), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    read (stateFile) self%resetSequence
    if (.not.self%resetSequence) call Pseudo_Random_Retrieve(self%pseudoSequenceObject,fgslStateFile)
   return
  end subroutine wetzel2010StateRestore

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
    use Hypergeometric_Functions
    implicit none
    double precision, intent(in   ) :: circularity

    wetzel2010CircularityCumulativeProbability=wetzel2010C0*(circularity**(wetzel2010CircularityGamma2+1.0d0))*Hypergeometric_2F1([-wetzel2010C1,1.0d0&
         &+wetzel2010CircularityGamma2],[2.0d0+wetzel2010CircularityGamma2],circularity)/(wetzel2010CircularityGamma2+1.0d0)
    return
  end function wetzel2010CircularityCumulativeProbability
