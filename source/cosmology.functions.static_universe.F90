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
  An implementation of the cosmological functions class for static universes. Intended for testing purposes.
  !!}

  !![
  <cosmologyFunctions name="cosmologyFunctionsStaticUniverse">
   <description>
    A cosmology functions class for static universes. Intended for testing purposes. Time is arbitrary (as there is no Big
    Bang), and expansion factor is fixed at $1$. Attempts to compute time from expansion factor will cause fatal
    errors. Expansion rates and the Hubble constant are set to zero. </description>
  </cosmologyFunctions>
  !!]
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  type, extends(cosmologyFunctionsClass) :: cosmologyFunctionsStaticUniverse
     !!{
     A cosmological functions class for cosmologies consisting of matter plus a cosmological constant.
     !!}
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
   contains
     final     ::                                  staticUniverseDestructor
     procedure :: epochValidate                 => staticUniverseEpochValidate
     procedure :: cosmicTime                    => staticUniverseCosmicTime
     procedure :: timeBigCrunch                 => staticUniverseTimeBigCrunch
     procedure :: expansionFactor               => staticUniverseExpansionFactor
     procedure :: expansionRate                 => staticUniverseExpansionRate
     procedure :: hubbleParameterEpochal        => staticUniverseHubbleParameterEpochal
     procedure :: hubbleParameterRateOfChange   => staticUniverseHubbleParameterRateOfChange
     procedure :: densityScalingEarlyTime       => staticUniverseDensityScalingEarlyTime
     procedure :: omegaMatterEpochal            => staticUniverseOmegaMatterEpochal
     procedure :: omegaMatterRateOfChange       => staticUniverseOmegaMatterRateOfChange
     procedure :: omegaDarkEnergyEpochal        => staticUniverseOmegaDarkEnergyEpochal
     procedure :: equationOfStateDarkEnergy     => staticUniverseEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy            => staticUniverseExponentDarkEnergy
     procedure :: equalityEpochMatterDarkEnergy => staticUniverseEqualityEpochMatterDarkEnergy
     procedure :: equalityEpochMatterCurvature  => staticUniverseEqualityEpochMatterCurvature
     procedure :: equalityEpochMatterRadiation  => staticUniverseEqualityEpochMatterRadiation
     procedure :: dominationEpochMatter         => staticUniverseDominationEpochMatter
     procedure :: temperatureCMBEpochal         => staticUniverseTemperatureCMBEpochal
     procedure :: distanceComoving              => staticUniverseDistanceComoving
     procedure :: distanceLuminosity            => staticUniverseDistanceLuminosity
     procedure :: distanceAngular               => staticUniverseDistanceAngular
     procedure :: timeAtDistanceComoving        => staticUniverseTimeAtDistanceComoving
     procedure :: distanceComovingConvert       => staticUniverseDistanceComovingConvert
     procedure :: matterDensityEpochal          => staticUniverseMatterDensityEpochal
  end type cosmologyFunctionsStaticUniverse

  interface cosmologyFunctionsStaticUniverse
     !!{
     Constructors for the matter plus cosmological constant cosmological functions class.
     !!}
     module procedure staticUniverseConstructorParameters
     module procedure staticUniverseConstructorInternal
  end interface cosmologyFunctionsStaticUniverse

contains

  function staticUniverseConstructorParameters(parameters) result(self)
    !!{
    Parameter-based constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (cosmologyFunctionsStaticUniverse)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(cosmologyParametersClass        ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=cosmologyFunctionsStaticUniverse(cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function staticUniverseConstructorParameters

  function staticUniverseConstructorInternal(cosmologyParameters_) result(self)
    !!{
    Constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    implicit none
    type (cosmologyFunctionsStaticUniverse)               , target :: self
    class(cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="*cosmologyParameters_"/>
    !!]

    return
  end function staticUniverseConstructorInternal

  subroutine staticUniverseDestructor(self)
    !!{
    Default constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    implicit none
    type(cosmologyFunctionsStaticUniverse), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine staticUniverseDestructor

  subroutine staticUniverseEpochValidate(self,timeIn,expansionFactorIn,collapsingIn,timeOut,expansionFactorOut,collapsingOut)
    !!{
    Validate a cosmic epoch, specified either by time or expansion factor, and optionally return time, expansion factor, and
    collapsing status.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactorIn , timeIn
    logical                                           , intent(in   ), optional :: collapsingIn
    double precision                                  , intent(  out), optional :: expansionFactorOut, timeOut
    logical                                           , intent(  out), optional :: collapsingOut
    !$GLC attributes unused :: self

    ! Check that we have a uniquely specified epoch.
    if (present(expansionFactorIn)) call Error_Report('time can not be determined from expansion factor'      //{introspection:location})
    if (.not.(present(timeIn).or. present(expansionFactorIn))) &
         & call Error_Report('one of "time" or "expansionFactor" should be specified'                         //{introspection:location})
    if (      present(timeIn).and.present(collapsingIn     ) ) &
         & call Error_Report('collapsing status of universe cannot be specified when epoch is defined by time'//{introspection:location})
    if (                          present(collapsingOut    ) ) then
       collapsingOut=.false.
       call Error_Report('collapsing status of universe is undefined'//{introspection:location})
    end if
    ! If we have a time, check that it is a valid, and compute outputs as required.
    if (present(timeIn)) then
        ! Set outputs.
       if (present(timeOut           )) timeOut           =timeIn
       if (present(expansionFactorOut)) expansionFactorOut=1.0d0
    end if
    ! If we have an expansion factor, check that it is valid, and compute outputs as required.
    if (present(expansionFactorIn)) then
       ! Validate.
       if ( expansionFactorIn /= 1.0d0) &
            & call Error_Report('expansion factor must be unity'//{introspection:location})
       ! Set outputs.
       if (present(timeOut           )) call Error_Report('time can not be determined from expansion factor'//{introspection:location})
       if (present(expansionFactorOut)) expansionFactorOut=expansionFactorIn
    end if
    return
  end subroutine staticUniverseEpochValidate

  double precision function staticUniverseCosmicTime(self,expansionFactor,collapsingPhase)
    !!{
    Return the cosmological time at a given expansion factor.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   )           :: expansionFactor
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, expansionFactor, collapsingPhase

    staticUniverseCosmicTime=0.0d0
    call Error_Report('time can not be determined from expansion factor'//{introspection:location})
    return
  end function staticUniverseCosmicTime

  double precision function staticUniverseTimeBigCrunch(self)
    !!{
    Return the time of the Big Crunch in a static cosmology. Since no Big Crunch occurs we return a negative value to indicate
    that.
    !!}
    implicit none
    class(cosmologyFunctionsStaticUniverse), intent(inout) :: self
    !$GLC attributes unused :: self

    staticUniverseTimeBigCrunch=-1.0d0
    return
  end function staticUniverseTimeBigCrunch

  double precision function staticUniverseExpansionFactor(self,time)
    !!{
    Returns the expansion factor at cosmological time {\normalfont \ttfamily time}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    staticUniverseExpansionFactor=1.0d0
    return
  end function staticUniverseExpansionFactor

  double precision function staticUniverseExpansionRate(self,expansionFactor)
    !!{
    Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: expansionFactor
    !$GLC attributes unused :: self, expansionFactor

    staticUniverseExpansionRate=0.0d0
    return
  end function staticUniverseExpansionRate

  double precision function staticUniverseHubbleParameterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the Hubble parameter at the request cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseHubbleParameterEpochal=0.0d0
    return
  end function staticUniverseHubbleParameterEpochal

  double precision function staticUniverseHubbleParameterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the rate of change of the Hubble parameter at the request cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseHubbleParameterRateOfChange=0.0d0
    return
  end function staticUniverseHubbleParameterRateOfChange

  double precision function staticUniverseOmegaMatterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseOmegaMatterEpochal=0.0d0
    call Error_Report('Omega_Matter is undefined in static universe'//{introspection:location})
    return
  end function staticUniverseOmegaMatterEpochal

  double precision function staticUniverseMatterDensityEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the matter density at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseMatterDensityEpochal=self%cosmologyParameters_%omegaMatter()*self%cosmologyParameters_%densityCritical()
    return
  end function staticUniverseMatterDensityEpochal

  double precision function staticUniverseOmegaMatterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the rate of change of the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseOmegaMatterRateOfChange=0.0d0
    call Error_Report('Omega_Matter is undefined in static universe'//{introspection:location})
    return
  end function staticUniverseOmegaMatterRateOfChange

  double precision function staticUniverseOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the dark energy density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseOmegaDarkEnergyEpochal=0.0d0
    call Error_Report('Omega_DarkEnergy is undefined in static universe'//{introspection:location})
    return
  end function staticUniverseOmegaDarkEnergyEpochal

  double precision function staticUniverseTemperatureCMBEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the temperature of the CMB at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor      , time
    logical                                           , intent(in   ), optional :: collapsingPhase
    !$GLC attributes unused :: self, time, expansionFactor, collapsingPhase

    staticUniverseTemperatureCMBEpochal=self%cosmologyParameters_%temperatureCMB()
    return
  end function staticUniverseTemperatureCMBEpochal

  subroutine staticUniverseDensityScalingEarlyTime(self,dominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   )           :: dominateFactor
    double precision                                  , intent(  out)           :: densityPower  , expansionFactorDominant
    double precision                                  , intent(  out), optional :: OmegaDominant
    !$GLC attributes unused :: self, dominateFactor, densityPower, expansionFactorDominant, OmegaDominant

    densityPower           =0.0d0
    expansionFactorDominant=0.0d0
    if (present(OmegaDominant)) OmegaDominant=0.0d0
    call Error_Report('early times are undefined in static universe'//{introspection:location})
    return
  end subroutine staticUniverseDensityScalingEarlyTime

  double precision function staticUniverseDominationEpochMatter(self,dominateFactor)
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: dominateFactor
    !$GLC attributes unused :: self, dominateFactor

    staticUniverseDominationEpochMatter=0.0d0
    call Error_Report('epochs are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseDominationEpochMatter

  double precision function staticUniverseEqualityEpochMatterDarkEnergy(self,requestType)
    !!{
    Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    integer                                  , intent(in   ), optional :: requestType
    !$GLC attributes unused :: self, requestType

    staticUniverseEqualityEpochMatterDarkEnergy=0.0d0
    call Error_Report('epochs are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseEqualityEpochMatterDarkEnergy

  double precision function staticUniverseEqualityEpochMatterCurvature(self,requestType)
    !!{
    Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    integer                                  , intent(in   ), optional :: requestType
    !$GLC attributes unused :: self, requestType

    staticUniverseEqualityEpochMatterCurvature=0.0d0
    call Error_Report('epochs are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseEqualityEpochMatterCurvature

  double precision function staticUniverseEqualityEpochMatterRadiation(self,requestType)
    !!{
    Return the epoch of matter-radiation magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    integer                                  , intent(in   ), optional :: requestType
    !$GLC attributes unused :: self, requestType

    staticUniverseEqualityEpochMatterRadiation=0.0d0
    call Error_Report('epochs are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseEqualityEpochMatterRadiation

  double precision function staticUniverseTimeAtDistanceComoving(self,comovingDistance)
    !!{
    Returns the cosmological time corresponding to given {\normalfont \ttfamily comovingDistance}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: comovingDistance
    !$GLC attributes unused :: self, comovingDistance

    staticUniverseTimeAtDistanceComoving=0.0d0
    call Error_Report('absolute times are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseTimeAtDistanceComoving

  double precision function staticUniverseDistanceComoving(self,time)
    !!{
    Returns the comoving distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    staticUniverseDistanceComoving=0.0d0
    call Error_Report('absolute times are undefined in static universe'//{introspection:location})
    return
  end function staticUniverseDistanceComoving

  double precision function staticUniverseDistanceLuminosity(self,time)
    !!{
    Returns the luminosity distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout) :: self
    double precision                                  , intent(in   ) :: time

    ! Compute the luminosity distance.
    staticUniverseDistanceLuminosity=self%distanceComoving(time)
    return
  end function staticUniverseDistanceLuminosity

  double precision function staticUniverseDistanceAngular(self,time,timeOrigin)
    !!{
    Returns the angular diameter distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   )           :: time
    double precision                                  , intent(in   ), optional :: timeOrigin

    staticUniverseDistanceAngular=self%distanceComoving(time)
    if (present(timeOrigin)) then
       if (timeOrigin < time) call Error_Report('expected timeOrigin ≥ time'//{introspection:location})
       staticUniverseDistanceAngular=+staticUniverseDistanceAngular     &
            &                        -self%distanceComoving(timeOrigin)
    end if
    return
  end function staticUniverseDistanceAngular

  double precision function staticUniverseDistanceComovingConvert(self,output,distanceLuminosity,distanceModulus,distanceModulusKCorrected,redshift)
    !!{
    Convert between different measures of distance.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    integer                                           , intent(in   )           :: output
    double precision                                  , intent(in   ), optional :: distanceModulus    , distanceModulusKCorrected, &
         &                                                                         redshift           , distanceLuminosity
    logical                                                                     :: gotComovingDistance
    double precision                                                            :: comovingDistance
    !$GLC attributes unused :: self

    ! Convert to comoving distance from whatever was supplied.
    gotComovingDistance=.false.
    comovingDistance   =-1.0d0
    if (present(distanceModulus)) then
       comovingDistance  =distanceLuminosity
       gotComovingDistance=.true.
    else if (present(distanceModulus)) then
       comovingDistance  =10.0d0**((distanceModulus          -25.0d0)/5.0d0)
       gotComovingDistance=.true.
    else if (present(distanceModulusKCorrected)) then
       comovingDistance  =10.0d0**((distanceModulusKCorrected-25.0d0)/5.0d0)
       gotComovingDistance=.true.
    end if
    if (present(redshift)) call Error_Report('zero redshift is undefined in static universe'//{introspection:location})
    if (.not.gotComovingDistance) &
         & call Error_Report('no distance measure provided'//{introspection:location})
    ! Convert to required distance measure.
    select case (output)
    case (distanceTypeComoving)
       staticUniverseDistanceComovingConvert=comovingDistance
    case default
       staticUniverseDistanceComovingConvert=-1.0d0
       call Error_Report('unrecognized output option'//{introspection:location})
    end select
    return
  end function staticUniverseDistanceComovingConvert

  double precision function staticUniverseEquationOfStateDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy equation of state.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    !$GLC attributes unused :: self, time, expansionFactor

    staticUniverseEquationOfStateDarkEnergy=-1.0d0
    return
  end function staticUniverseEquationOfStateDarkEnergy

  double precision function staticUniverseExponentDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy equation of state.
    !!}
    implicit none
    class           (cosmologyFunctionsStaticUniverse), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: expansionFactor, time
    !$GLC attributes unused :: self, time, expansionFactor

    staticUniverseExponentDarkEnergy=0.0d0
    return
  end function staticUniverseExponentDarkEnergy
