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
Contains a module which provides an object that implements cosmological functions.
!!}

module Cosmology_Functions
  !!{
  Provides an object that implements cosmological functions.
  !!}
  private

  !![
  <functionClass>
   <name>cosmologyFunctions</name>
   <descriptiveName>Cosmology Functions</descriptiveName>
   <description>
    A class that provides various cosmological functions. The background cosmology describes the evolution of an isotropic,
    homogeneous Universe within which our calculations are carried out. For the purposes of \glc, the background cosmology is
    used to relate expansion factor/redshift to cosmic time and to compute the density of various components (e.g. dark matter,
    dark energy, etc.) at different epochs.
   </description>
   <default>matterLambda</default>
   <method name="epochValidate" >
    <description>Check the given cosmic epoch is valid (aborting otherwise) and, optionally, return time or expansion factor associated with the epoch.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: timeIn            </argument>
    <argument>double precision, intent(in   ), optional :: expansionFactorIn </argument>
    <argument>logical         , intent(in   ), optional :: collapsingIn      </argument>
    <argument>double precision, intent(  out), optional :: timeOut           </argument>
    <argument>double precision, intent(  out), optional :: expansionFactorOut</argument>
    <argument>logical         , intent(  out), optional :: collapsingOut     </argument>
   </method>
   <method name="cosmicTime" >
    <description>Return the cosmological age at the given expansion factor.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="timeBigCrunch" >
    <description>Return the cosmological age at Big Crunch (or a negative value if no Big Crunch occurs).</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="expansionFactor" >
    <description>Returns the expansion factor at cosmological time {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="expansionRate" >
    <description>Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: expansionFactor</argument>
   </method>
   <method name="hubbleParameterEpochal" >
    <description>Returns the Hubble parameter at the requested cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time,expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="hubbleParameterRateOfChange" >
    <description>Returns the rate of change of the Hubble parameter at the requested cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time,expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="densityScalingEarlyTime" >
    <description>Compute the scaling of density with expansion factor at early times in the universe.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: dominateFactor</argument>
    <argument>double precision, intent(  out)           :: densityPower  , expansionFactorDominant</argument>
    <argument>double precision, intent(  out), optional :: OmegaDominant</argument>
   </method>
   <method name="omegaMatterEpochal" >
    <description>Return the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="omegaMatterRateOfChange" >
    <description>Return the rate of change of the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="omegaDarkEnergyEpochal" >
    <description>Return the dark energy density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="equationOfStateDarkEnergy" >
    <description>Return the equation of state paramerter, $w$, for the dark energy component.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time, expansionFactor</argument>
   </method>
   <method name="exponentDarkEnergy" >
    <description>Return the exponent of the dark energy density with expansion factor, i.e. $\gamma$ in $\rho_\mathrm{DE}(a) \propto a^\gamma$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time, expansionFactor</argument>
   </method>
   <method name="equalityEpochMatterDarkEnergy" >
    <description>Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ), optional :: requestType</argument>
   </method>
   <method name="equalityEpochMatterCurvature" >
    <description>Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ), optional :: requestType</argument>
   </method>
   <method name="equalityEpochMatterRadiation" >
    <description>Return the epoch of matter-radiation magnitude equality (either expansion factor or cosmic time).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ), optional :: requestType</argument>
   </method>
   <method name="dominationEpochMatter" >
    <description>Compute the epoch at which matter dominates over other forms of energy by a given factor.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: dominateFactor</argument>
   </method>
   <method name="temperatureCMBEpochal" >
    <description>Return the temperature of the cosmic microwave background at {\normalfont \ttfamily expansionFactor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
   <method name="distanceParticleHorizonComoving" >
    <description>Return the comoving distance to the particle horizon at the given {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="distanceComoving" >
    <description>Return the comoving distance to the given cosmic {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="distanceLuminosity" >
    <description>Return the luminosity distance to the given cosmic {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="distanceAngular" >
    <description>Return the angular diameter distance to the given cosmic {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time      </argument>
    <argument>double precision, intent(in   ), optional :: timeOrigin</argument>
   </method>
   <method name="timeAtDistanceComoving" >
    <description>Return the cosmic time corresponding to the given {\normalfont \ttfamily comovingDistance}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: comovingDistance</argument>
   </method>
   <method name="distanceComovingConvert" >
    <description>Convert between different measures of comoving distance.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer         , intent(in   )           :: output</argument>
    <argument>double precision, intent(in   ), optional :: distanceLuminosity, distanceModulus, distanceModulusKCorrected, redshift</argument>
   </method>
   <method name="redshiftFromExpansionFactor" >
    <description>Returns redshift for a given expansion factor.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: expansionFactor</argument>
    <code>
     !$GLC attributes unused :: self
     cosmologyFunctionsRedshiftFromExpansionFactor=1.0d0/expansionFactor-1.0d0
    </code>
   </method>
   <method name="expansionFactorFromRedshift" >
    <description>Returns expansion factor given a redshift.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: redshift</argument>
    <code>
     !$GLC attributes unused :: self
     cosmologyFunctionsExpansionFactorFromRedshift=1.0d0/(1.0d0+redshift)
    </code>
   </method>
   <method name="comovingVolumeElementRedshift" >
    <description>Returns the differential comoving volume element $\mathrm{d}V/\mathrm{d}z = r_\mathrm{c}^2(t) \mathrm{c} H^{-1}(t)$ (where $r_\mathrm{c}$ is the comoving distance to time $t$ and $H(t)$ is the Hubble parameter at that time) for unit solid angle at the specified {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Numerical_Constants_Physical Numerical_Constants_Prefixes</modules>
    <argument>double precision, intent(in   ) :: time</argument>
    <code>cosmologyFunctionsComovingVolumeElementRedshift=self%distanceComoving(time)**2*(speedLight/kilo)/self%hubbleParameterEpochal(time=time)</code>
   </method>
   <method name="comovingVolumeElementTime" >
    <description>Returns the differential comoving volume element $\mathrm{d}V/\mathrm{d}t = r_\mathrm{c}^2(t) \mathrm{c} a(t)$ (where $r_\mathrm{c}$ is the comoving distance to time $t$ and $a(t)$ is the expansion at that time) for unit solid angle at the specified {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Numerical_Constants_Astronomical Numerical_Constants_Physical</modules>
    <argument>double precision, intent(in   ) :: time</argument>
    <code>cosmologyFunctionsComovingVolumeElementTime=self%distanceComoving(time)**2*(gigaYear*speedLight/megaParsec)/self%expansionfactor(time)</code>
   </method>
   <method name="epochTime" >
    <description>Convenience function that returns the time corresponding to an epoch specified by time or expansion factor.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time           , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
    <modules>Error</modules>
    <code>
    if (present(time)) then
       cosmologyFunctionsEpochTime=time
    else if (present(expansionFactor)) then
       cosmologyFunctionsEpochTime=self%cosmicTime(expansionFactor,collapsingPhase)
    else
       cosmologyFunctionsEpochTime=-1.0d0
       call Error_Report('either "time" or "expansionFactor" must be given'//{introspection:location})
    end if
    </code>
   </method>
   <method name="matterDensityEpochal" >
    <description>Convenience function that returns the matter density at the specified epoch.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: time           , expansionFactor</argument>
    <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
   </method>
  </functionClass>
  !!]

  ! A recommended relative time tolerance to which other functions should approach the Big Crunch.
  double precision, parameter, public :: timeToleranceRelativeBigCrunch=1.0d-4

  ! Enumeration for different cosmological densities.
  !![
  <enumeration>
   <name>densityCosmological</name>
   <description>Enumeration of different cosmological densities.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <entry label="critical" />
   <entry label="mean"     />
  </enumeration>
  !!]

end module Cosmology_Functions
