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

!% Contains a module which provides an object that implements cosmological functions.

module Cosmology_Functions
  !% Provides an object that implements cosmological functions.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  !# <include directive="cosmologyFunctions" type="functionModules" >
  include 'cosmologyFunctions.functionModules.inc'
  !# </include>
  private

  !# <include directive="cosmologyFunctions" type="function" >
  !#  <descriptiveName>Cosmology Functions</descriptiveName>
  !#  <description>Object providing various cosmological functions.</description>
  !#  <default>matterLambda</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <method name="expansionFactorIsValid" >
  !#   <description>Returns true if the given expansion factor is valid one for this cosmology.</description>
  !#   <type>logical</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: expansionFactor</argument>
  !#  </method>
  !#  <method name="cosmicTimeIsValid" >
  !#   <description>Returns true if the given cosmic time is valid one for this cosmology.</description>
  !#   <type>logical</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="cosmicTime" >
  !#   <description>Return the cosmological age at the given expansion factor.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="expansionFactor" >
  !#   <description>Returns the expansion factor at cosmological time {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#   <bindC>true</bindC>
  !#  </method>
  !#  <method name="expansionRate" >
  !#   <description>Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: expansionFactor</argument>
  !#  </method>
  !#  <method name="hubbleParameterEpochal" >
  !#   <description>Returns the Hubble parameter at the requested cosmological time, {\tt time}, or expansion factor, {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time,expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="hubbleParameterRateOfChange" >
  !#   <description>Returns the rate of change of the Hubble parameter at the requested cosmological time, {\tt time}, or expansion factor, {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time,expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="densityScalingEarlyTime" >
  !#   <description>Compute the scaling of density with expansion factor at early times in the universe.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: dominateFactor</argument>
  !#   <argument>double precision, intent(  out)           :: densityPower  , expansionFactorDominant</argument>
  !#   <argument>double precision, intent(  out), optional :: OmegaDominant</argument>
  !#  </method>
  !#  <method name="omegaMatterEpochal" >
  !#   <description>Return the matter density parameter at expansion factor {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="omegaMatterRateOfChange" >
  !#   <description>Return the rate of change of the matter density parameter at expansion factor {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="omegaDarkEnergyEpochal" >
  !#   <description>Return the dark energy density parameter at expansion factor {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="equationOfStateDarkEnergy" >
  !#   <description></description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time, expansionFactor</argument>
  !#  </method>
  !#  <method name="exponentDarkEnergy" >
  !#   <description></description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time, expansionFactor</argument>
  !#  </method>
  !#  <method name="equalityEpochMatterDarkEnergy" >
  !#   <description>Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>integer, intent(in   ), optional :: requestType</argument>
  !#  </method>
  !#  <method name="equalityEpochMatterCurvature" >
  !#   <description>Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>integer, intent(in   ), optional :: requestType</argument>
  !#  </method>
  !#  <method name="dominationEpochMatter" >
  !#   <description>Compute the epoch at which matter dominates over other forms of energy by a given factor.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: dominateFactor</argument>
  !#  </method>
  !#  <method name="temperatureCMBEpochal" >
  !#   <description>Return the temperature of the cosmic microwave background at {\tt expansionFactor}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time  , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#  </method>
  !#  <method name="distanceComoving" >
  !#   <description>Return the comoving distance to the given cosmic {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="distanceLuminosity" >
  !#   <description>Return the luminosity distance to the given cosmic {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="distanceAngular" >
  !#   <description>Return the angular diameter distance to the given cosmic {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="timeAtDistanceComoving" >
  !#   <description>Return the cosmic time corresponding to the given {\tt comovingDistance}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: comovingDistance</argument>
  !#  </method>
  !#  <method name="distanceComovingConvert" >
  !#   <description>Convert between different measures of comoving distance.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>integer         , intent(in   )           :: output</argument>
  !#   <argument>double precision, intent(in   ), optional :: distanceModulus, redshift</argument>
  !#  </method>
  !#  <method name="redshiftFromExpansionFactor" >
  !#   <description>Returns redshift for a given expansion factor.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: expansionFactor</argument>
  !#   <code>redshiftFromExpansionFactor=1.0d0/expansionFactor-1.0d0</code>
  !#  </method>
  !#  <method name="expansionFactorFromRedshift" >
  !#   <description>Returns expansion factor given a redshift.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: redshift</argument>
  !#   <code>expansionFactorFromRedshift=1.0d0/(1.0d0+redshift)</code>
  !#  </method>
  !#  <method name="comovingVolumeElementRedshift" >
  !#   <description>Returns the differential comoving volume element ${\rm d}V/{\rm d}z = r_{\rm c}^2(t) {\rm c} H^{-1}(t)$ (where $r_{\rm c}$ is the comoving distance to time $t$ and $H(t)$ is the Hubble parameter at that time) for unit solid angle at the specified {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Numerical_Constants_Physical</modules>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#   <code>comovingVolumeElementRedshift=self%distanceComoving(time)**2*(speedLight/kilo)/self%hubbleParameterEpochal(time=time)</code>
  !#  </method>
  !#  <method name="comovingVolumeElementTime" >
  !#   <description>Returns the differential comoving volume element ${\rm d}V/{\rm d}t = r_{\rm c}^2(t) {\rm c} a(t)$ (where $r_{\rm c}$ is the comoving distance to time $t$ and $a(t)$ is the expansion at that time) for unit solid angle at the specified {\tt time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Numerical_Constants_Astronomical Numerical_Constants_Physical</modules>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#   <code>comovingVolumeElementTime=self%distanceComoving(time)**2*(gigaYear*speedLight/megaParsec)/self%expansionfactor(time)</code>
  !#  </method>
  !#  <method name="epochTime" >
  !#   <description>Convenience function that returns the time corresponding to an epoch specified by time or expansion factor.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time           , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsingPhase</argument>
  !#   <modules>Galacticus_Error</modules>
  !#   <code>
  !#   if (present(time)) then
  !#      epochTime=time
  !#   else if (present(expansionFactor)) then
  !#      epochTime=self%cosmicTime(expansionFactor,collapsingPhase)
  !#   else
  !#      call Galacticus_Error_Report('epochTime','either "time" or "expansionFactor" must be given')
  !#   end if
  !#   </code>
  !#  </method>
  include 'cosmologyFunctions.type.inc'
  !# </include>

end module Cosmology_Functions
