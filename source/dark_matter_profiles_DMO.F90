!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which provides an object that implements dark matter halo profiles.
!!}

module Dark_Matter_Profiles_DMO
  !!{
  Provides an object that implements dark matter halo profiles.
  !!}
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScale     , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_Generic, only : darkMatterProfileGeneric
  use :: Galacticus_Nodes            , only : treeNode
  private

  !![
  <functionClass>
   <name>darkMatterProfileDMO</name>
   <extends>darkMatterProfileGeneric</extends>
   <descriptiveName>Dark Matter Only Halo Profiles</descriptiveName>
   <description>
    Class providing dark matter-only halo profiles.
   </description>
   <default>NFW</default>
   <method name="density" >
    <description>Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="densityLogSlope" >
    <description>Returns the logarithmic slope of the density profile in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="radialMoment" >
    <description>Returns the {\normalfont \ttfamily m}$^\mathrm{th}$ radial moment of the dark matter profile of {\normalfont \ttfamily node} optionally between the given {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   )           :: moment</argument>
    <argument>double precision          , intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
   </method>
   <method name="energy" >
    <description>Return the total energy for the given {\normalfont \ttfamily node} in units of $M_\odot$ km$^2$ s$^{-2}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="rotationNormalization" >
    <description> Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for the given {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusFromSpecificAngularMomentum" >
    <description> Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} at which the specific angular momentum of a circular orbit equals {\normalfont \ttfamily specificAngularMomentum} (specified in units of km s$^{-1}$ Mpc.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: specificAngularMomentum</argument>
   </method>
   <method name="circularVelocity" >
    <description>Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="radiusCircularVelocityMaximum" >
    <description>Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="circularVelocityMaximum" >
    <description>Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radialVelocityDispersion" >
    <description>Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="potential" >
    <description>Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target   :: node</argument>
    <argument>double precision          , intent(in   )           :: radius</argument>
    <argument>integer                   , intent(  out), optional :: status</argument>
   </method>
   <method name="enclosedMass" >
    <description>Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="radiusEnclosingDensity" >
    <description>Returns the radius (in Mpc) enclosing a given density threshold (in $M_\odot \hbox{Mpc}^{-3}$) in the dark matter profile of {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: density</argument>
   </method>
   <method name="kSpace" >
    <description>Returns the normalized Fourier space density profile of the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily waveNumber} (given in units of Mpc$^{-1}$).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: wavenumber</argument>
   </method>
   <method name="freefallRadius" >
    <description>Returns the freefall radius (in Mpc) corresponding to the given {\normalfont \ttfamily time} (in Gyr) in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: time</argument>
   </method>
   <method name="freeFallRadiusIncreaseRate" >
    <description>Returns the rate of increase of the freefall radius (in Mpc/Gyr) corresponding to the given {\normalfont \ttfamily time} (in Gyr) in {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: time</argument>
   </method>
   <method name="radiusEnclosingMass" >
    <description>Returns the radius (in Mpc) enclosing a given mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type             (treeNode), intent(inout), target :: node</argument>
    <argument>double precision           , intent(in   )         :: mass</argument>
    <modules>Root_Finder Kind_Numbers</modules>
    <code>
     double precision                      :: radiusGuess
     type            (rootFinder), save    :: finder
     logical                     , save    :: finderConstructed=.false.
     double precision            , save    :: radiusPrevious   =-huge(0.0d0)
     integer         (kind_int8 ), save    :: uniqueIDPrevious =-1_kind_int8
     !$omp threadprivate(finder,finderConstructed,radiusPrevious,uniqueIDPrevious)
     if(mass &lt;= 0.0d0) then
        darkMatterProfileDMORadiusEnclosingMass=0.0d0
        return
     end if
     ! Initialize the root finder.
     if (.not.finderConstructed) then
        finder=rootFinder(                                                                    &amp;
             &amp;        rootFunction                 =darkMatterProfileDMOEnclosedMassRoot, &amp;
             &amp;        rangeExpandDownward          =0.5d0                               , &amp;
             &amp;        rangeExpandUpward            =2.0d0                               , &amp;
             &amp;        rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative       , &amp;
             &amp;        rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive       , &amp;
             &amp;        rangeExpandType              =rangeExpandMultiplicative           , &amp;
             &amp;        toleranceAbsolute            =0.0d0                               , &amp;
             &amp;        toleranceRelative            =1.0d-6                                &amp;
             &amp;       )
        finderConstructed=.true.
     end if
     if (node%uniqueID()     == uniqueIDPrevious) then
        radiusGuess     =radiusPrevious
     else
        radiusGuess     =self%darkMatterHaloScale_%radiusVirial(node)
        uniqueIDPrevious=node                     %uniqueID    (    )
     end if
     darkMatterProfileDMOSelf  => self
     darkMatterProfileDMONode  => node
     darkMatterProfileDMOMass_ =  mass
     radiusPrevious=finder%find(rootGuess=radiusGuess)
     darkMatterProfileDMORadiusEnclosingMass=radiusPrevious
     return
    </code>
   </method>
  </functionClass>
  !!]

  !![
  <functionClass>
   <name>darkMatterProfileHeating</name>
   <descriptiveName>Dark Matter Profile Heating</descriptiveName>
   <description>
    Class providing models of heating of dark matter profiles.
   </description>
   <default>null</default>
   <method name="specificEnergy" >
    <description>The specific energy of heating at the given {\normalfont \ttfamily radius} in the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode                 ), intent(inout) :: node</argument>
    <argument>double precision                           , intent(in   ) :: radius</argument>
    <argument>class           (darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_</argument>
   </method>
   <method name="specificEnergyGradient" >
    <description>The gradient of the specific energy of heating at the given {\normalfont \ttfamily radius} in the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode                 ), intent(inout) :: node</argument>
    <argument>double precision                           , intent(in   ) :: radius</argument>
    <argument>class           (darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_</argument>
   </method>
   <method name="specificEnergyIsEverywhereZero" >
    <description>Returns true if the specific energy is zero everywhere in the given {\normalfont \ttfamily node}.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type (treeNode                 ), intent(inout) :: node</argument>
    <argument>class(darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_</argument>
   </method>
  </functionClass>
  !!]

  ! Module-scope variables used in root finding.
  class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMOSelf   => null()
  type            (treeNode                 ), pointer :: darkMatterProfileDMONode   => null()
  double precision                                     :: darkMatterProfileDMOMass_
  !$omp threadprivate(darkMatterProfileDMOSelf,darkMatterProfileDMONode,darkMatterProfileDMOMass_)

contains

  double precision function darkMatterProfileDMOEnclosedMassRoot(radius)
    !!{
    Root function used in solving for the radius that encloses a given mass.
    !!}
    implicit none
    double precision,intent(in) :: radius

    darkMatterProfileDMOEnclosedMassRoot=darkMatterProfileDMOSelf%enclosedMass(darkMatterProfileDMONode,radius)-darkMatterProfileDMOMass_
    return
  end function darkMatterProfileDMOEnclosedMassRoot

end module Dark_Matter_Profiles_DMO
