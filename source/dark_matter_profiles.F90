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

!% Contains a module which provides an object that implements non-dark-matter-only dark matter halo profiles.

module Dark_Matter_Profiles
  !% Provides an object that implements non-dark-matter-only dark matter halo profiles.
  use :: Dark_Matter_Profiles_Generic, only : darkMatterProfileGeneric
  use :: Galacticus_Nodes            , only : treeNode
  use :: Kind_Numbers                , only : kind_int8
  private

  !# <functionClass>
  !#  <name>darkMatterProfile</name>
  !#  <extends>darkMatterProfileGeneric</extends>
  !#  <descriptiveName>Dark Matter Halo Profiles</descriptiveName>
  !#  <description>Object providing dark matter halo profiles.</description>
  !#  <default>adiabaticGnedin2004</default>
  !#  <method name="density" >
  !#   <description>Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="densityLogSlope" >
  !#   <description>Returns the logarithmic slope of the density profile in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="radialMoment" >
  !#   <description>Returns the {\normalfont \ttfamily m}$^\mathrm{th}$ radial moment of the dark matter profile of {\normalfont \ttfamily node} optionally between the given {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout)           :: node</argument>
  !#   <argument>double precision          , intent(in   )           :: moment</argument>
  !#   <argument>double precision          , intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
  !#  </method>
  !#  <method name="energy" >
  !#   <description>Return the total energy for the given {\normalfont \ttfamily node} in units of $M_\odot$ km$^2$ s$^{-2}$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="energyGrowthRate" >
  !#   <description>Return the rate of change of total energy for the given {\normalfont \ttfamily node} in units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="rotationNormalization" >
  !#   <description> Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for the given {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="radiusFromSpecificAngularMomentum" >
  !#   <description> Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} at which the specific angular momentum of a circular orbit equals {\normalfont \ttfamily specificAngularMomentum} (specified in units of km s$^{-1}$ Mpc.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: specificAngularMomentum</argument>
  !#  </method>
  !#  <method name="circularVelocity" >
  !#   <description>Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="radiusCircularVelocityMaximum" >
  !#   <description>Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="circularVelocityMaximum" >
  !#   <description>Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="radialVelocityDispersion" >
  !#   <description>Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="potential" >
  !#   <description>Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type            (treeNode), intent(inout), target   :: node</argument>
  !#   <argument>double precision          , intent(in   )           :: radius</argument>
  !#   <argument>integer                   , intent(  out), optional :: status</argument>
  !#  </method>
  !#  <method name="enclosedMass" >
  !#   <description>Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="radiusEnclosingDensity" >
  !#   <description>Returns the radius (in Mpc) enclosing a given density threshold (in $M_\odot \hbox{Mpc}^{-3}$) in the dark matter profile of {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type            (treeNode), intent(inout), target :: node</argument>
  !#   <argument>double precision          , intent(in   )         :: density</argument>
  !#  </method>
  !#  <method name="kSpace" >
  !#   <description>Returns the normalized Fourier space density profile of the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily waveNumber} (given in units of Mpc$^{-1}$).</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), target :: node</argument>
  !#   <argument>double precision          , intent(in   )         :: wavenumber</argument>
  !#  </method>
  !#  <method name="freefallRadius" >
  !#   <description>Returns the freefall radius (in Mpc) corresponding to the given {\normalfont \ttfamily time} (in Gyr) in {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="freeFallRadiusIncreaseRate" >
  !#   <description>Returns the rate of increase of the freefall radius (in Mpc/Gyr) corresponding to the given {\normalfont \ttfamily time} (in Gyr) in {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="radiusEnclosingMass" >
  !#   <description>Returns the radius (in Mpc) enclosing a given mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type             (treeNode), intent(inout), target :: node</argument>
  !#   <argument>double precision           , intent(in   )         :: mass</argument>
  !#   <modules>Root_Finder Kind_Numbers</modules>
  !#   <code>
  !#    double precision                      :: radiusGuess
  !#    type            (rootFinder), save    :: finder
  !#    logical                     , save    :: finderConstructed=.false.
  !#    double precision            , save    :: radiusPrevious   =-huge(0.0d0)
  !#    integer         (kind_int8 ), save    :: uniqueIDPrevious =-1_kind_int8
  !#    !$omp threadprivate(finder,finderConstructed,radiusPrevious,uniqueIDPrevious)
  !#    if(mass &lt;= 0) then
  !#       darkMatterProfileRadiusEnclosingMass=0.0d0
  !#       return
  !#    end if
  !#    ! Initialize the root finder.
  !#    if (.not.finderConstructed) then
  !#       finder=rootFinder(                                                             &amp;
  !#            &amp;        rootFunction                 =enclosedMassRoot             , &amp;
  !#            &amp;        rangeExpandDownward          =0.5d0                        , &amp;
  !#            &amp;        rangeExpandUpward            =2.0d0                        , &amp;
  !#            &amp;        rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
  !#            &amp;        rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &amp;
  !#            &amp;        rangeExpandType              =rangeExpandMultiplicative    , &amp;
  !#            &amp;        toleranceAbsolute            =0.0d+0                       , &amp;
  !#            &amp;        toleranceRelative            =1.0d-6                         &amp;
  !#            &amp;       )
  !#       finderConstructed=.true.
  !#    end if
  !#    if (node%uniqueID()     == uniqueIDPrevious) then
  !#       radiusGuess     =radiusPrevious
  !#    else
  !#       radiusGuess     =self%darkMatterHaloScale_%virialRadius(node)
  !#       uniqueIDPrevious=node                     %uniqueID    (    )
  !#    end if
  !#    radiusPrevious=finder%find(rootGuess=radiusGuess)
  !#    darkMatterProfileRadiusEnclosingMass=radiusPrevious
  !#    return
  !#    contains
  !#      double precision function enclosedMassRoot(radius)
  !#        !% Root function used in solving for the radius that encloses a given mass.
  !#        implicit none
  !#        double precision, intent(in) :: radius
  !#        enclosedMassRoot=self%enclosedMass(node,radius)-mass
  !#        return
  !#      end function enclosedMassRoot
  !#   </code>
  !#  </method>
  !# </functionClass>

end module Dark_Matter_Profiles
