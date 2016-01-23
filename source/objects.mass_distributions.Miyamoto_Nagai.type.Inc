!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

type, public, extends(massDistributionCylindrical) :: massDistributionMiyamotoNagai
   !% The Miyamoto-Nagai model \citep{miyamoto_three-dimensional_1975}.
   double precision                           :: a                        , b                  , &
        &                                        densityNormalization     , mass               , &
        &                                        shape                    , halfMassRadiusValue
   logical                                    :: surfaceDensityInitialized
   type            (table1DLogarithmicLinear) :: surfaceDensityTable
   logical                                    :: massEnclosedInitialized
   type            (table1DLogarithmicLinear) :: massEnclosedTable
 contains
   !@ <objectMethods>
   !@   <object>massDistributionMiyamotoNagai</object>
   !@   <objectMethod>
   !@     <method>initialize</method>
   !@     <description>Initialize the mass distribution.</description>
   !@     <type>\void</type>
   !@     <arguments>\doublezero a\argin, \doublezero b\argin, \doublezero mass\argin, \logicalzero\ isDimensionless\argin</arguments>
   !@   </objectMethod>
   !@   <objectMethod>
   !@     <method>surfaceDensityInitialize</method>
   !@     <description>Initialize the surface density tabulation.</description>
   !@     <type>\void</type>
   !@     <arguments></arguments>
   !@   </objectMethod>
   !@   <objectMethod>
   !@     <method>massEnclosedInitialize</method>
   !@     <description>Initialize the enclosed mass tabulation.</description>
   !@     <type>\void</type>
   !@     <arguments></arguments>
   !@   </objectMethod>
   !@ </objectMethods>
   procedure :: initialize                 => Mass_Distribution_Miyamoto_Nagai_Initialize
   procedure :: density                    => Mass_Distribution_Miyamoto_Nagai_Density
   procedure :: surfaceDensity             => Mass_Distribution_Miyamoto_Nagai_Surface_Density
   procedure :: surfaceDensityInitialize   => Mass_Distribution_Miyamoto_Nagai_Surface_Density_Initialize
   procedure :: surfaceDensityRadialMoment => Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt
   procedure :: massEnclosedBySphere       => Mass_Distribution_Miyamoto_Nagai_Mass_Enc_By_Sphere
   procedure :: massEnclosedInitialize     => Mass_Distribution_Miyamoto_Nagai_Mass_Enclosed_Initialize
   procedure :: potential                  => Mass_Distribution_Miyamoto_Nagai_Potential
   procedure :: rotationCurve              => Mass_Distribution_Miyamoto_Nagai_Rotation_Curve
   procedure :: rotationCurveGradient      => Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient
   procedure :: halfMassRadius             => Mass_Distribution_Miyamoto_Nagai_Half_Mass_Radius
end type massDistributionMiyamotoNagai

abstract interface
   double precision function integrandWeight(x)
     double precision, intent(in   ) :: x
   end function integrandWeight
end interface
