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

type, public, extends(massDistributionCylindrical) :: massDistributionExponentialDisk
   !% The exponential disk profile: $\rho(r,z)=\rho_0 \exp(-r/r_{\rm s}) \hbox{sech}^2(z/z_{\rm s})$.
   double precision :: scaleRadius         , scaleHeight                , &
        &              densityNormalization, surfaceDensityNormalization, &
        &              mass
   ! Tables used for rotation curves.
   logical                                    :: scaleLengthFactorSet                  =.false., rotationCurveInitialized              =.false.   , &
        &                                        rotationCurveGradientInitialized      =.false.
   double precision                           :: scaleLengthFactor
   double precision                           :: rotationCurveHalfRadiusMinimum                , rotationCurveHalfRadiusMaximum
   double precision                           :: rotationCurveGradientHalfRadiusMinimum        , rotationCurveGradientHalfRadiusMaximum
   type            (table1DLogarithmicLinear) :: rotationCurveTable                            , rotationCurveGradientTable
 contains
   !@ <objectMethods>
   !@   <object>massDistributionExponentialDisk</object>
   !@   <objectMethod>
   !@     <method>initialize</method>
   !@     <description>Initialize the mass distribution.</description>
   !@     <type>\void</type>
   !@     <arguments>\doublezero scaleRadius\argin, \doublezero scaleHeight\argin, \doublezero mass\argin, \logicalzero\ isDimensionless\argin</arguments>
   !@   </objectMethod>
   !@ </objectMethods>
   procedure :: initialize                 => Mass_Distribution_Exponential_Disk_Initialize
   procedure :: density                    => Mass_Distribution_Exponential_Disk_Density
   procedure :: surfaceDensity             => Mass_Distribution_Exponential_Disk_Surface_Density
   procedure :: massEnclosedBySphere       => Mass_Distribution_Exponential_Disk_Mass_Enc_By_Sphere
   procedure :: potential                  => Mass_Distribution_Exponential_Disk_Potential
   procedure :: rotationCurve              => Mass_Distribution_Exponential_Disk_Rotation_Curve
   procedure :: rotationCurveGradient      => Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient
   procedure :: halfMassRadius             => Mass_Distribution_Exponential_Disk_Half_Mass_Radius
   procedure :: surfaceDensityRadialMoment => Mass_Distribution_Exponential_Disk_Srfc_Dnsty_Rdl_Mmnt
   procedure :: stateStore                 => Mass_Distribution_Exponential_Disk_State_Store
   procedure :: stateRestore               => Mass_Distribution_Exponential_Disk_State_Restore
end type massDistributionExponentialDisk

! The radius (in units of the disk scale length) beyond which the disk is treated as a point mass for the purposes of computing
! rotation curves.
double precision, parameter :: radiusMaximum=30.0d0
