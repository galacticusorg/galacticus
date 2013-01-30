!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

type, public, extends(massDistributionSpherical) :: massDistributionSersic
   !% The S\'ersic density profile.
   double precision  :: radiusHalfMass,densityNormalization,mass
 contains
   procedure :: initialize              => Mass_Distribution_Sersic_Initialize
   procedure :: density                 => Mass_Distribution_Sersic_Density
   procedure :: densityRadialMoment     => Mass_Distribution_Sersic_Density_Radial_Moment
   procedure :: massEnclosedBySphere    => Mass_Distribution_Sersic_Mass_Enc_By_Sphere
   procedure :: potential               => Mass_Distribution_Sersic_Potential
   procedure :: halfMassRadius          => Mass_Distribution_Sersic_Half_Mass_Radius
   procedure :: halfMassRadiusProjected => Mass_Distribution_Sersic_Half_Mass_Radius_Projected
end type massDistributionSersic

! Public functions.
public :: Sersic_Profile_Tabulate_State_Retrieve,Sersic_Profile_Tabulate_State_Store

! Tabulation of the Sérsic profile.
logical                                                        :: sersicProfileTabulated=.false.
double precision                                               :: sersicProfileIndex,sersicCoefficient,sersicRadiusStart
logical                                                        :: sersicTableInitialized=.false.
integer                                                        :: sersicTableCount
double precision                                               :: sersicTableRadiusMinimum   =1.0d-3, &
     &                                                            sersicTableRadiusMaximum   =1.0d+3
double precision                                               :: sersicTable3dHalfMassRadius=1.0d+0
double precision                                               :: sersicTable2dHalfMassRadius
integer                            , parameter                 :: sersicTablePointsPerDecade=1000
double precision                   , allocatable, dimension(:) :: sersicTableRadius,sersicTableDensity,sersicTableEnclosedMass,sersicTablePotential
type            (fgsl_interp      )                            :: sersicTableInterpolationObject
type            (fgsl_interp_accel)                            :: sersicTableInterpolationAccelerator
logical                                                        :: sersicTableInterpolationReset=.true.
