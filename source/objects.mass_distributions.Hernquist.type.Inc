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

type, public, extends(massDistributionSpherical) :: massDistributionHernquist
   !% The Hernquist \citep{hernquist_analytical_1990} density profile.
   double precision  :: scaleLength,densityNormalization,mass
 contains
   procedure :: initialize           => Mass_Distribution_Hernquist_Initialize
   procedure :: density              => Mass_Distribution_Hernquist_Density
   procedure :: densityRadialMoment  => Mass_Distribution_Hernquist_Density_Radial_Moment
   procedure :: massEnclosedBySphere => Mass_Distribution_Hernquist_Mass_Enc_By_Sphere
   procedure :: potential            => Mass_Distribution_Hernquist_Potential
   procedure :: halfMassRadius       => Mass_Distribution_Hernquist_Half_Mass_Radius
end type massDistributionHernquist
  
