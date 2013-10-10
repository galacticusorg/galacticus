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

type, public, extends(massDistributionSpherical) :: massDistributionBetaProfile
   !% The $\beta$-profile: $\rho(r)=\rho_0/[1+(r/r_{\rm core})^2]^{3\beta/2}$
   double precision :: beta, coreRadius, densityNormalization
 contains
   !@ <objectMethods>
   !@   <object>massDistributionBetaProfile</object>
   !@   <objectMethod>
   !@     <method>initialize</method>
   !@     <description>Initialize the mass distribution.</description>
   !@     <type>\void</type>
   !@     <arguments>\doublezero beta\argin, \doublezero coreRadius\argin, \doublezero densityNormalization\argin, \doublezero mass\argin, \doublezero outerRadius\argin, \logicalzero\ isDimensionless\argin</arguments>
   !@   </objectMethod>
   !@ </objectMethods>
   procedure :: initialize=>Mass_Distribution_Beta_Profile_Initialize
   procedure :: density   =>Mass_Distribution_Beta_Profile_Density
end type massDistributionBetaProfile
