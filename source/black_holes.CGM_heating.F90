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
Contains a module which implements a class for black hole heating of the \gls{cgm}.
!!}

module Black_Hole_CGM_Heating
  !!{
  Implements a class for black hole heating of the \gls{cgm}.
  !!}
  use :: Galacticus_Nodes, only : nodeComponentBlackHole
  implicit none
  private

  !![
  <functionClass>
   <name>blackHoleCGMHeating</name>
   <descriptiveName>Black Hole CGM Heating</descriptiveName>
   <description>Class providing models of the heating rate (in $\mathrm{M}_\odot$ (km/s)$^2$ Gyr$^{-1}$) delivered by an active galactic nucleus
    to the \gls{cgm} of its host halo. AGN heating of the \gls{cgm} is a key mechanism for quenching star
    formation in massive galaxies (``radio-mode'' or ``maintenance-mode'' feedback), counterbalancing radiative
    cooling of the hot gas atmosphere. Implementations may couple the jet power, the wind power, or a fixed
    fraction of the AGN bolometric luminosity to the surrounding hot gas.</description>
   <default>jetPower</default>
   <method name="heatingRate" >
    <description>Compute the heating rate of the CGM due to the given \mono{blackHole}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(nodeComponentBlackHole), intent(inout) :: blackHole</argument>
   </method>
  </functionClass>
  !!]

end module Black_Hole_CGM_Heating
