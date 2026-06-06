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
Contains a module which implements a class for calculations of Population III supernovae.
!!}

module Supernovae_Population_III
  !!{
  Implements a class for calculations of Population III supernovae.
  !!}
  implicit none
  private

  !![
  <functionClass>
   <name>supernovaePopulationIII</name>
   <descriptiveName>Population III Supernovae</descriptiveName>
   <description>Class providing models of supernovae from Population III (metal-free) stars---the first generation
    of stars formed from primordial gas at high redshift. Very massive Population III stars ($M \gtrsim 100
    \mathrm{M}_\odot$) may explode as pair-instability supernovae with enormous energy output ($\sim 10^{53}$ erg), while
    lower-mass Population III stars undergo core-collapse. The cumulative energy input from Population III
    supernovae contributes to the reionization and early metal enrichment of the \gls{igm}.</description>
   <default>hegerWoosley2002</default>
   <method name="energyCumulative" >
    <description> Return the cumulative energy input from Population III supernovae from stars of given \mono{initialMass}, \mono{age} and \mono{metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: initialMass, age, metallicity</argument>
   </method>
  </functionClass>
  !!]

end module Supernovae_Population_III
