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
Contains a module which provides a class implementing solvers for collapse of spherical perturbations.
!!}

module Spherical_Collapse_Solvers
  !!{
  Provides a class implementing solvers for collapse of spherical perturbations.
  !!}
  use :: Tables, only : table1D, table2DLinLinLin
  private

  !![
  <functionClass>
   <name>sphericalCollapseSolver</name>
   <descriptiveName>Spherical Collapse Solvers</descriptiveName>
   <description>Class providing numerical solvers for the collapse of spherical density perturbations in a given
    cosmology---computing the critical linear overdensity $\delta_\mathrm{c}$, the virial density
    contrast $\Delta_\mathrm{vir}$, and the turnaround-to-virial radius ratio as functions of cosmic
    time by integrating the spherical collapse equations. These tabulated quantities are consumed by
    \refClass{criticalOverdensityClass} and \refClass{virialDensityContrastClass} implementations,
    allowing the halo collapse threshold and virial radius definition to vary with redshift and
    cosmological model.</description>
   <default>cllsnlssMttrCsmlgclCnstnt</default>
   <method name="criticalOverdensity" >
    <description>Returns a tabulation of the critical overdensity for collapse.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                      , intent(in   ) :: time                </argument>
    <argument>logical                               , intent(in   ) :: tableStore          </argument>
    <argument>class           (table1D), allocatable, intent(inout) :: criticalOverdensity_</argument>
   </method>
   <method name="virialDensityContrast" >
    <description>Returns a tabulation of the virial density contrast.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                      , intent(in   ) :: time                  </argument>
    <argument>logical                               , intent(in   ) :: tableStore            </argument>
    <argument>class           (table1D), allocatable, intent(inout) :: virialDensityContrast_</argument>
   </method>
   <method name="radiusTurnaround" >
    <description>Returns a tabulation of the ratio of turnaround to virial radii.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                      , intent(in   ) :: time               </argument>
    <argument>logical                               , intent(in   ) :: tableStore         </argument>
    <argument>class           (table1D), allocatable, intent(inout) :: radiusTurnaround_  </argument>
   </method>
   <method name="linearNonlinearMap" >
    <description>Returns a mapping of linear to nonlinear density overdensity.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                               , intent(in   ) :: time               </argument>
    <argument>class           (table2DLinLinLin), allocatable, intent(inout) :: linearNonlinearMap_</argument>
   </method>
  </functionClass>
  !!]

end module Spherical_Collapse_Solvers
