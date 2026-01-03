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
Contains a module which implements a class for stellar initial mass functions.
!!}

module Stellar_Populations_Initial_Mass_Functions
  !!{
  Implements a class for stellar initial mass functions.
  !!}
  use :: Tables, only : table1D
  implicit none
  private

  !![
  <functionClass>
   <name>initialMassFunction</name>
   <descriptiveName>Initial Mass Functions</descriptiveName>
   <description>
    Class providing stellar initial mass functions. All IMFs are assumed to be continuous in $M$, unless otherwise noted and
    normalized to unit mass.
   </description>
   <default>chabrier2001</default>
   <method name="massMinimum" >
    <description>Return the minimum mass in the initial mass function.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="massMaximum" >
    <description>Return the maximum mass in the initial mass function.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="phi" >
    <description>Return the initial mass function, $\phi(M)=\mathrm{d}N/\mathrm{d}M$, at the given mass $M=${\normalfont \ttfamily massInitial}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massInitial</argument>
   </method>
   <method name="numberCumulative" >
    <description>Return the integral of the initial mass function, $\int_{m_\mathrm{lower}}^{m_\mathrm{upper}} \phi(M) \mathrm{d}M$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massLower, massUpper</argument>
   </method>
   <method name="tabulate" >
    <description>Return the initial mass function, $\phi(M)=\mathrm{d}N/\mathrm{d}M$, at the given mass $M=${\normalfont \ttfamily initialMass}.</description>
    <type>void</type>
    <argument>class(table1D), allocatable, intent(inout) :: imfTable</argument>
    <pass>yes</pass>
   </method>
   <method name="label" >
    <description>Return the label for this \gls{imf}.</description>
    <type>type(varying_string)</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Stellar_Populations_Initial_Mass_Functions
