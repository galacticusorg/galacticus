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

!% Contains a module that provides a class implementing cooling functions.

module Cooling_Functions
  !% Provides a class implementing cooling functions.
  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Radiation_Fields             , only : radiationFieldClass
  private

  !# <functionClass>
  !#  <name>coolingFunction</name>
  !#  <descriptiveName>Cooling Function</descriptiveName>
  !#  <description>
  !#   Class providing cooling functions of gas, $\Lambda(\rho,T,\mathbf{Z})$.
  !#  </description>
  !#  <default>atomicCIECloudy</default>
  !#  <method name="coolingFunction" >
  !#   <description>Return the cooling function at the given temperature and hydrogen density for the specified set of abundances and radiation field. Units of the returned cooling function are the traditional ergs cm$^-3$ s$^{-1}$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
  !#   <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
  !#   <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities</argument>
  !#   <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
  !#  </method>
  !#  <method name="coolingFunctionDensityLogSlope" >
  !#   <description>Return $\d\ln\Lambda/\d\ln\rho$ for a cooling function at the given temperature and hydrogen density for the specified set of abundances and radiation field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
  !#   <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
  !#   <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities</argument>
  !#   <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
  !#  </method>
  !#  <method name="coolingFunctionTemperatureLogSlope" >
  !#   <description>Return $\d\ln\Lambda/\d\ln T$ for a cooling function at the given temperature and hydrogen density for the specified set of abundances and radiation field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision                     , intent(in   ) :: numberDensityHydrogen, temperature</argument>
  !#   <argument>type            (abundances         ), intent(in   ) :: gasAbundances</argument>
  !#   <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities</argument>
  !#   <argument>class           (radiationFieldClass), intent(inout) :: radiation</argument>
  !#  </method>
  !# </functionClass>

end module Cooling_Functions
