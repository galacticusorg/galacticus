!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that provides an object that implements incompleteness calculations for observed mass functions.

module Mass_Function_Incompletenesses
  !% Provides an object that implements incompleteness calculations for observed mass functions.
  use ISO_Varying_String
  !# <include directive="massFunctionIncompleteness" type="functionModules" >
  include 'massFunctionIncompleteness.functionModules.inc'
  !# </include>
  private

  !# <include directive="massFunctionIncompleteness" type="function" >
  !#  <descriptiveName>Mass Function Incompletenesses</descriptiveName>
  !#  <description>Object providing incompleteness calculations for observed mass functions.</description>
  !#  <default>complete</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="completeness" >
  !#   <description>Return the completeness of the observational sample at the given mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#  </method>
  include 'massFunctionIncompleteness.type.inc'
  !# </include>

end module Mass_Function_Incompletenesses
