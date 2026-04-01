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
Contains a module that implements a class of parameter mapping functions.
!!}

module Math_Operators_Unary
  !!{
  Implements a class of unary operators.
  !!}
  private

  !![
  <functionClass>
   <name>operatorUnary</name>
   <descriptiveName>Unary Operators</descriptiveName>
   <description>Class providing unary operators---invertible scalar mappings $f: \mathbb{R} \to \mathbb{R}$
    that transform a value and can be reversed. Common examples include the identity, logarithm,
    and various monotonic reparametrizations. These operators are used in parameter estimation and
    output analysis to transform model parameters or property values before comparison with
    observations (e.g.\ converting between linear and logarithmic scales).</description>
   <default>identity</default>
   <method name="operate" >
     <description>Apply the unary operator to the scalar input \mono{x}, returning the transformed value $f(x)$.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
   </method>
   <method name="unoperate" >
     <description>Reverse the unary operation by applying the inverse mapping $f^{-1}$ to the scalar input \mono{f}, returning the original value $x$ such that operate$(x) = f$.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: f</argument>
   </method>
  </functionClass>
  !!]

end module Math_Operators_Unary
