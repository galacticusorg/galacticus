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
Contains a module which provides a class that implements random number generators.
!!}

module Numerical_Random_Numbers
  !!{
  Provides a class that implements random number generators.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_long
  private

  !![
  <functionClass>
   <name>randomNumberGenerator</name>
   <descriptiveName>Random Number Generators</descriptiveName>
   <description>Class providing random number generators.</description>
   <default>GSL</default>
   <method name="mpiIndependent" >
    <description>Return true if this random number generator produces independent sequences per MPI process (when using the same seed and offsetting is requested).</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="openMPIndependent" >
    <description>Return true if this random number generator produces independent sequences per OpenMP thread (when using the same seed and offsetting is requested).</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="rangeMinimum" >
    <description>Return the smallest integer in the range for this generator.</description>
    <type>integer(c_long)</type>
    <pass>yes</pass>
   </method>
   <method name="rangeMaximum" >
    <description>Return the largest integer in the range for this generator.</description>
    <type>integer(c_long)</type>
    <pass>yes</pass>
   </method>
   <method name="sample" >
    <description>Return a random integer. If the optional argument {\normalfont \ttfamily n} is supplied the random integer will lie in the range 0 to  {\normalfont \ttfamily n}-1 inclusive (with all integers being equally likely). Otherwise, the integer will be drawn from the full range provided by the generator.</description>
    <type>integer(c_long)</type>
    <pass>yes</pass>
    <argument>integer(c_long), intent(in   ), optional :: n</argument>
   </method>
   <method name="uniformSample" >
    <description>Return a random number drawn from a uniform distribution on [0,1).</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="poissonSample" >
    <description>Return a random number drawn from a Poisson distribution with the given {\normalfont \ttfamily mean}.</description>
    <type>integer</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mean</argument>
   </method>
   <method name="standardNormalSample" >
    <description>Return a random number drawn from a standard normal distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="seedSet" >
    <description>Reset the seed for this random number generator.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer(c_long), intent(in   ) :: seed</argument>
    <argument>logical        , intent(in   ) :: offset</argument>
   </method>
   <method name="seed" >
    <description>Return the seed for this random number generator.</description>
    <type>integer(c_long)</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]
  
end module Numerical_Random_Numbers
