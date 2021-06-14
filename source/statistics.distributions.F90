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

!% Contains a module that implements a class of distributions.

module Statistics_Distributions
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  private

  !# <functionClass>
  !#  <name>distributionFunction1D</name>
  !#  <descriptiveName>One-dimensional Distribution Functions</descriptiveName>
  !#  <description>Class providing distribution functions of a single variable.</description>
  !#  <default>uniform</default>
  !#  <data>class(randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()</data>
  !#  <destructor>
  !#   <code>
  !#    call distributionFunction1DFinalize(self)
  !#    return
  !#   </code>
  !#  </destructor>
  !#  <method name="density" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>double precision, intent(in   ) :: x</argument>
  !#    <description>Return the probability density at {\normalfont \ttfamily x}.</description>
  !#  </method>
  !#  <method name="cumulative" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>double precision, intent(in   ) :: x</argument>
  !#    <description>Return the cumulative probability at {\normalfont \ttfamily x}.</description>
  !#  </method>
  !#  <method name="inverse" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <selfTarget>yes</selfTarget>
  !#    <argument>double precision, intent(in   ) :: p</argument>
  !#    <description>Return the value of the independent variable corresponding to cumulative probability {\normalfont \ttfamily p}.</description>
  !#  </method>
  !#  <method name="sample" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>class(randomNumberGeneratorClass), intent(inout), optional :: randomNumberGenerator_</argument>
  !#    <description>Return a random deviate from the distribution.</description>
  !#    <modules>Galacticus_Error</modules>
  !#    <code>
  !#     double precision :: uniformRandom
  !#     ! Draw a random number uniformly from 0 to 1 and use the inverse of our self to get the
  !#     ! corresponding random variate.
  !#     if (present(randomNumberGenerator_)) then
  !#        uniformRandom=     randomNumberGenerator_%uniformSample()
  !#     else
  !#        if (.not.associated(self%randomNumberGenerator_)) call Galacticus_Error_Report('no random number generator supplied'//{introspection:location})
  !#        uniformRandom=self%randomNumberGenerator_%uniformSample()
  !#     end if
  !#     distributionFunction1DSample=self%inverse(uniformRandom)
  !#     return
  !#    </code>
  !#  </method>
  !#  <method name="minimum" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <description>Returns the minimum possible value in the distribution.</description>
  !#  </method>
  !#  <method name="maximum" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <description>Returns the maximum possible value in the distribution.</description>
  !#  </method>
  !# </functionClass>

  ! Define a list of distributions.
  type, public :: distributionFunction1DList
     class(distributionFunction1DClass), pointer :: distributionFunction1D_
  end type distributionFunction1DList

contains
  
  subroutine distributionFunction1DFinalize(self)
    !% Finalizer for {\normalfont \ttfamily distributionFunction1D} objects.
    type(distributionFunction1DClass), intent(inout) :: self

    !# <objectDestructor name="self%randomNumberGenerator_"/>
    return
  end subroutine distributionFunction1DFinalize
  
end module Statistics_Distributions
