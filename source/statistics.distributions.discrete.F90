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
Contains a module that implements a class of discrete distributions.
!!}

module Statistics_Distributions_Discrete
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  private

  !![
  <functionClass>
   <name>distributionFunctionDiscrete1D</name>
   <descriptiveName>One-dimensional Discrete Distribution Functions</descriptiveName>
   <description>Class providing discrete probability distribution functions of a single integer variable---the probability mass
    function $p(x)$ (and its logarithm), the cumulative distribution function $P(x) = \sum_{x' \le x} p(x')$, and the quantile
    function $x(P)$. These distributions model count data such as the number of galaxies in a halo or the number of star formation
    events, and are used for drawing random variates and computing Poisson or binomial likelihoods in galaxy statistics and N-body
    halo occupation analyses.</description>
   <default>binomial</default>
   <destructor>
    <code>
     call distributionFunctionDiscrete1DFinalize(self)
     return
    </code>
   </destructor>
   <method name="mass" >
     <description>Return the probability mass function $p(x)$, giving the probability that the discrete random variable takes the integer value \mono{x}.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>integer, intent(in   ) :: x</argument>
   </method>
   <method name="massLogarithmic" >
     <description>Return the natural logarithm of the probability mass function $\ln p(x)$ evaluated at integer \mono{x}, which is more numerically stable for extremely small probabilities than computing $p(x)$ directly.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>integer, intent(in   ) :: x</argument>
   </method>
   <method name="cumulative" >
     <description>Return the cumulative distribution function $P(x) = \sum_{x' \le x} p(x')$, giving the probability that the discrete random variable takes a value less than or equal to integer \mono{x}.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>integer, intent(in   ) :: x</argument>
   </method>
   <method name="inverse" >
     <description>Return the value of the independent variable corresponding to cumulative probability \mono{p}.</description>
     <type>integer</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: p</argument>
   </method>
   <method name="sample" >
     <description>Return a random integer deviate drawn from this discrete probability distribution, using the inverse CDF method by default (drawing a uniform random number and applying the quantile function).</description>
     <type>integer</type>
     <pass>yes</pass>
     <modules>Error</modules>
     <argument>class(randomNumberGeneratorClass), intent(inout), optional :: randomNumberGenerator_</argument>
     <code>
      double precision :: uniformRandom
      ! Draw a random number uniformly from 0 to 1 and use the inverse of our self to get the
      ! corresponding random variate.
      if (present(randomNumberGenerator_)) then
         uniformRandom=     randomNumberGenerator_%uniformSample()
      else
         if (.not.associated(self%randomNumberGenerator_)) call Error_Report('no random number generator supplied'//{introspection:location})
         uniformRandom=self%randomNumberGenerator_%uniformSample()
      end if
      distributionFunctionDiscrete1DSample=self%inverse(uniformRandom)
      return
     </code>
   </method>
   <method name="minimum" >
     <description>Returns the minimum possible integer value in the support of this discrete distribution, i.e., the smallest integer $x$ for which the probability mass is non-zero.</description>
     <type>integer</type>
     <pass>yes</pass>
   </method>
   <method name="maximum" >
     <description>Returns the maximum possible integer value in the support of this discrete distribution, i.e., the largest integer $x$ for which the probability mass is non-zero.</description>
     <type>integer</type>
     <pass>yes</pass>
   </method>
   <data>class(randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()</data>
  </functionClass>
  !!]

contains
  
  subroutine distributionFunctionDiscrete1DFinalize(self)
    !!{
    Finalizer for \mono{distributionFunctionDiscrete1D} objects.
    !!}
    type(distributionFunctionDiscrete1DClass), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine distributionFunctionDiscrete1DFinalize

end module Statistics_Distributions_Discrete
