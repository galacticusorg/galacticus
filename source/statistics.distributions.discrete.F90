!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module that implements a class of discrete distributions.

module Statistics_Distributions_Discrete
  use Pseudo_Random
  private

  !# <functionClass>
  !#  <name>distributionFunctionDiscrete1D</name>
  !#  <descriptiveName>One-dimensional Discrete Distribution Functions</descriptiveName>
  !#  <description>Class providing discrete distribution functions of a single variable.</description>
  !#  <default>binomial</default>
  !#  <data>type(pseudoRandom) :: randomNumberGenerator</data>
  !#  <method name="mass" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>integer, intent(in   ) :: x</argument>
  !#    <description>Return the probability mass at {\normalfont \ttfamily x}.</description>
  !#  </method>
  !#  <method name="massLogarithmic" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>integer, intent(in   ) :: x</argument>
  !#    <description>Return the logarithm of the probability mass at {\normalfont \ttfamily x}.</description>
  !#  </method>
  !#  <method name="cumulative" >
  !#    <type>double precision</type>
  !#    <pass>yes</pass>
  !#    <argument>integer, intent(in   ) :: x</argument>
  !#    <description>Return the cumulative probability at {\normalfont \ttfamily x}.</description>
  !#  </method>
  !#  <method name="inverse" >
  !#    <type>integer</type>
  !#    <pass>yes</pass>
  !#    <argument>double precision, intent(in   ) :: p</argument>
  !#    <description>Return the value of the independent variable corresponding to cumulative probability {\normalfont \ttfamily p}.</description>
  !#  </method>
  !#  <method name="sample" >
  !#    <type>integer</type>
  !#    <pass>yes</pass>
  !#    <argument>integer                         , intent(in   ), optional :: incrementSeed</argument>
  !#    <argument>logical                         , intent(in   ), optional :: ompThreadOffset, mpiRankOffset</argument>
  !#    <argument>type            (pseudoRandom  ), intent(inout), optional :: randomNumberGenerator</argument>
  !#    <description>Return a random deviate from the distribution.</description>
  !#    <code>
  !#     logical          :: ompThreadOffset_, mpiRankOffset_
  !#     double precision :: uniformRandom    
  !#     ompThreadOffset_=.true.
  !#     mpiRankOffset_  =.true.
  !#     if (present(ompThreadOffset)) ompThreadOffset_=ompThreadOffset
  !#     if (present(mpiRankOffset  )) mpiRankOffset_  =mpiRankOffset
  !#     ! Draw a random number uniformly from 0 to 1 and use the inverse of our self to get the
  !#     ! corresponding random variate.
  !#     if (present(randomNumberGenerator)) then
  !#        uniformRandom=     randomNumberGenerator%uniformSample(                                  &amp;
  !#             &amp;                                            )
  !#     else
  !#        uniformRandom=self%randomNumberGenerator%uniformSample(                                  &amp;
  !#             &amp;                                             ompThreadOffset=ompThreadOffset_, &amp;
  !#             &amp;                                             mpiRankOffset  =mpiRankOffset_  , &amp;
  !#             &amp;                                             incrementSeed  =incrementSeed     &amp;
  !#             &amp;                                            )
  !#     end if
  !#     distributionFunctionDiscrete1DSample=self%inverse(uniformRandom)
  !#     return
  !#    </code>
  !#  </method>
  !#  <method name="samplerReset" >
  !#    <type>void</type>
  !#    <pass>yes</pass>
  !#    <description>Reset the sampler for the distribution.</description>
  !#    <code>
  !#     self%randomNumberGenerator=pseudoRandom()
  !#    </code>
  !#  </method>
  !#  <method name="minimum" >
  !#    <type>integer</type>
  !#    <pass>yes</pass>
  !#    <description>Returns the minimum possible value in the distribution.</description>
  !#  </method>
  !#  <method name="maximum" >
  !#    <type>integer</type>
  !#    <pass>yes</pass>
  !#    <description>Returns the maximum possible value in the distribution.</description>
  !#  </method>
  !# </functionClass>
  
end module Statistics_Distributions_Discrete
