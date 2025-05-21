!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module that implements a class of distributions.
!!}

module Statistics_Distributions
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  private

  !![
  <functionClass>
   <name>distributionFunction1D</name>
   <descriptiveName>One-dimensional Distribution Functions</descriptiveName>
   <description>Class providing distribution functions of a single variable.</description>
   <default>uniform</default>
   <data>class(randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()</data>
   <destructor>
    <code>
     call distributionFunction1DFinalize(self)
     return
    </code>
   </destructor>
   <method name="density" >
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Return the probability density at {\normalfont \ttfamily x}.</description>
   </method>
   <method name="cumulative" >
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Return the cumulative probability at {\normalfont \ttfamily x}.</description>
   </method>
   <method name="inverse" >
     <type>double precision</type>
     <pass>yes</pass>
     <selfTarget>yes</selfTarget>
     <argument>double precision, intent(in   ) :: p</argument>
     <description>Return the value of the independent variable corresponding to cumulative probability {\normalfont \ttfamily p}.</description>
     <modules>Root_Finder Error</modules>
     <code>
       ! Numerically solve for the inverse.
       type            (rootFinder), save      :: finder
       logical                     , save      :: initialized              =.false.
       double precision            , parameter :: fractionalStep           =1.0d-02
       double precision            , parameter :: toleranceRelative        =1.0d-06
       double precision            , parameter :: toleranceAbsoluteDefault =1.0d-16
       double precision            , parameter :: toleranceAbsoluteRelative=1.0d-06
       double precision                        :: pdfStep                          , pdfMidpoint, &amp;
         &amp;                                    toleranceAbsolute
       !$omp threadprivate(initialized,finder)
       
       if (p &lt; 0.0d0 .or. p &gt; 1.0d0) call Error_Report('probability out of range'//{introspection:location})
       if (.not.initialized) then
	 finder     =rootFinder(inverseRoot)
	 initialized=.true.
       end if
       self_      => self
       p_         =  p
       pdfStep    =  fractionalStep*self%maximum()-fractionalStep*self%minimum()
       pdfMidpoint=  0.5d0*(self%maximum()+self%minimum())
       if (self%minimum() == -huge(0.0d0) .or. self%maximum() == +huge(0.0d0)) then
          toleranceAbsolute=toleranceAbsoluteDefault
       else
          toleranceAbsolute=toleranceAbsoluteRelative*pdfStep
       end if
       call finder%tolerance  (                                                             &amp;
          &amp;                toleranceRelative            =toleranceRelative            , &amp;
          &amp;                toleranceAbsolute            =toleranceAbsolute              &amp;
          &amp;               )
       call finder%rangeExpand(                                                             &amp;
          &amp;                rangeExpandType              =rangeExpandAdditive          , &amp;
          &amp;                rangeExpandUpward            =pdfStep                      , &amp;
          &amp;                rangeExpandDownward          =pdfStep                      , &amp;
          &amp;                rangeUpwardLimit             =self%maximum()               , &amp;
          &amp;                rangeDownwardLimit           =self%minimum()               , &amp;
          &amp;                rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &amp;
          &amp;                rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &amp;
          &amp;               )
	 distributionFunction1DInverse=finder%find(rootGuess=pdfMidpoint)
     </code>
   </method>
   <method name="sample" >
     <type>double precision</type>
     <pass>yes</pass>
     <argument>class(randomNumberGeneratorClass), intent(inout), optional :: randomNumberGenerator_</argument>
     <description>Return a random deviate from the distribution.</description>
     <modules>Error</modules>
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
      distributionFunction1DSample=self%inverse(uniformRandom)
      return
     </code>
   </method>
   <method name="minimum" >
     <type>double precision</type>
     <pass>yes</pass>
     <description>Returns the minimum possible value in the distribution.</description>
     <code>
       distributionFunction1DMinimum=-huge(1.0d0)
     </code>
   </method>
   <method name="maximum" >
     <type>double precision</type>
     <pass>yes</pass>
     <description>Returns the maximum possible value in the distribution.</description>
     <code>
       distributionFunction1DMaximum=+huge(1.0d0)
     </code>
   </method>
  </functionClass>
  !!]

  ! Module-scope variables used in root-finding.
  class(distributionFunction1DClass), pointer :: self_
  double precision :: p_
  !$omp threadprivate(self_,p_)
  
  ! Define a list of distributions.
  type, public :: distributionFunction1DList
     class(distributionFunction1DClass), pointer :: distributionFunction1D_ => null()
  end type distributionFunction1DList

contains
  
  subroutine distributionFunction1DFinalize(self)
    !!{
    Destructor for \refClass{distributionFunction1DClass} objects.
    !!}
    type(distributionFunction1DClass), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine distributionFunction1DFinalize

  double precision function inverseRoot(x)
    !!{
    Root function used in numerically inverting cumulative distribution functions.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    inverseRoot=+self_%cumulative(x) &
         &      -p_
    return
  end function inverseRoot
  
end module Statistics_Distributions
