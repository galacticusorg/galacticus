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
  Implementation of a normal 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DNormal">
   <description>
    A normal distribution, optionally with lower and upper limits:
    \begin{equation}
     P(x) \propto \left\{ \begin{array}{ll} \exp[-(x-\mu)^2/2S] &amp; \hbox{ if } x_\mathrm{l} \leq x \leq x_\mathrm{u} \\ 0 &amp; \hbox{ otherwise.}  \end{array} \right.
    \end{equation}
    Specified using:
    \begin{description}
    \item[{\normalfont \ttfamily [mean]}] The mean, $\mu$;
    \item[{\normalfont \ttfamily [variance]}] The variance, $S$;
    \item[{\normalfont \ttfamily [minimum]}] The lower limit of the range, $x_\mathrm{l}$;
    \item[{\normalfont \ttfamily [maximum]}] The upper limit of the range, $x_\mathrm{u}$.
    \end{description}
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DNormal
     !!{
     Implementation of a normal 1D distribution function.
     !!}
     private
     logical          :: limitLowerExists, limitUpperExists
     double precision :: limitLower      , limitUpper      , &
          &              cdfAtLowerLimit , cdfAtUpperLimit , &
          &              mean            , variance
   contains
     procedure :: density    => normalDensity
     procedure :: cumulative => normalCumulative
     procedure :: inverse    => normalInverse
     procedure :: minimum    => normalMinimum
     procedure :: maximum    => normalMaximum
  end type distributionFunction1DNormal

  interface distributionFunction1DNormal
     !!{
     Constructors for the \refClass{distributionFunction1DNormal} 1D distribution function class.
     !!}
     module procedure normalConstructorParameters
     module procedure normalConstructorInternal
  end interface distributionFunction1DNormal

contains

  function normalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNormal} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DNormal)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    double precision                                              :: mean                  , variance  , &
         &                                                           limitLower            , limitUpper

    !![
    <inputParameter>
      <name>mean</name>
      <description>The mean of the normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>variance</name>
      <description>The variance of the normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('limitLower')) then
       !![
       <inputParameter>
         <name>limitLower</name>
         <description>The lower limit of the normal distribution.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    if (parameters%isPresent('limitUpper')) then
       !![
       <inputParameter>
         <name>limitUpper</name>
         <description>The upper limit of the normal distribution.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <conditionalCall>
      <call>self=distributionFunction1DNormal(mean,variance,randomNumberGenerator_=randomNumberGenerator_{conditions})</call>
      <argument name="limitLower" value="limitLower" parameterPresent="parameters"/>
      <argument name="limitUpper" value="limitUpper" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function normalConstructorParameters

  function normalConstructorInternal(mean,variance,limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNormal} 1D distribution function class.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    type            (distributionFunction1DNormal)                                  :: self
    double precision                              , intent(in   )                   :: mean                  , variance
    double precision                              , intent(in   ), optional         :: limitLower            , limitUpper
    class           (randomNumberGeneratorClass  ), intent(in   ), optional, target :: randomNumberGenerator_
    !![
    <constructorAssign variables="mean, variance, limitLower, limitUpper, *randomNumberGenerator_"/>
    !!]

    self%limitLowerExists=present(limitLower)
    self%limitUpperExists=present(limitUpper)
    if (self%limitLowerExists) then
       self%limitLower     =limitLower
       self%cdfAtLowerLimit=0.5d0*(1.0d0+Error_Function((limitLower-mean)/sqrt(2.0d0*variance)))
    else
       self%limitUpper     =-huge(0.0d0)
       self%cdfAtLowerLimit=      0.0d0
    end if
    if (self%limitUpperExists) then
       self%limitUpper     =limitUpper
       self%cdfAtUpperLimit=0.5d0*(1.0d0+Error_Function((limitUpper-mean)/sqrt(2.0d0*variance)))
    else
       self%limitUpper     =+huge(0.0d0)
       self%cdfAtUpperLimit=      1.0d0
    end if
    return
  end function normalConstructorInternal

  double precision function normalMinimum(self)
    !!{
    Return the minimum possible value of a uniform distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DNormal), intent(inout) :: self

    if (self%limitLowerExists) then
       normalMinimum=self%limitLower
    else
       normalMinimum=0.0d0
       call Error_Report('no minimum exists'//{introspection:location})
    end if
    return
  end function normalMinimum

  double precision function normalMaximum(self)
    !!{
    Return the maximum possible value of a uniform distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DNormal), intent(inout) :: self

    if (self%limitUpperExists) then
       normalMaximum=self%limitUpper
    else
       normalMaximum=0.0d0
       call Error_Report('no maximum exists'//{introspection:location})
    end if
    return
  end function normalMaximum

  double precision function normalDensity(self,x)
    !!{
    Return the density of a normal distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DNormal), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    if     (                                                   &
         &   (self%limitLowerExists .and. x < self%limitLower) &
         &  .or.                                               &
         &   (self%limitUpperExists .and. x > self%limitUpper) &
         & ) then
       normalDensity=0.0d0
    else
       normalDensity=     &
            & exp(                    &
            &     -0.5d0              &
            &     *(x-self%mean)**2   &
            &     /self%variance      &
            &    )                    &
            & /sqrt(                  &
            &        2.0d0            &
            &       *Pi               &
            &       *self%variance    &
            &      )                  &
            & /(                      &
            &   +self%cdfAtUpperLimit &
            &   -self%cdfAtLowerLimit &
            &  )
    end if
    return
  end function normalDensity

  double precision function normalCumulative(self,x)
    !!{
    Return the cumulative probability of a normal distribution.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    class           (distributionFunction1DNormal), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    if      (self%limitLowerExists .and. x < self%limitLower) then
       normalCumulative=0.0d0
    else if (self%limitUpperExists .and. x > self%limitUpper) then
       normalCumulative=1.0d0
    else
       normalCumulative=+(                                                                       &
            &             +0.5d0*(1.0d0+Error_Function((x-self%mean)/sqrt(2.0d0*self%variance))) &
            &             -self%cdfAtLowerLimit                                                  &
            &            )                                                                       &
            &           /(                                                                       &
            &             +self%cdfAtUpperLimit                                                  &
            &             -self%cdfAtLowerLimit                                                  &
            &            )
    end if
    return
  end function normalCumulative

  double precision function normalInverse(self,p)
    !!{
    Return the inverse of a normal distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DNormal), intent(inout), target :: self
    double precision                              , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    normalInverse=+self%mean                                     &
         &        +sqrt(self%variance)                           &
         &        *normalStandardInverse(                        &
         &                                p                      &
         &                               *(                      &
         &                                 +self%cdfAtUpperLimit &
         &                                 -self%cdfAtLowerLimit &
         &                               )                       &
         &                               +  self%cdfAtLowerLimit &
         &                              )
    return
  end function normalInverse

  double precision function normalStandardInverse(p)
    !!{
    Evaluates the inverse of the standard normal cumulative distribution function. Based on the Fortran90 version by John
    Burkardt (itself based on the original Fortran 77 version by Michael Wichura), using the algorithm of
    \cite{wichura_percentage_1988}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   )                :: p
    double precision, parameter    , dimension (8) :: a=[                            &
         &                                               3.38713287279636660800d+00, &
         &                                               1.33141667891784377450d+02, &
         &                                               1.97159095030655144270d+03, &
         &                                               1.37316937655094611250d+04, &
         &                                               4.59219539315498714570d+04, &
         &                                               6.72657709270087008530d+04, &
         &                                               3.34305755835881281050d+04, &
         &                                               2.50908092873012267270d+03  &
         &                                              ]
    double precision, parameter    , dimension (8) :: b=[                            &
         &                                               1.00000000000000000000d+00, &
         &                                               4.23133307016009112520d+01, &
         &                                               6.87187007492057908300d+02, &
         &                                               5.39419602142475110770d+03, &
         &                                               2.12137943015865958670d+04, &
         &                                               3.93078958000927106100d+04, &
         &                                               2.87290857357219426740d+04, &
         &                                               5.22649527885285456100d+03  &
         &                                              ]
    double precision, parameter    , dimension (8) :: c=[                            &
         &                                               1.42343711074968357734d+00, &
         &                                               4.63033784615654529590d+00, &
         &                                               5.76949722146069140550d+00, &
         &                                               3.64784832476320460504d+00, &
         &                                               1.27045825245236838258d+00, &
         &                                               2.41780725177450611770d-01, &
         &                                               2.27238449892691845833d-02, &
         &                                               7.74545014278341407640d-04  &
         &                                              ]
    double precision, parameter    , dimension (8) :: d=[                            &
         &                                               1.00000000000000000000d+00, &
         &                                               2.05319162663775882187d+00, &
         &                                               1.67638483018380384940d+00, &
         &                                               6.89767334985100004550d-01, &
         &                                               1.48103976427480074590d-01, &
         &                                               1.51986665636164571966d-02, &
         &                                               5.47593808499534494600d-04, &
         &                                               1.05075007164441684324d-09  &
         &                                              ]
    double precision, parameter    , dimension (8) :: e=[                            &
         &                                               6.65790464350110377720d+00, &
         &                                               5.46378491116411436990d+00, &
         &                                               1.78482653991729133580d+00, &
         &                                               2.96560571828504891230d-01, &
         &                                               2.65321895265761230930d-02, &
         &                                               1.24266094738807843860d-03, &
         &                                               2.71155556874348757815d-05, &
         &                                               2.01033439929228813265d-07  &
         &                                              ]
    double precision, parameter    , dimension (8) :: f=[                            &
         &                                               1.00000000000000000000d+00, &
         &                                               5.99832206555887937690d-01, &
         &                                               1.36929880922735805310d-01, &
         &                                               1.48753612908506148525d-02, &
         &                                               7.86869131145613259100d-04, &
         &                                               1.84631831751005468180d-05, &
         &                                               1.42151175831644588870d-07, &
         &                                               2.04426310338993978564d-15  &
         &                                              ]
    double precision, parameter                    :: const1=0.180625d0
    double precision, parameter                    :: const2=1.6d0
    double precision, parameter                    :: split1=0.425d0
    double precision, parameter                    :: split2=5.0d0
    double precision                               :: q, r

    if (p <= 0.0d0) then
       normalStandardInverse=-huge(p)
       return
    end if
    if (1.0d0 <= p) then
       normalStandardInverse=+huge(p)
       return
    end if
    q=p-0.5d0
    if (abs(q) <= split1) then
       r=const1-q*q
       normalStandardInverse=+q*normalPolynomialEvaluate(8,a,r) &
            &                  /normalPolynomialEvaluate(8,b,r)
    else
       if (q < 0.0d0) then
          r=p
       else
          r=1.0d0-p
       end if
       if (r <= 0.0d0) then
          normalStandardInverse=-1.0d0
          call Error_Report('out of range - this should not happen'//{introspection:location})
       end if
       r=sqrt(-log(r))
       if (r <= split2) then
          r=r-const2
          normalStandardInverse=+normalPolynomialEvaluate(8,c,r) &
               &                /normalPolynomialEvaluate(8,d,r)
       else
          r=r-split2
          normalStandardInverse=+normalPolynomialEvaluate(8,e,r) &
               &                /normalPolynomialEvaluate(8,f,r)
       end if
       if (q < 0.0d0) normalStandardInverse=-normalStandardInverse
    end if
    return
  end function normalStandardInverse

  double precision function normalPolynomialEvaluate(n,a,x)
    !!{
    Evaluates a polynomial based on the implementation by John Burkardt. For sanity's sake, the value of N indicates the
    \emph{number} of coefficients, or more precisely, the \emph{order} of the polynomial, rather than the \emph{degree} of the
    polynomial. The two quantities differ by 1, but cause a great deal of confusion. Given {\normalfont \ttfamily n} and {\normalfont \ttfamily a}, the form of the
    polynomial is:
    \begin{equation}
    p(x) = a(1) + a(2) * x + \ldots + a(n-1) * x^{n-2} + a(n) * x^{n-1}.
    \end{equation}
    !!}
    implicit none
    integer         , intent(in   )               :: n
    double precision, intent(in   ), dimension(n) :: a
    double precision, intent(in   )               :: x
    integer                                       :: i

    normalPolynomialEvaluate=0.0d0
    do i=n,1,-1
       normalPolynomialEvaluate=normalPolynomialEvaluate*x+a(i)
    end do
    return
  end function normalPolynomialEvaluate
