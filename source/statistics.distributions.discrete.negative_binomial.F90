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

  !% Implementation of a negative binomial 1D discrete distibution function.
  
  !# <distributionFunctionDiscrete1D name="distributionFunctionDiscrete1DNegativeBinomial">
  !#  <description>A negative binomial 1D discrete distribution function class.</description>
  !# </distributionFunctionDiscrete1D>
  type, extends(distributionFunctionDiscrete1DClass) :: distributionFunctionDiscrete1DNegativeBinomial
     !% Implementation of a negativeBinomial 1D discrete distibution function.
     private
     double precision :: probabilitySuccess, countFailures
   contains
     procedure :: mass            => negativeBinomialMass
     procedure :: massLogarithmic => negativeBinomialMassLogarithmic
     procedure :: cumulative      => negativeBinomialCumulative
     procedure :: inverse         => negativeBinomialInverse
     procedure :: minimum         => negativeBinomialMinimum
     procedure :: maximum         => negativeBinomialMaximum
  end type distributionFunctionDiscrete1DNegativeBinomial

  interface distributionFunctionDiscrete1DNegativeBinomial
     !% Constructors for the {\normalfont \ttfamily negativeBinomial} 1D discrete distribution function class.
     module procedure negativeBinomialConstructorParameters
     module procedure negativeBinomialConstructorInternal
  end interface distributionFunctionDiscrete1DNegativeBinomial

contains

  function negativeBinomialConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily negativeBinomial} 1D discrete distribution function class which builds
    !% the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (distributionFunctionDiscrete1DNegativeBinomial)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    double precision                                                                :: probabilitySuccess, countFailures
    
    !# <inputParameter>
    !#   <name>probabilitySuccess</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The probability of success for a single trial.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>countFailures</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of failures.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=distributionFunctionDiscrete1DNegativeBinomial(probabilitySuccess,countFailures)
    !# <inputParametersValidate source="parameters"/>
    return
  end function negativeBinomialConstructorParameters

  function negativeBinomialConstructorInternal(probabilitySuccess,countFailures) result(self)
    !% Constructor for ``negativeBinomial'' 1D distribution function class.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (distributionFunctionDiscrete1DNegativeBinomial)                :: self
    double precision                                                , intent(in   ) :: probabilitySuccess, countFailures
    !# <constructorAssign variables="probabilitySuccess, countFailures"/>

    if (probabilitySuccess <  0.0d0 .or. probabilitySuccess > 1.0d0) call Galacticus_Error_Report('p ∈ [0,1]'//{introspection:location})
    if (countFailures      <= 0.0d0                                ) call Galacticus_Error_Report('r ∈ (0,∞]'//{introspection:location})
    return
  end function negativeBinomialConstructorInternal

  double precision function negativeBinomialMass(self,x)
    !% Return the mass of a negative binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    use Gamma_Functions , only : Gamma_Function_Logarithmic
    implicit none
    class  (distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    integer                                                , intent(in   ) :: x

    if (x < 0) call Galacticus_Error_Report('k ∈ [0,n]'//{introspection:location})
    negativeBinomialMass=+exp(                                                        &
         &                    +Gamma_Function_Logarithmic(dble(x)+self%countFailures) &
         &                    -Gamma_Function_Logarithmic(       +self%countFailures) &
         &                    -Gamma_Function_Logarithmic(dble(x)+     1.0d0        ) &
         &                   )                                                        &
         &               *       self%probabilitySuccess **self%countFailures         &
         &               *(1.0d0-self%probabilitySuccess)**     x
    return
  end function negativeBinomialMass

  double precision function negativeBinomialMassLogarithmic(self,x)
    !% Return the logarithmic mass of a negative binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    use Gamma_Functions , only : Gamma_Function_Logarithmic
    implicit none
    class  (distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    integer                                                , intent(in   ) :: x

    if (x < 0) call Galacticus_Error_Report('k ∈ [0,n]'//{introspection:location})
    negativeBinomialMassLogarithmic=+Gamma_Function_Logarithmic(dble(x)+self%countFailures) &
         &                          -Gamma_Function_Logarithmic(       +self%countFailures) &
         &                          -Gamma_Function_Logarithmic(dble(x)+     1.0d0        ) &
         &                          +self%countFailures*log(     +self%probabilitySuccess)  &
         &                          +dble(x)           *log(1.0d0-self%probabilitySuccess)
    return
  end function negativeBinomialMassLogarithmic

  double precision function negativeBinomialCumulative(self,x)
    !% Return the cumulative probability of a negativeBinomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    integer                                                , intent(in   ) :: x
    !GCC$ attributes unused :: self, x
    
    negativeBinomialCumulative=0.0d0
    call Galacticus_Error_Report('cumulative distribution function is not implemented'//{introspection:location})
    return
  end function negativeBinomialCumulative

  integer function negativeBinomialInverse(self,p)
    !% Return the inverse of a negative binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    double precision                                                , intent(in   ) :: p
    !GCC$ attributes unused :: self, p
    
    negativeBinomialInverse=0
    call Galacticus_Error_Report('inverse function is not implemented'//{introspection:location})
    return
  end function negativeBinomialInverse

  integer function negativeBinomialMinimum(self)
    !% Return the minimum possible value in a negative binomial discrete distribution.
    implicit none
    class(distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    !GCC$ attributes unused :: self

    negativeBinomialMinimum=0
    return
  end function negativeBinomialMinimum

  integer function negativeBinomialMaximum(self)
    !% Return the maximum possible value in a negative binomial discrete distribution.
    implicit none
    class(distributionFunctionDiscrete1DNegativeBinomial), intent(inout) :: self
    !GCC$ attributes unused :: self

    negativeBinomialMaximum=huge(0)
    return
  end function negativeBinomialMaximum
