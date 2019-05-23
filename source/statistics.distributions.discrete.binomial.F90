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

  !% Implementation of a binomial 1D discrete distibution function.
  
  !# <distributionFunctionDiscrete1D name="distributionFunctionDiscrete1DBinomial">
  !#  <description>A binomial 1D discrete distribution function class.</description>
  !# </distributionFunctionDiscrete1D>
  type, extends(distributionFunctionDiscrete1DClass) :: distributionFunctionDiscrete1DBinomial
     !% Implementation of a binomial 1D discrete distibution function.
     private
     double precision                            :: probabilitySuccess
     integer                                     :: countTrials
     double precision, dimension(:), allocatable :: probabilityCumulative
   contains
     procedure :: density    => binomialMass
     procedure :: cumulative => binomialCumulative
     procedure :: inverse    => binomialInverse
     procedure :: minimum    => binomialMinimum
     procedure :: maximum    => binomialMaximum
  end type distributionFunctionDiscrete1DBinomial

  interface distributionFunctionDiscrete1DBinomial
     !% Constructors for the {\normalfont \ttfamily binomial} 1D discrete distribution function class.
     module procedure binomialConstructorParameters
     module procedure binomialConstructorInternal
  end interface distributionFunctionDiscrete1DBinomial

contains

  function binomialConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily binomial} 1D discrete distribution function class which builds
    !% the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (distributionFunctionDiscrete1DBinomial)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: probabilitySuccess
    integer                                                                 :: countTrials
    
    !# <inputParameter>
    !#   <name>probabilitySuccess</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The probability of success for a single trial.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>countTrials</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of trials.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    self=distributionFunctionDiscrete1DBinomial(probabilitySuccess,countTrials)
    !# <inputParametersValidate source="parameters"/>
    return
  end function binomialConstructorParameters

  function binomialConstructorInternal(probabilitySuccess,countTrials) result(self)
    !% Constructor for ``binomial'' 1D distribution function class.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (distributionFunctionDiscrete1DBinomial)                :: self
    double precision                                        , intent(in   ) :: probabilitySuccess
    integer                                                 , intent(in   ) :: countTrials
    double precision                                        , parameter     :: tolerance         =1.0d-4
    integer                                                                 :: i
    !# <constructorAssign variables="probabilitySuccess, countTrials"/>

    if (probabilitySuccess <  0.0d0 .or. probabilitySuccess > 1.0d0) call Galacticus_Error_Report('p ∈ [0,1]'//{introspection:location})
    if (countTrials        <= 0                                    ) call Galacticus_Error_Report('n ∈ [1,∞]'//{introspection:location})
    ! Build the cumulative distribution.
    allocate(self%probabilityCumulative(0:countTrials))
    do i=0,countTrials
       self%probabilityCumulative(i)=self%mass(i)
       if (i > 0) self%probabilityCumulative(i)=+self%probabilityCumulative(i  ) &
            &                                   +self%probabilityCumulative(i-1)
    end do
    if (self%probabilityCumulative(countTrials) < 1.0-tolerance) call Galacticus_Error_Report('CDF(n) < 1'//{introspection:location})
    self%probabilityCumulative(countTrials)=1.0d0
    return
  end function binomialConstructorInternal

  double precision function binomialMass(self,x)
    !% Return the density of a binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    use Factorials      , only : Factorial
    implicit none
    class  (distributionFunctionDiscrete1DBinomial), intent(inout) :: self
    integer                                        , intent(in   ) :: x

    if (x < 0 .or. x > self%countTrials) call Galacticus_Error_Report('k∈[0,n]'//{introspection:location})
    binomialMass=+            Factorial            (self%countTrials  ) &
         &       /            Factorial            (                 x) &
         &       /            Factorial            (self%countTrials-x) &
         &       *       self%probabilitySuccess **                  x  &
         &       *(1.0d0-self%probabilitySuccess)**(self%countTrials-x)
    return
  end function binomialMass

  double precision function binomialCumulative(self,x)
    !% Return the cumulative probability of a binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (distributionFunctionDiscrete1DBinomial), intent(inout) :: self
    integer                                        , intent(in   ) :: x

    if (x < 0 .or. x > self%countTrials) call Galacticus_Error_Report('k∈[0,n]'//{introspection:location})
    binomialCumulative=self%probabilityCumulative(x)
    return
  end function binomialCumulative

  integer function binomialInverse(self,p)
    !% Return the inverse of a binomial discrete distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (distributionFunctionDiscrete1DBinomial), intent(inout) :: self
    double precision                                        , intent(in   ) :: p
    
    if (p < 0.0d0 .or. p > 1.0d0) call Galacticus_Error_Report('p∈[0,1]'//{introspection:location})
    binomialInverse=0
    do while (self%probabilityCumulative(binomialInverse) < p)
       binomialInverse=binomialInverse+1
    end do
    return
  end function binomialInverse

  integer function binomialMinimum(self)
    !% Return the minimum possible value in a binomial discrete distribution.
    implicit none
    class(distributionFunctionDiscrete1DBinomial), intent(inout) :: self
    !GCC$ attributes unused :: self

    binomialMinimum=0
    return
  end function binomialMinimum

  integer function binomialMaximum(self)
    !% Return the maximum possible value in a binomial discrete distribution.
    implicit none
    class(distributionFunctionDiscrete1DBinomial), intent(inout) :: self
    !GCC$ attributes unused :: self

    binomialMaximum=huge(0)
    return
  end function binomialMaximum
