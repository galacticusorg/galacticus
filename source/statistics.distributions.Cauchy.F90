!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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
  
  !% Implementation of a 1D Cauchy distibution function.
  
  !# <distributionFunction1D name="distributionFunction1DCauchy">
  !#  <description>A 1D Cauchy distibution function.</description>
  !# </distributionFunction1D>
  type, extends(distributionFunction1DClass) :: distributionFunction1DCauchy
     !% Implementation of a 1D Cauchy distibution function.
     private
     double precision :: median, scale
   contains
     procedure :: density    => cauchyDensity
     procedure :: cumulative => cauchyCumulative
     procedure :: inverse    => cauchyInverse
  end type distributionFunction1DCauchy

  interface distributionFunction1DCauchy
     !% Constructors for the {\normalfont \ttfamily cauchy} 1D distribution function class.
     module procedure cauchyConstructorParameters
     module procedure cauchyConstructorInternal
     module procedure cauchyConstructorProbability
  end interface distributionFunction1DCauchy

contains

  function cauchyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily cauchy} 1D distribution function class which builds
    !% the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (distributionFunction1DCauchy)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: median    , scale

    !# <inputParameter>
    !#   <name>median</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The median of the Cauchy distribution function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scale</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The scale parameter of the Cauchy distribution function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=distributionFunction1DCauchy(median,scale)
    !# <inputParametersValidate source="parameters"/>
    return
  end function cauchyConstructorParameters
  
  function cauchyConstructorProbability(median,limit,probability) result(self)
    !% Constructor for ``cauchy'' 1D distribution function class.
    use Numerical_Constants_Math
    type            (distributionFunction1DCauchy)                :: self
    double precision                              , intent(in   ) :: median     , limit, &
         &                                                           probability

    self=distributionFunction1DCauchy(median,limit/tan(0.5d0*Pi*(1.0d0-probability)))
    return
  end function cauchyConstructorProbability

  function cauchyConstructorInternal(median,scale) result(self)
    !% Constructor for ``cauchy'' 1D distribution function class.
    type            (distributionFunction1DCauchy)                :: self
    double precision                              , intent(in   ) :: median, scale
    !# <constructorAssign variables="median, scale"/>
    
    return
  end function cauchyConstructorInternal

  double precision function cauchyDensity(self,x)
    !% Return the density of a Cauchy distribution.
    use Numerical_Constants_Math
    implicit none
    class           (distributionFunction1DCauchy), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    cauchyDensity=1.0d0/Pi/self%scale/(1.0d0+((x-self%median)/self%scale)**2)
    return
  end function cauchyDensity

  double precision function cauchyCumulative(self,x)
    !% Return the cumulative probability of a Cauchy distribution.
    use Numerical_Constants_Math
    implicit none
    class           (distributionFunction1DCauchy), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    cauchyCumulative=0.5d0+atan((x-self%median)/self%scale)/Pi
    return
  end function cauchyCumulative

  double precision function cauchyInverse(self,p)
    !% Return the inverse of a Cauchy distribution.
    use Galacticus_Error
    use Numerical_Constants_Math
    implicit none
    class           (distributionFunction1DCauchy), intent(inout) :: self
    double precision                              , intent(in   ) :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Galacticus_Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    cauchyInverse=self%median+self%scale*tan(Pi*(p-0.5d0))
    return
  end function cauchyInverse
