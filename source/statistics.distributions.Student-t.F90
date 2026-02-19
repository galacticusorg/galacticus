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
  Implementation of a 1D Student-t distribution function.
  !!}

  ! Add dependency on GSL library.
  !; gsl
  
  use, intrinsic :: ISO_C_Binding, only : c_double
  
  !![
  <distributionFunction1D name="distributionFunction1DStudentT">
   <description>
    Student's t-distribution:
    \begin{equation}
     P(x) \propto \left(1 + {x^2\over \nu}\right)^{-(\nu+1)/2}
    \end{equation}
    Specified using:
    \begin{description}
    \item[{\normalfont \ttfamily degreesOfFreedom}] The number of degrees of freedom, $\nu$.
    \end{description}
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DStudentT
     !!{
     Implementation of a 1D Student-t distribution function.
     !!}
     private
     double precision :: degreesOfFreedom
   contains
     !![
     <methods>
       <method description="The upper-tail cumulative distribution function." method="cumulativeUpper" />
       <method description="The upper-tail inverse cumulative distribution function." method="inverseUpper" />
     </methods>
     !!]
     procedure :: density         => studentTDensity
     procedure :: cumulative      => studentTCumulative
     procedure :: inverse         => studentTInverse
     procedure :: cumulativeUpper => studentTCumulativeUpper
     procedure :: inverseUpper    => studentTInverseUpper
  end type distributionFunction1DStudentT

  interface distributionFunction1DStudentT
     !!{
     Constructors for the \refClass{distributionFunction1DStudentT} 1D distribution function class.
     !!}
     module procedure studentTConstructorParameters
     module procedure studentTConstructorInternal
  end interface distributionFunction1DStudentT

  interface
     function gsl_ran_tdist_pdf(x,nu) bind(c,name='gsl_ran_tdist_pdf')
       !!{
       Template for the GSL Student t-distribution probability density.
       !!}
       import
       real(c_double)        :: gsl_ran_tdist_pdf
       real(c_double), value :: x                , nu
     end function gsl_ran_tdist_pdf

     function gsl_cdf_tdist_P(x,nu) bind(c,name='gsl_cdf_tdist_P')
       !!{
       Template for the GSL Student t-distribution cumulative probability function.
       !!}
       import
       real(c_double)        :: gsl_cdf_tdist_P
       real(c_double), value :: x              , nu
     end function gsl_cdf_tdist_P

     function gsl_cdf_tdist_Pinv(P,nu) bind(c,name='gsl_cdf_tdist_Pinv')
       !!{
       Template for the GSL Student t-distribution inverse cumulative probability function.
       !!}
       import
       real(c_double)        :: gsl_cdf_tdist_Pinv
       real(c_double), value :: P              , nu
     end function gsl_cdf_tdist_Pinv

     function gsl_cdf_tdist_Q(x,nu) bind(c,name='gsl_cdf_tdist_Q')
       !!{
       Template for the GSL Student t-distribution upper-tail cumulative probability function.
       !!}
       import
       real(c_double)        :: gsl_cdf_tdist_Q
       real(c_double), value :: x              , nu
     end function gsl_cdf_tdist_Q

     function gsl_cdf_tdist_Qinv(Q,nu) bind(c,name='gsl_cdf_tdist_Qinv')
       !!{
       Template for the GSL Student t-distribution inverse upper-tail cumulative probability function.
       !!}
       import
       real(c_double)        :: gsl_cdf_tdist_Qinv
       real(c_double), value :: Q              , nu
     end function gsl_cdf_tdist_Qinv
  end interface
  
contains

  function studentTConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DStudentT} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DStudentT)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass    ), pointer       :: randomNumberGenerator_
    double precision                                                :: degreesOfFreedom

    !![
    <inputParameter>
      <name>degreesOfFreedom</name>
      <description>The degrees of freedom of the Student-t distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DStudentT(degreesOfFreedom,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function studentTConstructorParameters

  function studentTConstructorInternal(degreesOfFreedom,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DStudentT} 1D distribution function class.
    !!}
    type            (distributionFunction1DStudentT)                                  :: self
    double precision                                , intent(in   )                   :: degreesOfFreedom
    class           (randomNumberGeneratorClass    ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="degreesOfFreedom, *randomNumberGenerator_"/>
    !!]

    return
  end function studentTConstructorInternal

  double precision function studentTDensity(self,x)
    !!{
    Return the density of a Student-t distribution.
    !!}
    implicit none
    class           (distributionFunction1DStudentT), intent(inout) :: self
    double precision                                , intent(in   ) :: x

    studentTDensity=GSL_Ran_tDist_PDF(x,self%degreesOfFreedom)
    return
  end function studentTDensity

  double precision function studentTCumulative(self,x)
    !!{
    Return the cumulative probability of a Student-t distribution.
    !!}
    implicit none
    class           (distributionFunction1DStudentT), intent(inout) :: self
    double precision                                , intent(in   ) :: x

    studentTCumulative=GSL_CDF_tDist_P(x,self%degreesOfFreedom)
    return
  end function studentTCumulative

  double precision function studentTInverse(self,p)
    !!{
    Return the inverse of the cumulative probability of a Student-t distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DStudentT), intent(inout), target :: self
    double precision                                , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    studentTInverse=GSL_CDF_tDist_Pinv(p,self%degreesOfFreedom)
    return
  end function studentTInverse

  double precision function studentTCumulativeUpper(self,x)
    !!{
    Return the upper-tail cumulative probability of a Student-t distribution.
    !!}
    implicit none
    class           (distributionFunction1DStudentT), intent(inout) :: self
    double precision                                , intent(in   ) :: x

    studentTCumulativeUpper=GSL_CDF_tDist_Q(x,self%degreesOfFreedom)
    return
  end function studentTCumulativeUpper

  double precision function studentTInverseUpper(self,q)
    !!{
    Return the inverse of the upper-tail cumulative probability of a Student-t distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DStudentT), intent(inout), target :: self
    double precision                                , intent(in   )         :: q

    if (q < 0.0d0 .or. q > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    studentTInverseUpper=GSL_CDF_tDist_Qinv(q,self%degreesOfFreedom)
    return
  end function studentTInverseUpper
