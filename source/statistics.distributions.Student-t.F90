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

  !% Implementation of a 1D Student-t distribution function.
  
  !# <distributionFunction1D name="distributionFunction1DStudentT">
  !#  <description>A 1D Student-t distibution function.</description>
  !# </distributionFunction1D>
  type, extends(distributionFunction1DClass) :: distributionFunction1DStudentT
     !% Implementation of a 1D Student-t distribution function.
     private
     double precision :: degreesOfFreedom
   contains
     procedure :: density    => studentTDensity
     procedure :: cumulative => studentTCumulative
     procedure :: inverse    => studentTInverse
  end type distributionFunction1DStudentT

  interface distributionFunction1DStudentT
     !% Constructors for the {\normalfont \ttfamily studentT} 1D distribution function class.
     module procedure studentTConstructorParameters
     module procedure studentTConstructorInternal
  end interface distributionFunction1DStudentT

contains

  function studentTConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily studentT} 1D distribution function class which builds
    !% the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (distributionFunction1DStudentT)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: degreesOfFreedom

    !# <inputParameter>
    !#   <name>degreesOfFreedom</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The degrees of freedom of the Student-t distribution function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=distributionFunction1DStudentT(degreesOfFreedom)
    !# <inputParametersValidate source="parameters"/>
    return
  end function studentTConstructorParameters
  
  function studentTConstructorInternal(degreesOfFreedom) result(self)
    !% Constructor for ``studentT'' 1D distribution function class.
    type            (distributionFunction1DStudentT)                :: self
    double precision                                , intent(in   ) :: degreesOfFreedom
    !# <constructorAssign variables="degreesOfFreedom"/>
    
    self%randomNumberGenerator=pseudoRandom()
    return
  end function studentTConstructorInternal

  double precision function studentTDensity(self,x)
    !% Return the density of a Student-t distribution.
    use FGSL, only : FGSL_Ran_tDist_PDF
    implicit none
    class           (distributionFunction1DStudentT), intent(inout) :: self
    double precision                                , intent(in   ) :: x

    studentTDensity=FGSL_Ran_tDist_PDF(x,self%degreesOfFreedom)
    return
  end function studentTDensity

  double precision function studentTCumulative(self,x)
    !% Return the cumulative probability of a Student-t distribution.
    use FGSL, only : FGSL_CDF_tDist_P
    implicit none
    class           (distributionFunction1DStudentT), intent(inout) :: self
    double precision                                , intent(in   ) :: x

    studentTCumulative=FGSL_CDF_tDist_P(x,self%degreesOfFreedom)
    return
  end function studentTCumulative

  double precision function studentTInverse(self,p)
    !% Return the cumulative probability of a Student-t distribution.
    use Galacticus_Error, only : Galacticus_Error_Report
    use FGSL            , only : FGSL_CDF_tDist_Pinv
    implicit none
    class           (distributionFunction1DStudentT), intent(inout), target :: self
    double precision                                , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Galacticus_Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    studentTInverse=FGSL_CDF_tDist_Pinv(p,self%degreesOfFreedom)
    return
  end function studentTInverse
