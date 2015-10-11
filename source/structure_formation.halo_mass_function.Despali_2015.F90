!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a \cite{despali_universality_2015} dark matter halo mass function class.
  use Cosmological_Mass_Variance
  use Critical_Overdensities

  !# <haloMassFunction name="haloMassFunctionDespali2015">
  !#  <description>The halo mass function is computed from the function given by \cite{despali_universality_2015}.</description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionShethTormen) :: haloMassFunctionDespali2015
     !% A halo mass function class using the fitting function of \cite{despali_universality_2015}.
     private
     class(virialDensityContrastClass                        ), pointer :: virialDensityContrast_
     type (virialDensityContrastSphericalCollapseMatterLambda)          :: referenceDensityContrast
    contains
     !@ <objectMethods>
     !@   <object>haloMassFunctionDespali2015</object>
     !@   <objectMethod>
     !@     <method>x</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ time\argin, \doublezero\ mass\argin</arguments>
     !@     <description>Return the parameter $x$ in the \cite{despali_universality_2015} halo mass function fit.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                  despali2015Destructor
     procedure :: x             => despali2015X
     procedure :: a             => despali2015A
     procedure :: p             => despali2015P
     procedure :: normalization => despali2015Normalization
  end type haloMassFunctionDespali2015

  interface haloMassFunctionDespali2015
     !% Constructors for the {\normalfont \ttfamily despali2015} halo mass function class.
     module procedure despali2015ConstructorParameters
     module procedure despali2015ConstructorInternal
  end interface haloMassFunctionDespali2015

contains

  function despali2015ConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily despali2015} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionDespali2015)                :: despali2015ConstructorParameters
    type(inputParameters            ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="cosmologyParameters"      name="despali2015ConstructorParameters%cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="despali2015ConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="despali2015ConstructorParameters%criticalOverdensity_"      source="parameters"/>
    despali2015ConstructorParameters%virialDensityContrast_   => virialDensityContrast                             ()
    despali2015ConstructorParameters%referenceDensityContrast =  virialDensityContrastSphericalCollapseMatterLambda()
   return
  end function despali2015ConstructorParameters

  function despali2015ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,virialDensityContrast_)
    !% Internal constructor for the {\normalfont \ttfamily despali2015} halo mass function class.
    implicit none
    type (haloMassFunctionDespali2015  )                        :: despali2015ConstructorInternal
    class(cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_    
    class(cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class(criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    class(virialDensityContrastClass   ), target, intent(in   ) :: virialDensityContrast_

    despali2015ConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    despali2015ConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    despali2015ConstructorInternal%criticalOverdensity_      => criticalOverdensity_
    despali2015ConstructorInternal%virialDensityContrast_    => virialDensityContrast_
    despali2015ConstructorInternal%referenceDensityContrast  =  virialDensityContrastSphericalCollapseMatterLambda()
    return
  end function despali2015ConstructorInternal

  subroutine despali2015Destructor(self)
    !% Destructor for the {\normalfont \ttfamily despali2015} halo mass function class.
    implicit none
    type(haloMassFunctionDespali2015), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_" />
    !# <objectDestructor name="self%criticalOverdensity_"      />
    return
  end subroutine despali2015Destructor

  double precision function despali2015A(self,time,mass)
    !% Return the parameter $a$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    
    double precision                                             :: x
    
    x=self%x(time,mass)
    despali2015A=+0.4332d0*x**2 &
         &       +0.2263d0*x    &
         &       +0.7665d0
    return
  end function despali2015A
 
  double precision function despali2015P(self,time,mass)
    !% Return the parameter $p$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    
    double precision                                             :: x
    
    x=self%x(time,mass)
    despali2015P=-0.1151d0*x**2 &
         &       +0.2554d0*x    &
         &       +0.2488d0
    return
  end function despali2015P
  
  double precision function despali2015Normalization(self,time,mass)
    !% Return the normalization, $A$, in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    
    double precision                                             :: x
    
    x=self%x(time,mass)
    despali2015Normalization=-0.1362d0*x    &
         &                   +0.3292d0
    return
  end function despali2015Normalization
 
  double precision function despali2015X(self,time,mass)
    !% Return the parameter $x$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    

    despali2015X=log10(                                                                    &
         &             +self%virialDensityContrast_  %densityContrast(mass=mass,time=time) &
         &             /self%referenceDensityContrast%densityContrast(mass=mass,time=time) &
         &            )
    return
  end function despali2015X
 
