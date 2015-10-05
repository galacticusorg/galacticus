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

!% Contains a module which implements a \cite{sheth_ellipsoidal_2001} dark matter halo mass function class.
  use Cosmological_Mass_Variance
  use Critical_Overdensities

  !# <haloMassFunction name="haloMassFunctionShethTormen">
  !#  <description>The halo mass function is computed from the function given by \cite{sheth_ellipsoidal_2001}.</description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionShethTormen
     !% A halo mass function class using the fitting function of \cite{sheth_ellipsoidal_2001}.
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_
     double precision                                         :: aValue                   , pValue, &
          &                                                      normalizationValue
   contains
     !@ <objectMethods>
     !@   <object>haloMassFunctionShethTormen</object>
     !@   <objectMethod>
     !@     <method>a</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ time\argin, \doublezero\ mass\argin</arguments>
     !@     <description>Return the parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>p</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ time\argin, \doublezero\ mass\argin</arguments>
     !@     <description>Return the parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>normalization</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ time\argin, \doublezero\ mass\argin</arguments>
     !@     <description>Return the parameter $A$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                  shethTormenDestructor
     procedure :: differential  => shethTormenDifferential
     procedure :: a             => shethTormenA
     procedure :: p             => shethTormenP
     procedure :: normalization => shethTormenNormalization
  end type haloMassFunctionShethTormen

  interface haloMassFunctionShethTormen
     !% Constructors for the {\normalfont \ttfamily shethTormen} halo mass function class.
     module procedure shethTormenConstructorParameters
     module procedure shethTormenConstructorInternal
  end interface haloMassFunctionShethTormen

contains

  function shethTormenConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily shethTormen} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionShethTormen)                :: shethTormenConstructorParameters
    type(inputParameters            ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="cosmologyParameters"      name="shethTormenConstructorParameters%cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="shethTormenConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="shethTormenConstructorParameters%criticalOverdensity_"      source="parameters"/>
    !# <inputParameter>
    !#   <name>a</name>
    !#   <source>parameters</source>
    !#   <variable>shethTormenConstructorParameters%aValue</variable>
    !#   <defaultValue>0.707d0</defaultValue>
    !#   <description>The parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>p</name>
    !#   <source>parameters</source>
    !#   <variable>shethTormenConstructorParameters%pValue</variable>
    !#   <defaultValue>0.3d0</defaultValue>
    !#   <description>The parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <source>parameters</source>
    !#   <variable>shethTormenConstructorParameters%normalizationValue</variable>
    !#   <defaultValue>0.3221836349d0</defaultValue>
    !#   <description>The normalization parameter $A$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
   return
  end function shethTormenConstructorParameters

  function shethTormenConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,normalization)
    !% Internal constructor for the {\normalfont \ttfamily shethTormen} halo mass function class.
    implicit none
    type            (haloMassFunctionShethTormen  )                        :: shethTormenConstructorInternal
    class           (cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_    
    class           (cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class           (criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    double precision                                       , intent(in   ) :: a                             , p, &
         &                                                                    normalization

    shethTormenConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    shethTormenConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    shethTormenConstructorInternal%criticalOverdensity_      => criticalOverdensity_
    shethTormenConstructorInternal%aValue                    =  a
    shethTormenConstructorInternal%pValue                    =  p
    shethTormenConstructorInternal%normalizationValue        =  normalization
   return
  end function shethTormenConstructorInternal

  subroutine shethTormenDestructor(self)
    !% Destructor for the {\normalfont \ttfamily shethTormen} halo mass function class.
    implicit none
    type(haloMassFunctionShethTormen), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_" />
    !# <objectDestructor name="self%criticalOverdensity_"      />
    return
  end subroutine shethTormenDestructor

  double precision function shethTormenDifferential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time   , mass    
    double precision                                             :: alpha  , nu  , &
         &                                                          nuPrime

    ! Compute the mass function.
    nu                     =+(                                                                  &
         &                    +self%criticalOverdensity_     %value       (time=time,mass=mass) &
         &                    /self%cosmologicalMassVariance_%rootVariance(               mass) &
         &                   )**2
    nuPrime                =+self%a(time,mass)                                                  &
         &                  *nu
    alpha                  =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass))
    shethTormenDifferential=+self%cosmologyParameters_%OmegaMatter    () &
         &                  *self%cosmologyParameters_%densityCritical() &
         &                  /mass**2                                     &
         &                  *alpha                                       &
         &                  *sqrt(                                       &
         &                        +2.0d0                                 &
         &                        *nuPrime                               &
         &                        /Pi                                    &
         &                       )                                       &
         &                  *self%normalization(time,mass)               &
         &                  *(                                           &
         &                    +1.0d0                                     &
         &                    +1.0d0                                     &
         &                    /nuPrime**self%p(time,mass)                &
         &                  )                                            &
         &                  *exp(                                        &
         &                       -0.5d0                                  &
         &                       *nuPrime                                &
         &                  )
    return
  end function shethTormenDifferential

  double precision function shethTormenA(self,time,mass)
    !% Return the parameter $a$ in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    

    shethTormenA=self%aValue
    return
  end function shethTormenA
 
  double precision function shethTormenP(self,time,mass)
    !% Return the parameter $p$ in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    

    shethTormenP=self%pValue
    return
  end function shethTormenP
  
  double precision function shethTormenNormalization(self,time,mass)
    !% Return the normalization, $A$, in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass    

    shethTormenNormalization=self%normalizationValue
    return
  end function shethTormenNormalization
 
