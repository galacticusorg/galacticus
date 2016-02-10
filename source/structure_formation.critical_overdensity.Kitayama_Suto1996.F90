!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a critical overdensity class based on the fitting functions of
!% \cite{kitayama_semianalytic_1996}.
  use Linear_Growth
  use Cosmology_Functions

  !# <criticalOverdensity name="criticalOverdensityKitayamaSuto1996" defaultThreadPrivate="yes">
  !#  <description>Provides a critical overdensity class based on the fitting functions of \cite{kitayama_semianalytic_1996}, and is therefore valid only for flat cosmological models.</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensityKitayamaSuto1996
     !% A critical overdensity class based on the fitting functions of \cite{kitayama_semianalytic_1996}.
     private
     class(linearGrowthClass), pointer :: linearGrowth_
    contains
     final     ::                   kitayamaSuto1996Destructor
     procedure :: value          => kitayamaSuto1996Value
     procedure :: gradientTime   => kitayamaSuto1996GradientTime
     procedure :: gradientMass   => kitayamaSuto1996GradientMass
  end type criticalOverdensityKitayamaSuto1996

  interface criticalOverdensityKitayamaSuto1996
     !% Constructors for the {\normalfont \ttfamily kitayamaSuto1996} critical overdensity class.
     module procedure kitayamaSuto1996ConstructorParameters
     module procedure kitayamaSuto1996ConstructorInternal
  end interface criticalOverdensityKitayamaSuto1996

contains

  function kitayamaSuto1996ConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily kitayamaSuto1996} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(criticalOverdensityKitayamaSuto1996)                :: kitayamaSuto1996ConstructorParameters
    type(inputParameters                    ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    !# <objectBuilder class="linearGrowth"             name="kitayamaSuto1996ConstructorParameters%linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="kitayamaSuto1996ConstructorParameters%cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="kitayamaSuto1996ConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
   return
  end function kitayamaSuto1996ConstructorParameters

  function kitayamaSuto1996ConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_)
    !% Internal constructor for the {\normalfont \ttfamily kitayamaSuto1996} critical overdensity class.
    implicit none
    type (criticalOverdensityKitayamaSuto1996)                        :: kitayamaSuto1996ConstructorInternal
    class(cosmologyFunctionsClass            ), target, intent(in   ) :: cosmologyFunctions_
    class(linearGrowthClass                  ), target, intent(in   ) :: linearGrowth_
    class(cosmologicalMassVarianceClass      ), target, intent(in   ) :: cosmologicalMassVariance_

    kitayamaSuto1996ConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    kitayamaSuto1996ConstructorInternal%linearGrowth_             => linearGrowth_
    kitayamaSuto1996ConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    return
  end function kitayamaSuto1996ConstructorInternal

  subroutine kitayamaSuto1996Destructor(self)
    !% Destructor for the {\normalfont \ttfamily kitayamaSuto1996} critical overdensity class.
    implicit none
    type(criticalOverdensityKitayamaSuto1996), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    !# <objectDestructor name="self%linearGrowth_"      />
    return
  end subroutine kitayamaSuto1996Destructor

  double precision function kitayamaSuto1996Value(self,time,expansionFactor,collapsing,mass)
    !% Return the critical overdensity at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time               , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass
    double precision                                                               :: time_
    
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    kitayamaSuto1996Value=+(3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)                                  &
         &                *(1.0d0+0.0123d0*log10(self%cosmologyFunctions_%omegaMatterEpochal(time_))) &
         &                /                      self%linearGrowth_      %value             (time_)
    return
  end function kitayamaSuto1996Value

  double precision function kitayamaSuto1996GradientTime(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to time of critical overdensity at the given time and mass.
    use Numerical_Constants_Math
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time               , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass
    double precision                                                               :: time_              , expansionFactor_
    
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    kitayamaSuto1996GradientTime=+(3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)                                                                 &
         &                       *(                                                                                                         &
         &                                +0.0123d0*      self%cosmologyFunctions_%omegaMatterRateOfChange             (time_           )   &
         &                         /                      self%cosmologyFunctions_%omegaMatterEpochal                  (time_           )   &
         &                         /                log  (10.0d0                                                                         )  &
         &                         -(1.0d0+0.0123d0*log10(self%cosmologyFunctions_%omegaMatterEpochal                  (time_           ))) &
         &                         *                      self%linearGrowth_      %logarithmicDerivativeExpansionFactor(time_           )   &
         &                         *                      self%cosmologyFunctions_%expansionRate                       (expansionFactor_)   &
         &                        )                                                                                                         &
         &                       /                        self%linearGrowth_      %value                               (time_           )
    return
  end function kitayamaSuto1996GradientTime

  double precision function kitayamaSuto1996GradientMass(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time      , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass

    kitayamaSuto1996GradientMass=0.0d0
    return
  end function kitayamaSuto1996GradientMass
