!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements an fixed critical overdensity class.
  use Linear_Growth
  use Cosmology_Functions

  !# <criticalOverdensity name="criticalOverdensityFixed">
  !#  <description>The critical overdensity is set to a fixed number divided by the linear growth factor.</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensityFixed
     !% A fixed critical overdensity class.
     private
     double precision                             :: criticalOverdensity
     class           (linearGrowthClass), pointer :: linearGrowth_
    contains
     final     ::                   fixedDestructor
     procedure :: value          => fixedValue
     procedure :: gradientTime   => fixedGradientTime
     procedure :: gradientMass   => fixedGradientMass
  end type criticalOverdensityFixed

  interface criticalOverdensityFixed
     !% Constructors for the fixed critical overdensity class.
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface criticalOverdensityFixed

contains

  function fixedConstructorParameters(parameters)
    !% Constructor for the fixed critical overdensity class which takes a parameter set as input.
    use Input_Parameters2
    use Numerical_Constants_Math
    implicit none
    type(criticalOverdensityFixed)                :: fixedConstructorParameters
    type(inputParameters         ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>criticalOverdensity</name>
    !#   <source>parameters</source>
    !#   <variable>fixedConstructorParameters%criticalOverdensity</variable>
    !#   <defaultValue>(3.0d0/20.0d0)*(12.0d0*Pi)**(2.0d0/3.0d0)</defaultValue>
    !#   <description>The value to use for the critical overdensity for collapse of dark matter halos when using a fixed value.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="linearGrowth"             name="fixedConstructorParameters%linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="fixedConstructorParameters%cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="fixedConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
   return
  end function fixedConstructorParameters

  function fixedConstructorInternal(criticalOverdensity,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_)
    !% Internal constructor for the fixed critical overdensity class.
    implicit none
    type            (criticalOverdensityFixed)                        :: fixedConstructorInternal
    double precision                                  , intent(in   ) :: criticalOverdensity
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_    
    class           (linearGrowthClass       ), target, intent(in   ) :: linearGrowth_    
    class(cosmologicalMassVarianceClass      ), target, intent(in   ) :: cosmologicalMassVariance_

    fixedConstructorInternal%criticalOverdensity       =  criticalOverdensity
    fixedConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    fixedConstructorInternal%linearGrowth_             => linearGrowth_
    fixedConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !% Destructor for the fixed critical overdensity class.
    implicit none
    type(criticalOverdensityFixed), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    !# <objectDestructor name="self%linearGrowth_"      />
    return
  end subroutine fixedDestructor

  double precision function fixedValue(self,time,expansionFactor,collapsing,mass)
    !% Return the critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time               , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    double precision                                                    :: time_
    !GCC$ attributes unused :: mass
    
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    fixedValue=+self%criticalOverdensity              &
         &     /self%linearGrowth_      %value(time_)
    return
  end function fixedValue

  double precision function fixedGradientTime(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to time of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time               , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    double precision                                                    :: time_              , expansionFactor_
    !GCC$ attributes unused :: mass
    
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    fixedGradientTime=-self%criticalOverdensity                                                        &
         &            *self%linearGrowth_      %logarithmicDerivativeExpansionFactor(time_           ) &
         &            *self%cosmologyFunctions_%expansionRate                       (expansionFactor_) &
         &            /self%linearGrowth_      %value                               (time_           )    
    return
  end function fixedGradientTime

  double precision function fixedGradientMass(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time      , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    !GCC$ attributes unused :: self, time, expansionFactor, collapsing, mass

    fixedGradientMass=0.0d0
    return
  end function fixedGradientMass
