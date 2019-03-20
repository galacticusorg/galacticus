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

!% Contains a module which implements an peak-background split critical overdensity class.

  !# <criticalOverdensity name="criticalOverdensityPeakBackgroundSplit">
  !#  <description>The ciritical overdensity is given by some other critical overdensity class offset by the halo envinronmental overdensity.</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensityPeakBackgroundSplit
     !% A peak-background split critical overdensity class.
     private
     class(criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     class(haloEnvironmentClass    ), pointer :: haloEnvironment_ => null()
    contains
     final     ::                    peakBackgroundSplitDestructor
     procedure :: value           => peakBackgroundSplitValue
     procedure :: gradientTime    => peakBackgroundSplitGradientTime
     procedure :: gradientMass    => peakBackgroundSplitGradientMass
     procedure :: isMassDependent => peakBackgroundSplitIsMassDependent
  end type criticalOverdensityPeakBackgroundSplit

  interface criticalOverdensityPeakBackgroundSplit
     !% Constructors for the {\normalfont \ttfamily peakBackgroundSplit} critical overdensity class.
     module procedure peakBackgroundSplitConstructorParameters
     module procedure peakBackgroundSplitConstructorInternal
  end interface criticalOverdensityPeakBackgroundSplit

contains
  
  function peakBackgroundSplitConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily peakBackgroundSplit} critical overdensity class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (criticalOverdensityPeakBackgroundSplit)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(criticalOverdensityClass              ), pointer       :: criticalOverdensity_
    class(haloEnvironmentClass                  ), pointer       :: haloEnvironment_
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassVariance_

    ! Check and read parameters.
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !# <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    self=criticalOverdensityPeakBackgroundSplit(criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,cosmologicalMassVariance_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="criticalOverdensity_"     />
    !# <objectDestructor name="haloEnvironment_"         />
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    return
  end function peakBackgroundSplitConstructorParameters

  function peakBackgroundSplitConstructorInternal(criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,cosmologicalMassVariance_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily peakBackgroundSplit} critical overdensity class.
    implicit none
    type (criticalOverdensityPeakBackgroundSplit)                        :: self
    class(criticalOverdensityClass              ), target, intent(in   ) :: criticalOverdensity_    
    class(haloEnvironmentClass                  ), target, intent(in   ) :: haloEnvironment_    
    class(cosmologyFunctionsClass               ), target, intent(in   ) :: cosmologyFunctions_    
    class(cosmologicalMassVarianceClass         ), target, intent(in   ) :: cosmologicalMassVariance_
    !# <constructorAssign variables="*criticalOverdensity_, *haloEnvironment_, *cosmologyFunctions_, *cosmologicalMassVariance_"/>
    
    return
  end function peakBackgroundSplitConstructorInternal

  subroutine peakBackgroundSplitDestructor(self)
    !% Destructor for the {\normalfont \ttfamily peakBackgroundSplit} critical overdensity class.
    implicit none
    type(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self

    !# <objectDestructor name="self%criticalOverdensity_"     />
    !# <objectDestructor name="self%haloEnvironment_"         />
    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    return
  end subroutine peakBackgroundSplitDestructor

  double precision function peakBackgroundSplitValue(self,time,expansionFactor,collapsing,mass,node)
    !% Return the critical overdensity for collapse at the given time and mass.
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    ! Get the critical overdensity at zero environmental overdensity.
    peakBackgroundSplitValue=+self%criticalOverdensity_%value(time,expansionFactor,collapsing,mass,node)
    ! Offset the critical overdensity by the environmental overdensity. Note that the convention used with Galacticus is that the
    ! critical overdensity is scaled by 1/D(t), where D(t) is the growth factor, while σ(M) remains constant at its z=0 value. The
    ! environmental overdensity function provides the overdensity at z=0. This scales with redshift as δₑ(t)=δₑ(t[z=0])
    ! D(t). Applying the 1/D(t) scaling used in our definition of critical overdensity then means that the environmental
    ! contribution to the critical overdensity should always use the z=0 value.
    if (present(node))                                                                                  &
         & peakBackgroundSplitValue=+peakBackgroundSplitValue                                           &
         &                          -self%haloEnvironment_   %overdensityLinear(node,presentDay=.true.)       
    return
  end function peakBackgroundSplitValue

  double precision function peakBackgroundSplitGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !% Return the gradient with respect to time of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    peakBackgroundSplitGradientTime=self%criticalOverdensity_%gradientTime(time,expansionFactor,collapsing,mass,node)
    return
  end function peakBackgroundSplitGradientTime

  double precision function peakBackgroundSplitGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    peakBackgroundSplitGradientMass=self%criticalOverdensity_%gradientMass(time,expansionFactor,collapsing,mass,node)
    return
  end function peakBackgroundSplitGradientMass

  logical function peakBackgroundSplitIsMassDependent(self)
    !% Return whether the critical overdensity is mass dependent.
    implicit none
    class(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self

    peakBackgroundSplitIsMassDependent=self%criticalOverdensity_%isMassDependent()
    return
  end function peakBackgroundSplitIsMassDependent
