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

  !% An implementation of cosmological density field mass variance which modifies another member of the class by offsetting for
  !% the peak-background split.

  use Cosmology_Parameters
  
  !# <cosmologicalMassVariance name="cosmologicalMassVariancePeakBackgroundSplit">
  !#  <description>
  !#   The cosmological mass variance is computed by taking the variance from some other mass variance class, $\sigma^2(M)$, and
  !#   offsetting it by the variance of the background in the peak-background split model, $\sigma^2(M_\mathrm{e})$, where
  !#   $M_\mathrm{e}$ is the mass contained within the region defined as the background.
  !#  </description>
  !# </cosmologicalMassVariance>
  
  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVariancePeakBackgroundSplit
     !% A cosmological mass variance class computing variance from a filtered power spectrum.
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
     class           (haloEnvironmentClass         ), pointer :: haloEnvironment_
     double precision                                         :: varianceBackground       , massBackground
   contains
     final     ::                                       variancePeakBackgroundSplitDestructor
     procedure :: sigma8                             => variancePeakBackgroundSplitSigma8
     procedure :: powerNormalization                 => variancePeakBackgroundSplitPowerNormalization
     procedure :: rootVariance                       => variancePeakBackgroundSplitRootVariance
     procedure :: rootVarianceLogarithmicGradient    => variancePeakBackgroundSplitRootVarianceLogarithmicGradient
     procedure :: rootVarianceAndLogarithmicGradient => variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient
     procedure :: mass                               => variancePeakBackgroundSplitMass
  end type cosmologicalMassVariancePeakBackgroundSplit

  interface cosmologicalMassVariancePeakBackgroundSplit
     !% Constructors for the {\normalfont \ttfamily peakBackgroundSplit} cosmological mass variance class.
     module procedure variancePeakBackgroundSplitConstructorParameters
     module procedure variancePeakBackgroundSplitConstructorInternal
  end interface cosmologicalMassVariancePeakBackgroundSplit

contains

  function variancePeakBackgroundSplitConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily peakBackgroundSplit} cosmological mass variance class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (cosmologicalMassVariancePeakBackgroundSplit)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(haloEnvironmentClass                       ), pointer       :: haloEnvironment_
    class(cosmologicalMassVarianceClass              ), pointer       :: cosmologicalMassVariance_
    class(cosmologyParametersClass                   ), pointer       :: cosmologyParameters_

    !# <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>    
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>    
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>    
    ! Construct the instance.
    self=cosmologicalMassVariancePeakBackgroundSplit(haloEnvironment_,cosmologicalMassVariance_,cosmologyParameters_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="haloEnvironment_"         />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    !# <objectDestructor name="cosmologyParameters_"     />
    return
  end function variancePeakBackgroundSplitConstructorParameters

  function variancePeakBackgroundSplitConstructorInternal(haloEnvironment_,cosmologicalMassVariance_,cosmologyParameters_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily peakBackgroundSplit} linear growth class.
    implicit none
    type (cosmologicalMassVariancePeakBackgroundSplit)                        :: self
    class(haloEnvironmentClass                       ), intent(in   ), target :: haloEnvironment_
    class(cosmologicalMassVarianceClass              ), intent(in   ), target :: cosmologicalMassVariance_
    class(cosmologyParametersClass                   ), intent(in   ), target :: cosmologyParameters_
    !# <constructorAssign variables="*haloEnvironment_, *cosmologicalMassVariance_, *cosmologyParameters_"/>

    ! Evaluate and store the variance of the background.
    self%massBackground    =self%haloEnvironment_         %environmentMass(                   )
    self%varianceBackground=self%cosmologicalMassVariance_%rootVariance   (self%massBackground)**2
    return
  end function variancePeakBackgroundSplitConstructorInternal

  subroutine variancePeakBackgroundSplitDestructor(self)
    !% Destructor for the {\normalfont \ttfamily peakBackgroundSplit} linear growth class.
    implicit none
    type(cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    
    !# <objectDestructor name="self%haloEnvironment_"         />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%cosmologyParameters_"     />
    return
  end subroutine variancePeakBackgroundSplitDestructor

  double precision function variancePeakBackgroundSplitPowerNormalization(self)
    !% Return the normalization of the power spectrum.
    use Galacticus_Error
    implicit none
    class(cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    !GCC$ attributes unused :: self

    variancePeakBackgroundSplitPowerNormalization=0.0d0
    call Galacticus_Error_Report('power spectrum normalization is not well-defined in peak-background split model'//{introspection:location})
    return
  end function variancePeakBackgroundSplitPowerNormalization
  
  double precision function variancePeakBackgroundSplitSigma8(self)
    !% Return the value of $\sigma_8$.
    use Numerical_Constants_Math
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    ! Radius for σ(M) normalization in Mpc/h.
    double precision                                             , parameter     :: radiusNormalization=8.0d0

    variancePeakBackgroundSplitSigma8=self%rootVariance(                                                                 &
         &                                              +(                                                               &
         &                                                +4.0d0                                                         &
         &                                                /3.0d0                                                         &
         &                                                *Pi                                                            &
         &                                               )                                                               &
         &                                              *  self%cosmologyParameters_%OmegaMatter    (                  ) &
         &                                              *  self%cosmologyParameters_%densityCritical(                  ) &
         &                                              *(                                                               &
         &                                                +radiusNormalization                                           &
         &                                                /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
         &                                               )**3                                                            &
         &                                             )
    return
  end function variancePeakBackgroundSplitSigma8
  
  double precision function variancePeakBackgroundSplitRootVariance(self,mass)
    !% Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    !% \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass
    double precision                                                             :: varianceTotal

    varianceTotal=self%cosmologicalMassVariance_%rootVariance(mass)**2
    if (varianceTotal > self%varianceBackground) then
       variancePeakBackgroundSplitRootVariance=+sqrt(                         &
            &                                              varianceTotal      &
            &                                        -self%varianceBackground &
            &                                       )
    else
       variancePeakBackgroundSplitRootVariance=+0.0d0
    end if
    return
  end function variancePeakBackgroundSplitRootVariance
  
  double precision function variancePeakBackgroundSplitRootVarianceLogarithmicGradient(self,mass)
    !% Return the logairhtmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    !% region containing the given {\normalfont \ttfamily mass} on average.
    use Numerical_Constants_Math
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass
    double precision                                                             :: varianceTotal

    varianceTotal=self%cosmologicalMassVariance_%rootVariance(mass)**2
    if (varianceTotal > self%varianceBackground) then
       variancePeakBackgroundSplitRootVarianceLogarithmicGradient=+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass)    &
            &                                                     *                               varianceTotal                            &
            &                                                     /self                          %rootVariance                   (mass)**2
    else
       variancePeakBackgroundSplitRootVarianceLogarithmicGradient=+0.0d0
    end if
    return
  end function variancePeakBackgroundSplitRootVarianceLogarithmicGradient
  
  subroutine variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient(self,mass,rootVariance,rootVarianceLogarithmicGradient)
    !% Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    !% spherical region containing the given {\normalfont \ttfamily mass} on average.
    use Numerical_Constants_Math
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass
    double precision                                             , intent(  out) :: rootVariance , rootVarianceLogarithmicGradient
    double precision                                                             :: varianceTotal

    varianceTotal=self%cosmologicalMassVariance_%rootVariance(mass)**2
    if (varianceTotal > self%varianceBackground) then
       rootVariance                   =+sqrt(                                                                         &
            &                                +     varianceTotal                                                      &
            &                                -self%varianceBackground                                                 &
            &                               )
       rootVarianceLogarithmicGradient=+      self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass)    &
            &                          *                                     varianceTotal                            &
            &                          /                                     rootVariance                         **2
    else
       rootVariance                   =+0.0d0
       rootVarianceLogarithmicGradient=+0.0d0
    end if
    return
  end subroutine variancePeakBackgroundSplitRootVarianceAndLogarithmicGradient

  double precision function variancePeakBackgroundSplitMass(self,rootVariance)
    !% Return the mass corrresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    implicit none
    class           (cosmologicalMassVariancePeakBackgroundSplit), intent(inout) :: self
    double precision                                             , intent(in   ) :: rootVariance

    variancePeakBackgroundSplitMass=self%cosmologicalMassVariance_%mass(                                  &
         &                                                              +sqrt(                            &
         &                                                                    +     rootVariance      **2 &
         &                                                                    +self%varianceBackground    &
         &                                                                   )                            &
         &                                                             )
    return
  end function variancePeakBackgroundSplitMass
 
