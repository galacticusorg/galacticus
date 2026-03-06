!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  An implementation of cosmological density field mass variance which adds variance to mimic that associated with the formation of
  artificial halos.
  !!}
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Linear_Growth      , only : linearGrowthClass

  !![
  <cosmologicalMassVariance name="cosmologicalMassVarianceArtificialHalos">
   <description>
    The mass variance of cosmological density fields is computed by adding to the variance of another class to mimic that
    associated with the formation of artificial halos. Specifically,
    \begin{equation}
     \sigma^2(M,t) \rightarrow \sigma^2(M,t) + S_0 \left(\frac{M}{M_0}\right)^\alpha D^2(t),
    \end{equation}    
    where $S_0=${\normalfont \ttfamily [normalization]}, $M_0=${\normalfont \ttfamily [massZeroPoint]}, and $\alpha=${\normalfont
    \ttfamily [exponent]}.
   </description>
  </cosmologicalMassVariance>
  !!]
  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceArtificialHalos
     !!{
     A cosmological mass variance class that scales the variance from another class.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (linearGrowthClass            ), pointer :: linearGrowth_             => null()
     double precision                                         :: normalization                      , massZeroPoint, &
          &                                                      exponent
   contains
     final     ::                                        artificialHalosDestructor
     procedure :: sigma8                              => artificialHalosSigma8
     procedure :: powerNormalization                  => artificialHalosPowerNormalization
     procedure :: rootVariance                        => artificialHalosRootVariance
     procedure :: rootVarianceLogarithmicGradient     => artificialHalosRootVarianceLogarithmicGradient
     procedure :: rootVarianceLogarithmicGradientTime => artificialHalosRootVarianceLogarithmicGradientTime
     procedure :: rootVarianceAndLogarithmicGradient  => artificialHalosRootVarianceAndLogarithmicGradient
     procedure :: mass                                => artificialHalosMass   
     procedure :: growthIsMassDependent               => artificialHalosGrowthIsMassDependent
     procedure :: descriptorNormalizationOnly         => artificialHalosDescriptorNormalizationOnly
  end type cosmologicalMassVarianceArtificialHalos

  interface cosmologicalMassVarianceArtificialHalos
     !!{
     Constructors for the {\normalfont \ttfamily artificialHalos} cosmological mass variance class.
     !!}
     module procedure artificialHalosConstructorParameters
     module procedure artificialHalosConstructorInternal
  end interface cosmologicalMassVarianceArtificialHalos

contains

  function artificialHalosConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily artificialHalos} cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologicalMassVarianceArtificialHalos)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                      ), pointer       :: linearGrowth_
    double precision                                                         :: normalization            , massZeroPoint, &
         &                                                                      exponent
    
    !![
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <description></description>
    </inputParameter>
    <inputParameter>
      <name>massZeroPoint</name>
      <source>parameters</source>
      <description></description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <description></description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    self=cosmologicalMassVarianceArtificialHalos(normalization, massZeroPoint, exponent,cosmologicalMassVariance_,cosmologyFunctions_,linearGrowth_)
    !![
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="linearGrowth_"            />
    !!]    
    return
  end function artificialHalosConstructorParameters

  function artificialHalosConstructorInternal(normalization, massZeroPoint, exponent,cosmologicalMassVariance_,cosmologyFunctions_,linearGrowth_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily artificialHalos} cosmological mass variance class.
    !!}
    implicit none
    type            (cosmologicalMassVarianceArtificialHalos)                        :: self
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                      ), intent(in   ), target :: linearGrowth_
    double precision                                         , intent(in   )         :: normalization            , massZeroPoint, &
         &                                                                              exponent
    !![
    <constructorAssign variables="normalization, massZeroPoint, exponent, *cosmologicalMassVariance_, *cosmologyFunctions_, *linearGrowth_"/>
    !!]
    
    return
  end function artificialHalosConstructorInternal

  subroutine artificialHalosDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily artificialHalos} cosmological mass variance class.
    !!}
    implicit none
    type   (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyFunctions_"      />
    !!]
    return
  end subroutine artificialHalosDestructor

  double precision function artificialHalosPowerNormalization(self)
    !!{
    Return the normalization of the power spectrum. This is left unchanged by the presence of artificial halo.
    !!}
    implicit none
    class(cosmologicalMassVarianceArtificialHalos), intent(inout) :: self

    artificialHalosPowerNormalization=self%cosmologicalMassVariance_%powerNormalization()
    return
  end function artificialHalosPowerNormalization

  double precision function artificialHalosSigma8(self)
    !!{
    Return the value of $\sigma_8$. This is left unchanged by the presence of artificial halo.
    !!}
    implicit none
    class(cosmologicalMassVarianceArtificialHalos), intent(inout) :: self

    artificialHalosSigma8=self%cosmologicalMassVariance_%sigma8()
    return
  end function artificialHalosSigma8

  double precision function artificialHalosRootVariance(self,mass,time)
    !!{
    Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self
    double precision                                         , intent(in   ) :: mass, time

    artificialHalosRootVariance=+sqrt(                                                           &
         &                            +self%cosmologicalMassVariance_%rootVariance(mass,time)**2 &
         &                            +self%normalization                                        &
         &                            *(                                                         &
         &                              +     mass                                               &
         &                              /self%massZeroPoint                                      &
         &                             )**self%exponent                                          &
         &                            *self%linearGrowth_            %value       (     time)**2 &
         &                            )
    return
  end function artificialHalosRootVariance

  double precision function artificialHalosRootVarianceLogarithmicGradient(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self
    double precision                                         , intent(in   ) :: mass            , time
    double precision                                                         :: varianceOriginal, varianceArtificial

    varianceOriginal                              =+self%cosmologicalMassVariance_%rootVariance(mass,time)**2
    varianceArtificial                            =+self%normalization                                        &
         &                                         *(                                                         &
         &                                           +     mass                                               &
         &                                           /self%massZeroPoint                                      &
         &                                          )**self%exponent                                          &
         &                                         *self%linearGrowth_            %value       (     time)**2
    artificialHalosRootVarianceLogarithmicGradient=+(                                                                                                    &
         &                                           +varianceOriginal  *self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time)       &
         &                                           +varianceArtificial*self%exponent                                                            /2.0d0 &
         &                                          )                                                                                                    &
         &                                         /(                                                                                                    &
         &                                           +varianceOriginal                                                                                   &
         &                                           +varianceArtificial                                                                                 &
         &                                          )
    return
  end function artificialHalosRootVarianceLogarithmicGradient

  double precision function artificialHalosRootVarianceLogarithmicGradientTime(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to time of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self
    double precision                                         , intent(in   ) :: mass            , time
    double precision                                                         :: varianceOriginal, varianceArtificial

    varianceOriginal                              =+self%cosmologicalMassVariance_%rootVariance(mass,time)**2
    varianceArtificial                            =+self%normalization                                        &
         &                                         *(                                                         &
         &                                           +     mass                                               &
         &                                           /self%massZeroPoint                                      &
         &                                          )**self%exponent                                          &
         &                                         *self%linearGrowth_            %value       (     time)**2
    artificialHalosRootVarianceLogarithmicGradientTime=+(                                                                                                                                              &
         &                                               +varianceOriginal  *self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime (mass,                                          time)  &
         &                                               +varianceArtificial*self%linearGrowth_            %logarithmicDerivativeExpansionFactor(                                               time)  &
         &                                               *                   self%cosmologyFunctions_      %expansionRate                       (     self%cosmologyFunctions_ %expansionFactor(time)) &
         &                                               *                                                                                                                                      time   &
         &                                              )                                                                                                                                              &
         &                                             /(                                                                                                                                              &
         &                                               +varianceOriginal                                                                                                                             &
         &                                               +varianceArtificial                                                                                                                           &
         &                                              )
   return
  end function artificialHalosRootVarianceLogarithmicGradientTime

  subroutine artificialHalosRootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient)
    !!{
    Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    spherical region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self
    double precision                                         , intent(in   ) :: mass            , time
    double precision                                         , intent(  out) :: rootVariance    , rootVarianceLogarithmicGradient
    double precision                                                         :: varianceOriginal, varianceArtificial

    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,rootVarianceLogarithmicGradient)
    varianceOriginal              =+self%cosmologicalMassVariance_%rootVariance(mass,time)**2
    varianceArtificial            =+self%normalization                                        &
         &                         *(                                                         &
         &                           +     mass                                               &
         &                           /self%massZeroPoint                                      &
         &                          )**self%exponent                                          &
         &                         *self%linearGrowth_            %value       (     time)**2
    rootVarianceLogarithmicGradient=+(                                                                                                    &
         &                            +varianceOriginal  *self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time)       &
         &                            +varianceArtificial*self%exponent                                                            /2.0d0 &
         &                           )                                                                                                    &
         &                          /(                                                                                                    &
         &                            +varianceOriginal                                                                                   &
         &                            +varianceArtificial                                                                                 &
         &                           )
    rootVariance                   =+sqrt(                    &
         &                                +varianceOriginal   &
         &                                +varianceArtificial &
         &                               )
    return
  end subroutine artificialHalosRootVarianceAndLogarithmicGradient

  double precision function artificialHalosMass(self,rootVariance,time)
    !!{
    Return the mass corresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    !!}
    implicit none
    class           (cosmologicalMassVarianceArtificialHalos), intent(inout) :: self
    double precision                                         , intent(in   ) :: rootVariance, time

    artificialHalosMass=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function artificialHalosMass

  logical function artificialHalosGrowthIsMassDependent(self)
    !!{
    Return true if the growth rate of the variance is mass-dependent.
    !!}
    implicit none
    class(cosmologicalMassVarianceArtificialHalos), intent(inout) :: self

    artificialHalosGrowthIsMassDependent=self%cosmologicalMassVariance_%growthIsMassDependent()
    return
  end function artificialHalosGrowthIsMassDependent

  subroutine artificialHalosDescriptorNormalizationOnly(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object, for power spectrum normalization usage
    only (i.e. we exclude the window function).
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (cosmologicalMassVarianceArtificialHalos), intent(inout)           :: self
    type     (inputParameters                        ), intent(inout)           :: descriptor
    logical                                           , intent(in   ), optional :: includeClass, includeFileModificationTimes
    type     (inputParameters                        )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('cosmologicalMassVariance','filteredPower')
    parameters=descriptor%subparameters('cosmologicalMassVariance')
    if (associated(self%cosmologyFunctions_      )) &
         & call self%cosmologyFunctions_      %descriptor                 (parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    if (associated(self%cosmologicalMassVariance_)) &
         & call self%cosmologicalMassVariance_%descriptorNormalizationOnly(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    if (associated(self%linearGrowth_            )) &
         & call self%linearGrowth_            %descriptor                 (parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    return
  end subroutine artificialHalosDescriptorNormalizationOnly
