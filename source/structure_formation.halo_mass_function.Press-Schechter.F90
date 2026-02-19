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
Implements a \cite{press_formation_1974} dark matter halo mass function class.
!!}
  use :: Cosmological_Density_Field    , only : cosmologicalMassVarianceClass
  use :: Excursion_Sets_First_Crossings, only : excursionSetFirstCrossingClass

  !![
  <haloMassFunction name="haloMassFunctionPressSchechter">
   <description>
    A dark matter halo mass function class using the function given by \cite{press_formation_1974}. Specifically,
    \begin{equation}
    n(M,t) = 2 {\Omega_\mathrm{M} \rho_\mathrm{crit} \over M^2} \alpha \sigma^2(M) f[S(M,t)],
    \end{equation}
    where $\alpha = \mathrm{d}\ln\sigma/\mathrm{d}\ln M$ and $f[S]$ is the excursion set barrier first crossing distribution
    for variance $S(M)=\sigma^2(M)$, computed using the selected \refClass{excursionSetFirstCrossingClass}.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionPressSchechter
     !!{
     A halo mass function class using the model of \cite{press_formation_1974}.
     !!}
     private
     class(cosmologicalMassVarianceClass ), pointer :: cosmologicalMassVariance_  => null()
     class(excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_ => null()
    contains
     final     ::                 pressSchechterDestructor
     procedure :: differential => pressSchechterDifferential
  end type haloMassFunctionPressSchechter

  interface haloMassFunctionPressSchechter
     !!{
     Constructors for the \refClass{haloMassFunctionPressSchechter} halo mass function class.
     !!}
     module procedure pressSchechterConstructorParameters
     module procedure pressSchechterConstructorInternal
  end interface haloMassFunctionPressSchechter

contains

  function pressSchechterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionPressSchechter} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionPressSchechter)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(cosmologicalMassVarianceClass ), pointer       :: cosmologicalMassVariance_
    class(excursionSetFirstCrossingClass), pointer       :: excursionSetFirstCrossing_
    class(cosmologyParametersClass      ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyParameters"       name="cosmologyParameters_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"  name="cosmologicalMassVariance_"  source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing" name="excursionSetFirstCrossing_" source="parameters"/>
    !!]
    self=haloMassFunctionPressSchechter(cosmologyParameters_,cosmologicalMassVariance_,excursionSetFirstCrossing_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"      />
    <objectDestructor name="cosmologicalMassVariance_" />
    <objectDestructor name="excursionSetFirstCrossing_"/>
    !!]
    return
  end function pressSchechterConstructorParameters

  function pressSchechterConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,excursionSetFirstCrossing_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionPressSchechter} halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionPressSchechter)                        :: self
    class(cosmologyParametersClass      ), target, intent(in   ) :: cosmologyParameters_
    class(cosmologicalMassVarianceClass ), target, intent(in   ) :: cosmologicalMassVariance_
    class(excursionSetFirstCrossingClass), target, intent(in   ) :: excursionSetFirstCrossing_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *excursionSetFirstCrossing_"/>
    !!]

    return
  end function pressSchechterConstructorInternal

  subroutine pressSchechterDestructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionPressSchechter} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionPressSchechter), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"       />
    <objectDestructor name="self%cosmologicalMassVariance_"  />
    <objectDestructor name="self%excursionSetFirstCrossing_" />
    !!]
    return
  end subroutine pressSchechterDestructor

  double precision function pressSchechterDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloMassFunctionPressSchechter), intent(inout), target   :: self
    double precision                                , intent(in   )           :: time , mass
    type            (treeNode                      ), intent(inout), optional :: node
    double precision                                                          :: alpha, variance

    if (.not.present(node)) call Error_Report('"node" must be present'//{introspection:location})
    alpha                     =abs(self%cosmologicalMassVariance_ %rootVarianceLogarithmicGradient(mass,time))
    variance                  =    self%cosmologicalMassVariance_ %rootVariance                   (mass,time) **2
    if (variance > 0.0d0) then
       pressSchechterDifferential=+2.0d0                                                                  &
            &                     *   self%cosmologyParameters_      %OmegaMatter    (                  ) &
            &                     *   self%cosmologyParameters_      %densityCritical(                  ) &
            &                     *   self%excursionSetFirstCrossing_%probability    (variance,time,node) &
            &                     /mass**2                                                                &
            &                     *alpha                                                                  &
            &                     *variance
    else
       pressSchechterDifferential=+0.0d0
    end if
    return
  end function pressSchechterDifferential
