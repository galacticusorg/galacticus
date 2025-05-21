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
  An implementation of \cite{bryan_statistical_1998} dark matter halo virial density contrasts.
  !!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  ! Enumeration for different fitting function types.
  !![
  <enumeration>
   <name>bryanNorman1998Fit</name>
   <description>Specifies fit type for \cite{bryan_statistical_1998} virial density contrast.</description>
   <entry label="flatUniverse" />
   <entry label="zeroLambda"   />
  </enumeration>
  !!]

  !![
  <virialDensityContrast name="virialDensityContrastBryanNorman1998">
   <description>
    A dark matter halo virial density contrast class using the fitting functions given by \cite{bryan_statistical_1998}. As such,
    it is valid only for $\Omega_\Lambda=0$ or $\Omega_\mathrm{M}+\Omega_\Lambda=1$ cosmologies and will either abort on other
    cosmologies (if {\normalfont \ttfamily [allowUnsupportedCosmology]=false}), or revert to a numerical solution from top-hat
    collapse (if {\normalfont \ttfamily [allowUnsupportedCosmology]=true}).
   </description>
   <deepCopy>
     <functionClass variables="virialDensityContrastTopHat_"/>
   </deepCopy>
   <stateStorable>
     <functionClass variables="virialDensityContrastTopHat_"/>
   </stateStorable>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastBryanNorman1998
     !!{
     A dark matter halo virial density contrast class using the fitting functions of \cite{bryan_statistical_1998}.
     !!}
     private
     class  (cosmologyParametersClass                                      ), pointer :: cosmologyParameters_         => null()
     class  (cosmologyFunctionsClass                                       ), pointer :: cosmologyFunctions_          => null()
     type   (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: virialDensityContrastTopHat_ => null()
     type   (enumerationBryanNorman1998FitType                             )          :: fitType
     logical                                                                          :: useFittingFunction           =  .true., allowUnsupportedCosmology
   contains
     final     ::                                bryanNorman1998Destructor
     procedure :: densityContrast             => bryanNorman1998DensityContrast
     procedure :: densityContrastRateOfChange => bryanNorman1998DensityContrastRateOfChange
     procedure :: turnAroundOverVirialRadii   => bryanNorman1998TurnAroundOverVirialRadii
  end type virialDensityContrastBryanNorman1998

  interface virialDensityContrastBryanNorman1998
     !!{
     Constructors for the \refClass{virialDensityContrastBryanNorman1998} dark matter halo virial density contrast class.
     !!}
     module procedure bryanNorman1998ConstructorParameters
     module procedure bryanNorman1998ConstructorInternal
  end interface virialDensityContrastBryanNorman1998

contains

  function bryanNorman1998ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastBryanNorman1998} dark matter halo virial density contrast class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (virialDensityContrastBryanNorman1998)                :: self
    type   (inputParameters                     ), intent(inout) :: parameters
    class  (cosmologyParametersClass            ), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    logical                                                      :: allowUnsupportedCosmology

    !![
    <inputParameter>
      <name>allowUnsupportedCosmology</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, unsupported cosmologies revert to using a numerical solution from the top-hat collapse model. Otherwise, unsupported cosmologies result in an error.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=virialDensityContrastBryanNorman1998(allowUnsupportedCosmology,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function bryanNorman1998ConstructorParameters

  function bryanNorman1998ConstructorInternal(allowUnsupportedCosmology,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{virialDensityContrastBryanNorman1998} dark matter halo virial density contrast class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type   (virialDensityContrastBryanNorman1998)                        :: self
    class  (cosmologyParametersClass            ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass             ), intent(in   ), target :: cosmologyFunctions_
    logical                                      , intent(in   )         :: allowUnsupportedCosmology
    !![
    <constructorAssign variables="allowUnsupportedCosmology, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Check that fitting formulae are applicable to this cosmology.
    self%useFittingFunction=.true.
    if (self%cosmologyParameters_%OmegaDarkEnergy() == 0.0d0) then
       self%fitType=bryanNorman1998FitZeroLambda
    else if (.not.Values_Differ(self%cosmologyParameters_%OmegaMatter()+self%cosmologyParameters_%OmegaDarkEnergy(),1.0d0,absTol=1.0d-6)) then
       self%fitType=bryanNorman1998FitFlatUniverse
    else if (self%allowUnsupportedCosmology) then
       self%useFittingFunction=.false.
       allocate(self%virialDensityContrastTopHat_)
       !![
       <referenceConstruct isResult="yes" owner="self" object="virialDensityContrastTopHat_">
	 <constructor>
	   virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                              &amp;
	    &amp;                                                         tableStore         =.true.                  , &amp;
	    &amp;                                                         cosmologyFunctions_=self%cosmologyFunctions_  &amp;
	    &amp;                                                        )
	 </constructor>
       </referenceConstruct>
       !!]
    else
       call Error_Report('no fitting formula available for this cosmology'//{introspection:location})
    end if
    return
  end function bryanNorman1998ConstructorInternal

  subroutine bryanNorman1998Destructor(self)
    !!{
    Destructor for the \refClass{virialDensityContrastBryanNorman1998} virial density contrast class.
    !!}
    implicit none
    type(virialDensityContrastBryanNorman1998), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%cosmologyFunctions_"  />
    !!]
    if (.not.self%useFittingFunction) then
       !![
       <objectDestructor name="self%virialDensityContrastTopHat_" />
       !!]
    end if
    return
  end subroutine bryanNorman1998Destructor

  double precision function bryanNorman1998DensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming the fitting function of \cite{bryan_statistical_1998}.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (virialDensityContrastBryanNorman1998), intent(inout)           :: self
    double precision                                      , intent(in   )           :: mass
    double precision                                      , intent(in   ), optional :: time      , expansionFactor
    logical                                               , intent(in   ), optional :: collapsing
    double precision                                                                :: x

    if (self%useFittingFunction) then
       x=self%cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0
       select case (self%fitType%ID)
       case (bryanNorman1998FitZeroLambda  %ID)
          bryanNorman1998DensityContrast=(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)/(x+1.0d0)
       case (bryanNorman1998FitFlatUniverse%ID)
          bryanNorman1998DensityContrast=(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/(x+1.0d0)
       case default
          bryanNorman1998DensityContrast=0.0d0
          call Error_Report('invalid fit type'//{introspection:location})
       end select
    else
       bryanNorman1998DensityContrast=self%virialDensityContrastTopHat_%densityContrast(mass,time,expansionFactor,collapsing)
    end if
    return
  end function bryanNorman1998DensityContrast

  double precision function bryanNorman1998DensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming the fitting function of \cite{bryan_statistical_1998}.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (virialDensityContrastBryanNorman1998), intent(inout)           :: self
    double precision                                      , intent(in   )           :: mass
    double precision                                      , intent(in   ), optional :: time      , expansionFactor
    logical                                               , intent(in   ), optional :: collapsing
    double precision                                                                :: x

    if (self%useFittingFunction) then
       x=self%cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0
       select case (self%fitType%ID)
       case (bryanNorman1998FitZeroLambda  %ID)
          bryanNorman1998DensityContrastRateOfChange=                                               &
               & (                                                                                  &
               &  +(            +60.0d0  -64.0d0*x   )                                              &
               &  -(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)                                              &
               &  /(x+1.0d0)                                                                        &
               & )                                                                                  &
               & *self%cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing) &
               & / (x+1.0d0)
       case (bryanNorman1998FitFlatUniverse%ID)
          bryanNorman1998DensityContrastRateOfChange=                                               &
               & (                                                                                  &
               &  +(            +82.0d0  -78.0d0*x   )                                              &
               &  -(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)                                              &
               &  /(x+1.0d0)                                                                        &
               & )                                                                                  &
               & *self%cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing) &
               & / (x+1.0d0)
       case default
          bryanNorman1998DensityContrastRateOfChange=0.0d0
          call Error_Report('invalid fit type'//{introspection:location})
       end select
    else
       bryanNorman1998DensityContrastRateOfChange=self%virialDensityContrastTopHat_%densityContrastRateOfChange(mass,time,expansionFactor,collapsing)
    end if
    return
  end function bryanNorman1998DensityContrastRateOfChange

  double precision function bryanNorman1998TurnAroundOverVirialRadii(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    constant universe.
    !!}
    implicit none
    class           (virialDensityContrastBryanNorman1998), intent(inout)           :: self
    double precision                                      , intent(in   )           :: mass
    double precision                                      , intent(in   ), optional :: time      , expansionFactor
    logical                                               , intent(in   ), optional :: collapsing

    if (self%useFittingFunction) then
       ! In simple cosmological constant dark energy universes, this ratio is always precisely 2 (e.g. Percival 2005;
       ! http://adsabs.harvard.edu/abs/2005A%26A...443..819P)
       bryanNorman1998TurnAroundOverVirialRadii=2.0d0
    else
       bryanNorman1998TurnAroundOverVirialRadii=self%turnAroundOverVirialRadii(mass,time,expansionFactor,collapsing)
    end if
    return
  end function bryanNorman1998TurnAroundOverVirialRadii
