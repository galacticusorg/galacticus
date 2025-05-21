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
  Implements a stellar feedback class which performs a simple calculation of energy feedback from stellar populations.
  !!}

  use :: Stellar_Astrophysics      , only : stellarAstrophysicsClass
  use :: Stellar_Astrophysics_Winds, only : stellarWindsClass
  use :: Supernovae_Population_III , only : supernovaePopulationIIIClass
  use :: Supernovae_Type_Ia        , only : supernovaeTypeIaClass

  !![
  <stellarFeedback name="stellarFeedbackStandard">
   <description>
    A stellar feedback class which assumes that the cumulative energy input from a stellar population is equal to the total
    number of (Type II and Type Ia) supernovae multiplied by {\normalfont \ttfamily [supernovaEnergy]} (specified in ergs) plus
    any Population III-specific supernovae energy plus the integrated energy input from stellar winds. The minimum mass of a
    star required to form a Type II supernova is specified (in $M_\odot$) via the {\normalfont \ttfamily
    [initialMassForSupernovaeTypeII]} parameter.
   </description>
  </stellarFeedback>
  !!]
  type, extends(stellarFeedbackClass) :: stellarFeedbackStandard
     !!{
     A stellar feedback class which performs a simple calculation of energy feedback from stellar populations.
     !!}
     private
     class           (supernovaeTypeIaClass       ), pointer :: supernovaeTypeIa_              => null()
     class           (supernovaePopulationIIIClass), pointer :: supernovaePopulationIII_       => null()
     class           (stellarWindsClass           ), pointer :: stellarWinds_                  => null()
     class           (stellarAstrophysicsClass    ), pointer :: stellarAstrophysics_           => null()
     double precision                                        :: initialMassForSupernovaeTypeII          , supernovaEnergy
   contains
     final     ::                          standardDestructor
     procedure :: energyInputCumulative => standardEnergyInputCumulative
  end type stellarFeedbackStandard

  interface stellarFeedbackStandard
     !!{
     Constructors for the \refClass{stellarFeedbackStandard} stellar feedback class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface stellarFeedbackStandard

  ! Module-scope variables used in integrands.
  class           (stellarFeedbackStandard), pointer :: self_
  double precision                                   :: massInitial_, metallicity_
  !$omp threadprivate(self_,massInitial_,metallicity_)

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarFeedbackStandard} stellar feedback class which takes a parameter list as input.
    !!}
    use :: Input_Parameters                , only : inputParameter, inputParameters
    use :: Numerical_Constants_Astronomical, only : massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    type            (stellarFeedbackStandard     )                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (supernovaeTypeIaClass       ), pointer       :: supernovaeTypeIa_
    class           (supernovaePopulationIIIClass), pointer       :: supernovaePopulationIII_
    class           (stellarWindsClass           ), pointer       :: stellarWinds_
    class           (stellarAstrophysicsClass    ), pointer       :: stellarAstrophysics_
    double precision                                              :: initialMassForSupernovaeTypeII, supernovaEnergy

    !![
    <inputParameter>
      <name>initialMassForSupernovaeTypeII</name>
      <defaultValue>8.0d0</defaultValue>
      <description>The minimum mass that a star must have in order that is result in a Type II supernova.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>supernovaEnergy</name>
      <defaultValue>1.0d51</defaultValue>
      <description>The energy produced by a supernova (in ergs).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    ! Convert energy to M☉ (km/s)².
    supernovaEnergy=supernovaEnergy*ergs/massSolar/kilo**2
    !![
    <objectBuilder class="supernovaeTypeIa"        name="supernovaeTypeIa_"        source="parameters"/>
    <objectBuilder class="supernovaePopulationIII" name="supernovaePopulationIII_" source="parameters"/>
    <objectBuilder class="stellarWinds"            name="stellarWinds_"            source="parameters"/>
    <objectBuilder class="stellarAstrophysics"     name="stellarAstrophysics_"     source="parameters"/>
    !!]
    self=stellarFeedbackStandard(initialMassForSupernovaeTypeII,supernovaEnergy,supernovaeTypeIa_,supernovaePopulationIII_,stellarWinds_,stellarAstrophysics_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="supernovaeTypeIa_"       />
    <objectDestructor name="supernovaePopulationIII_"/>
    <objectDestructor name="stellarWinds_"           />
    <objectDestructor name="stellarAstrophysics_"    />
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(initialMassForSupernovaeTypeII,supernovaEnergy,supernovaeTypeIa_,supernovaePopulationIII_,stellarWinds_,stellarAstrophysics_) result(self)
    !!{
    Constructor for the \refClass{stellarFeedbackStandard} stellar feedback class which takes a parameter list as input.
    !!}
    implicit none
    type            (stellarFeedbackStandard     )                        :: self
    class           (supernovaeTypeIaClass       ), intent(in   ), target :: supernovaeTypeIa_
    class           (supernovaePopulationIIIClass), intent(in   ), target :: supernovaePopulationIII_
    class           (stellarWindsClass           ), intent(in   ), target :: stellarWinds_
    class           (stellarAstrophysicsClass    ), intent(in   ), target :: stellarAstrophysics_
    double precision                              , intent(in   )         :: initialMassForSupernovaeTypeII, supernovaEnergy
    !![
    <constructorAssign variables="initialMassForSupernovaeTypeII, supernovaEnergy, *supernovaeTypeIa_, *supernovaePopulationIII_, *stellarWinds_, *stellarAstrophysics_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
   !!{
   Destructor for the \refClass{stellarFeedbackStandard} stellar feedback class.
   !!}
    implicit none
    type(stellarFeedbackStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%supernovaeTypeIa_"       />
    <objectDestructor name="self%supernovaePopulationIII_"/>
    <objectDestructor name="self%stellarWinds_"           />
    <objectDestructor name="self%stellarAstrophysics_"    />
    !!]
    return
  end subroutine standardDestructor

  double precision function standardEnergyInputCumulative(self,initialMassFunction_,initialMass,age,metallicity)
    !!{
    Compute the cumulative energy input from a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age} and {\normalfont \ttfamily metallicity}.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    use :: Numerical_Integration           , only : integrator
    implicit none
    class           (stellarFeedbackStandard ), intent(inout), target :: self
    double precision                          , intent(in   )         :: age                                                    , initialMass, metallicity
    class           (initialMassFunctionClass), intent(inout)         :: initialMassFunction_
    double precision                          , parameter             :: populationIIIMaximumMetallicity=1.0d-4*metallicitySolar
    double precision                                                  :: energySNe                                              , energyWinds, lifetime
    type            (integrator              )                        :: integrator_
    
    ! Begin with zero energy input.
    standardEnergyInputCumulative=0.0d0
    ! Check if the star is sufficiently massive to result in a Type II supernova.
    if (initialMass >= self%initialMassForSupernovaeTypeII) then
       ! Get the lifetime of the star.
       lifetime=self%stellarAstrophysics_%lifetime(initialMass,metallicity)
       ! If lifetime is exceeded, assume a SNe has occurred.
       if (age >= lifetime) then
          energySNe=0.0d0
          ! Check for pair instability supernovae.
          if (metallicity <= populationIIIMaximumMetallicity) energySNe=energySNe+self%supernovaePopulationIII_%energyCumulative(initialMass,age,metallicity)
          ! Population II star - normal supernova.
          energySNe=energySNe+self%supernovaEnergy
          ! Add the supernova energy.
          standardEnergyInputCumulative=standardEnergyInputCumulative+energySNe
       end if
    end if
    ! Add in contribution from Type Ia supernovae.
    standardEnergyInputCumulative=+standardEnergyInputCumulative                                                            &
         &                        +self%supernovaeTypeIa_%number         (initialMassFunction_,initialMass,age,metallicity) &
         &                        *self                  %supernovaEnergy
    ! Add in the contribution from stellar winds.
    self_                         =>  self
    massInitial_                  =   initialMass
    metallicity_                  =   metallicity
    integrator_                   =   integrator           (standardWindEnergyIntegrand,toleranceAbsolute=1.0d-3*standardEnergyInputCumulative,toleranceRelative=1.0d-3)
    energyWinds                   =   integrator_%integrate(0.0d0,age)
    standardEnergyInputCumulative =  +standardEnergyInputCumulative &
         &                           +energyWinds
    return
  end function standardEnergyInputCumulative

  double precision function standardWindEnergyIntegrand(age)
    !!{
    Integrand used in evaluating cumulative energy input from winds.
    !!}
    implicit none
    double precision, intent(in   ) :: age

    standardWindEnergyIntegrand=+0.5d0                                                                  &
         &                      *self_%stellarWinds_%rateMassLoss    (massInitial_,age,metallicity_)    &
         &                      *self_%stellarWinds_%velocityTerminal(massInitial_,age,metallicity_)**2
    return
  end function standardWindEnergyIntegrand
