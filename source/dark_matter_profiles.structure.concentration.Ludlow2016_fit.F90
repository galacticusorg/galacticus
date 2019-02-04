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

  !% An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2016}
  !% fitting function.
  
  use Cosmology_Functions
  use Cosmology_Parameters
  use Cosmological_Density_Field

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationLudlow2016Fit">
  !#  <description>Dark matter halo concentrations are computed using the fitting function of \cite{ludlow_mass-concentration-redshift_2016}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationLudlow2016Fit
     !% A dark matter halo profile concentration class implementing the fitting function of
     !% \cite{ludlow_mass-concentration-redshift_2016}.
     private
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class(cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     type (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type (darkMatterProfileEinasto     ), pointer :: darkMatterProfileDefinition_     => null()
   contains
     final     ::                                ludlow2016FitDestructor
     procedure :: concentration               => ludlow2016FitConcentration
     procedure :: densityContrastDefinition   => ludlow2016FitDensityContrastDefinition
     procedure :: darkMatterProfileDefinition => ludlow2016FitDarkMatterProfileDefinition
  end type darkMatterProfileConcentrationLudlow2016Fit

  interface darkMatterProfileConcentrationLudlow2016Fit
     !% Constructors for the {\normalfont \ttfamily ludlow2016Fit} dark matter halo profile concentration class.
     module procedure ludlow2016FitConstructorParameters
     module procedure ludlow2016FitConstructorInternal
  end interface darkMatterProfileConcentrationLudlow2016Fit

contains

  function ludlow2016FitConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily ludlow2016Fit} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type (darkMatterProfileConcentrationLudlow2016Fit)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                    ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass                   ), pointer       :: cosmologyParameters_     
    class(cosmologicalMassVarianceClass              ), pointer       :: cosmologicalMassVariance_

    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    self=darkMatterProfileConcentrationLudlow2016Fit(cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="cosmologyParameters_"     />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    return
  end function ludlow2016FitConstructorParameters
  
  function ludlow2016FitConstructorInternal(cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_) result(self)
    !% Constructor for the {\normalfont \ttfamily ludlow2016Fit} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type (darkMatterProfileConcentrationLudlow2016Fit)                        :: self
    class(cosmologyFunctionsClass                    ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass                   ), intent(in   ), target :: cosmologyParameters_     
    class(cosmologicalMassVarianceClass              ), intent(in   ), target :: cosmologicalMassVariance_
    !# <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *cosmologicalMassVariance_"/>

    return
  end function ludlow2016FitConstructorInternal

  subroutine ludlow2016FitDestructor(self)
    !% Destructor for the {\normalfont \ttfamily ludlow2016Fit} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationLudlow2016Fit), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"             />
    !# <objectDestructor name="self%cosmologyParameters_"            />
    !# <objectDestructor name="self%cosmologicalMassVariance_"       />
    !# <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !# <objectDestructor name="self%darkMatterProfileDefinition_"    />
    return
  end subroutine ludlow2016FitDestructor

  double precision function ludlow2016FitConcentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    !% \cite{ludlow_mass-concentration-redshift_2016} fitting function.
    use Galacticus_Error
    use Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileConcentrationLudlow2016Fit), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    class           (nodeComponentBasic                         ), pointer                :: basic
    double precision                                             , parameter              :: criticalOverdensitySphericalCollapse=1.686d0
    double precision                                                                      :: peakHeight                                  , expansionFactor, &
         &                                                                                   c0                                          , beta           , &
         &                                                                                   gamma1                                      , gamma2         , &
         &                                                                                   nu0
    
    basic           =>  node                          %basic          (                                            )
    peakHeight      =  +                                                     criticalOverdensitySphericalCollapse    &
         &             /self%cosmologicalMassVariance_%rootVariance   (basic%mass                                ())
    expansionFactor =   self%cosmologyFunctions_      %expansionFactor(basic%time                                ())
    if (expansionFactor < 0.1d0) call Galacticus_Error_Report('redshift out of range of fitting function'//{introspection:location})
    c0                        =+3.395d0*expansionFactor**0.215d0
    beta                      =+0.307d0/expansionFactor**0.540d0
    gamma1                    =+0.628d0*expansionFactor**0.047d0
    gamma2                    =+0.317d0*expansionFactor**0.893d0
    nu0                       =+(                              &
         &                       +4.13500d0                    &
         &                       -0.56400d0/expansionFactor    &
         &                       -0.21000d0/expansionFactor**2 &
         &                       +0.05570d0/expansionFactor**3 &
         &                       -0.00348d0/expansionFactor**4 &
         &                      )
    ludlow2016FitConcentration=+c0                                                             &
         &                     /       (peakHeight/nu0)               **              gamma1   &
         &                     /(1.0d0+(peakHeight/nu0)**(1.0d0/beta))**(beta*(gamma2-gamma1))
    return
  end function ludlow2016FitConcentration

  function ludlow2016FitDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    implicit none
    class(virialDensityContrastClass                 ), pointer       :: ludlow2016FitDensityContrastDefinition
    class(darkMatterProfileConcentrationLudlow2016Fit), intent(inout) :: self
    
    if (.not.associated(self%virialDensityContrastDefinition_)) then
       allocate(self%virialDensityContrastDefinition_)
       !# <referenceConstruct owner="self" object="virialDensityContrastDefinition_" constructor="virialDensityContrastFixed(200.0d0,fixedDensityTypeCritical,self%cosmologyFunctions_)"/>
    end if
    ludlow2016FitDensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function ludlow2016FitDensityContrastDefinition

  function ludlow2016FitDarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{ludlow_mass-concentration-redshift_2016} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: ludlow2016FitDarkMatterProfileDefinition
    class(darkMatterProfileConcentrationLudlow2016Fit       ), intent(inout) :: self
    type (darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition_
    class(virialDensityContrastClass                        ), pointer       :: virialDensityContrastDefinition_

    if (.not.associated(self%darkMatterProfileDefinition_)) then
       allocate(self%darkMatterProfileDefinition_  )
       allocate(     darkMatterHaloScaleDefinition_)
       !# <referenceAcquire                target="virialDensityContrastDefinition_" source     ="self%densityContrastDefinition()"/>
       !# <referenceConstruct              object="darkMatterHaloScaleDefinition_"   constructor="darkMatterHaloScaleVirialDensityContrastDefinition(self%cosmologyParameters_,self%cosmologyFunctions_,virialDensityContrastDefinition_)"/>
       !# <referenceConstruct owner="self" object="darkMatterProfileDefinition_"     constructor="darkMatterProfileEinasto                          (darkMatterHaloScaleDefinition_                                                     )"/>
       !# <objectDestructor name="darkMatterHaloScaleDefinition_"  />
       !# <objectDestructor name="virialDensityContrastDefinition_"/>
    end if
    ludlow2016FitDarkMatterProfileDefinition => self%darkMatterProfileDefinition_
    return
  end function ludlow2016FitDarkMatterProfileDefinition
