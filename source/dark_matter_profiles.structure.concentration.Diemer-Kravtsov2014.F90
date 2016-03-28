!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of dark matter halo profile concentrations using the
  !% \cite{diemer_universal_2014} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationDiemerKravtsov2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{diemer_universal_2014}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDiemerKravtsov2014
     !% A dark matter halo profile concentration class implementing the algorithm of
     !% \cite{diemer_universal_2014}.
     private
     double precision :: kappa, phi0, phi1, eta0, eta1, alpha, beta
   contains
     final     ::                                diemerKravtsov2014Destructor
     procedure :: concentration               => diemerKravtsov2014Concentration
     procedure :: densityContrastDefinition   => diemerKravtsov2014DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => diemerKravtsov2014DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDiemerKravtsov2014

  interface darkMatterProfileConcentrationDiemerKravtsov2014
     !% Constructors for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile concentration class.
     module procedure diemerKravtsov2014ConstructorParameters
     module procedure diemerKravtsov2014ConstructorInternal
  end interface darkMatterProfileConcentrationDiemerKravtsov2014

contains

  function diemerKravtsov2014ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo
    !% profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014)                :: diemerKravtsov2014ConstructorParameters
    type(inputParameters                                 ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>kappa</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%kappa</variable>
    !#   <defaultValue>0.69d0</defaultValue>
    !#   <description>The parameter $\kappa$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>phi0</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%phi0</variable>
    !#   <defaultValue>6.58d0</defaultValue>
    !#   <description>The parameter $\phi_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>phi1</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%phi1</variable>
    !#   <defaultValue>1.37d0</defaultValue>
    !#   <description>The parameter $\phi_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta0</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%eta0</variable>
    !#   <defaultValue>6.82d0</defaultValue>
    !#   <description>The parameter $\eta_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta1</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%eta1</variable>
    !#   <defaultValue>1.42d0</defaultValue>
    !#   <description>The parameter $\eta_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%alpha</variable>
    !#   <defaultValue>1.12d0</defaultValue>
    !#   <attachedTo>module</attachedTo>
    !#   <description>The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%beta</variable>
    !#   <defaultValue>1.69d0</defaultValue>
    !#   <description>The parameter $\beta$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    return
  end function diemerKravtsov2014ConstructorParameters

  function diemerKravtsov2014ConstructorInternal(kappa,phi0,phi1,eta0,eta1,alpha,beta)
    !% Constructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile
    !% concentration class.
    use Galacticus_Error
    implicit none
    type            (darkMatterProfileConcentrationDiemerKravtsov2014)                :: diemerKravtsov2014ConstructorInternal
    double precision                                                  , intent(in   ) :: kappa, phi0, phi1, eta0, eta1, alpha, beta
    
    diemerKravtsov2014ConstructorInternal%kappa=kappa
    diemerKravtsov2014ConstructorInternal%phi0 =phi0
    diemerKravtsov2014ConstructorInternal%phi1 =phi1
    diemerKravtsov2014ConstructorInternal%eta0 =eta0
    diemerKravtsov2014ConstructorInternal%eta1 =eta1
    diemerKravtsov2014ConstructorInternal%alpha=alpha
    diemerKravtsov2014ConstructorInternal%beta =beta
    return
  end function diemerKravtsov2014ConstructorInternal
  
  subroutine diemerKravtsov2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile
    !% concentration class.
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine diemerKravtsov2014Destructor

  double precision function diemerKravtsov2014Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    !% using the \cite{diemer_universal_2014} algorithm.
    use Numerical_Constants_Math
    use Power_Spectra
    use Critical_Overdensity
    use Cosmology_Parameters
    implicit none
    class           (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)          :: self
    type            (treeNode                                        ), intent(inout), pointer :: node
    class           (nodeComponentBasic                              )               , pointer :: basic
    class           (cosmologyParametersClass                        )               , pointer :: cosmologyParameters_
    double precision                                                                           :: radiusHaloLagrangian, peakHeight        , &
         &                                                                                        wavenumber          , powerSpectrumSlope, &
         &                                                                                        concentrationMinimum, peakHeightMinimum
    
    cosmologyParameters_           => cosmologyParameters      ()
    basic                          => node               %basic()
    radiusHaloLagrangian           =  +(                                                &
         &                              +3.0d0                                          &
         &                              *basic%mass()                                   &
         &                              /4.0d0                                          &
         &                              /Pi                                             &
         &                              /cosmologyParameters_%densityCritical()         &
         &                              /cosmologyParameters_%OmegaMatter    ()         &
         &                             )**(1.0d0/3.0d0)
    peakHeight                     =                                                    &
         &                            +Critical_Overdensity_for_Collapse(basic%time())  &
         &                            /Cosmological_Mass_Root_Variance  (basic%mass())
    wavenumber                     =                                                    &
         &                            +self%kappa                                       &
         &                            *2.0d0                                            &
         &                            *Pi                                               &
         &                            /radiusHaloLagrangian
    powerSpectrumSlope             = +Power_Spectrum_Logarithmic_Derivative(wavenumber)
    concentrationMinimum           =                                                    &
         &                           +self%phi0                                         &
         &                           +self%phi1                                         &
         &                           *powerSpectrumSlope
    peakHeightMinimum              =                                                    &
         &                           +self%eta0                                         &
         &                           +self%eta1                                         &
         &                           *powerSpectrumSlope
    diemerKravtsov2014Concentration=                                                    &
         &                           +0.5d0                                             &
         &                           *concentrationMinimum                              &
         &                           *(                                                 &
         &                             +(peakHeight/peakHeightMinimum)**(-self%alpha)   &
         &                             +(peakHeight/peakHeightMinimum)**(+self%beta )   &
         &                            )
    return
  end function diemerKravtsov2014Concentration

  function diemerKravtsov2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    implicit none
    class(virialDensityContrastClass                      ), pointer       :: diemerKravtsov2014DensityContrastDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: diemerKravtsov2014DensityContrastDefinition)
    select type (diemerKravtsov2014DensityContrastDefinition)
    type is (virialDensityContrastFixed)
       diemerKravtsov2014DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select    
    return
  end function diemerKravtsov2014DensityContrastDefinition

  function diemerKravtsov2014DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: diemerKravtsov2014DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014  ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: diemerKravtsov2014DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition                )
    select type (diemerKravtsov2014DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition                =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          diemerKravtsov2014DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function diemerKravtsov2014DarkMatterProfileDefinition
