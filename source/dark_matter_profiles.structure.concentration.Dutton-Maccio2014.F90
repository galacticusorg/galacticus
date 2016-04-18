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

  !% An implementation of dark matter halo profile concentrations using the \cite{dutton_cold_2014} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationDuttonMaccio2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{dutton_cold_2014}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDuttonMaccio2014
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{dutton_cold_2014}.
     private
     double precision :: a1                   , a2                  , &
          &              a3                   , a4                  , &
          &              b1                   , b2
     integer          :: densityContrastMethod, densityProfileMethod
   contains
     final     ::                                duttonMaccio2014Destructor
     procedure :: concentration               => duttonMaccio2014Concentration
     procedure :: densityContrastDefinition   => duttonMaccio2014DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => duttonMaccio2014DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDuttonMaccio2014

  interface darkMatterProfileConcentrationDuttonMaccio2014
     !% Constructors for the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.
     module procedure duttonMaccio2014ConstructorParameters
     module procedure duttonMaccio2014ConstructorInternalType
     module procedure duttonMaccio2014ConstructorInternalDefined
  end interface darkMatterProfileConcentrationDuttonMaccio2014
  
  ! Density contrast methods.
  !# <enumeration>
  !#  <name>duttonMaccio2014DensityContrastMethod</name>
  !#  <description>Enumeration of density contrast methods available in the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.</description>
  !#  <visibility>private</visibility>
  !#  <entry label="virial" />
  !#  <entry label="mean200"/>
  !# </enumeration>

  ! Density profile methods.
  !# <enumeration>
  !#  <name>duttonMaccio2014DensityProfileMethod</name>
  !#  <description>Enumeration of density profile methods available in the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.</description>
  !#  <visibility>private</visibility>
  !#  <entry label="NFW"    />
  !#  <entry label="einasto"/>
  !# </enumeration>
  
contains

  function duttonMaccio2014ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationDuttonMaccio2014)                :: duttonMaccio2014ConstructorParameters
    type            (inputParameters                               ), intent(inout) :: parameters
    type            (varying_string                                )                :: fitType
    double precision                                                                :: a1                                   , a2, &
         &                                                                             a3                                   , a4, &
         &                                                                             b1                                   , b2
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>fitType</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('nfwVirial')</defaultValue>
    !#   <description>The type of halo definition for which the concentration-mass relation should be computed. Allowed values are {\normalfont \ttfamily nfwVirial}, {\normalfont \ttfamily nfwMean200}, {\normalfont \ttfamily einastoMean200}, and {\normalfont \ttfamily userDefined}.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    if (fitType == "userDefined") then
       !# <inputParameter>
       !#   <name>a1</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $a_1$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>a2</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $a_2$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>a3</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $a_3$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>a4</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $a_4$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>b1</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $b_1$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>b2</name>
       !#   <source>parameters</source>
       !#   <description>Parameter $b_2$ in the \cite{dutton_cold_2014} halo concentration--mass relation.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       duttonMaccio2014ConstructorParameters=duttonMaccio2014ConstructorInternalDefined(a1,a2,a3,a4,b1,b2)
    else
       duttonMaccio2014ConstructorParameters=duttonMaccio2014ConstructorInternalType   (char(fitType))
    end if
    return
  end function duttonMaccio2014ConstructorParameters

  function duttonMaccio2014ConstructorInternalType(fitType)
    !% Constructor for the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type     (darkMatterProfileConcentrationDuttonMaccio2014)                :: duttonMaccio2014ConstructorInternalType
    character(len=*                                         ), intent(in   ) :: fitType

    select case (fitType)
    case ('nfwVirial'     )
       duttonMaccio2014ConstructorInternalType%a1                   =+0.537d0
       duttonMaccio2014ConstructorInternalType%a2                   =+1.025d0
       duttonMaccio2014ConstructorInternalType%a3                   =-0.718d0
       duttonMaccio2014ConstructorInternalType%a4                   =+1.080d0
       duttonMaccio2014ConstructorInternalType%b1                   =-0.097d0
       duttonMaccio2014ConstructorInternalType%b2                   =+0.024d0
       duttonMaccio2014ConstructorInternalType%densityContrastMethod=duttonMaccio2014DensityContrastMethodVirial
       duttonMaccio2014ConstructorInternalType%densityProfileMethod =duttonMaccio2014DensityProfileMethodNFW
    case ('nfwMean200'    )
       duttonMaccio2014ConstructorInternalType%a1                   =+0.520d0
       duttonMaccio2014ConstructorInternalType%a2                   =+0.905d0
       duttonMaccio2014ConstructorInternalType%a3                   =-0.617d0
       duttonMaccio2014ConstructorInternalType%a4                   =+1.210d0
       duttonMaccio2014ConstructorInternalType%b1                   =-0.101d0
       duttonMaccio2014ConstructorInternalType%b2                   =+0.026d0
       duttonMaccio2014ConstructorInternalType%densityContrastMethod=duttonMaccio2014DensityContrastMethodMean200
       duttonMaccio2014ConstructorInternalType%densityProfileMethod =duttonMaccio2014DensityProfileMethodNFW
    case ('einastoMean200')
       duttonMaccio2014ConstructorInternalType%a1                   =+0.459d0
       duttonMaccio2014ConstructorInternalType%a2                   =+0.977d0
       duttonMaccio2014ConstructorInternalType%a3                   =-0.490d0
       duttonMaccio2014ConstructorInternalType%a4                   =+1.303d0
       duttonMaccio2014ConstructorInternalType%b1                   =-0.130d0
       duttonMaccio2014ConstructorInternalType%b2                   =+0.029d0
       duttonMaccio2014ConstructorInternalType%densityContrastMethod=duttonMaccio2014DensityContrastMethodMean200
       duttonMaccio2014ConstructorInternalType%densityProfileMethod =duttonMaccio2014DensityProfileMethodEinasto
    case default
       call Galacticus_Error_Report('duttonMaccio2014ConstructorInternalType','unrecognized fit type [available types are: nfwVirial, nfwMean200, einastoMean200]')
    end select
    return
  end function duttonMaccio2014ConstructorInternalType
  
  function duttonMaccio2014ConstructorInternalDefined(a1,a2,a3,a4,b1,b2)
    !% Constructor for the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class with user defined
    !% parameters.
    use Galacticus_Error
    implicit none
    type           (darkMatterProfileConcentrationDuttonMaccio2014)                :: duttonMaccio2014ConstructorInternalDefined
    double precision                                               , intent(in   ) :: a1, a2, a3, a4, b1, b2

    duttonMaccio2014ConstructorInternalDefined%a1=a1
    duttonMaccio2014ConstructorInternalDefined%a2=a2
    duttonMaccio2014ConstructorInternalDefined%a3=a3
    duttonMaccio2014ConstructorInternalDefined%a4=a4
    duttonMaccio2014ConstructorInternalDefined%b1=b1
    duttonMaccio2014ConstructorInternalDefined%b2=b2
    return
  end function duttonMaccio2014ConstructorInternalDefined

  subroutine duttonMaccio2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily duttonMaccio2014} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationDuttonMaccio2014), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine duttonMaccio2014Destructor

  double precision function duttonMaccio2014Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the \cite{dutton_cold_2014}
    !% algorithm.
    use Cosmology_Functions
    implicit none
    class           (darkMatterProfileConcentrationDuttonMaccio2014), intent(inout)          :: self
    type            (treeNode                                      ), intent(inout), pointer :: node
    class           (nodeComponentBasic                            )               , pointer :: basic
    class           (cosmologyFunctionsClass                       )               , pointer :: cosmologyFunctions_
    double precision                                                , parameter              :: littleHubbleConstantDuttonMaccio2014= 0.671d0
    double precision                                                , parameter              :: massNormalization                   =12.000d0
    double precision                                                                         :: redshift                                     , logarithmHaloMass, &
         &                                                                                      parameterA                                   , parameterB

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Get the basic component.
    basic               => node%basic()
    ! Compute the concentration.
    logarithmHaloMass            =+log10(                                                       &
         &                                littleHubbleConstantDuttonMaccio2014                  &
         &                               *basic%mass()                                          &
         &                              )                                                       &
         &                        -massNormalization
    redshift                     =max(                                                                &
         &                            0.0d0                                                         , &
         &                            cosmologyFunctions_ %redshiftFromExpansionFactor(               &
         &                             cosmologyFunctions_%expansionFactor             (              &
         &                                                                              basic%time()  &
         &                                                                             )              &
         &                                                                            )               &
         &                           )
    parameterA                   =self%a1+(self%a2-self%a1)*exp(self%a3*redshift**self%a4)
    parameterB                   =self%b1+self%b2*redshift
    duttonMaccio2014Concentration=10.0d0**(parameterA+parameterB*logarithmHaloMass)
    return
  end function duttonMaccio2014Concentration

  function duttonMaccio2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{dutton_cold_2014} algorithm.
    implicit none
    class(virialDensityContrastClass                    ), pointer       :: duttonMaccio2014DensityContrastDefinition
    class(darkMatterProfileConcentrationDuttonMaccio2014), intent(inout) :: self
    
    select case (self%densityContrastMethod)
    case (duttonMaccio2014DensityContrastMethodMean200   )
       allocate(virialDensityContrastFixed                         :: duttonMaccio2014DensityContrastDefinition)
       select type (duttonMaccio2014DensityContrastDefinition)
       type is (virialDensityContrastFixed)
          duttonMaccio2014DensityContrastDefinition=virialDensityContrastFixed                        (200.0d0,virialDensityContrastFixedDensityTypeMean)
       end select
    case (duttonMaccio2014DensityContrastMethodVirial)
       allocate(virialDensityContrastSphericalCollapseMatterLambda :: duttonMaccio2014DensityContrastDefinition)
       select type (duttonMaccio2014DensityContrastDefinition)
       type is (virialDensityContrastSphericalCollapseMatterLambda)
          duttonMaccio2014DensityContrastDefinition=virialDensityContrastSphericalCollapseMatterLambda(                                                 )
       end select
    end select
    return
  end function duttonMaccio2014DensityContrastDefinition

  function duttonMaccio2014DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: duttonMaccio2014DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDuttonMaccio2014    ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition)
    select type (darkMatterHaloScaleDefinition)
    type is (darkMatterHaloScaleVirialDensityContrastDefinition)
       select case (self%densityProfileMethod)
       case (duttonMaccio2014DensityProfileMethodNFW    )
          allocate(darkMatterProfileNFW     :: duttonMaccio2014DarkMatterProfileDefinition)
          select type (duttonMaccio2014DarkMatterProfileDefinition)
          type is (darkMatterProfileNFW    )
             darkMatterHaloScaleDefinition              =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             duttonMaccio2014DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
          end select
       case (duttonMaccio2014DensityProfileMethodEinasto)
          allocate(darkMatterProfileEinasto :: duttonMaccio2014DarkMatterProfileDefinition)
          select type (duttonMaccio2014DarkMatterProfileDefinition)
          type is (darkMatterProfileEinasto)
             darkMatterHaloScaleDefinition              =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             duttonMaccio2014DarkMatterProfileDefinition=darkMatterProfileEinasto                          (darkMatterHaloScaleDefinition   )
          end select
       end select
    end select
    return
  end function duttonMaccio2014DarkMatterProfileDefinition
