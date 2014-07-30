!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo profile concentrations using the \cite{munoz-cuartas_redshift_2011} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationMunozCuartas2011">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{munoz-cuartas_redshift_2011}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationMunozCuartas2011
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{munoz-cuartas_redshift_2011}.
     private
   contains
     procedure :: concentration               => munozCuartas2011Concentration
     procedure :: densityContrastDefinition   => munozCuartas2011DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => munozCuartas2011DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationMunozCuartas2011

  interface darkMatterProfileConcentrationMunozCuartas2011
     !% Constructors for the {\tt munozCuartas2011} dark matter halo profile concentration class.
     module procedure munozCuartas2011DefaultConstructor
  end interface darkMatterProfileConcentrationMunozCuartas2011

contains

  function munozCuartas2011DefaultConstructor()
    !% Default constructor for the {\tt munozCuartas2011} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationMunozCuartas2011), target  :: munozCuartas2011DefaultConstructor

    return
  end function munozCuartas2011DefaultConstructor

  double precision function munozCuartas2011Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{munoz-cuartas_redshift_2011} algorithm.
    use Cosmology_Functions
    use Cosmology_Parameters
    implicit none
    class           (darkMatterProfileConcentrationMunozCuartas2011), intent(inout)          :: self
    type            (treeNode                                      ), intent(inout), pointer :: node
    class           (nodeComponentBasic                            )               , pointer :: basic
    class           (cosmologyFunctionsClass                       )               , pointer :: cosmologyFunctions_
    class           (cosmologyParametersClass                      )               , pointer :: cosmologyParameters_
    double precision                                                , parameter              :: alpha               =-110.001d0, beta               =2469.720d0, &
         &                                                                                      gamma               =16.885d0  , m                  =0.097d0   , &
         &                                                                                      w                   =0.029d0
    double precision                                                                         :: a                                           , b                             , &
         &                                                                                      concentrationLogarithmic                    , haloMassLogarithmic           , &
         &                                                                                      redshift

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    ! Get the default cosmology functions object.
    cosmologyFunctions_  => cosmologyFunctions ()
    ! Get the basic component.
    basic                => node%basic         ()
    ! Compute the concentration.
    redshift                     =cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(basic%time()))
    a                            =w*redshift-m
    b                            =alpha/(redshift+gamma)+beta/(redshift+gamma)**2
    haloMassLogarithmic          =log10(basic%mass()*cosmologyParameters_%HubbleConstant(unitsLittleH))
    concentrationLogarithmic     =a*haloMassLogarithmic+b
    munozCuartas2011Concentration=10.0d0**concentrationLogarithmic
    return
  end function munozCuartas2011Concentration

  function munozCuartas2011DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{munoz-cuartas_redshift_2011} algorithm.
    implicit none
    class(virialDensityContrastClass                    ), pointer       :: munozCuartas2011DensityContrastDefinition
    class(darkMatterProfileConcentrationMunozCuartas2011), intent(inout) :: self
    
    allocate(virialDensityContrastBryanNorman1998 :: munozCuartas2011DensityContrastDefinition)
    select type (munozCuartas2011DensityContrastDefinition)
    type is (virialDensityContrastBryanNorman1998)
      munozCuartas2011DensityContrastDefinition=virialDensityContrastBryanNorman1998()
    end select
    return
  end function munozCuartas2011DensityContrastDefinition
  
  function munozCuartas2011DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{munoz-cuartas_redshift_2011} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: munozCuartas2011DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationMunozCuartas2011    ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: munozCuartas2011DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition              )
    select type (munozCuartas2011DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition              =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          munozCuartas2011DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function munozCuartas2011DarkMatterProfileDefinition

