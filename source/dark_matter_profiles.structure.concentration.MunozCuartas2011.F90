!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the \cite{munoz-cuartas_redshift_2011} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_MunozCuartas2011
  !% Implements the \cite{munoz-cuartas_redshift_2011} NFW halo concentration algorithm.
  implicit none
  private
  public :: Dark_Matter_Concentrations_MunozCuartas2011_Initialize

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_MunozCuartas2011_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_MunozCuartas2011_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``Munoz-Cuartas2011'' halo concentration module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                    ), intent(in   )          :: darkMatterConcentrationMethod
    procedure(Dark_Matter_Profile_Concentration_MunozCuartas2011), intent(inout), pointer :: Dark_Matter_Profile_Concentration_Get

    if (darkMatterConcentrationMethod == 'Munoz-Cuartas2011') Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_MunozCuartas2011

    return
  end subroutine Dark_Matter_Concentrations_MunozCuartas2011_Initialize

  double precision function Dark_Matter_Profile_Concentration_MunozCuartas2011(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{munoz-cuartas_redshift_2011}.
    use Galacticus_Nodes
    use Cosmology_Functions
    use Cosmological_Parameters
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    double precision                    , parameter              :: alpha                   =-110.001d0, beta               =2469.720d0, &
         &                                                          gamma                   =16.885d0  , m                  =0.097d0   , &
         &                                                          w                       =0.029d0
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                                             :: a                                  , b                             , &
         &                                                          concentrationLogarithmic           , haloMassLogarithmic           , &
         &                                                          redshift

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    ! Compute the concentration.
    redshift                =Redshift_from_Expansion_Factor(Expansion_Factor(thisBasicComponent%time()))
    a                       =w*redshift-m
    b                       =alpha/(redshift+gamma)+beta/(redshift+gamma)**2
    haloMassLogarithmic     =log10(thisBasicComponent%mass()*Little_H_0())
    concentrationLogarithmic=a*haloMassLogarithmic+b
    Dark_Matter_Profile_Concentration_MunozCuartas2011=10.0d0**concentrationLogarithmic
    return
  end function Dark_Matter_Profile_Concentration_MunozCuartas2011

end module Dark_Matter_Profiles_Concentrations_MunozCuartas2011
