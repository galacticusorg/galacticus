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

!% Contains a module which implements the \cite{gao_redshift_2008} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_Gao2008
  !% Implements the \cite{gao_redshift_2008} NFW halo concentration algorithm.
  implicit none
  private
  public :: Dark_Matter_Concentrations_Gao2008_Initialize

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_Gao2008_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_Gao2008_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``Gao2008'' halo concentration module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                           ), intent(in   )          :: darkMatterConcentrationMethod
    procedure(Dark_Matter_Profile_Concentration_Gao2008), intent(inout), pointer :: Dark_Matter_Profile_Concentration_Get

    if (darkMatterConcentrationMethod == 'Gao2008') Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_Gao2008
    return
  end subroutine Dark_Matter_Concentrations_Gao2008_Initialize

  double precision function Dark_Matter_Profile_Concentration_Gao2008(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{gao_redshift_2008}. More
    !% specifically, we fit the redshift dependence of the parameters $A$ and $B$ in the fitting formula of
    !% \cite{gao_redshift_2008} and use these fits to find $A$ and $B$ at any given redshift, from which we then compute the
    !% concentration. Note that the fits of \cite{gao_redshift_2008} were computed using Einasto profile fits and utilizing
    !% $M_{200}$ and $c_{200}$.
    use Galacticus_Nodes
    use Cosmology_Functions
    implicit none
    type            (treeNode               ), intent(inout), pointer :: thisNode
    double precision                         , parameter              :: littleHubbleConstantGao2008=0.73d0
    class           (nodeComponentBasic     )               , pointer :: thisBasicComponent
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault
    double precision                                                  :: logarithmExpansionFactor          , logarithmHaloMass, &
         &                                                               parameterA                        , parameterB

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    ! Compute the concentration.
    logarithmHaloMass                        =log10(littleHubbleConstantGao2008*thisBasicComponent%mass() )
    logarithmExpansionFactor                 =log10(cosmologyFunctionsDefault%expansionFactor(           thisBasicComponent%time()))
    parameterA                               =-0.140d0*exp(-((logarithmExpansionFactor+0.05d0)/0.35d0)**2)
    parameterB                               = 2.646d0*exp(-((logarithmExpansionFactor+0.00d0)/0.50d0)**2)
    Dark_Matter_Profile_Concentration_Gao2008=10.0d0**(parameterA*logarithmHaloMass+parameterB)
    return
  end function Dark_Matter_Profile_Concentration_Gao2008

end module Dark_Matter_Profiles_Concentrations_Gao2008
