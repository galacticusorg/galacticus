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

  !% An implementation of dark matter halo profile concentrations using the \cite{gao_redshift_2008} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationGao2008">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{gao_redshift_2008}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationGao2008
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{gao_redshift_2008}.
     private
   contains
     procedure :: concentration => gao2008Concentration
  end type darkMatterProfileConcentrationGao2008

  interface darkMatterProfileConcentrationGao2008
     !% Constructors for the {\tt gao2008} dark matter halo profile concentration class.
     module procedure gao2008DefaultConstructor
  end interface darkMatterProfileConcentrationGao2008

contains

  function gao2008DefaultConstructor()
    !% Default constructor for the {\tt gao2008} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationGao2008), target  :: gao2008DefaultConstructor

    return
  end function gao2008DefaultConstructor

  double precision function gao2008Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{gao_redshift_2008} algorithm.
    use Cosmology_Functions
    implicit none
    class           (darkMatterProfileConcentrationGao2008), intent(inout)          :: self
    type            (treeNode                             ), intent(inout), pointer :: node
    class           (nodeComponentBasic                   )               , pointer :: basic
    class           (cosmologyFunctionsClass              )               , pointer :: cosmologyFunctions_
    double precision                                       , parameter              :: littleHubbleConstantGao2008=0.73d0
    double precision                                                                :: logarithmExpansionFactor, logarithmHaloMass, &
         &                                                                             parameterA              , parameterB

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Get the basic component.
    basic               => node%basic()
    ! Compute the concentration.
    logarithmHaloMass       =log10(littleHubbleConstantGao2008        *basic%mass())
    logarithmExpansionFactor=log10(cosmologyFunctions_%expansionFactor(basic%time()))
    parameterA              =-0.140d0*exp(-((logarithmExpansionFactor+0.05d0)/0.35d0)**2)
    parameterB              = 2.646d0*exp(-((logarithmExpansionFactor+0.00d0)/0.50d0)**2)
    gao2008Concentration    =10.0d0**(parameterA*logarithmHaloMass+parameterB)
    return
  end function gao2008Concentration
