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

  !% An implementation of dark matter halo profile concentrations using the \cite{zhao_accurate_2009} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationZhao2009">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{zhao_accurate_2009}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationZhao2009
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{zhao_accurate_2009}.
     private
   contains
     procedure :: concentration => zhao2009Concentration
  end type darkMatterProfileConcentrationZhao2009

  interface darkMatterProfileConcentrationZhao2009
     !% Constructors for the {\tt zhao2009} dark matter halo profile concentration class.
     module procedure zhao2009DefaultConstructor
  end interface darkMatterProfileConcentrationZhao2009

contains

  function zhao2009DefaultConstructor()
    !% Default constructor for the {\tt zhao2009} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationZhao2009), target  :: zhao2009DefaultConstructor
    return
  end function zhao2009DefaultConstructor

  double precision function zhao2009Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{zhao_accurate_2009} algorithm.
    use Dark_Matter_Halo_Formation_Times
    implicit none
    class           (darkMatterProfileConcentrationZhao2009), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    class           (nodeComponentBasic                    )               , pointer :: basic
    double precision                                        , parameter              :: concentrationMinimum =4.00d0
    double precision                                        , parameter              :: formationMassFraction=0.04d0
    double precision                                                                 :: timeFormation               , timeNode

    ! Get the basic component.
    basic => node%basic()
    ! Compute the concentration.
    timeNode     =basic%time()
    timeFormation=Dark_Matter_Halo_Formation_Time(node,formationMassFraction)
    ! Compute the concentration from the formation time using the Zhao et al. (2009) fitting formula.
    zhao2009Concentration=concentrationMinimum*(1.0d0+(timeNode/3.75d0/timeFormation)**8.4d0)**0.125d0
   return
  end function zhao2009Concentration
