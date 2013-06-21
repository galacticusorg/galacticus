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

!% Contains a module which implements the \cite{zhao_accurate_2009} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_Zhao2009
  !% Implements the \cite{zhao_accurate_2009} NFW halo concentration algorithm.
  use Galacticus_Nodes
  implicit none
  private
  public :: Dark_Matter_Concentrations_Zhao2009_Initialize

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_Zhao2009_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_Zhao2009_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``Zhao2009'' halo concentration module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                            ), intent(in   )          :: darkMatterConcentrationMethod
    procedure(Dark_Matter_Profile_Concentration_Zhao2009), intent(inout), pointer :: Dark_Matter_Profile_Concentration_Get

    if (darkMatterConcentrationMethod == 'Zhao2009')                                          &
         & Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_Zhao2009

    return
  end subroutine Dark_Matter_Concentrations_Zhao2009_Initialize

  double precision function Dark_Matter_Profile_Concentration_Zhao2009(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{zhao_accurate_2009}.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Formation_Times
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    double precision                    , parameter              :: concentrationMinimum =4.00d0
    double precision                    , parameter              :: formationMassFraction=0.04d0
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                                             :: timeFormation               , timeNode

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Compute the concentration.
    timeNode     =thisBasicComponent%time()
    timeFormation=Dark_Matter_Halo_Formation_Time(thisNode,formationMassFraction)

    ! Compute the concentration from the formation time using the Zhao et al. (2009) fitting formula.
    Dark_Matter_Profile_Concentration_Zhao2009=concentrationMinimum*(1.0d0+(timeNode/3.75d0/timeFormation)**8.4d0)**0.125d0

    return
  end function Dark_Matter_Profile_Concentration_Zhao2009

end module Dark_Matter_Profiles_Concentrations_Zhao2009
