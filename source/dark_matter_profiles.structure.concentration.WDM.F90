!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !% An implementation of warm dark matter halo profile concentrations using the
  !% \cite{schneider_non-linear_2012} modifier.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationWDM">
  !#  <description>Dark matter halo concentrations are computed using the modifier of \cite{schneider_non-linear_2012}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationWDM
     !% A dark matter halo profile concentration class implementing the modifier of
     !% \cite{schneider_non-linear_2012}.
     private
     class(darkMatterProfileConcentrationClass), pointer :: cdmConcentration
   contains
     final     ::                                wdmDestructor
     procedure :: concentration               => wdmConcentration
     procedure :: densityContrastDefinition   => wdmDensityContrastDefinition
     procedure :: darkMatterProfileDefinition => wdmDarkMatterProfileDefinition
  end type darkMatterProfileConcentrationWDM

  interface darkMatterProfileConcentrationWDM
     !% Constructors for the {\normalfont \ttfamily WDM} dark matter halo profile concentration
     !% class.
     module procedure wdmConstructorParameters
     module procedure wdmConstructorInternal
  end interface darkMatterProfileConcentrationWDM

contains

  function wdmConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily wdm} dark matter halo profile concentration class.
    use Input_Parameters2
    implicit none
    type(darkMatterProfileConcentrationWDM)                :: wdmConstructorParameters
    type(inputParameters                  ), intent(inout) :: parameters

    !# <objectBuilder class="darkMatterProfileConcentration" name="wdmConstructorParameters%cdmConcentration" source="parameters"/>
    return
  end function wdmConstructorParameters

  function wdmConstructorInternal(cdmMethod)
    !% Generic constructor for the \gls{wdm} dark matter halo concentration class.
    implicit none
    type (darkMatterProfileConcentrationWDM  )                        :: wdmConstructorInternal
    class(darkMatterProfileConcentrationClass), intent(in   ), target :: cdmMethod

    ! Construct the object.
    wdmConstructorInternal%cdmConcentration => cdmMethod
    return
  end function wdmConstructorInternal

  subroutine wdmDestructor(self)
    !% Destructor for the {\normalfont \ttfamily wdm} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationWDM), intent(inout) :: self

    !# <objectDestructor name="self%cdmConcentration"/>
    return
  end subroutine wdmDestructor

  double precision function wdmConcentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    !% using the warm dark matter modifier of \cite{schneider_non-linear_2012}.
    use Transfer_Functions
    implicit none
    class           (darkMatterProfileConcentrationWDM), intent(inout)          :: self
    type            (treeNode                         ), intent(inout), pointer :: node
    class           (nodeComponentBasic               )               , pointer :: basic
    class           (transferFunctionClass            )               , pointer :: transferFunction_
    ! Parameters of Schneider et al. (2012)'s fitting formula.
    double precision                                   , parameter              :: gamma1           =15.0d0, gamma2=0.3d0

    ! Get required objects.
    basic             => node            %basic()
    transferFunction_ => transferFunction      ()
    ! Get WDM concentration
    wdmConcentration=                                         &
         &  self%cdmConcentration%concentration(node)         &
         & /(                                                 &
         &    1.0d0                                           &
         &   +gamma1                                          &
         &   *(transferFunction_%halfModeMass()/basic%mass()) &
         &  )**gamma2
    return
  end function wdmConcentration

  function wdmDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of
    !% concentration in the warm dark matter modifier of \cite{schneider_non-linear_2012}.
    implicit none
    class(virialDensityContrastClass       ), pointer       :: wdmDensityContrastDefinition
    class(darkMatterProfileConcentrationWDM), intent(inout) :: self

    wdmDensityContrastDefinition => self%cdmConcentration%densityContrastDefinition()
    return
  end function wdmDensityContrastDefinition

  function wdmDarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% warm dark matter modifier of \cite{schneider_non-linear_2012}.
    implicit none
    class(darkMatterProfileClass           ), pointer       :: wdmDarkMatterProfileDefinition
    class(darkMatterProfileConcentrationWDM), intent(inout) :: self

    wdmDarkMatterProfileDefinition => self%cdmConcentration%darkMatterProfileDefinition()
    return
  end function wdmDarkMatterProfileDefinition
