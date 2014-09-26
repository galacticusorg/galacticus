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

  !% An implementation of \cite{bryan_statistical_1998} dark matter halo virial density contrasts.

  !# <virialDensityContrast name="virialDensityContrastBryanNorman1998">
  !#  <description>\cite{bryan_statistical_1998} dark matter halo virial density contrasts.</description>
  !# </virialDensityContrast>

  type, extends(virialDensityContrastClass) :: virialDensityContrastBryanNorman1998
     !% A dark matter halo virial density contrast class using the fitting functions of \cite{bryan_statistical_1998}.
     private
     integer   :: fitType
   contains
     procedure :: densityContrast             => bryanNorman1998DensityContrast
     procedure :: densityContrastRateOfChange => bryanNorman1998DensityContrastRateOfChange
  end type virialDensityContrastBryanNorman1998

  interface virialDensityContrastBryanNorman1998
     !% Constructors for the {\tt bryanNorman1998} dark matter halo virial density contrast class.
     module procedure bryanNorman1998DefaultConstructor
  end interface virialDensityContrastBryanNorman1998

  ! Labels for different fitting function types.
  integer, parameter :: bryanNorman1998FitTypeFlatUniverse=1, bryanNorman1998FitTypeZeroLambda=0

contains

  function bryanNorman1998DefaultConstructor()
    !% Default constructor for the {\tt bryanNorman1998} dark matter halo virial density contrast class.
    use Cosmology_Parameters
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    type (virialDensityContrastBryanNorman1998), target  :: bryanNorman1998DefaultConstructor
    class(cosmologyParametersClass            ), pointer :: cosmologyParameters_

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    ! Check that fitting formulae are applicable to this cosmology.
    if (cosmologyParameters_%OmegaDarkEnergy() == 0.0d0) then
       bryanNorman1998DefaultConstructor%fitType=bryanNorman1998FitTypeZeroLambda
    else if (.not.Values_Differ(cosmologyParameters_%OmegaMatter()+cosmologyParameters_%OmegaDarkEnergy(),1.0d0,absTol=1.0d-6)) then
       bryanNorman1998DefaultConstructor%fitType=bryanNorman1998FitTypeFlatUniverse
    else
       call Galacticus_Error_Report('bryanNorman1998DefaultConstructor','no fitting formula available for this cosmology')
    end if
    return
  end function bryanNorman1998DefaultConstructor

  double precision function bryanNorman1998DensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming the fitting function of \cite{bryan_statistical_1998}.
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastBryanNorman1998), intent(inout)           :: self
    double precision                                      , intent(in   )           :: mass
    double precision                                      , intent(in   ), optional :: time               , expansionFactor
    logical                                               , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass             ), pointer                 :: cosmologyFunctions_
    double precision                                                                :: x

    cosmologyFunctions_ => cosmologyFunctions()
    x=cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0
    select case (self%fitType)
    case (bryanNorman1998FitTypeZeroLambda)
       bryanNorman1998DensityContrast=(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)/cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)
    case (bryanNorman1998FitTypeFlatUniverse)
       bryanNorman1998DensityContrast=(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)
    end select
    return
  end function bryanNorman1998DensityContrast

  double precision function bryanNorman1998DensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming the fitting function of \cite{bryan_statistical_1998}.
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastBryanNorman1998), intent(inout)           :: self
    double precision                                      , intent(in   )           :: mass
    double precision                                      , intent(in   ), optional :: time      , expansionFactor
    logical                                               , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass             ), pointer                 :: cosmologyFunctions_
    double precision                                                                :: x

    cosmologyFunctions_ => cosmologyFunctions()
    x=cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0
    select case (self%fitType)
    case (bryanNorman1998FitTypeZeroLambda)
       bryanNorman1998DensityContrastRateOfChange=                                           &
            & (                                                                              &
            &  +(            +60.0d0  -64.0d0*x   )                                          &
            &  -(18.0d0*Pi**2+60.0d0*x-32.0d0*x**2)                                          &
            &  /cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing) &
            & )                                                                              &
            & * cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing) &
            & / cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing)
    case (bryanNorman1998FitTypeFlatUniverse)
     bryanNorman1998DensityContrastRateOfChange=                                             &
            & (                                                                              &
            &  +(            +82.0d0  -78.0d0*x   )                                          &
            &  -(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)                                          &
            &  /cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing) &
            & )                                                                              &
            & * cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing) &
            & / cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing)
    end select
    return
  end function bryanNorman1998DensityContrastRateOfChange
