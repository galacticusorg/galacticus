!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements the dust extinction curve of :cite:t:`cardelli_relationship_1989`.
  !!}

  !![
  <dustExtinctionCurve name="dustExtinctionCurveCardelli1989" docformat="rst">
   <description>
   The Milky Way dust extinction curve of :cite:t:`cardelli_relationship_1989`. Their equation (1) gives
   :math:`A(\lambda)/A(\mathrm{V}) = a(x) + b(x)/R_\mathrm{V}` with :math:`x=1/\lambda` in inverse microns, and
   :math:`a` and :math:`b` the fitting functions of their equations (2)--(5). The parameter ``Rv`` selects the
   total-to-selective extinction ratio; :math:`R_\mathrm{V}=3.1` is the canonical diffuse interstellar medium value.

   The fit spans :math:`0.3\,\mu\mathrm{m}^{-1} \le x &lt; 8\,\mu\mathrm{m}^{-1}`, i.e. :math:`0.125\,\mu\mathrm{m} &lt;
   \lambda \le 3.33\,\mu\mathrm{m}$. Outside that range this implementation returns zero, i.e. *no* attenuation ---
   the behavior of the original ``stellarSpectraDustAttenuation`` implementation, retained here for consistency. Use
   ``wavelengthRange`` to detect the boundary.
   </description>
  </dustExtinctionCurve>
  !!]
  type, extends(dustExtinctionCurveClass) :: dustExtinctionCurveCardelli1989
     !!{RST
     The dust extinction curve of :cite:t:`cardelli_relationship_1989`.
     !!}
     private
     double precision :: Rv
   contains
     !![
     <methods docformat="rst">
       <method method="a" description="Return the fitting function :math:`a(x)` of :cite:t:`cardelli_relationship_1989`."/>
       <method method="b" description="Return the fitting function :math:`b(x)` of :cite:t:`cardelli_relationship_1989`."/>
     </methods>
     !!]
     procedure :: attenuationRelative => cardelli1989AttenuationRelative
     procedure :: wavelengthRange     => cardelli1989WavelengthRange
     procedure :: a                   => cardelli1989A
     procedure :: b                   => cardelli1989B
  end type dustExtinctionCurveCardelli1989

  interface dustExtinctionCurveCardelli1989
     !!{RST
     Constructors for the :galacticus-class:`dustExtinctionCurveCardelli1989` dust extinction curve class.
     !!}
     module procedure cardelli1989ConstructorParameters
     module procedure cardelli1989ConstructorInternal
  end interface dustExtinctionCurveCardelli1989

  ! Range of validity of the fit, in inverse microns (Cardelli et al. 1989).
  double precision, parameter :: cardelli1989XMinimum=0.3d0
  double precision, parameter :: cardelli1989XMaximum=8.0d0

contains

  function cardelli1989ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustExtinctionCurveCardelli1989` dust extinction curve class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustExtinctionCurveCardelli1989)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: Rv

    !![
    <inputParameter docformat="rst">
      <name>Rv</name>
      <defaultValue>3.1d0</defaultValue>
      <description>
      The total-to-selective extinction ratio, :math:`R_\mathrm{V}=A(\mathrm{V})/E(B-V)`, in the
      :cite:t:`cardelli_relationship_1989` extinction curve.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustExtinctionCurveCardelli1989(Rv)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cardelli1989ConstructorParameters

  function cardelli1989ConstructorInternal(Rv) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustExtinctionCurveCardelli1989` dust extinction curve class.
    !!}
    implicit none
    type            (dustExtinctionCurveCardelli1989)                :: self
    double precision                                 , intent(in   ) :: Rv
    !![
    <constructorAssign variables="Rv"/>
    !!]

    return
  end function cardelli1989ConstructorInternal

  double precision function cardelli1989AttenuationRelative(self,wavelength) result(attenuationRelative)
    !!{RST
    Return the relative dust opacity for the extinction curve of :cite:t:`cardelli_relationship_1989`.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveCardelli1989), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavelength
    double precision                                                 :: x

    ! Eqn. (1) of Cardelli et al. (1989).
    x                  =+1.0d0                &
         &              /(                    &
         &                +wavelength         &
         &                /micronsToAngstroms &
         &               )
    attenuationRelative=+self%a(x)            &
         &              +self%b(x)            &
         &              /self%Rv
    return
  end function cardelli1989AttenuationRelative

  subroutine cardelli1989WavelengthRange(self,wavelengthMinimum,wavelengthMaximum)
    !!{RST
    Return the range of wavelengths over which the :cite:t:`cardelli_relationship_1989` extinction curve is defined.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveCardelli1989), intent(inout) :: self
    double precision                                 , intent(  out) :: wavelengthMinimum, wavelengthMaximum
    !$GLC attributes unused :: self

    wavelengthMinimum=micronsToAngstroms/cardelli1989XMaximum
    wavelengthMaximum=micronsToAngstroms/cardelli1989XMinimum
    return
  end subroutine cardelli1989WavelengthRange

  double precision function cardelli1989A(self,x) result(a)
    !!{RST
    Return the fitting function :math:`a(x)` for the extinction curve of :cite:t:`cardelli_relationship_1989`.
    !!}
    implicit none
    class           (dustExtinctionCurveCardelli1989), intent(inout) :: self
    double precision                                 , intent(in   ) :: x
    double precision                                                 :: y   , Fa
    !$GLC attributes unused :: self

    a=0.0d0
    if      (x >= cardelli1989XMinimum .and. x < 1.1d0               ) then
       a=+0.57400d0*x**1.61d0
    else if (x >= 1.1d0                .and. x < 3.3d0               ) then
       y      =x-1.82d0
       a      =+1.00000d0      &
            &  +0.17699d0*y    &
            &  -0.50447d0*y**2 &
            &  -0.02427d0*y**3 &
            &  +0.72085d0*y**4 &
            &  +0.01979d0*y**5 &
            &  -0.77530d0*y**6 &
            &  +0.32999d0*y**7
    else if (x >= 3.3d0                .and. x < cardelli1989XMaximum) then
       if (x < 5.9d0) then
          Fa     =+0.0d0
       else
          Fa     =-0.044730d0*(x-5.9d0)**2 &
               &  -0.009779d0*(x-5.9d0)**3
       end if
       a         =+1.752d0                         &
            &     -0.316d0*  x                     &
            &     -0.104d0/((x-4.67d0)**2+0.341d0) &
            &     +Fa
    end if
    return
  end function cardelli1989A

  double precision function cardelli1989B(self,x) result(b)
    !!{RST
    Return the fitting function :math:`b(x)` for the extinction curve of :cite:t:`cardelli_relationship_1989`.
    !!}
    implicit none
    class           (dustExtinctionCurveCardelli1989), intent(inout) :: self
    double precision                                 , intent(in   ) :: x
    double precision                                                 :: y   , Fb
    !$GLC attributes unused :: self

    b=0.0d0
    if      (x >= cardelli1989XMinimum .and. x < 1.1d0               ) then
       b      =-0.52700d0*x**1.61d0
    else if (x >= 1.1d0                .and. x < 3.3d0               ) then
       y      =x-1.82d0
       b      =+1.41338d0*y    &
            &  +2.28305d0*y**2 &
            &  +1.07233d0*y**3 &
            &  -5.38434d0*y**4 &
            &  -0.62251d0*y**5 &
            &  +5.30260d0*y**6 &
            &  -2.09002d0*y**7
    else if (x >= 3.3d0                .and. x < cardelli1989XMaximum) then
       if (x < 5.9d0) then
          Fb     =+0.0d0
       else
          Fb     =+0.2130d0*(x-5.9d0)**2 &
               &  +0.1207d0*(x-5.9d0)**3
       end if
       b         =-3.090d0                         &
            &     +1.825d0*  x                     &
            &     +1.206d0/((x-4.62d0)**2+0.263d0) &
            &     +Fb
    end if
    return
  end function cardelli1989B
