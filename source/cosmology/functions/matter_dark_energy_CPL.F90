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
  An implementation of the cosmological functions class for cosmologies consisting of collisionless matter and dark energy with an equation of state of the form: :math:`P=\rho^w` with :math:`w(a)=w_0+w_\mathrm{a}(1-a)`.
  !!}

  !![
  <cosmologyFunctions name="cosmologyFunctionsMatterDarkEnergyCPL" docformat="rst">
   <description>
   Cosmological relations are computed assuming a universe that contains only collisionless matter and dark energy with an
   equation of state :math:`w(a)=w_0+w_\mathrm{a}(1-a)` :cite:p:`chevallier_accelerating_2001,linder_exploring_2003`, with
   :math:`w_0=`\ ``[darkEnergyEquationOfStateW0]``, and :math:`w_\mathrm{a}=`\ ``[darkEnergyEquationOfStateWA]``. This
   parameterization is commonly referred to as the "CPL" (Chevallier-Polarski-Linder) form. All functionality is inherited from
   the :galacticus-class:`cosmologyFunctionsMatterDarkEnergy` class---only the dark energy equation of state, and the
   corresponding dark energy density exponent, are changed.
   </description>
  </cosmologyFunctions>
  !!]
  type, extends(cosmologyFunctionsMatterDarkEnergy) :: cosmologyFunctionsMatterDarkEnergyCPL
     !!{RST
     A cosmological functions class for cosmologies consisting of matter plus dark energy with equation of state :math:`w(a)=w_0+w_\mathrm{a}(1-a)`.
     !!}
     private
     double precision :: darkEnergyEquationOfStateWA
   contains
     procedure :: equationOfStateDarkEnergy    => matterDarkEnergyCPLEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy           => matterDarkEnergyCPLExponentDarkEnergy
     procedure :: exponentDarkEnergyDerivative => matterDarkEnergyCPLExponentDarkEnergyDerivative
  end type cosmologyFunctionsMatterDarkEnergyCPL

  interface cosmologyFunctionsMatterDarkEnergyCPL
     !!{RST
     Constructors for the :galacticus-class:`cosmologyFunctionsMatterDarkEnergyCPL` cosmological functions class.
     !!}
     module procedure matterDarkEnergyCPLConstructorParameters
     module procedure matterDarkEnergyCPLConstructorInternal
  end interface cosmologyFunctionsMatterDarkEnergyCPL

  ! Maximum value of |ln(a)| for which series expansions are used in evaluating the dark energy exponent and its derivative. Both
  ! expressions contain removable singularities at a=1, and so must be evaluated via a series expansion sufficiently close to that
  ! point to avoid catastrophic cancellation.
  double precision, parameter :: matterDarkEnergyCPLExpansionFactorLogarithmicSmall=1.0d-3

contains

  function matterDarkEnergyCPLConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`cosmologyFunctionsMatterDarkEnergyCPL` cosmological functions class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergyCPL)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    double precision                                                       :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateWA

    !![
    <inputParameter docformat="rst">
      <name>darkEnergyEquationOfStateW0</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      The equation of state parameter for dark energy, :math:`w_0`, defined such that :math:`P=\rho^w` with
      :math:`w(a)=w_0+w_\mathrm{a}(1-a)`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>darkEnergyEquationOfStateWA</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The equation of state parameter for dark energy, :math:`w_\mathrm{a}`, defined such that :math:`P=\rho^w` with
      :math:`w(a)=w_0+w_\mathrm{a}(1-a)`.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    ! Use it to construct a matter plus dark energy cosmological functions class.
    self=cosmologyFunctionsMatterDarkEnergyCPL(                             &
         &                                     cosmologyParameters_       , &
         &                                     darkEnergyEquationOfStateW0, &
         &                                     darkEnergyEquationOfStateWA  &
         &                                    )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function matterDarkEnergyCPLConstructorParameters

  function matterDarkEnergyCPLConstructorInternal(cosmologyParameters_,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateWA) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`cosmologyFunctionsMatterDarkEnergyCPL` cosmological functions class.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergyCPL)                        :: self
    class           (cosmologyParametersClass             ), intent(in   ), target :: cosmologyParameters_
    double precision                                       , intent(in   )         :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateWA
    !![
    <constructorAssign variables="*cosmologyParameters_, darkEnergyEquationOfStateW0, darkEnergyEquationOfStateWA"/>
    !!]

    ! The equation of state parameter of our parent class corresponds to a different parameterization of the dark energy equation
    ! of state and is unused by this class - set it to zero so that it is always well-defined.
    self%darkEnergyEquationOfStateW1=0.0d0
    call self%initialize()
    return
  end function matterDarkEnergyCPLConstructorInternal

  double precision function matterDarkEnergyCPLEquationOfStateDarkEnergy(self,time,expansionFactor)
    !!{RST
    Return the dark energy equation of state, :math:`w(a)=w_0+w_\mathrm{a}(1-a)`.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyCPL), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: expansionFactor      , time
    double precision                                                                 :: expansionFactorActual

    matterDarkEnergyCPLEquationOfStateDarkEnergy=self%darkEnergyEquationOfStateW0
    if (self%darkEnergyEquationOfStateWA /= 0.0d0) then
       if      (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time           )) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       matterDarkEnergyCPLEquationOfStateDarkEnergy=+matterDarkEnergyCPLEquationOfStateDarkEnergy &
            &                                       +self%darkEnergyEquationOfStateWA             &
            &                                       *(1.0d0-expansionFactorActual)
    end if
    return
  end function matterDarkEnergyCPLEquationOfStateDarkEnergy

  double precision function matterDarkEnergyCPLExponentDarkEnergy(self,time,expansionFactor)
    !!{RST
    Return the dark energy exponent, :math:`\epsilon(a) \equiv \ln [\rho_\mathrm{DE}(a)/\rho_\mathrm{DE,0}] / \ln a`, such that
    :math:`\rho_\mathrm{DE}(a) = \rho_\mathrm{DE,0} a^{\epsilon(a)}`. For the CPL equation of state, :math:`w(a) =
    w_0+w_\mathrm{a}(1-a)`, integrating :math:`\d \ln \rho_\mathrm{DE} / \d \ln a = -3[1+w(a)]` gives

    .. math::

       \epsilon(a) = -3(1+w_0+w_\mathrm{a}) + 3 w_\mathrm{a} \frac{a-1}{\ln a}.

    The final term has a removable singularity at :math:`a=1`, where :math:`(a-1)/\ln a \rightarrow 1` and, therefore,
    :math:`\epsilon(1) = -3(1+w_0)` as expected. Writing :math:`x = \ln a`, that term is :math:`(\mathrm{e}^x-1)/x = \sum_{n \ge
    1} x^{n-1}/n!`, which is used for small :math:`|x|` to avoid catastrophic cancellation.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyCPL), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: expansionFactor      , time
    double precision                                                                 :: expansionFactorActual, expansionFactorLogarithmic

    matterDarkEnergyCPLExponentDarkEnergy=-3.0d0*(1.0d0+self%darkEnergyEquationOfStateW0)
    if (self%darkEnergyEquationOfStateWA /= 0.0d0) then
       if      (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time           )) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       expansionFactorLogarithmic=log(expansionFactorActual)
       if (abs(expansionFactorLogarithmic) < matterDarkEnergyCPLExpansionFactorLogarithmicSmall) then
          ! Series expansion of 3wₐ[(a-1)/ln a - 1] in x=ln a, valid close to a=1.
          matterDarkEnergyCPLExponentDarkEnergy=+matterDarkEnergyCPLExponentDarkEnergy &
               &                                +3.0d0                                 &
               &                                *self%darkEnergyEquationOfStateWA      &
               &                                *expansionFactorLogarithmic            &
               &                                *(                                     &
               &                                  +1.0d0/2.0d0                         &
               &                                  +      expansionFactorLogarithmic    &
               &                                  /      6.0d0                         &
               &                                  +      expansionFactorLogarithmic**2 &
               &                                  /     24.0d0                         &
               &                                 )
       else
          matterDarkEnergyCPLExponentDarkEnergy=+matterDarkEnergyCPLExponentDarkEnergy &
               &                                -3.0d0                                 &
               &                                *self%darkEnergyEquationOfStateWA      &
               &                                +3.0d0                                 &
               &                                *self%darkEnergyEquationOfStateWA      &
               &                                *(expansionFactorActual-1.0d0)         &
               &                                / expansionFactorLogarithmic
       end if
    end if
    return
  end function matterDarkEnergyCPLExponentDarkEnergy

  double precision function matterDarkEnergyCPLExponentDarkEnergyDerivative(self,time,expansionFactor)
    !!{RST
    Return the derivative of the dark energy exponent with respect to expansion factor,

    .. math::

       \frac{\d \epsilon}{\d a} = 3 w_\mathrm{a} \frac{\ln a - (a-1)/a}{\ln^2 a}.

    As for :math:`\epsilon(a)` itself this has a removable singularity at :math:`a=1`, where it tends to :math:`3w_\mathrm{a}/2`.
    Writing :math:`x = \ln a` the ratio is :math:`(x-1+\mathrm{e}^{-x})/x^2 = \sum_{n \ge 2} (-1)^n x^{n-2}/n!`, which is used for
    small :math:`|x|`.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyCPL), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: expansionFactor      , time
    double precision                                                                 :: expansionFactorActual, expansionFactorLogarithmic

    if      (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time           )) then
       expansionFactorActual=self%expansionFactor(time)
    else
       if (self%darkEnergyEquationOfStateWA /= 0.0d0) call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
       expansionFactorActual=1.0d0
    end if
    expansionFactorLogarithmic=log(expansionFactorActual)
    if (abs(expansionFactorLogarithmic) < matterDarkEnergyCPLExpansionFactorLogarithmicSmall) then
       matterDarkEnergyCPLExponentDarkEnergyDerivative=+3.0d0                                 &
            &                                          *self%darkEnergyEquationOfStateWA      &
            &                                          *(                                     &
            &                                            +1.0d0/  2.0d0                       &
            &                                            -      expansionFactorLogarithmic    &
            &                                            /      6.0d0                         &
            &                                            +      expansionFactorLogarithmic**2 &
            &                                            /     24.0d0                         &
            &                                           )
    else
       matterDarkEnergyCPLExponentDarkEnergyDerivative=+3.0d0                                 &
            &                                          *self%darkEnergyEquationOfStateWA      &
            &                                          *(                                     &
            &                                            +expansionFactorLogarithmic          &
            &                                            -(expansionFactorActual-1.0d0)       &
            &                                            / expansionFactorActual              &
            &                                           )                                     &
            &                                          / expansionFactorLogarithmic**2
    end if
    return
  end function matterDarkEnergyCPLExponentDarkEnergyDerivative
