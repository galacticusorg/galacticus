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

  !!{
  Implementation of a power-law stellar feedback model.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsPowerLaw">
   <description>
    A stellar feedback outflow class in which the outflow rate is a power-law in circular velocity. Specifically the outflow
    rate is given by:
    \begin{equation}
     \dot{M}_\mathrm{outflow} = \left({V_\mathrm{outflow} \over V}\right)^{\alpha_\mathrm{outflow}} {\dot{E} \over
     E_\mathrm{canonical}},
    \end{equation}
    where $V_\mathrm{outflow}=${\normalfont \ttfamily [velocityCharacteristic]} (in km/s) and
    $\alpha_\mathrm{outflow}=${\normalfont \ttfamily [exponent]} are input parameters, $V$ is the characteristic velocity of
    the component, $\dot{E}$ is the rate of energy input from stellar populations and $E_\mathrm{canonical}$ is the total
    energy input by a canonical stellar population normalized to $1 M_\odot$ after infinite time.
   </description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsPowerLaw
     !!{
     Implementation of a power-law stellar feedback model.
     !!}
     private
     double precision :: velocityCharacteristic_, exponent
   contains
     !![
     <methods>
       <method description="Return the characteristic velocity for power law feedback models." method="velocityCharacteristic" />
     </methods>
     !!]
     procedure :: outflowRate            => powerLawOutflowRate
     procedure :: velocityCharacteristic => powerLawVelocityCharacteristic
  end type stellarFeedbackOutflowsPowerLaw

  interface stellarFeedbackOutflowsPowerLaw
     !!{
     Constructors for the power-law stellar feedback class.
     !!}
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface stellarFeedbackOutflowsPowerLaw

contains

  function powerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the power-law stellar feedback class which takes a parameter set as input.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsPowerLaw)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: velocityCharacteristic, exponent

    !![
    <inputParameter>
      <name>velocityCharacteristic</name>
      <source>parameters</source>
      <defaultValue>250.0d0</defaultValue>
      <description>The velocity scale at which the \gls{sne}-driven outflow rate equals the star formation rate in disks.</description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <defaultValue>3.5d0</defaultValue>
      <description>The velocity scaling of the \gls{sne}-driven outflow rate in disks.</description>
    </inputParameter>
    !!]
    self=stellarFeedbackOutflowsPowerLaw(velocityCharacteristic,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(velocityCharacteristic_,exponent) result(self)
    !!{
    Internal constructor for the power-law stellar feedback class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (stellarFeedbackOutflowsPowerLaw)                :: self
    double precision                                 , intent(in   ) :: velocityCharacteristic_, exponent
    character       (len=13                         )                :: label
    !![
    <constructorAssign variables="velocityCharacteristic_, exponent"/>
    !!]

    if (velocityCharacteristic_ < 0.0d0) then
       write (label,'(e13.6)') velocityCharacteristic_
       call Error_Report('characteristic velocity must be non-negative ['//trim(adjustl(label))//' < 0]'//{introspection:location})
    end if
    return
  end function powerLawConstructorInternal

  subroutine powerLawOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the given {\normalfont \ttfamily
    component}. The outflow rate is given by \begin{equation} \dot{M}_\mathrm{outflow} = \left({V_\mathrm{outflow} \over
    V}\right)^{\alpha_\mathrm{outflow}}, \end{equation} where $V_\mathrm{outflow}$(={\normalfont
    \ttfamily velocityCharacteristic}) is the velocity scale at which outflow rate equals star formation rate and $\alpha_{\mathrm
    outflow}$(={\normalfont \ttfamily exponent}) controls the scaling with velocity. Note that the velocity
    $V$ is whatever characteristic value returned by the disk/spheroid component. This scaling is functionally similar to
    that adopted by \cite{cole_hierarchical_2000}, except that they specifically used the circular velocity at half-mass radius.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk                     , nodeComponentSpheroid
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (stellarFeedbackOutflowsPowerLaw), intent(inout) :: self
    class           (nodeComponent                  ), intent(inout) :: component
    double precision                                 , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                 , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    double precision                                                 :: velocity
    !$GLC attributes unused :: rateStarFormation

    ! Get velocity.
    select type (component)
    class is (nodeComponentDisk    )
       velocity=component%velocity()
    class is (nodeComponentSpheroid)
       velocity=component%velocity()
    class default
       velocity=0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Check for zero velocity.
    if (velocity <= 0.0d0) then
       ! No well defined answer in this case.
       rateOutflowEjective=+0.0d0
    else
       rateOutflowEjective=+(                                        &
            &                +self%velocityCharacteristic(component) &
            &                /     velocity                          &
            &               )**self%exponent                         &
            &              *rateEnergyInput                          &
            &              /feedbackEnergyInputAtInfinityCanonical
    end if
    rateOutflowExpulsive  =+0.0d0
    return
  end subroutine powerLawOutflowRate

  double precision function powerLawVelocityCharacteristic(self,component)
    !!{
    Return the characteristic velocity for power-law feedback models. In this case the characteristic velocity is a
    constant.
    !!}
    implicit none
    class(stellarFeedbackOutflowsPowerLaw), intent(inout) :: self
    class(nodeComponent                  ), intent(inout) :: component
    !$GLC attributes unused :: component

    powerLawVelocityCharacteristic=self%velocityCharacteristic_
    return
  end function powerLawVelocityCharacteristic
