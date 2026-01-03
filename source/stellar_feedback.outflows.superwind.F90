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
  Implementation of a ``superwind'' stellar feedback model.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsSuperWind">
   <description>
    A stellar feedback outflow class which implements a ``superwind''. The outflow rate is given by:
    \begin{equation}
     \dot{M}_\mathrm{outflow} = \beta_\mathrm{superwind} {\dot{E} \over E_\mathrm{canonical}} \left\{ \begin{array}{ll} \left(
     V_\mathrm{superwind}/V\right)^2 &amp; \hbox{ if } V &gt; V_\mathrm{superwind} \\ 1 &amp; \hbox{ otherwise,} \end{array}
     \right.
    \end{equation}
    where $V_\mathrm{superwind}=${\normalfont \ttfamily [velocityCharacteristic]} (in km/s) and
    $\beta_\mathrm{superwind}=${\normalfont \ttfamily [massLoading]} are input parameters, $V$ is the characteristic velocity
    of the component, $\dot{E}$ is the rate of energy input from stellar populations and $E_\mathrm{canonical}$ is the total
    energy input by a canonical stellar population normalized to $1 M_\odot$ after infinite time.
   </description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsSuperWind
     !!{
     Implementation of a ``superwind'' stellar feedback model.
     !!}
     private
     double precision :: velocityCharacteristic, massLoading
   contains
     procedure :: outflowRate => superWindOutflowRate
  end type stellarFeedbackOutflowsSuperWind

  interface stellarFeedbackOutflowsSuperWind
     !!{
     Constructors for the superwind stellar feedback class.
     !!}
     module procedure superWindConstructorParameters
     module procedure superWindConstructorInternal
  end interface stellarFeedbackOutflowsSuperWind

contains

  function superWindConstructorParameters(parameters) result(self)
    !!{
    Constructor for the superwind stellar feedback class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarFeedbackOutflowsSuperWind)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: velocityCharacteristic, massLoading

    !![
    <inputParameter>
      <name>velocityCharacteristic</name>
      <source>parameters</source>
      <defaultValue>200.0d0</defaultValue>
      <description>The velocity scale at which the \gls{sne}-driven superwind outflow rate transitions to a constant.</description>
    </inputParameter>
    <inputParameter>
      <name>massLoading</name>
      <source>parameters</source>
      <defaultValue>2.0d0</defaultValue>
      <description>The mass-loading of ``superwind'' outflows.</description>
    </inputParameter>
    !!]
    self=stellarFeedbackOutflowsSuperWind(velocityCharacteristic,massLoading)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function superWindConstructorParameters

  function superWindConstructorInternal(velocityCharacteristic,massLoading) result(self)
    !!{
    Internal constructor for the superwind stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsSuperWind)                :: self
    double precision                                  , intent(in   ) :: velocityCharacteristic, massLoading
    !![
    <constructorAssign variables="velocityCharacteristic, massLoading"/>
    !!]

    return
  end function superWindConstructorInternal

  subroutine superWindOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the given {\normalfont \ttfamily component}. The outflow
    rate is given by 
    \begin{equation}
    \dot{M}_\mathrm{outflow} = f_\mathrm{SW,0} \left\{ \begin{array}{ll} 1 & \hbox{ if } V_\mathrm{disk} < V_\mathrm{disk,SW} \\ (V_\mathrm{disk,SW}/V_\mathrm{disk})^2 &  \hbox{ if } V_\mathrm{disk} \ge V_\mathrm{disk,SW} \end{array} \right. ,
    \end{equation}
    where $V_\mathrm{disk,SW}=${\normalfont \ttfamily [velocityCharacteristic]} and $f_\mathrm{SW,0}=${\normalfont \ttfamily
    [massLoading]}. Note that the velocity $V_\mathrm{ disk}$ is whatever characteristic value returned by the disk
    method. This scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000} and \cite{baugh_can_2005},
    except that they specifically used the circular velocity at half-mass radius.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk                     , nodeComponentSpheroid
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (stellarFeedbackOutflowsSuperWind), intent(inout) :: self
    class           (nodeComponent                   ), intent(inout) :: component
    double precision                                  , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                  , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    double precision                                                  :: velocity           , outflowRateToStarFormationRate
    !$GLC attributes unused :: rateStarFormation
    
    select type (component)
    class is (nodeComponentDisk    )
       velocity=component%velocity()
    class is (nodeComponentSpheroid)
       velocity=component%velocity()
    class default
       velocity=0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Check for zero velocity disk.
    if (velocity <= 0.0d0) then
       rateOutflowExpulsive=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=+      self%massLoading                 &
            &                         *(                                      &
            &                           +    self%velocityCharacteristic      &
            &                           /max(                                 &
            &                                     velocity              ,     &
            &                                self%velocityCharacteristic      &
            &                               )                                 &
            &                         )**2
       rateOutflowExpulsive          =+outflowRateToStarFormationRate         &
            &                         *rateEnergyInput                        &
            &                         /feedbackEnergyInputAtInfinityCanonical
    end if
    rateOutflowEjective              =+0.0d0
    return
  end subroutine superWindOutflowRate
