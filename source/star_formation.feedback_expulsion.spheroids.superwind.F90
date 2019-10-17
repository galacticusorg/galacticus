!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a ``superwind'' expulsive outflow rate due to star formation feedback in galactic spheroids.

  !# <starFormationExpulsiveFeedbackSpheroids name="starFormationExpulsiveFeedbackSpheroidsSuperWind">
  !#  <description>A superwind expulsive outflow rate due to star formation feedback in galactic spheroids.</description>
  !# </starFormationExpulsiveFeedbackSpheroids>
  type, extends(starFormationExpulsiveFeedbackSpheroidsClass) :: starFormationExpulsiveFeedbackSpheroidsSuperWind
     !% Implementation of a superwind expulsive outflow rate due to star formation feedback in galactic spheroids.
     private
     double precision :: velocityCharacteristic, massLoading
   contains
     procedure :: outflowRate => superWindOutflowRate
  end type starFormationExpulsiveFeedbackSpheroidsSuperWind

  interface starFormationExpulsiveFeedbackSpheroidsSuperWind
     !% Constructors for the superwind expulsive star formation feedback in spheroids class.
     module procedure superWindConstructorParameters
     module procedure superWindConstructorInternal
  end interface starFormationExpulsiveFeedbackSpheroidsSuperWind

contains

  function superWindConstructorParameters(parameters) result(self)
    !% Constructor for the superwind expulsive star formation feedback in spheroids class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationExpulsiveFeedbackSpheroidsSuperWind)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    double precision                                                                  :: velocityCharacteristic, massLoading

    !# <inputParameter>
    !#   <name>velocityCharacteristic</name>
    !#   <source>parameters</source>
    !#   <defaultValue>200.0d0</defaultValue>
    !#   <description>The velocity scale at which the \gls{sne}-driven superwind outflow rate transitions to a constant in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massLoading</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>The mass-loading of ``superwind'' expulsive outflows in spheroids..</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=starFormationExpulsiveFeedbackSpheroidsSuperWind(velocityCharacteristic,massLoading)
    !# <inputParametersValidate source="parameters"/>
    return
  end function superWindConstructorParameters

  function superWindConstructorInternal(velocityCharacteristic,massLoading) result(self)
    !% Internal constructor for the superwind expulsive star formation feedback from spheroids class.
    implicit none
    type            (starFormationExpulsiveFeedbackSpheroidsSuperWind)                :: self
    double precision                                                  , intent(in   ) :: velocityCharacteristic, massLoading
    !# <constructorAssign variables="velocityCharacteristic, massLoading"/>

    return
  end function superWindConstructorInternal

  double precision function superWindOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns the expulsive outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic spheroid of {\normalfont \ttfamily node}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_\mathrm{outflow} = f_\mathrm{SW,0} \left\{ \begin{array}{ll} 1 & \hbox{ if } V_\mathrm{spheroid} < V_\mathrm{spheroid,SW} \\ (V_\mathrm{spheroid,SW}/V_\mathrm{spheroid})^2 &  \hbox{ if } V_\mathrm{spheroid} \ge V_\mathrm{spheroid,SW} \end{array} \right. ,
    !% \end{equation}
    !% where $V_\mathrm{spheroid,SW}=${\normalfont \ttfamily [velocityCharacteristic]} and $f_\mathrm{SW,0}=${\normalfont \ttfamily
    !% [massLoadnig]}. Note that the velocity $V_\mathrm{ spheroid}$ is whatever characteristic value returned by the spheroid
    !% method. This scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000} and \cite{baugh_can_2005},
    !% except that they specifically used the circular velocity at half-mass radius.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid                 , treeNode
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (starFormationExpulsiveFeedbackSpheroidsSuperWind), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    class           (nodeComponentSpheroid                           ), pointer       :: spheroid
    double precision                                                  , intent(in   ) :: rateEnergyInput , rateStarFormation
    double precision                                                                  :: velocitySpheroid, outflowRateToStarFormationRate
    !GCC$ attributes unused :: rateStarFormation

    spheroid         => node%spheroid    ()
    velocitySpheroid =  spheroid%velocity()
    ! Check for zero velocity spheroid.
    if (velocitySpheroid <= 0.0d0) then
       superWindOutflowRate=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=+      self%massLoading                 &
            &                         *(                                      &
            &                           +    self%velocityCharacteristic      &
            &                           /max(                                 &
            &                                     velocitySpheroid      ,     &
            &                                self%velocityCharacteristic      &
            &                               )                                 &
            &                         )**2
       superWindOutflowRate          =+outflowRateToStarFormationRate         &
            &                         *rateEnergyInput                        &
            &                         /feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function superWindOutflowRate
