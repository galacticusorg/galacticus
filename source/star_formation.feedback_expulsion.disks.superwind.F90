!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a ``superwind'' expulsive outflow rate due to star formation feedback in galactic disks.

  !# <starFormationExpulsiveFeedbackDisks name="starFormationExpulsiveFeedbackDisksSuperWind">
  !#  <description>A superwind expulsive outflow rate due to star formation feedback in galactic disks.</description>
  !# </starFormationExpulsiveFeedbackDisks>
  type, extends(starFormationExpulsiveFeedbackDisksClass) :: starFormationExpulsiveFeedbackDisksSuperWind
     !% Implementation of a superwind expulsive outflow rate due to star formation feedback in galactic disks.
     private
     double precision :: velocityCharacteristic, massLoading
   contains
     procedure :: outflowRate => superWindOutflowRate
  end type starFormationExpulsiveFeedbackDisksSuperWind

  interface starFormationExpulsiveFeedbackDisksSuperWind
     !% Constructors for the superwind expulsive star formation feedback in disks class.
     module procedure superWindConstructorParameters
     module procedure superWindConstructorInternal
  end interface starFormationExpulsiveFeedbackDisksSuperWind

contains

  function superWindConstructorParameters(parameters) result(self)
    !% Constructor for the superwind expulsive star formation feedback in disks class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationExpulsiveFeedbackDisksSuperWind)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: velocityCharacteristic, massLoading

    !# <inputParameter>
    !#   <name>velocityCharacteristic</name>
    !#   <source>parameters</source>
    !#   <defaultValue>200.0d0</defaultValue>
    !#   <description>The velocity scale at which the \gls{sne}-driven superwind outflow rate transitions to a constant in disks.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massLoading</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>The mass-loading of ``superwind'' expulsive outflows in disks..</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=starFormationExpulsiveFeedbackDisksSuperWind(velocityCharacteristic,massLoading)
    !# <inputParametersValidate source="parameters"/>
    return
  end function superWindConstructorParameters

  function superWindConstructorInternal(velocityCharacteristic,massLoading) result(self)
    !% Internal constructor for the superwind expulsive star formation feedback from disks class.
    implicit none
    type            (starFormationExpulsiveFeedbackDisksSuperWind)                :: self
    double precision                                              , intent(in   ) :: velocityCharacteristic, massLoading
    !# <constructorAssign variables="velocityCharacteristic, massLoading"/>

    return
  end function superWindConstructorInternal

  double precision function superWindOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns the expulsive outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\normalfont \ttfamily node}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_\mathrm{outflow} = f_\mathrm{SW,0} \left\{ \begin{array}{ll} 1 & \hbox{ if } V_\mathrm{disk} < V_\mathrm{disk,SW} \\ (V_\mathrm{disk,SW}/V_\mathrm{disk})^2 &  \hbox{ if } V_\mathrm{disk} \ge V_\mathrm{disk,SW} \end{array} \right. ,
    !% \end{equation}
    !% where $V_\mathrm{disk,SW}=${\normalfont \ttfamily [velocityCharacteristic]} and $f_\mathrm{SW,0}=${\normalfont \ttfamily
    !% [massLoadnig]}. Note that the velocity $V_\mathrm{ disk}$ is whatever characteristic value returned by the disk
    !% method. This scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000} and \cite{baugh_can_2005},
    !% except that they specifically used the circular velocity at half-mass radius.
    use :: Galacticus_Nodes, only : nodeComponentDisk                     , treeNode
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (starFormationExpulsiveFeedbackDisksSuperWind), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    class           (nodeComponentDisk                           ), pointer       :: disk
    double precision                                              , intent(in   ) :: rateEnergyInput, rateStarFormation
    double precision                                                              :: velocityDisk   , outflowRateToStarFormationRate
    !$GLC attributes unused :: rateStarFormation

    disk         => node%disk    ()
    velocityDisk =  disk%velocity()
    ! Check for zero velocity disk.
    if (velocityDisk <= 0.0d0) then
       superWindOutflowRate=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=+      self%massLoading                 &
            &                         *(                                      &
            &                           +    self%velocityCharacteristic      &
            &                           /max(                                 &
            &                                     velocityDisk          ,     &
            &                                self%velocityCharacteristic      &
            &                               )                                 &
            &                         )**2
       superWindOutflowRate          =+outflowRateToStarFormationRate         &
            &                         *rateEnergyInput                        &
            &                         /feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function superWindOutflowRate
