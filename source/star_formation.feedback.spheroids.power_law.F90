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

  !% Implementation of a power-law outflow rate due to star formation feedback in galactic spheroids.
  
  !# <starFormationFeedbackSpheroids name="starFormationFeedbackSpheroidsPowerLaw">
  !#  <description>A power-law outflow rate due to star formation feedback in galactic spheroids.</description>
  !# </starFormationFeedbackSpheroids>
  type, extends(starFormationFeedbackSpheroidsClass) :: starFormationFeedbackSpheroidsPowerLaw
     !% Implementation of a power-law outflow rate due to star formation feedback in galactic spheroids.
     private
     double precision :: velocityCharacteristic_, exponent
   contains
     !@ <objectMethods>
     !@   <object>starFormationFeedbackSpheroidsPowerLaw</object>
     !@   <objectMethod>
     !@     <method>velocityCharacteristic</method>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} node\arginout</arguments>
     !@     <type>\doublezero</type>
     !@     <description>Return the characteristic velocity for power law spheroid feedback models.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: outflowRate            => powerLawOutflowRate
     procedure :: velocityCharacteristic => powerLawVelocityCharacteristic
  end type starFormationFeedbackSpheroidsPowerLaw

  interface starFormationFeedbackSpheroidsPowerLaw
     !% Constructors for the power-law star formation feedback in spheroids class.
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface starFormationFeedbackSpheroidsPowerLaw

contains

  function powerLawConstructorParameters(parameters) result(self)
    !% Constructor for the power-law star formation feedback in spheroids class which takes a parameter set as input.
    use Galacticus_Error
    implicit none
    type            (starFormationFeedbackSpheroidsPowerLaw)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: velocityCharacteristic, exponent

    !# <inputParameter>
    !#   <name>velocityCharacteristic</name>
    !#   <source>parameters</source>
    !#   <defaultValue>250.0d0</defaultValue>
    !#   <description>The velocity scale at which the \gls{sne}-driven outflow rate equals the star formation rate in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponent</name>
    !#   <source>parameters</source>
    !#   <defaultValue>3.5d0</defaultValue>
    !#   <description>The velocity scaling of the \gls{sne}-driven outflow rate in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=starFormationFeedbackSpheroidsPowerLaw(velocityCharacteristic,exponent)
    !# <inputParametersValidate source="parameters"/>
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(velocityCharacteristic_,exponent) result(self)
    !% Internal constructor for the power-law star formation feedback from spheroids class.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (starFormationFeedbackSpheroidsPowerLaw)                :: self
    double precision                                        , intent(in   ) :: velocityCharacteristic_, exponent
    character       (len=13                                )                :: label
    !# <constructorAssign variables="velocityCharacteristic_, exponent"/>    

    if (velocityCharacteristic_ < 0.0d0) then
       write (label,'(e13.6)') velocityCharacteristic_
       call Galacticus_Error_Report('characteristic velocity must be non-negative ['//trim(adjustl(label))//' < 0]'//{introspection:location})
    end if
    return
  end function powerLawConstructorInternal

  double precision function powerLawOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic spheroid of {\normalfont \ttfamily
    !% thisNode}. The outflow rate is given by \begin{equation} \dot{M}_\mathrm{outflow} = \left({V_\mathrm{spheroid,outflow} \over
    !% V_\mathrm{spheroid}}\right)^{\alpha_\mathrm{spheroid,outflow}}, \end{equation} where $V_\mathrm{spheroid,outflow}$(={\normalfont
    !% \ttfamily velocityCharacteristic}) is the velocity scale at which outflow rate equals star formation rate and $\alpha_{\mathrm
    !% spheroid,outflow}$(={\normalfont \ttfamily exponent}) controls the scaling with velocity. Note that the velocity
    !% $V_\mathrm{spheroid}$ is whatever characteristic value returned by the spheroid method. This scaling is functionally similar to
    !% that adopted by \cite{cole_hierarchical_2000}, but that they specifically used the circular velocity at half-mass radius.
    use Stellar_Feedback
    use Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class           (starFormationFeedbackSpheroidsPowerLaw), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: rateEnergyInput , rateStarFormation
    class           (nodeComponentSpheroid                 ), pointer       :: spheroid
    double precision                                                        :: velocitySpheroid
    !GCC$ attributes unused :: rateStarFormation

    ! Get spheroid circular velocity.
    spheroid         => node    %spheroid()
    velocitySpheroid =  spheroid%velocity()
    ! Check for zero velocity spheroid.
    if (velocitySpheroid <= 0.0d0) then
       ! No well defined answer in this case.
       powerLawOutflowRate=+0.0d0
    else
       powerLawOutflowRate=+(                                      &
            &                +self%velocityCharacteristic(node)    &
            &                /     velocitySpheroid                &
            &               )**self%exponent                       &
            &              *rateEnergyInput                        &
            &              /feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function powerLawOutflowRate

  double precision function powerLawVelocityCharacteristic(self,node)
    !% Return the characteristic velocity for power-law feedback models in spheroids. In this case the characteristic velocity is a
    !% constant.
    implicit none
    class(starFormationFeedbackSpheroidsPowerLaw), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    !GCC$ attributes unused :: node

    powerLawVelocityCharacteristic=self%velocityCharacteristic_
    return
  end function powerLawVelocityCharacteristic
