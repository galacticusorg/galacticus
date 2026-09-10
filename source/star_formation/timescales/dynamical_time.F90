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

  !!{RST
  Implementation of a timescale for star formation which scales with the component dynamical time.
  !!}

  !![
  <starFormationTimescale name="starFormationTimescaleDynamicalTime" docformat="rst">
   <description>
   A star formation timescale class in which the star formation timescale scales with the component dynamical time. Specifically:

   .. math::

      \tau_\star = \epsilon_\star^{-1} \tau_\mathrm{dynamical} \left( {V \over 200\hbox{km/s}} \right)^{\alpha_\star},

   where :math:`\epsilon_\star=`\ ``[efficiency]`` and :math:`\alpha_\star=`\ ``[exponentVelocity]`` are input parameters, :math:`\tau_\mathrm{dynamical}\equiv r/V` is the dynamical timescale of the :term:`component` and :math:`r` and :math:`V` are the characteristic radius and velocity respectively of the component. The timescale is not allowed to fall below a minimum value specified by ``[timescaleMinimum]`` (in Gyr).
   </description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleDynamicalTime
     !!{RST
     Implementation of a timescale for star formation which scales with the dynamical time.
     !!}
     private
     double precision :: efficiency                 , exponentVelocity , &
          &              timescaleMinimum
     logical          :: diskSupported              , spheroidSupported, &
          &              nuclearStarClusterSupported
   contains
     procedure :: timescale => dynamicalTimeTimescale
  end type starFormationTimescaleDynamicalTime

  interface starFormationTimescaleDynamicalTime
     !!{RST
     Constructors for the :galacticus-class:`starFormationTimescaleDynamicalTime` timescale for star formation class.
     !!}
     module procedure dynamicalTimeConstructorParameters
     module procedure dynamicalTimeConstructorInternal
  end interface starFormationTimescaleDynamicalTime

contains

  function dynamicalTimeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`starFormationTimescaleDynamicalTime` timescale for star formation class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleDynamicalTime)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: efficiency      , exponentVelocity, &
         &                                                                  timescaleMinimum

    !![
    <inputParameter docformat="rst">
      <name>efficiency</name>
      <defaultValue>0.01d0</defaultValue>
      <description>
      The efficiency of star formation for the dynamical time method.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentVelocity</name>
      <defaultValue>-1.50d0</defaultValue>
      <description>
      The velocity exponent for star formation for the dynamical time method.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>timescaleMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>
      The minimum allowed timescale for star formation (in Gyr) in the dynamical time prescription, preventing unphysically short formation timescales in high-density or high-velocity systems.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationTimescaleDynamicalTime(efficiency,exponentVelocity,timescaleMinimum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function dynamicalTimeConstructorParameters

  function dynamicalTimeConstructorInternal(efficiency,exponentVelocity,timescaleMinimum) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`starFormationTimescaleDynamicalTime` timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleDynamicalTime)                :: self
    double precision                                     , intent(in   ) :: efficiency      , exponentVelocity, &
         &                                                                  timescaleMinimum
    !![
    <constructorAssign variables="efficiency, exponentVelocity, timescaleMinimum"/>
    <componentPropertyAssert class="disk"     properties="velocity radius" require="gettable" assignTo="self%diskSupported"              />
    <componentPropertyAssert class="spheroid" properties="velocity radius" require="gettable" assignTo="self%spheroidSupported"          />
    <componentPropertyAssert class="NSC"      properties="velocity radius" require="gettable" assignTo="self%nuclearStarClusterSupported"/>
    !!]
    return
  end function dynamicalTimeConstructorInternal

  double precision function dynamicalTimeTimescale(self,component)
    !!{RST
    Returns the timescale (in Gyr) for star formation in the given ``component``. The timescale is given by

    .. math::

       \tau_\star = \epsilon_\star^{-1} \tau_\mathrm{dynamical} \left( {V \over 200\hbox{km/s}} \right)^{\alpha_\star},

    where :math:`\epsilon_\star`\ (=\ ``efficiency``) is a star formation efficiency and :math:`\alpha_\star`\ (=\ ``exponentVelocity``) controls the scaling with velocity. Note that :math:`\tau_\mathrm{dynamical}=R/V` where the radius and velocity are whatever characteristic values returned by the component. This scaling is functionally similar to that adopted by :cite:t:`cole_hierarchical_2000`, but they specifically used the half-mass radius and circular velocity at that radius.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponent    , nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (starFormationTimescaleDynamicalTime), intent(inout) :: self
    class           (nodeComponent                      ), intent(inout) :: component
    double precision                                     , parameter     :: velocityZeroPoint=200.0d0            !   (km/s)
    double precision                                                     :: velocity                 , radius, &
         &                                                                  timeDynamical

    select type (component)
    class is (nodeComponentDisk              )
       !![
       <componentPropertyAssert class="disk"     properties="velocity radius" require="gettable" condition="self%diskSupported"              />
       !!]
       velocity=component%velocity()
       radius  =component%radius  ()
    class is (nodeComponentSpheroid          )
       !![
       <componentPropertyAssert class="spheroid" properties="velocity radius" require="gettable" condition="self%spheroidSupported"          />
       !!]
       velocity=component%velocity()
       radius  =component%radius  ()
    class is (nodeComponentNSC)
       !![
       <componentPropertyAssert class="NSC"      properties="velocity radius" require="gettable" condition="self%nuclearStarClusterSupported"/>
       !!]
       velocity=component%velocity()
       radius  =component%radius  ()
    class default
       velocity=0.0d0
       radius  =0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Check for zero velocity.
    if (velocity <= 0.0d0) then
       ! No well defined answer in this case.
       dynamicalTimeTimescale=0.0d0
    else if (self%efficiency == 0.0d0) then
       ! No star formation occurs if the efficiency is zero.
       dynamicalTimeTimescale=0.0d0
    else
       ! Get the dynamical time in Gyr.
       timeDynamical=+MpcPerKmPerSToGyr &
            &        *radius            &
            &        /velocity
       ! Compute the star formation timescale using a simple scaling factor.
       dynamicalTimeTimescale=max(                           &
            &                     +timeDynamical             &
            &                     *(                         &
            &                       +velocity                &
            &                       /velocityZeroPoint       &
            &                      )**self%exponentVelocity  &
            &                     /self%efficiency         , &
            &                     +self%timescaleMinimum     &
            &                    )
    end if
    return
  end function dynamicalTimeTimescale
