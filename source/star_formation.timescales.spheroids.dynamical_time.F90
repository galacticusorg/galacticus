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

  !% Implementation of a timescale for star formation in galactic spheroids which scales with the spheroid dynamical time.
  
  !# <starFormationTimescaleSpheroids name="starFormationTimescaleSpheroidsDynamicalTime">
  !#  <description>A timescale for star formation in galactic spheroids which scales with the spheroid dynamical time.</description>
  !# </starFormationTimescaleSpheroids>
  type, extends(starFormationTimescaleSpheroidsClass) :: starFormationTimescaleSpheroidsDynamicalTime
     !% Implementation of a timescale for star formation in galactic spheroids which scales with the dynamical time.
     private
     double precision :: efficiency      , exponentVelocity, &
          &              timescaleMinimum
   contains
     procedure :: timescale => dynamicalTimeTimescale
  end type starFormationTimescaleSpheroidsDynamicalTime

  interface starFormationTimescaleSpheroidsDynamicalTime
     !% Constructors for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in spheroids class.
     module procedure dynamicalTimeConstructorParameters
     module procedure dynamicalTimeConstructorInternal
  end interface starFormationTimescaleSpheroidsDynamicalTime

contains

  function dynamicalTimeConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in spheroids class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleSpheroidsDynamicalTime)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: efficiency      , exponentVelocity, &
         &                                                                           timescaleMinimum
    
    !# <inputParameter>
    !#   <name>efficiency</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The efficiency of star formation in spheroids for the dynamical time method.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-1.500d0</defaultValue>
    !#   <description>The velocity exponent for star formation in spheroids for the dynamical time method.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timescaleMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>The minimum timescale for star formation in spheroids.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=starFormationTimescaleSpheroidsDynamicalTime(efficiency,exponentVelocity,timescaleMinimum)
    !# <inputParametersValidate source="parameters"/>
    return
  end function dynamicalTimeConstructorParameters

  function dynamicalTimeConstructorInternal(efficiency,exponentVelocity,timescaleMinimum) result(self)
    !% Internal constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in spheroids class.
    use Galacticus_Error, only : Galacticus_Error_Report , Galacticus_Component_List
    use Array_Utilities
    use Galacticus_Nodes, only : defaultSpheroidComponent
    implicit none
    type            (starFormationTimescaleSpheroidsDynamicalTime)                :: self
    double precision                                              , intent(in   ) :: efficiency      , exponentVelocity, &
         &                                                                           timescaleMinimum
    !# <constructorAssign variables="efficiency, exponentVelocity, timescaleMinimum"/>

    if     (                                                                                                            &
         &  .not.(                                                                                                      &
         &         defaultSpheroidComponent%velocityIsGettable()                                                        &
         &        .and.                                                                                                 &
         &         defaultSpheroidComponent%  radiusIsGettable()                                                        &
         &       )                                                                                                      &
         & ) call Galacticus_Error_Report                                                                               &
         &        (                                                                                                     &
         &         'spheroid component must have gettable radius and velocity properties.'//                            &
         &         Galacticus_Component_List(                                                                           &
         &                                   'spheroid'                                                              ,  &
         &                                    defaultSpheroidComponent%velocityAttributeMatch(requireGettable=.true.)   &
         &                                   .intersection.                                                             &
         &                                    defaultSpheroidComponent%  radiusAttributeMatch(requireGettable=.true.)   &
         &                                  )                                                                        // &
         &         {introspection:location}                                                                             &
         &        )
    return
  end function dynamicalTimeConstructorInternal

  double precision function dynamicalTimeTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic spheroid of {\normalfont \ttfamily node}. The timescale is given by
    !% \begin{equation}
    !% \tau_\star = \epsilon_\star^{-1} \tau_\mathrm{dynamical, spheroid} \left( {V_\mathrm{spheroid} \over 200\hbox{km/s}} \right)^{\alpha_\star},
    !% \end{equation}
    !% where $\epsilon_\star$(={\normalfont \ttfamily efficiency}) is a star formation efficiency and $\alpha_\star$(={\normalfont \ttfamily
    !% exponentVelocity}) controls the scaling with velocity. Note that $\tau_\mathrm{dynamical,spheroid}=R_\mathrm{
    !% spheroid}/V_\mathrm{spheroid}$ where the radius and velocity are whatever characteristic values returned by the spheroid method. This
    !% scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000}, but that they specifically used the
    !% half-mass radius and circular velocity at that radius.
    use Numerical_Constants_Astronomical
    use Galacticus_Nodes                , only : nodeComponentSpheroid
    implicit none
    class           (starFormationTimescaleSpheroidsDynamicalTime), intent(inout)         :: self
    type            (treeNode                                    ), intent(inout), target :: node
    class           (nodeComponentSpheroid                       ), pointer               :: spheroid
    double precision                                              , parameter             :: velocityZeroPoint=200.0d+0               !   (km/s)
    double precision                                              , parameter             :: velocityTiny     =  1.0d-3               !   (km/s)
    double precision                                                                      :: velocitySpheroid         , timeDynamical

    spheroid         => node%spheroid    ()
    velocitySpheroid =  spheroid%velocity()
    ! Check for zero velocity spheroid.
    if (velocitySpheroid <= velocityTiny) then
       ! No well defined answer in this case.
       dynamicalTimeTimescale=0.0d0
    else if (self%efficiency == 0.0d0) then
       ! No star formation occurs if the efficiency is zero.
       dynamicalTimeTimescale=0.0d0
    else
       ! Get the dynamical time in Gyr.
       timeDynamical=+Mpc_per_km_per_s_To_Gyr &
            &        *spheroid%radius()       &
            &        /velocitySpheroid
       ! Compute the star formation timescale using a simple scaling factor.
       dynamicalTimeTimescale=max(                           &
            &                     +timeDynamical             &
            &                     *(                         &
            &                       +velocitySpheroid        &
            &                       /velocityZeroPoint       &
            &                      )**self%exponentVelocity  &
            &                     /self%efficiency         , &
            &                     +self%timescaleMinimum     &
            &                    )
    end if
    return
  end function dynamicalTimeTimescale
