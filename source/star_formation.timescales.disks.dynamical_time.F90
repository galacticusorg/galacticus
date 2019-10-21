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

  !% Implementation of a timescale for star formation in galactic disks which scales with the disk dynamical time.

  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksDynamicalTime">
  !#  <description>A timescale for star formation in galactic disks which scales with the disk dynamical time.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksDynamicalTime
     !% Implementation of a timescale for star formation in galactic disks which scales with the dynamical time.
     private
     double precision :: efficiency      , exponentVelocity, &
          &              timescaleMinimum
   contains
     procedure :: timescale => dynamicalTimeTimescale
  end type starFormationTimescaleDisksDynamicalTime

  interface starFormationTimescaleDisksDynamicalTime
     !% Constructors for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in disks class.
     module procedure dynamicalTimeConstructorParameters
     module procedure dynamicalTimeConstructorInternal
  end interface starFormationTimescaleDisksDynamicalTime

contains

  function dynamicalTimeConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in disks class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleDisksDynamicalTime)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    double precision                                                          :: efficiency      , exponentVelocity, &
         &                                                                       timescaleMinimum

    !# <inputParameter>
    !#   <name>efficiency</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The efficiency of star formation in disks for the dynamical time method.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-1.500d0</defaultValue>
    !#   <description>The velocity exponent for star formation in disks for the dynamical time method.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timescaleMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>The minimum timescale for star formation in disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=starFormationTimescaleDisksDynamicalTime(efficiency,exponentVelocity,timescaleMinimum)
    !# <inputParametersValidate source="parameters"/>
    return
  end function dynamicalTimeConstructorParameters

  function dynamicalTimeConstructorInternal(efficiency,exponentVelocity,timescaleMinimum) result(self)
    !% Internal constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation in disks class.
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Galacticus_Error, only : Galacticus_Component_List, Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultDiskComponent
    implicit none
    type            (starFormationTimescaleDisksDynamicalTime)                :: self
    double precision                                          , intent(in   ) :: efficiency      , exponentVelocity, &
         &                                                                       timescaleMinimum
    !# <constructorAssign variables="efficiency, exponentVelocity, timescaleMinimum"/>

    if     (                                                                                                        &
         &  .not.(                                                                                                  &
         &         defaultDiskComponent%velocityIsGettable()                                                        &
         &        .and.                                                                                             &
         &         defaultDiskComponent%  radiusIsGettable()                                                        &
         &       )                                                                                                  &
         & ) call Galacticus_Error_Report                                                                           &
         &        (                                                                                                 &
         &         'disk component must have gettable radius and velocity properties.'//                            &
         &         Galacticus_Component_List(                                                                       &
         &                                   'disk'                                                              ,  &
         &                                    defaultDiskComponent%velocityAttributeMatch(requireGettable=.true.)   &
         &                                   .intersection.                                                         &
         &                                    defaultDiskComponent%  radiusAttributeMatch(requireGettable=.true.)   &
         &                                  )                                                                    // &
         &         {introspection:location}                                                                         &
         &        )
    return
  end function dynamicalTimeConstructorInternal

  double precision function dynamicalTimeTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node}. The timescale is given by
    !% \begin{equation}
    !% \tau_\star = \epsilon_\star^{-1} \tau_\mathrm{dynamical, disk} \left( {V_\mathrm{disk} \over 200\hbox{km/s}} \right)^{\alpha_\star},
    !% \end{equation}
    !% where $\epsilon_\star$(={\normalfont \ttfamily effiency}) is a star formation efficiency and $\alpha_\star$(={\normalfont \ttfamily
    !% exponentVelocity}) controls the scaling with velocity. Note that $\tau_\mathrm{dynamical,disk}=R_\mathrm{
    !% disk}/V_\mathrm{disk}$ where the radius and velocity are whatever characteristic values returned by the disk method. This
    !% scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000}, but that they specifically used the
    !% half-mass radius and circular velocity at that radius.
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , treeNode
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (starFormationTimescaleDisksDynamicalTime), intent(inout), target :: self
    type            (treeNode                                ), intent(inout), target :: node
    class           (nodeComponentDisk                       ), pointer               :: disk
    double precision                                          , parameter             :: velocityZeroPoint=200.0d0                !   (km/s)
    double precision                                                                  :: velocityDisk             , timeDynamical

    disk         => node%disk    ()
    velocityDisk =  disk%velocity()
    ! Check for zero velocity disk.
    if (velocityDisk <= 0.0d0) then
       ! No well defined answer in this case.
       dynamicalTimeTimescale=0.0d0
    else if (self%efficiency == 0.0d0) then
       ! No star formation occurs if the efficiency is zero.
       dynamicalTimeTimescale=0.0d0
    else
       ! Get the dynamical time in Gyr.
       timeDynamical=+Mpc_per_km_per_s_To_Gyr &
            &        *disk%radius()           &
            &        /velocityDisk
       ! Compute the star formation timescale using a simple scaling factor.
       dynamicalTimeTimescale=max(                           &
            &                   +timeDynamical             &
            &                   *(                         &
            &                     +velocityDisk            &
            &                     /velocityZeroPoint       &
            &                    )**self%exponentVelocity  &
            &                   /self%efficiency         , &
            &                   +self%timescaleMinimum     &
            &                  )
    end if
    return
  end function dynamicalTimeTimescale
