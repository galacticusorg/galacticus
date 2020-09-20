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

  !% Implementation of a simple tidal stripping of spheroids class.

  use :: Satellites_Tidal_Fields, only : satelliteTidalField, satelliteTidalFieldClass

  !# <tidalStrippingSpheroids name="tidalStrippingSpheroidsSimple">
  !#  <description>A simple model of tidal stripping in galactic spheroids.</description>
  !# </tidalStrippingSpheroids>
  type, extends(tidalStrippingSpheroidsClass) :: tidalStrippingSpheroidsSimple
     !% Implementation of a simple model of tidal stripping of galactic spheroids.
     private
     class           (satelliteTidalFieldClass), pointer :: satelliteTidalField_ => null()
     double precision                                    :: rateFractionalMaximum
   contains
     final     ::                 simpleDestructor
     procedure :: rateMassLoss => simpleRateMassLoss
  end type tidalStrippingSpheroidsSimple

  interface tidalStrippingSpheroidsSimple
     !% Constructors for the {\normalfont \ttfamily simple} model of tidal stripping of spheroids class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface tidalStrippingSpheroidsSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily simple} timescale for star formation feedback in spheroids class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (tidalStrippingSpheroidsSimple)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (satelliteTidalFieldClass     ), pointer       :: satelliteTidalField_
    double precision                                               :: rateFractionalMaximum

    !# <inputParameter>
    !#   <name>rateFractionalMaximum</name>
    !#   <defaultValue>10.0d0</defaultValue>
    !#   <description>The maximum fractional mass loss rate per dynamical time in the simple model of mass loss from spheroids due to tidal stripping.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="satelliteTidalField" name="satelliteTidalField_" source="parameters"/>
    self=tidalStrippingSpheroidsSimple(rateFractionalMaximum,satelliteTidalField_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="satelliteTidalField_"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(rateFractionalMaximum,satelliteTidalField_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily simple} model of tidal stripping of spheroids class.
    implicit none
    type            (tidalStrippingSpheroidsSimple)                        :: self
    double precision                               , intent(in   )         :: rateFractionalMaximum
    class           (satelliteTidalFieldClass     ), intent(in   ), target :: satelliteTidalField_
    !# <constructorAssign variables="rateFractionalMaximum, *satelliteTidalField_"/>

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily simple} model of tidal stripping of spheroids class.
    implicit none
    type(tidalStrippingSpheroidsSimple), intent(inout) :: self

    !# <objectDestructor name="self%satelliteTidalField_"/>
    return
  end subroutine simpleDestructor

  double precision function simpleRateMassLoss(self,node)
    !% Computes the mass loss rate from spheroids due to tidal stripping assuming a simple model. Specifically, the mass loss
    !% rate is
    !% \begin{equation}
    !% \dot{M} = -\alpha M/\tau_\mathrm{spheroid},
    !% \end{equation}
    !% where
    !% \begin{equation}
    !% \alpha = F_\mathrm{tidal}/F_\mathrm{gravity},
    !% \end{equation}
    !% $F_\mathrm{tidal}=\mathcal{F}_\mathrm{tidal} r_{1/2}$, $\mathcal{F}_\mathrm{tidal}$ is the tidal field from the host halo (see \href{https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Development.pdf\#methods.satelliteTidalField}{here}), and
    !% \begin{equation}
    !% F_\mathrm{gravity} = V_{1/2}^2(r_{1/2})/r_{1/2}
    !% \end{equation}
    !% is the gravitational restoring force in the spheroid at the half-mass radius, $r_\mathrm{1/2}$.
    use :: Galactic_Structure_Rotation_Curves, only : Galactic_Structure_Rotation_Curve
    use :: Galacticus_Nodes                  , only : nodeComponentSpheroid            , treeNode
    use :: Numerical_Constants_Astronomical  , only : gigaYear                         , megaParsec
    use :: Numerical_Constants_Prefixes      , only : kilo
    implicit none
    class           (tidalStrippingSpheroidsSimple), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (nodeComponentSpheroid        ), pointer       :: spheroid
    double precision                                               :: forceGravitational    , forceTidal    , &
         &                                                            rateMassLossFractional, radiusHalfMass, &
         &                                                            tidalTensorRadial     , timeDynamical , &
         &                                                            velocityRotation

    ! Assume no mass loss rate due to tidal by default.
    simpleRateMassLoss=0.0d0
    ! Return immediately if this is not a satellite.
    if (.not.node%isSatellite()) return
    ! Get the spheroid component.
    spheroid          => node    %spheroid                              (    )
    ! Get the tidal field due to the host halo.
    tidalTensorRadial =  self    %satelliteTidalField_%tidalTensorRadial(node)
    ! Get the spheroid half-mass radius.
    radiusHalfMass    =  spheroid%halfMassRadius                        (    )
    ! Get the tidal force exerted at the half-mass radius.
    forceTidal        =  +tidalTensorRadial &
         &               *radiusHalfMass
    ! Return if the tidal field is compressive.
    if (forceTidal <= 0.0d0) return
    ! Compute the spheroid rotation curve.
    velocityRotation=Galactic_Structure_Rotation_Curve(                &
         &                                             node          , &
         &                                             radiusHalfMass  &
         &                                            )
    ! Compute the gravitational restoring force in the spheroid midplane.
    forceGravitational=+velocityRotation**2 &
         &             /radiusHalfMass
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (forceTidal < self%rateFractionalMaximum*forceGravitational) then
       rateMassLossFractional=+forceTidal         &
            &                 /forceGravitational
    else
       rateMassLossFractional=+self%rateFractionalMaximum
    end if
    ! Compute the dynamical time.
    timeDynamical     =+megaParsec               &
         &             /kilo                     &
         &             /gigaYear                 &
         &             *  spheroid%radius     () &
         &             /  spheroid%velocity   ()
    ! Compute the mass loss rate.
    simpleRateMassLoss=+rateMassLossFractional   &
         &             *(                        &
         &               +spheroid%massGas    () &
         &               +spheroid%massStellar() &
         &              )                        &
         &             /timeDynamical
    return
  end function simpleRateMassLoss
