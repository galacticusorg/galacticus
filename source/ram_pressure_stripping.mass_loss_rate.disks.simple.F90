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

  !% Implementation of a simple ram pressure stripping of disks class.

  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForceClass

  !# <ramPressureStrippingDisks name="ramPressureStrippingDisksSimple">
  !#  <description>A simple model of ram pressure stripping in galactic disks.</description>
  !# </ramPressureStrippingDisks>
  type, extends(ramPressureStrippingDisksClass) :: ramPressureStrippingDisksSimple
     !% Implementation of a simple model of ram pressure stripping of galactic disks.
     private
     class           (hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_ => null()
     double precision                                        :: rateFractionalMaximum
   contains
     final     ::                 simpleDestructor
     procedure :: rateMassLoss => simpleRateMassLoss
  end type ramPressureStrippingDisksSimple

  interface ramPressureStrippingDisksSimple
     !% Constructors for the {\normalfont \ttfamily simple} model of ram pressure stripping of disks class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface ramPressureStrippingDisksSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily simple} timescale for star formation feedback in disks class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (ramPressureStrippingDisksSimple)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (hotHaloRamPressureForceClass   ), pointer       :: hotHaloRamPressureForce_
    double precision                                                 :: rateFractionalMaximum

    !# <inputParameter>
    !#   <name>rateFractionalMaximum</name>
    !#   <defaultValue>10.0d0</defaultValue>
    !#   <description>The maximum fractional mass loss rate per dynamical time in the simple model of mass loss from disks due to ram pressure stripping.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    self=ramPressureStrippingDisksSimple(rateFractionalMaximum,hotHaloRamPressureForce_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="hotHaloRamPressureForce_"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(rateFractionalMaximum,hotHaloRamPressureForce_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily simple} model of ram pressure stripping of disks class.
    implicit none
    type            (ramPressureStrippingDisksSimple)                        :: self
    double precision                                 , intent(in   )         :: rateFractionalMaximum
    class           (hotHaloRamPressureForceClass   ), intent(in   ), target :: hotHaloRamPressureForce_
    !# <constructorAssign variables="rateFractionalMaximum, *hotHaloRamPressureForce_"/>

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily simple} model of ram pressure stripping of disks class.
    implicit none
    type(ramPressureStrippingDisksSimple), intent(inout) :: self

    !# <objectDestructor name="self%hotHaloRamPressureForce_"/>
    return
  end subroutine simpleDestructor

  double precision function simpleRateMassLoss(self,node)
    !% Computes the mass loss rate from disks due to ram pressure stripping assuming a simple model. Specifically, the mass loss
    !% rate is
    !% \begin{equation}
    !% \dot{M}_\mathrm{gas} = -\alpha M_\mathrm{gas}/\tau_\mathrm{disk},
    !% \end{equation}
    !% where
    !% \begin{equation}
    !% \alpha = F_\mathrm{ram}/F_\mathrm{gravity},
    !% \end{equation}
    !% $F_\mathrm{ram}$ is the ram pressure force from the hot halo (see \href{https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Development.pdf\#methods.hotHaloRamPressureForce}{here}), and
    !% \begin{equation}
    !% F_\mathrm{gravity} = 2 \pi \mathrm{G} \Sigma_\mathrm{gas}(r_{1/2}) \Sigma_\mathrm{total}(r_{1/2})
    !% \end{equation}
    !% is the gravitational restoring force in the disk at the half-mass radius, $r_\mathrm{1/2}$.
    use :: Galactic_Structure_Options          , only : componentTypeDisk                 , coordinateSystemCylindrical, massTypeAll, massTypeGaseous
    use :: Galactic_Structure_Surface_Densities, only : Galactic_Structure_Surface_Density
    use :: Galacticus_Nodes                    , only : nodeComponentDisk                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gigaYear                          , megaParsec
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Constants_Astronomical        , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes        , only : kilo
    implicit none
    class           (ramPressureStrippingDisksSimple), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    class           (nodeComponentDisk              ), pointer       :: disk
    double precision                                                 :: forceGravitational    , forceRamPressure   , &
         &                                                              rateMassLossFractional, radiusHalfMass     , &
         &                                                              surfaceDensityGas     , surfaceDensityTotal, &
         &                                                              timeDynamical

    ! Assume no mass loss rate due to ram pressure by default.
    simpleRateMassLoss=0.0d0
    ! Return immediately if this is not a satellite.
    if (.not.node%isSatellite()) return
    ! Get the disk component.
    disk                => node%disk                          (    )
    ! Get the ram pressure force due to the hot halo.
    forceRamPressure    =  self%hotHaloRamPressureForce_%force(node)
    ! Get the disk half-mass radius.
    radiusHalfMass      =  disk%halfMassRadius                (    )
    ! Compute the disk densities at the half mass radius.
    surfaceDensityGas   =  Galactic_Structure_Surface_Density(                                              &
         &                                                    node                                        , &
         &                                                    [radiusHalfMass,0.0d0,0.0d0]                , &
         &                                                    coordinateSystem=coordinateSystemCylindrical, &
         &                                                    massType        =massTypeGaseous            , &
         &                                                    componentType   =componentTypeDisk            &
         &                                                   )
    surfaceDensityTotal =  Galactic_Structure_Surface_Density(                                              &
         &                                                    node                                        , &
         &                                                    [radiusHalfMass,0.0d0,0.0d0]                , &
         &                                                    coordinateSystem=coordinateSystemCylindrical, &
         &                                                    massType        =massTypeAll                , &
         &                                                    componentType   =componentTypeDisk            &
         &                                                   )
    ! Compute the gravitational restoring force in the disk midplane.
    forceGravitational  =  +2.0d0                           &
         &                 *Pi                              &
         &                 *gravitationalConstantGalacticus &
         &                 *surfaceDensityGas               &
         &                 *surfaceDensityTotal
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (forceRamPressure < self%rateFractionalMaximum*forceGravitational) then
       rateMassLossFractional=forceRamPressure/forceGravitational
    else
       rateMassLossFractional=self%rateFractionalMaximum
    end if
    ! Compute the dynamical time.
    timeDynamical     =+megaParsec             &
         &             /kilo                   &
         &             /gigaYear               &
         &             *disk%radius  ()        &
         &             /disk%velocity()
    ! Compute the mass loss rate.
    simpleRateMassLoss=+rateMassLossFractional &
         &             *disk%massGas ()        &
         &             /timeDynamical
    return
  end function simpleRateMassLoss
