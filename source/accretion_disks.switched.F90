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
  Implementation of a switched (ADAF/thin) accretion disk.
  !!}

  !![
  <accretionDisks name="accretionDisksSwitched" docformat="rst">
   <description>
   A circumnuclear accretion disk class in which accretion switches between radiatively efficient (Shakura-Sunyaev) and inefficient (:term:`ADAF`) modes. This is intended to crudely model the fact that accretion disks are expected to be radiatively inefficient at both high accretion rates (i.e. as they approach the Eddington luminosity the radiation pressure from a radiatively efficient flow would begin to disrupt the flow itself by definition), while accretion flows at low accretion rates (:math:`\dot{M}_\mathrm{acc} &lt; \alpha^2 \dot{M}_\mathrm{Edd}`, where :math:`\alpha\sim 0.1` is the usual parameter controlling the rate of angular momentum transport in a :cite:t:`shakura_black_1973` accretion disk) are also radiatively inefficient as radiative processes are too inefficient at the associated low densities to radiate energy at the rate it is being liberated. A more thorough discussion is given by :cite:t:`begelman_accreting_2014`.

   The properties of the switched disk (e.g. radiative efficiency, jet power), are a linear combination of those of the Shakura-Sunyaev and :term:`ADAF` modes, with the :term:`ADAF` fraction being given by:

   .. math::

      f_\mathrm{ADAF} = [1+\exp(y_\mathrm{min})]^{-1} + [1+\exp(y_\mathrm{max})]^{-1},
      \label{eq:SwitchedDiskADAFFraction}

   where

   .. math::

      y_\mathrm{min} &amp; = +\log(x/x_\mathrm{min})/\Delta_{x}, \\
      y_\mathrm{max} &amp; = -\log(x/x_\mathrm{max})/\Delta_{x},

   and,

   .. math::

      x = \dot{M} / \dot{M}_\mathrm{Eddington}.

   Here, :math:`x_\mathrm{min}=`\ ``[accretionRateThinDiskMinimum]``, :math:`x_\mathrm{max}=`\ ``[accretionRateThinDiskMaximum]``, and :math:`\Delta_{x}`\ =\ ``[accretionRateTransitionWidth]``. If either ``[accretionRateThinDiskMinimum]`` or ``[accretionRateThinDiskMaximum]`` is set to "``none``" then the corresponding term in eqn. () is excluded.

   Additionally, if ``[scaleADAFRadiativeEfficiency]`` is set to ``true`` then the radiative efficiency of the :term:`ADAF` component is reduced by a factor :math:`x/x_\mathrm{min}` when :math:`x&lt;x_\mathrm{min}`.
   </description>
  </accretionDisks>
  !!]
  type, extends(accretionDisksClass) :: accretionDisksSwitched
     !!{RST
     Implementation of an accretion disk class in which accretion switches between thin-disk and ADAF modes.
     !!}
     private
     ! Parameters controlling the range of accretion rates over which the accretion disk will be an ADAF.
     class           (accretionDisksClass), pointer :: accretionDisksADAF_                     => null(), accretionDisksShakuraSunyaev_           => null()
     double precision                               :: accretionRateThinDiskMaximum                     , accretionRateThinDiskMinimum                     , &
          &                                            accretionRateTransitionWidth
     double precision                               :: accretionRateThinDiskMaximumLogarithmic          , accretionRateThinDiskMinimumLogarithmic
     logical                                        :: accretionRateThinDiskMaximumExists               , accretionRateThinDiskMinimumExists
     ! Option controlling ADAF radiative efficiency.
     logical                                        :: scaleADAFRadiativeEfficiency
   contains
     !![
     <methods>
       <method description="Return the fraction of the accretion flow to be represented as an ADAF."                        method="fractionADAF"                  />
       <method description="Return the scaling of radiative efficiency of the ADAF component in a switched accretion disk." method="efficiencyRadiativeScalingADAF"/>
     </methods>
     !!]
     final     ::                                   switchedDestructor
     procedure :: efficiencyRadiative            => switchedEfficiencyRadiative
     procedure :: powerJet                       => switchedPowerJet
     procedure :: rateSpinUp                     => switchedRateSpinUp
     procedure :: fractionADAF                   => switchedFractionADAF
     procedure :: efficiencyRadiativeScalingADAF => switchedEfficiencyRadiativeScalingADAF
  end type accretionDisksSwitched

  interface accretionDisksSwitched
     !!{RST
     Constructors for the switched accretion disk class.
     !!}
     module procedure switchedConstructorParameters
     module procedure switchedConstructorInternal
  end interface accretionDisksSwitched

contains

  function switchedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the switched accretion disk class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionDisksSwitched)                :: self
    type            (inputParameters       ), intent(inout) :: parameters
    class           (accretionDisksClass   ), pointer       :: accretionDisksADAF_          , accretionDisksShakuraSunyaev_
    character       (len=128               )                :: accretionRateThinDiskMinimum , accretionRateThinDiskMaximum
    double precision                                        :: accretionRateTransitionWidth , accretionRateThinDiskMinimum_, &
         &                                                     accretionRateThinDiskMaximum_
    logical                                                 :: scaleADAFRadiativeEfficiency

    !![
    <inputParameter docformat="rst">
      <name>accretionRateThinDiskMinimum</name>
      <source>parameters</source>
      <defaultValue>'0.01d0'</defaultValue>
      <description>
      The accretion rate (in Eddington units) below which a switched accretion disk becomes an ADAF.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>accretionRateThinDiskMaximum</name>
      <source>parameters</source>
      <defaultValue>'0.3d0'</defaultValue>
      <description>
      The accretion rate (in Eddington units) above which a switched accretion disk becomes an ADAF.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>accretionRateTransitionWidth</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>
      The width (in :math:`\ln[\dot{M}/\dot{M}_\mathrm{Eddington}]`) over which transitions between accretion disk states occur.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scaleADAFRadiativeEfficiency</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      Specifies whether the radiative efficiency of the ADAF component in a switched accretion disk scales with accretion rate.
      </description>
    </inputParameter>
    <objectBuilder class="accretionDisks" parameterName="accretionDisksADAF"           name="accretionDisksADAF_"           source="parameters"/>
    <objectBuilder class="accretionDisks" parameterName="accretionDisksShakuraSunyaev" name="accretionDisksShakuraSunyaev_" source="parameters"/>
    !!]
    ! If minimum or maximum accretion rate for thin disk does not exist set suitable values.
    if (accretionRateThinDiskMinimum == "none") then
       accretionRateThinDiskMinimum_=-huge(0.0d0)
    else
       read (accretionRateThinDiskMinimum,*) accretionRateThinDiskMinimum_
    end if
    if (accretionRateThinDiskMaximum == "none") then
       accretionRateThinDiskMaximum_=+huge(0.0d0)
    else
       read (accretionRateThinDiskMaximum,*) accretionRateThinDiskMaximum_
    end if
    ! Build the object.
    self=accretionDisksSwitched(accretionDisksADAF_,accretionDisksShakuraSunyaev_,accretionRateThinDiskMinimum_,accretionRateThinDiskMaximum_,accretionRateTransitionWidth,scaleADAFRadiativeEfficiency)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionDisksADAF_"          />
    <objectDestructor name="accretionDisksShakuraSunyaev_"/>
    !!]
    return
  end function switchedConstructorParameters

  function switchedConstructorInternal(accretionDisksADAF_,accretionDisksShakuraSunyaev_,accretionRateThinDiskMinimum,accretionRateThinDiskMaximum,accretionRateTransitionWidth,scaleADAFRadiativeEfficiency) result(self)
    !!{RST
    Internal constructor for the switched accretion disk class.
    !!}
    implicit none
    type            (accretionDisksSwitched)                        :: self
    class           (accretionDisksClass   ), intent(in   ), target :: accretionDisksADAF_         , accretionDisksShakuraSunyaev_
    double precision                        , intent(in   )         :: accretionRateThinDiskMinimum, accretionRateThinDiskMaximum , &
         &                                                             accretionRateTransitionWidth
    logical                                 , intent(in   )         :: scaleADAFRadiativeEfficiency
    !![
    <constructorAssign variables="*accretionDisksADAF_, *accretionDisksShakuraSunyaev_, accretionRateThinDiskMinimum, accretionRateThinDiskMaximum, accretionRateTransitionWidth, scaleADAFRadiativeEfficiency"/>
    !!]

    self%accretionRateThinDiskMinimumExists=accretionRateThinDiskMinimum >      0.0d0
    self%accretionRateThinDiskMaximumExists=accretionRateThinDiskMaximum < huge(0.0d0)
    if (self%accretionRateThinDiskMinimumExists) self%accretionRateThinDiskMinimumLogarithmic=log(self%accretionRateThinDiskMinimum)
    if (self%accretionRateThinDiskMaximumExists) self%accretionRateThinDiskMaximumLogarithmic=log(self%accretionRateThinDiskMaximum)
    return
  end function switchedConstructorInternal

  subroutine switchedDestructor(self)
    !!{RST
    Destructor for the switched accretion disk class.
    !!}
    implicit none
    type(accretionDisksSwitched), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionDisksADAF_"           />
    <objectDestructor name="self%accretionDisksShakuraSunyaev_" />
    !!]
    return
  end subroutine switchedDestructor

  double precision function switchedEfficiencyRadiative(self,blackHole,accretionRateMass,accretionDiskType) result(efficiencyRadiative)
    !!{RST
    Return the radiative efficiency of a switched (ADAF/thin) accretion disk.
    !!}
    implicit none
    class           (accretionDisksSwitched          ), intent(inout)           :: self
    class           (nodeComponentBlackHole          ), intent(inout)           :: blackHole
    double precision                                  , intent(in   )           :: accretionRateMass
    type            (enumerationAccretionDiskTypeType), intent(in   ), optional :: accretionDiskType
    double precision                                                            :: fractionADAF               , efficiencyRadiativeADAF, &
         &                                                                         efficiencyRadiativeThinDisk
    !![
    <optionalArgument name="accretionDiskType" defaultsTo="accretionDiskTypeAny"/>
    !!]

    fractionADAF               =self                              %fractionADAF       (blackHole,accretionRateMass)
    efficiencyRadiativeThinDisk=self%accretionDisksShakuraSunyaev_%efficiencyRadiative(blackHole,accretionRateMass)
    efficiencyRadiativeADAF    =self%accretionDisksADAF_          %efficiencyRadiative(blackHole,accretionRateMass)
    if (self%scaleADAFRadiativeEfficiency                     )                                      &
         & efficiencyRadiativeADAF=+efficiencyRadiativeADAF                                          &
         &                         *self%efficiencyRadiativeScalingADAF(blackHole,accretionRateMass)
    efficiencyRadiative=+0.0d0
    if (accretionDiskType_ == accretionDiskTypeAny .or. accretionDiskType_ == accretionDiskTypeADAF) &
         & efficiencyRadiative=+                      efficiencyRadiative                            &
         &                     +        fractionADAF *efficiencyRadiativeADAF
    if (accretionDiskType_ == accretionDiskTypeAny .or. accretionDiskType_ == accretionDiskTypeThin) &
         & efficiencyRadiative=+                      efficiencyRadiative                            &
         &                     +(+1.0d0-fractionADAF)*efficiencyRadiativeThinDisk
    return
  end function switchedEfficiencyRadiative

  double precision function switchedPowerJet(self,blackHole,accretionRateMass)
    !!{RST
    Return the jet power of a switched (ADAF/thin) accretion disk.
    !!}
    implicit none
    class           (accretionDisksSwitched), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                                        :: fractionADAF

    fractionADAF    =+                      self                              %fractionADAF(blackHole,accretionRateMass)
    switchedPowerJet=+(+1.0d0-fractionADAF)*self%accretionDisksShakuraSunyaev_%powerJet    (blackHole,accretionRateMass) &
         &           +        fractionADAF *self%accretionDisksADAF_          %powerJet    (blackHole,accretionRateMass)
    return
  end function switchedPowerJet

  double precision function switchedRateSpinUp(self,blackHole,accretionRateMass)
    !!{RST
    Computes the spin up rate of the given ``blackHole`` due to accretion from a switched (ADAF/thin) accretion disk.
    !!}
    implicit none
    class           (accretionDisksSwitched), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                                        :: fractionADAF

    fractionADAF      =+                     self                              %fractionADAF(blackHole,accretionRateMass)
    switchedRateSpinUp=+(1.0d0-fractionADAF)*self%accretionDisksShakuraSunyaev_%rateSpinUp  (blackHole,accretionRateMass) &
         &             +       fractionADAF *self%accretionDisksADAF_          %rateSpinUp  (blackHole,accretionRateMass)
    return
  end function switchedRateSpinUp

  double precision function switchedFractionADAF(self,blackHole,accretionRateMass)
    !!{RST
    Decide which type of accretion disk to use.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Eddington_Accretion_Rate
    implicit none
    class           (accretionDisksSwitched), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                        , parameter     :: exponentialArgumentMaximum    =60.0d0
    double precision                                        :: accretionRateLogarithmic             , fractionADAF          , &
         &                                                     argument                             , accretionRateEddington, &
         &                                                     accretionRateMassDimensionless

    ! Get the Eddington accretion rate.
    accretionRateEddington=Black_Hole_Eddington_Accretion_Rate(blackHole)
    ! Check that a black hole is present.
    if (accretionRateEddington > 0.0d0 .and. accretionRateMass > 0.0d0) then
       ! Compute the accretion rate in Eddington units.
       accretionRateMassDimensionless=     +accretionRateMass              &
            &                              /accretionRateEddington
       accretionRateLogarithmic      =+log(                                &
            &                              +accretionRateMassDimensionless &
            &                             )
       ! Compute the ADAF fraction.
       fractionADAF=0.0d0
       if (self%accretionRateThinDiskMinimumExists) then
          argument    =min(+(accretionRateLogarithmic-self%accretionRateThinDiskMinimumLogarithmic)/self%accretionRateTransitionWidth,exponentialArgumentMaximum)
          fractionADAF=fractionADAF+1.0d0/(1.0d0+exp(argument))
       end if
       if (self%accretionRateThinDiskMaximumExists) then
          argument    =min(-(accretionRateLogarithmic-self%accretionRateThinDiskMaximumLogarithmic)/self%accretionRateTransitionWidth,exponentialArgumentMaximum)
          fractionADAF=fractionADAF+1.0d0/(1.0d0+exp(argument))
       end if
       switchedFractionADAF=fractionADAF
    else
       ! No black hole present: assume a thin disk.
       switchedFractionADAF=0.0d0
    end if
    return
  end function switchedFractionADAF

  double precision function switchedEfficiencyRadiativeScalingADAF(self,blackHole,accretionRateMass)
    !!{RST
    Determine the scaling of radiative efficiency of the ADAF component in a switched accretion disk.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_Eddington_Accretion_Rate
    implicit none
    class           (accretionDisksSwitched), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: accretionRateMass
    double precision                                        :: accretionRateEddington, accretionRateMassDimensionless

    ! Get the Eddington accretion rate.
    accretionRateEddington=Black_Hole_Eddington_Accretion_Rate(blackHole)
    ! Check that a black hole is present.
    if (accretionRateEddington > 0.0d0 .and. accretionRateMass > 0.0d0) then
       ! Compute the accretion rate in Eddington units.
       accretionRateMassDimensionless=+accretionRateMass      &
            &                         /accretionRateEddington
       ! If below the critical accretion rate for transition to a thin disk, reduce the radiative efficiency by a factor
       ! proportional to the accretion rate.
       if (self%accretionRateThinDiskMinimumExists .and. accretionRateMassDimensionless < self%accretionRateThinDiskMinimum) then
          switchedEfficiencyRadiativeScalingADAF=accretionRateMassDimensionless/self%accretionRateThinDiskMinimum
       else
          switchedEfficiencyRadiativeScalingADAF=1.0d0
       end if
    else
       ! No black hole present: return unit scaling.
       switchedEfficiencyRadiativeScalingADAF=1.0d0
    end if
    return
  end function switchedEfficiencyRadiativeScalingADAF
