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
  Implements a transfer function accelerator class which tabulates a transfer function for rapid interpolation.
  !!}

  use :: Numerical_Ranges, only : rangeLattice
  use :: Tables          , only : table1DLogarithmicLinear

  !![
  <transferFunction name="transferFunctionAccelerator" docformat="rst">
   <description>
   A transfer function class which accelerates calculations of another transfer function class by pre-tabulating the transfer function over a grid of wavenumbers and then using rapid interpolation for subsequent evaluations. The density of the tabulation grid in wavenumber is set by ``[wavenumberPerDecade]``.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionAccelerator
     !!{RST
     A transfer function class which accelerates calculations of another transfer function class by tabulation for rapid interpolation.
     !!}
     private
     type            (table1DLogarithmicLinear)          :: transferTable
     class           (transferFunctionClass   ), pointer :: transferFunction_       => null()
     ! Lattice to which the tabulated wavenumbers are pinned. This is the source of truth for the extent of the tabulation: the
     ! bounds below are derived from it, and are retained because they are what the test for a sufficient tabulation reads.
     type            (rangeLattice            )          :: latticeWavenumber
     double precision                                    :: wavenumberMinimum                  , wavenumberMaximum
     integer                                             :: tablePointsPerDecade
     logical                                             :: tableInitialized        =  .false.
   contains
     final     ::                          acceleratorDestructor
     procedure :: value                 => acceleratorValue
     procedure :: logarithmicDerivative => acceleratorLogarithmicDerivative
     procedure :: halfModeMass          => acceleratorHalfModeMass
     procedure :: quarterModeMass       => acceleratorQuarterModeMass
     procedure :: fractionModeMass      => acceleratorFractionModeMass
     procedure :: epochTime             => acceleratorEpochTime
  end type transferFunctionAccelerator

  interface transferFunctionAccelerator
     !!{RST
     Constructors for the accelerator transfer function class.
     !!}
     module procedure acceleratorConstructorParameters
     module procedure acceleratorConstructorInternal
  end interface transferFunctionAccelerator

  ! Seed range of wavenumbers to tabulate. Any two tabulations therefore contain this range, and so always overlap - which is
  ! what makes their shared lattice points reusable no matter which wavenumbers each was asked for. Note that these are exact
  ! powers of ten, so on a lattice anchored to whole decades they are themselves anchor points: a tabulation which is never
  ! asked for a wavenumber outside them spans exactly this range, as it did before being pinned.
  double precision, parameter :: wavenumberTableSeedMinimum=1.0d-6, wavenumberTableSeedMaximum=1.0d+6

contains

  function acceleratorConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the accelerator transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    type  (transferFunctionAccelerator)                :: self
    type  (inputParameters            ), intent(inout) :: parameters
    class (transferFunctionClass      ), pointer       :: transferFunction_
    class (cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    integer                                            :: tablePointsPerDecade

    !![
    <inputParameter docformat="rst">
      <name>tablePointsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>
      The number of points per decade of wavenumber at which to tabulate the transfer function.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunction_"    source="parameters"/>
    !!]
    self=transferFunctionAccelerator(transferFunction_,cosmologyParameters_,tablePointsPerDecade)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="transferFunction_"/>
    !!]
    return
  end function acceleratorConstructorParameters

  function acceleratorConstructorInternal(transferFunction_,cosmologyParameters_,tablePointsPerDecade) result(self)
    !!{RST
    Internal constructor for the accelerator transfer function class.
    !!}
    implicit none
    type   (transferFunctionAccelerator)                        :: self
    class  (cosmologyParametersClass   ), intent(in   ), target :: cosmologyParameters_
    class  (transferFunctionClass      ), intent(in   ), target :: transferFunction_
    integer                             , intent(in   )         :: tablePointsPerDecade
    !![
    <constructorAssign variables="*cosmologyParameters_, *transferFunction_, tablePointsPerDecade"/>
    !!]

    self%tableInitialized            =.false.
    self%wavenumberMinimum           =wavenumberTableSeedMinimum
    self%wavenumberMaximum           =wavenumberTableSeedMaximum
    return
  end function acceleratorConstructorInternal

  subroutine acceleratorDestructor(self)
    !!{RST
    Destructor for the accelerator transfer function class.
    !!}
    implicit none
    type(transferFunctionAccelerator), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%transferFunction_"   />
    !!]
    if (self%tableInitialized) call self%transferTable%destroy()
    return
  end subroutine acceleratorDestructor

  double precision function acceleratorValue(self,wavenumber)
    !!{RST
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber

    call acceleratorTabulate(self,wavenumber)
    acceleratorValue=exp(self%transferTable%interpolate(wavenumber))
    return
  end function acceleratorValue

  double precision function acceleratorLogarithmicDerivative(self,wavenumber)
    !!{RST
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber

    call acceleratorTabulate(self,wavenumber)
    ! The tabulation holds the logarithm of the transfer function against wavenumber, so its gradient is d(ln T)/dk; the
    ! logarithmic derivative sought here is with respect to the logarithm of the wavenumber.
    acceleratorLogarithmicDerivative=+self%transferTable%interpolateGradient(wavenumber) &
         &                           *                                       wavenumber
    return
  end function acceleratorLogarithmicDerivative

  double precision function acceleratorEpochTime(self)
    !!{RST
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionAccelerator), intent(inout) :: self

    acceleratorEpochTime=self%transferFunction_%epochTime()
    return
  end function acceleratorEpochTime

  double precision function acceleratorHalfModeMass(self,status)
    !!{RST
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative to a :term:`CDM` transfer function
    !!}
    implicit none
    class  (transferFunctionAccelerator), intent(inout), target   :: self
    integer                             , intent(  out), optional :: status

    acceleratorHalfModeMass=self%transferFunction_%halfModeMass(status)
    return
  end function acceleratorHalfModeMass

  double precision function acceleratorQuarterModeMass(self,status)
    !!{RST
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative to a :term:`CDM` transfer function
    !!}
    implicit none
    class  (transferFunctionAccelerator), intent(inout), target   :: self
    integer                             , intent(  out), optional :: status

    acceleratorQuarterModeMass=self%transferFunction_%quarterModeMass(status)
    return
  end function acceleratorQuarterModeMass

  double precision function acceleratorFractionModeMass(self,fraction,status)
    !!{RST
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative to a :term:`CDM` transfer function
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout), target   :: self
    double precision                             , intent(in   )           :: fraction
    integer                                      , intent(  out), optional :: status

    acceleratorFractionModeMass=self%transferFunction_%fractionModeMass(fraction,status)
    return
  end function acceleratorFractionModeMass

  subroutine acceleratorTabulate(self,wavenumber)
    !!{RST
    Tabulate the transfer function for rapid interpolation.
    !!}
    use :: Numerical_Ranges, only : Range_Pinned, gridSchemePerDecade
    implicit none
    class           (transferFunctionAccelerator), intent(inout)                            :: self
    double precision                             , intent(in   )                            :: wavenumber
    logical                                                                                 :: makeTable
    integer                                                                                 :: i
    type            (rangeLattice                )                                          :: latticeWavenumber
    logical                                                     , allocatable, dimension(:) :: isComputed

    makeTable=.not.self%tableInitialized
    if (.not.makeTable)                                        &
         & makeTable= wavenumber < self%wavenumberMinimum      &
         &           .or.                                      &
         &            wavenumber > self%wavenumberMaximum
    if (makeTable) then
       ! Find the range of wavenumbers to tabulate, pinning it to an absolute lattice so that the wavenumbers evaluated - and
       ! therefore every value interpolated between them - depend only on which lattice points are spanned, and not on the
       ! sequence of wavenumbers which happened to be requested. The request is passed as the target and the range already
       ! tabulated is unioned in through `latticeCurrent`; folding the latter into the target instead - as the `min`/`max`
       ! against the current bounds formerly did - would apply the safety margin to an already margined bound and so ratchet the
       ! range outward on every retabulation. The margin is a factor of e, which is the margin of one unit in the natural
       ! logarithm of the wavenumber that was applied before.
       latticeWavenumber=Range_Pinned(                                                                 &
            &                                        [wavenumber]                                    , &
            &                                         self%tablePointsPerDecade                      , &
            &                                         gridSchemePerDecade                            , &
            &                         marginFactor  = exp(1.0d0)                                     , &
            &                         rangeCurrent  =[wavenumberTableSeedMinimum,wavenumberTableSeedMaximum], &
            &                         latticeCurrent= self%latticeWavenumber                           &
            &                        )
       ! Extend the tabulation onto the new lattice, preserving the values already computed - the transfer function being
       ! tabulated is expensive to evaluate, which is the entire purpose of this class. The abscissae come from the lattice, so
       ! a given wavenumber is bit-identical between one tabulation and another however many each spans.
       call self%transferTable%extend(latticeWavenumber,isComputed)
       do i=1,latticeWavenumber%count
          if (isComputed(i)) cycle
          call self%transferTable%populate(log(self%transferFunction_%value(self%transferTable%x(i))),i)
       end do
       self%latticeWavenumber=latticeWavenumber
       self%wavenumberMinimum=latticeWavenumber%minimum()
       self%wavenumberMaximum=latticeWavenumber%maximum()
       self%tableInitialized =.true.
    end if
    return
  end subroutine acceleratorTabulate

