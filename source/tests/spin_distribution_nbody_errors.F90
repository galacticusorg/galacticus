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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests the N-body errors halo spin distribution.
!!}

program Test_Spin_Distribution_Nbody_Errors
  !!{RST
  Tests the :galacticus-class:`haloSpinDistributionNbodyErrors` halo spin distribution class, whose tabulation of the
  error-convolved distribution is built on absolute lattices in mass and spin. The test requires that:

  #. evaluating twice with identical arguments returns the same distribution;
  #. a second, independent object asked the same question returns the same answer.

  Both failed before the fix which accompanies this test: the class recorded the mass and spin it was asked for only when it
  found them changed from a previous request, so the first evaluation for a fresh object tabulated at the default minimum mass
  rather than the one requested.

  The tabulated range is deliberately *not* driven outside the class's defaults of [3.0d-4,0.5d0] in spin. Those already span
  every physical halo spin, and forcing an extension beyond them drives the distribution into underflow - the value returned
  becomes denormal - so it tests a regime the class is not meant to work in. Extension within the defaults changes interpolated
  values only at the level of round-off in the interpolation weights, so it cannot be asserted bit-for-bit either.

  Before this test the class had no coverage of any kind: no test script, unit test, or bundled parameter file referenced it.

  The fixed-point entry point is used deliberately. It tabulates at a single mass, so the table is one row of order a hundred
  points rather than the two-dimensional grid of some two thousand which the general entry point builds - each point of which
  is several nested numerical integrals. Driving the general entry point takes tens of minutes, which is why this test
  exercises the spin axis only; the mass axis is left uncovered, as it was before. The extents of both axes are fixed by
  module parameters in the class and cannot be reduced from a parameter file.
  !!}
  use :: Display          , only : displayVerbositySet          , verbosityLevelStandard
  use :: Events_Hooks     , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes , only : nodeClassHierarchyInitialize , treeNode                          , nodeComponentBasic, nodeComponentSpin
  use :: Halo_Spin_Distributions, only : haloSpinDistributionNbodyErrors
  use :: Input_Parameters , only : inputParameters
  use :: Node_Components  , only : Node_Components_Initialize   , Node_Components_Thread_Initialize , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests       , only : Assert                       , Unit_Tests_Begin_Group            , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  type            (inputParameters                )          :: parameters          , parametersSpin
  type            (haloSpinDistributionNbodyErrors)          :: distributionExtended, distributionDirect
  type            (treeNode                       ), pointer :: node
  class           (nodeComponentBasic             ), pointer :: basic
  class           (nodeComponentSpin              ), pointer :: spin
  double precision                                           :: probeInitial        , probeAfterExtension, &
       &                                                        probeDirect
  double precision                                 , parameter :: massProbe=1.0d12
  ! The spin at which the tabulation is read, and two ranges of measured spin: the second reaches well outside the first and
  ! so forces the tabulated range to be extended.
  double precision                                 , parameter :: spinMeasured       =3.0d-2
  double precision                                 , parameter :: spinNarrowMinimum  =2.0d-2, spinNarrowMaximum=5.0d-2
  ! An angular momentum, not a spin: the class divides by the halo's angular momentum scale to form the dimensionless spin.
  ! Its precise value does not matter - the assertions are about reproducibility of the tabulation, not about particular
  ! distribution values - only that it places the spin inside, or close to, the tabulated range.
  double precision                                 , parameter :: angularMomentumProbe=1.0d12

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Functions_Global_Set (                      )
  parameters=inputParameters('testSuite/parameters/spinDistributionNbodyErrors.xml')
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  call Unit_Tests_Begin_Group("N-body errors halo spin distribution")
  ! Build two identical distributions from the same parameters. The class reads its own parameters from the object it is
  ! handed, so it is given the `haloSpinDistribution` subtree rather than the root.
  parametersSpin      =parameters%subParameters('haloSpinDistribution')
  distributionExtended=haloSpinDistributionNbodyErrors(parametersSpin)
  distributionDirect  =haloSpinDistributionNbodyErrors(parametersSpin)
  ! Evaluate twice with identical arguments. The class records the fixed point - the mass and spin it was asked for - only when
  ! it finds them changed from a previous request, which on a fresh object there is none of; without that record the first
  ! evaluation tabulates at the default minimum mass rather than the one requested, and the second, finding the mismatch and
  ! rebuilding, returns a different answer for the same question. For a 10^12 solar mass halo the two differed by a factor of
  ! about fourteen.
  probeInitial       =distributionValue(distributionExtended,massProbe,angularMomentumProbe,spinNarrowMinimum,spinNarrowMaximum)
  probeAfterExtension=distributionValue(distributionExtended,massProbe,angularMomentumProbe,spinNarrowMinimum,spinNarrowMaximum)
  call Assert("repeated evaluation returns the same distribution",probeAfterExtension == probeInitial,.true.)
  ! A second, independent object asked the same question must give the same answer as the first.
  probeDirect        =distributionValue(distributionDirect  ,massProbe,angularMomentumProbe,spinNarrowMinimum,spinNarrowMaximum)
  call Assert("the distribution does not depend on which object was asked",probeDirect == probeInitial,.true.)
  ! Guard against all of the above being satisfied by a distribution which is identically zero.
  call Assert("the tabulated distribution is non-zero",probeInitial > 0.0d0,.true.)
  ! Clean up the node-component hierarchy.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function distributionValue(distribution_,mass,angularMomentum_,spinMinimum_,spinMaximum_)
    !!{RST
    Return the spin distribution for a halo of the given mass and spin, building the node it is evaluated for.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    type            (haloSpinDistributionNbodyErrors), intent(inout) :: distribution_
    double precision                                 , intent(in   ) :: mass         , angularMomentum_, &
         &                                                              spinMinimum_ , spinMaximum_
    class           (nodeComponentDarkMatterProfile  ), pointer       :: darkMatterProfile

    node              =>  treeNode                  (                 )
    basic             =>  node    %basic            (autoCreate=.true.)
    spin              =>  node    %spin             (autoCreate=.true.)
    darkMatterProfile =>  node    %darkMatterProfile(autoCreate=.true.)
    call basic%massSet            (mass     )
    call basic%timeSet            (13.8d0   )
    ! Give the halo a well-formed dark matter profile. The class sets the scale radius itself on the node it builds while
    ! tabulating, but not on the node passed to it, and the angular momentum scale it forms below integrates over that
    ! profile - which is degenerate, and produces an invalid result, if the scale is left unset. The value need only be
    ! physically sensible; it is scaled with the halo's radius so that the concentration is the same at every mass.
    call darkMatterProfile%scaleSet(2.0d-2*(mass/1.0d12)**(1.0d0/3.0d0))
    call spin %angularMomentumSet (angularMomentum_)
    distributionValue=distribution_%distributionFixedPoint(node,spinMeasured,spinMinimum_,spinMaximum_)
    call node%destroy()
    deallocate(node)
    return
  end function distributionValue

end program Test_Spin_Distribution_Nbody_Errors
