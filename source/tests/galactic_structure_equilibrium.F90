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
Contains a program which tests the equilibrium galactic structure solver against independently computed reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Galactic_Structure_Equilibrium
  !!{RST
  Tests the equilibrium galactic structure solver against reference values computed independently by the
  ``galacticStructureEquilibrium.py`` script in the `galacticusDevTools
  &lt;https://github.com/galacticusorg/galacticusDevTools&gt;`_ repository.

  An exponential disk and a Hernquist spheroid, of given masses and angular momenta, are placed in an NFW halo which contracts
  adiabatically in response to them following :cite:t:`gnedin_response_2004`. The solver finds the radius at which each
  component's specific angular momentum matches that of a circular orbit in the combined potential,

  .. math::

     j = v(r) r, \qquad v^2(r) = \frac{\mathrm{G} M_\mathrm{DM}(&lt;r)}{r} + v^2_\mathrm{baryonic}(r),

  where :math:`M_\mathrm{DM}` is the *contracted* dark matter profile. The reference script solves the same system as a plain
  two-level root find, where the solver reaches the fixed point by iterating :math:`r \rightarrow \sqrt{j r / v}` with
  heuristics to break oscillations. Agreement is therefore a statement about the solution, not about the algorithm.

  The grid varies the halo concentration and the baryonic mass fraction, which set how centrally the dark matter sits and how
  strongly it contracts, together with the spin, which sets where the components land. It deliberately does *not* rest on
  varying the halo mass: with the concentration and the baryonic mass fractions fixed, and the angular momentum expressed in
  units of the halo's own, the problem is scale free in mass and every length scales as :math:`M^{1/3}`. Two halo masses are
  included at one fixed set of dimensionless parameters to exercise that scaling, and the remaining twelve models vary the
  quantities which change the solution.

  All physics classes are declared in the parameter file rather than constructed here. The dark matter profile component builds
  its ``darkMatterProfile`` through an object builder reading those parameters, so the contracted profile the node reports -
  which is the profile the solver equilibrates against - can only be set there.
  !!}
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Display                   , only : displayMessage                    , displayVerbositySet              , verbosityLevelStandard
  use :: Error                     , only : Error_Handler_Register
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galactic_Structure_Solvers, only : galacticStructureSolverClass
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize        , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                    nodeComponentDisk                 , nodeComponentSpheroid            , treeNode
  use :: Input_Parameters          , only : inputParameters
  use :: ISO_Varying_String        , only : varying_string                    , assignment(=)                    , var_str
  use :: Node_Components           , only : Node_Components_Initialize        , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                            , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  ! The reference models and their solved radii, emitted by `galacticStructureEquilibrium.py --fortran`.
  integer                                      , parameter :: countModels=14
  double precision, dimension(countModels)     , parameter :: massHalo               =[ 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+12, 1.000000000000000d+11, 1.000000000000000d+13]
  double precision, dimension(countModels)     , parameter :: concentration          =[ 5.000000000000000d+00, 5.000000000000000d+00, 5.000000000000000d+00, 5.000000000000000d+00, 1.000000000000000d+01, 1.000000000000000d+01, 1.000000000000000d+01, 1.000000000000000d+01, 2.000000000000000d+01, 2.000000000000000d+01, 2.000000000000000d+01, 2.000000000000000d+01, 1.000000000000000d+01, 1.000000000000000d+01]
  double precision, dimension(countModels)     , parameter :: massDisk               =[ 1.000000000000000d+10, 5.000000000000000d+10, 1.000000000000000d+10, 5.000000000000000d+10, 1.000000000000000d+10, 5.000000000000000d+10, 1.000000000000000d+10, 5.000000000000000d+10, 1.000000000000000d+10, 5.000000000000000d+10, 1.000000000000000d+10, 5.000000000000000d+10, 4.000000000000000d+09, 4.000000000000000d+11]
  double precision, dimension(countModels)     , parameter :: massSpheroid           =[ 2.500000000000000d+09, 1.250000000000000d+10, 2.500000000000000d+09, 1.250000000000000d+10, 2.500000000000000d+09, 1.250000000000000d+10, 2.500000000000000d+09, 1.250000000000000d+10, 2.500000000000000d+09, 1.250000000000000d+10, 2.500000000000000d+09, 1.250000000000000d+10, 1.000000000000000d+09, 1.000000000000000d+11]
  double precision, dimension(countModels)     , parameter :: angularMomentumDisk    =[ 1.034385050495032d+10, 5.171925252475159d+10, 2.585962626237580d+10, 1.292981313118790d+11, 1.034385050495032d+10, 5.171925252475159d+10, 2.585962626237580d+10, 1.292981313118790d+11, 1.034385050495032d+10, 5.171925252475159d+10, 2.585962626237580d+10, 1.292981313118790d+11, 8.914060142547512d+08, 1.920476040013492d+12]
  double precision, dimension(countModels)     , parameter :: angularMomentumSpheroid=[ 6.464906565593948d+08, 3.232453282796975d+09, 1.616226641398487d+09, 8.081133206992437d+09, 6.464906565593948d+08, 3.232453282796975d+09, 1.616226641398487d+09, 8.081133206992437d+09, 6.464906565593948d+08, 3.232453282796975d+09, 1.616226641398487d+09, 8.081133206992437d+09, 5.571287589092195d+07, 1.200297525008433d+11]
  double precision, dimension(countModels)     , parameter :: radiusDiskReference    =[ 5.866076444012095d-03, 2.376271619863644d-03, 1.463563271584771d-02, 9.391846908615399d-03, 4.808213137308325d-03, 2.204504141405828d-03, 1.155191947632160d-02, 8.064051202100745d-03, 3.784058239341330d-03, 1.982335180512968d-03, 8.894749682392230d-03, 6.689961994158126d-03, 1.187397875626599d-03, 5.511412720563722d-03]
  double precision, dimension(countModels)     , parameter :: radiusSpheroidReference=[ 2.010649877337007d-03, 7.505791935858934d-04, 5.235467947185564d-03, 3.125177265434999d-03, 1.664179501597806d-03, 7.012607684294530d-04, 4.075153800186249d-03, 2.693429437297989d-03, 1.315663495473915d-03, 6.361119265877454d-04, 3.060245598370788d-03, 2.226357100342361d-03, 3.812302710950055d-04, 1.769514169349703d-03]
  ! The assertion tolerance. The disk has Galacticus' default scale height (0.137 scale lengths), so its rotation curve is that
  ! of a thick disk, which Galacticus tabulates at 100 points per decade and interpolates with a cubic spline in ln V_c^2 vs
  ! ln R. `galacticStructureEquilibrium.py --tolerance` re-solves every model with that tabulation in place of its own (which
  ! is accurate to 1.4e-9) and reports a largest fractional shift of 3e-11, so the tabulation no longer limits the agreement.
  ! What does is the contraction's `toleranceRelative` of 1.0e-6 (see the note in the parameter file): the radii agree with
  ! the reference values to 5.5e-6. The tolerance below is set well above that, rather than tuned to it.
  double precision                                , parameter :: toleranceRelative       =1.0d-4
  class           (cosmologyFunctionsClass       ), pointer   :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass      ), pointer   :: darkMatterHaloScale_
  class           (galacticStructureSolverClass  ), pointer   :: galacticStructureSolver_
  type            (treeNode                      ), pointer   :: node
  class           (nodeComponentBasic            ), pointer   :: basic
  class           (nodeComponentDarkMatterProfile), pointer   :: dmProfile
  class           (nodeComponentDisk             ), pointer   :: disk
  class           (nodeComponentSpheroid         ), pointer   :: spheroid
  type            (inputParameters               )            :: parameters
  character       (len=128                       )            :: message
  integer                                                     :: iModel
  double precision                                            :: radiusVirial                  , time              , &
       &                                                         differenceDisk                , differenceSpheroid, &
       &                                                         differenceMaximum

  ! Set verbosity level, and register the error handler so that failures in the numerical libraries are reported rather than
  ! silently ignored.
  call displayVerbositySet              (verbosityLevelStandard)
  call Error_Handler_Register           (                      )
  ! Initialize event hooks, global functions and node components.
  parameters=inputParameters(var_str('testSuite/parameters/galacticStructureEquilibrium.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="cosmologyFunctions"      name="cosmologyFunctions_"      source="parameters"/>
  <objectBuilder class="darkMatterHaloScale"     name="darkMatterHaloScale_"     source="parameters"/>
  <objectBuilder class="galacticStructureSolver" name="galacticStructureSolver_" source="parameters"/>
  !!]
  call Unit_Tests_Begin_Group("Equilibrium galactic structure solver")
  time             =cosmologyFunctions_%cosmicTime(1.0d0)
  differenceMaximum=0.0d0
  do iModel=1,countModels
     ! Build a fresh node for every model. A node memoizes the mass distributions built from its components, so reusing one and
     ! resetting its properties would leave the previous model's mass distribution in place.
     node      => treeNode                  (                 )
     basic     => node    %basic            (autoCreate=.true.)
     dmProfile => node    %darkMatterProfile(autoCreate=.true.)
     disk      => node    %disk             (autoCreate=.true.)
     spheroid  => node    %spheroid         (autoCreate=.true.)
     call basic    %timeSet            (time                            )
     call basic    %timeLastIsolatedSet(time                            )
     call basic    %massSet            (massHalo           (iModel)     )
     ! Set the halo scale radius from the concentration. The virial radius follows from the fixed density contrast, so this is
     ! the only structural property the halo needs.
     radiusVirial=darkMatterHaloScale_%radiusVirial(node)
     call dmProfile%scaleSet           (radiusVirial/concentration(iModel))
     ! Set the baryonic masses and angular momenta. The radii are left unset, so that the solver makes its own initial estimate
     ! rather than being handed one.
     call disk     %massStellarSet     (massDisk               (iModel) )
     call disk     %angularMomentumSet (angularMomentumDisk    (iModel) )
     call spheroid %massStellarSet     (massSpheroid           (iModel) )
     call spheroid %angularMomentumSet (angularMomentumSpheroid(iModel) )
     ! Solve for the equilibrium structure.
     call galacticStructureSolver_%solve(node)
     ! Compare against the reference values. The fractional difference of each comparison is reported, not just whether it
     ! passed, so that a reader of the log can see how much of the tolerance is actually being used - a test which passes only
     ! because its tolerance is loose looks identical to one which passes sharply, unless the margin is shown.
     differenceDisk    =abs(disk    %radius()/radiusDiskReference    (iModel)-1.0d0)
     differenceSpheroid=abs(spheroid%radius()/radiusSpheroidReference(iModel)-1.0d0)
     differenceMaximum =max(differenceMaximum,differenceDisk,differenceSpheroid)
     write (message,'(a,i0,a,e12.5,a,e12.5)') 'model ',iModel,': fractional difference, disk ',differenceDisk,', spheroid ',differenceSpheroid
     call displayMessage(trim(message))
     write (message,'(a,i0)') 'disk radius, model ',iModel
     call Assert(trim(message),disk    %radius(),radiusDiskReference    (iModel),relTol=toleranceRelative)
     write (message,'(a,i0)') 'spheroid radius, model ',iModel
     call Assert(trim(message),spheroid%radius(),radiusSpheroidReference(iModel),relTol=toleranceRelative)
     call node%destroy()
     deallocate(node)
  end do
  write (message,'(a,e12.5)') 'largest fractional difference over all models: ',differenceMaximum
  call displayMessage(trim(message))
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  !![
  <objectDestructor name="cosmologyFunctions_"     />
  <objectDestructor name="darkMatterHaloScale_"    />
  <objectDestructor name="galacticStructureSolver_"/>
  !!]
end program Test_Galactic_Structure_Equilibrium
