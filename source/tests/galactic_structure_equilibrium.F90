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
  double precision, dimension(countModels)     , parameter :: radiusDiskReference    =[ 5.722708324603595d-03, 2.249875686443800d-03, 1.448789024853587d-02, 9.065523769848861d-03, 4.715863781207704d-03, 2.097389645936857d-03, 1.146094768275026d-02, 7.827402922198780d-03, 3.729079372572252d-03, 1.897719856957710d-03, 8.840832145499192d-03, 6.529352026972797d-03, 1.134187960632437d-03, 5.264434173289563d-03]
  double precision, dimension(countModels)     , parameter :: radiusSpheroidReference=[ 1.950213735459256d-03, 7.112390750122593d-04, 5.164556304050153d-03, 3.005577553793464d-03, 1.623010828022040d-03, 6.667408704051618d-04, 4.032736624498625d-03, 2.604216157573111d-03, 1.290105574755845d-03, 6.075901791636204d-04, 3.036686284409648d-03, 2.165248795921700d-03, 3.634672526447422d-04, 1.687065541259741d-03]
  ! The assertion tolerance. `galacticStructureEquilibrium.py --tolerance` re-solves every model with the Bessel factor of the
  ! disk rotation curve replaced by the tabulated-and-interpolated form Galacticus uses - 100 points per decade in a
  ! `table1DLogarithmicLinear`, interpolated linearly - and reports a largest fractional shift of 1.7e-5. That tabulation is
  ! the largest known difference between the two calculations, so the tolerance below is set a few times above it rather than
  ! tuned to whatever the comparison happens to give. The two solver tolerances which would otherwise dominate it are tightened
  ! in the parameter file; see the note there.
  double precision                              , parameter :: toleranceRelative       =1.0d-4
  class           (cosmologyFunctionsClass     ), pointer   :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass    ), pointer   :: darkMatterHaloScale_
  class           (galacticStructureSolverClass), pointer   :: galacticStructureSolver_
  type            (treeNode                    ), pointer   :: node
  class           (nodeComponentBasic          ), pointer   :: basic
  class           (nodeComponentDarkMatterProfile), pointer :: dmProfile
  class           (nodeComponentDisk           ), pointer   :: disk
  class           (nodeComponentSpheroid       ), pointer   :: spheroid
  type            (inputParameters             )            :: parameters
  character       (len=128                     )            :: message
  integer                                                   :: iModel
  double precision                                          :: radiusVirial                  , time              , &
       &                                                       differenceDisk                , differenceSpheroid, &
       &                                                       differenceMaximum

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
