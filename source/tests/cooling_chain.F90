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
Contains a program which tests the chain of calculations determining the rate at which gas cools out of a hot halo, against
independently computed reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Cooling_Chain
  !!{RST
  Tests the chain of calculations which determines the rate at which gas cools out of a hot halo, against reference values
  computed independently by the ``coolingChain.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository. The chain tested is

  .. math::

     M_\mathrm{halo}, z \rightarrow \bar{\rho} \rightarrow r_\mathrm{vir}, V_\mathrm{vir}, T_\mathrm{vir}, \tau_\mathrm{dyn}
     \rightarrow \rho_0 \rightarrow t_\mathrm{cool} \rightarrow r_\mathrm{cool} \rightarrow \dot{r}_\mathrm{cool} \rightarrow
     \dot{M}_\mathrm{cool}.

  The reference implementation is written from the definitions of the physics rather than transcribed from this code. The two
  tabulated Cloudy files are read directly by both---they are data, not code---but the interpolation scheme applied to them is
  reimplemented independently, since that scheme is itself part of what is being checked.

  The test halo is given a hot atmosphere containing one tenth of its mass, with an outer radius equal to the virial radius and
  primordial composition. Cases are chosen to exercise all three branches of
  :galacticus-class:`coolingRateWhiteFrenk1991`: the *saturated* branch, in which the cooling radius has reached the outer
  radius and infall is limited to the halo dynamical timescale; the *interior* branch, in which the cooling radius lies strictly
  inside the halo and the rate is :math:`4 \pi r^2 \rho(r) \dot{r}`; and the branch in which the cooling time exceeds the time
  available everywhere, so that nothing cools. The remaining branch, in which cooling is suppressed above
  ``[velocityCutOff]``, would require a virial velocity above :math:`10^4` km/s---a halo of some :math:`10^{17}
  \mathrm{M}_\odot`---and is better exercised by lowering that parameter.

  Note that the virial density contrast is pinned to a fixed value of 200 relative to the *mean* matter density (the default
  being spherical collapse in a cosmological constant cosmology), so that the virial radius is analytic and this test does not
  also depend on the spherical collapse solver, which has its own tests. Note also that
  :galacticus-class:`cosmologyFunctionsMatterLambda` contains no radiation term in its expansion rate, so radiation is excluded
  from the reference calculation also.

  Tolerances are:

  * halo scales (:math:`r_\mathrm{vir}`, :math:`V_\mathrm{vir}`, :math:`T_\mathrm{vir}`, :math:`\tau_\mathrm{dyn}`):
    :math:`5\times 10^{-4}` relative. Two effects contribute, and the larger is not the one that might be expected. The smaller
    is the difference between the GSL constants used here and the CODATA 2018/IAU 2015 values used by the reference, which
    shift :math:`\rho_\mathrm{crit}` by :math:`4.5\times 10^{-4}` and hence :math:`r_\mathrm{vir}` by
    :math:`2\times 10^{-5}`. The larger is that the mean density of a halo is not evaluated directly: for a density contrast
    which is not mass dependent it is interpolated linearly in :math:`\ln t` from a table built at 100 points per decade, and
    the expansion factor at a given cosmic time is itself interpolated from a further tabulation, built at 300 points per
    decade. The resulting error therefore
    depends on where the epoch happens to fall within an interpolation interval, rather than varying smoothly with redshift:
    for this cosmology :math:`z=0` falls at a fraction :math:`0.99` of the way through its interval, where the error in
    :math:`\bar{\rho}` is only :math:`1.6\times 10^{-5}`, while :math:`z=2` falls at :math:`0.51`---almost exactly the worst
    case---where it is :math:`2.7\times 10^{-4}`. Since :math:`r_\mathrm{vir} \propto \bar{\rho}^{-1/3}`,
    :math:`T_\mathrm{vir} \propto r_\mathrm{vir}^{-1}` and :math:`\tau_\mathrm{dyn} \propto r_\mathrm{vir}^{3/2}`, that
    propagates to :math:`9\times 10^{-5}`, :math:`9\times 10^{-5}` and :math:`1.4\times 10^{-4}` respectively. The differences
    actually seen are :math:`1.1\times 10^{-4}` and :math:`1.7\times 10^{-4}` at :math:`z=2`, some 25% larger, the remainder
    being contributed by the expansion factor tabulation which this estimate does not include;
  * the hot gas density at the cooling radius: :math:`10^{-3}` relative. In addition to the error in
    :math:`\rho_0 \propto \bar{\rho}` above, this inherits the error in the cooling radius below, amplified by the logarithmic
    slope of the :math:`\beta`-profile;
  * the cooling radius, its growth rate, and the mass cooling rate: :math:`10^{-3}` relative. The cooling radius is located by a
    root find with a relative tolerance of :math:`10^{-6}`, but that error is amplified by the steep density gradient of the
    :math:`\beta`-profile; and the cooling function and electron density are obtained by bilinear interpolation in tables
    tabulated every :math:`0.025` dex in temperature, which the two implementations index independently.
  !!}
  use            :: Abundances_Structure       , only : zeroAbundances
  use            :: Cooling_Radii              , only : coolingRadius                      , coolingRadiusClass
  use            :: Cooling_Rates              , only : coolingRate                        , coolingRateClass
  use            :: Coordinates                , only : coordinateSpherical                , assignment(=)
  use            :: Cosmology_Functions        , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use            :: Dark_Matter_Halo_Scales    , only : darkMatterHaloScale                , darkMatterHaloScaleClass
  use            :: Display                    , only : displayVerbositySet                , verbosityLevelStandard
  use            :: Error                      , only : Error_Handler_Register
  use            :: Events_Hooks               , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities , only : Functions_Global_Set
  use            :: Galactic_Structure_Options , only : componentTypeHotHalo               , massTypeGaseous
  use            :: Galacticus_Nodes           , only : nodeClassHierarchyInitialize       , nodeComponentBasic               , &
  &                                                     nodeComponentHotHalo               , treeNode
  use            :: IO_HDF5                    , only : ioHDF5AccessInitialize
  use            :: ISO_Varying_String         , only : assignment(=)                      , varying_string
  use            :: Input_Parameters           , only : inputParameters
  use            :: Mass_Distributions         , only : massDistributionClass
  use            :: Node_Components            , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, &
  &                                                     Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Unit_Tests                 , only : Assert                             , Unit_Tests_Begin_Group           , &
  &                                                     Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_
  class           (coolingRadiusClass      ), pointer                   :: coolingRadius_
  class           (coolingRateClass        ), pointer                   :: coolingRate_
  class           (massDistributionClass   ), pointer                   :: massDistribution_
  type            (treeNode                ), pointer                   :: node
  class           (nodeComponentBasic      ), pointer                   :: basic
  class           (nodeComponentHotHalo    ), pointer                   :: hotHalo
  type            (coordinateSpherical     )                            :: coordinates
  type            (varying_string          )                            :: parameterFile
  type            (inputParameters         )                            :: parameters
  ! The fraction of the halo mass placed in the hot atmosphere.
  double precision                          , parameter                 :: massHotFraction             =1.0000000000d-01
  ! Test cases: halo masses and redshifts chosen to span the three branches of the cooling rate.
  integer                                   , parameter                 :: countCases                  =8
  double precision                          , dimension(countCases)     :: massHalo                    =[3.1622776602d+10,1.0000000000d+11,1.0000000000d+11,1.7782794100d+11,1.7782794100d+11,3.1622776602d+11,3.1622776602d+11,1.0000000000d+12]
  double precision                          , dimension(countCases)     :: redshift                    =[0.0000000000d+00,0.0000000000d+00,1.0000000000d+00,0.0000000000d+00,2.0000000000d+00,0.0000000000d+00,1.0000000000d+00,0.0000000000d+00]
  ! Reference values, computed by coolingChain.py.
  double precision                          , dimension(countCases)     :: radiusVirialReference       =[9.8329092600d-02,1.4432737000d-01,7.2163685000d-02,1.7485660100d-01,5.8285533600d-02,2.1184360800d-01,1.0592180400d-01,3.1094389300d-01]
  double precision                          , dimension(countCases)     :: velocityVirialReference     =[3.7191705600d+01,5.4589958200d+01,7.7201859300d+01,6.6137244300d+01,1.1455306700d+02,8.0127100700d+01,1.1331683200d+02,1.1761050000d+02]
  double precision                          , dimension(countCases)     :: temperatureVirialReference  =[4.9587896600d+04,1.0683388500d+05,2.1366776900d+05,1.5681069800d+05,4.7043209300d+05,2.3016662700d+05,4.6033325400d+05,4.9587896600d+05]
  double precision                          , dimension(countCases)     :: timeDynamicalReference      =[2.5851838400d+00,2.5851838400d+00,9.1400051100d-01,2.5851838400d+00,4.9751886200d-01,2.5851838400d+00,9.1400051100d-01,2.5851838400d+00]
  double precision                          , dimension(countCases)     :: densityCoolingReference     =[3.9409158400d+11,3.9409158400d+11,7.0356612300d+12,1.1054644800d+12,6.3371035800d+13,2.9682564700d+12,3.3275707500d+13,0.0000000000d+00]
  double precision                          , dimension(countCases)     :: radiusCoolingReference      =[9.8329092600d-02,1.4432737000d-01,4.5551075700d-02,9.5545747200d-02,1.7776532100d-02,4.9554152500d-02,1.2203044200d-02,0.0000000000d+00]
  double precision                          , dimension(countCases)     :: radiusCoolingGrowthReference=[0.0000000000d+00,0.0000000000d+00,3.0547175400d-02,2.4049726400d-02,3.5150506100d-02,2.5348454400d-02,5.1941258600d-02,0.0000000000d+00]
  double precision                          , dimension(countCases)     :: rateCoolingReference        =[1.2232312500d+09,3.8681968600d+09,5.6038063300d+09,3.0499134800d+09,8.8455804400d+09,2.3217896900d+09,3.2343424200d+09,0.0000000000d+00]
  double precision                          , dimension(countCases)     :: radiusVirialComputed                                                                                                                                                 , &
       &                                                                   velocityVirialComputed                                                                                                                                               , &
       &                                                                   temperatureVirialComputed                                                                                                                                            , &
       &                                                                   timeDynamicalComputed                                                                                                                                                , &
       &                                                                   densityCoolingComputed                                                                                                                                               , &
       &                                                                   radiusCoolingComputed                                                                                                                                                , &
       &                                                                   radiusCoolingGrowthComputed                                                                                                                                          , &
       &                                                                   rateCoolingComputed
  double precision                                                      :: time
  integer                                                               :: i

  ! Establish error handlers, so that errors reported by the GSL are trapped rather than aborting.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize the lock used to serialize access to the HDF5 library. The cooling function and chemical state classes read
  ! tabulated data from HDF5 files, and would fail at their first access without this.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cooling chain")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/coolingChainComponents.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
  <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
  <objectBuilder class="coolingRadius"       name="coolingRadius_"       source="parameters"/>
  <objectBuilder class="coolingRate"         name="coolingRate_"         source="parameters"/>
  !!]

  ! Build a halo with a hot atmosphere at each test mass and redshift, and evaluate the chain.
  do i=1,countCases
     time    =  cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift(i)))
     node    => treeNode           (                 )
     basic   => node               %basic  (autoCreate=.true.)
     hotHalo => node               %hotHalo(autoCreate=.true.)
     call basic%massSet            (massHalo(i))
     call basic%timeSet            (time       )
     call basic%timeLastIsolatedSet(time       )
     ! Halo scales.
     radiusVirialComputed     (i)=darkMatterHaloScale_%radiusVirial      (node)
     velocityVirialComputed   (i)=darkMatterHaloScale_%velocityVirial    (node)
     temperatureVirialComputed(i)=darkMatterHaloScale_%temperatureVirial (node)
     timeDynamicalComputed    (i)=darkMatterHaloScale_%timescaleDynamical(node)
     ! Populate the hot atmosphere. The outer radius getter clamps the stored value to lie between one tenth of the virial
     ! radius and the virial radius, so storing the virial radius returns it unchanged.
     call hotHalo%massSet       (massHotFraction*massHalo(i))
     call hotHalo%outerRadiusSet(radiusVirialComputed    (i))
     call hotHalo%abundancesSet (zeroAbundances             )
     ! The cooling chain.
     radiusCoolingComputed      (i)=coolingRadius_%radius          (node)
     radiusCoolingGrowthComputed(i)=coolingRadius_%radiusGrowthRate(node)
     rateCoolingComputed        (i)=coolingRate_  %rate            (node)
     ! The hot gas density at the cooling radius. This exercises the beta-profile built by the hot halo component itself, rather
     ! than any object constructed by this test.
     if (radiusCoolingComputed(i) > 0.0d0) then
        coordinates                =  [radiusCoolingComputed(i),0.0d0,0.0d0]
        massDistribution_          => node             %massDistribution(componentTypeHotHalo,massTypeGaseous)
        densityCoolingComputed  (i)=  massDistribution_%density         (coordinates                         )
        !![
        <objectDestructor name="massDistribution_"/>
        !!]
     else
        densityCoolingComputed  (i)=  0.0d0
     end if
     call node%destroy()
     deallocate(node)
  end do

  ! Halo scales.
  call Unit_Tests_Begin_Group("Halo scales")
  call Assert("virial radius"                ,radiusVirialComputed       ,radiusVirialReference       ,relTol=5.0d-4)
  call Assert("virial velocity"              ,velocityVirialComputed     ,velocityVirialReference     ,relTol=5.0d-4)
  call Assert("virial temperature"           ,temperatureVirialComputed  ,temperatureVirialReference  ,relTol=5.0d-4)
  call Assert("dynamical time"               ,timeDynamicalComputed      ,timeDynamicalReference      ,relTol=5.0d-4)
  call Unit_Tests_End_Group()
  ! The hot atmosphere.
  call Unit_Tests_Begin_Group("Hot atmosphere")
  call Assert("density at the cooling radius",densityCoolingComputed     ,densityCoolingReference     ,relTol=1.0d-3)
  call Unit_Tests_End_Group()
  ! The cooling chain.
  call Unit_Tests_Begin_Group("Cooling radius and rate")
  call Assert("cooling radius"               ,radiusCoolingComputed      ,radiusCoolingReference      ,relTol=1.0d-3)
  call Assert("cooling radius growth rate"   ,radiusCoolingGrowthComputed,radiusCoolingGrowthReference,relTol=1.0d-3)
  call Assert("mass cooling rate"            ,rateCoolingComputed        ,rateCoolingReference        ,relTol=1.0d-3)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  !![
  <objectDestructor name="cosmologyFunctions_" />
  <objectDestructor name="darkMatterHaloScale_"/>
  <objectDestructor name="coolingRadius_"      />
  <objectDestructor name="coolingRate_"        />
  !!]
end program Test_Cooling_Chain
