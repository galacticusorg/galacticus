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
Contains a program which tests the satellite tidal heating rate against independent reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Satellite_Tidal_Heating_Rate
  !!{RST
  Tests the :cite:t:`gnedin_tidal_1999` satellite tidal heating rate against reference values computed independently by the
  ``satelliteTidalHeatingRate.py`` script in the `galacticusDevTools
  &lt;https://github.com/galacticusorg/galacticusDevTools&gt;`_ repository.

  The rate is compared pointwise, with the path-integrated tidal tensor set directly, rather than through an integrated orbit.
  The accumulation of that path integral along an orbit is a separate step, and a comparison of the accumulated heating alone
  could not distinguish a wrong rate from a correct rate wrongly accumulated.

  The satellite is placed off the coordinate axes, and every path integral has non-zero off-diagonal elements, so that both the
  orientation of the host's tidal tensor and the double contraction of the two tensors are exercised; on an axis the tidal
  tensor is diagonal and would test neither. One configuration on the x-axis is kept as the simplest case. The radius, speed,
  satellite concentration and bound mass move the adiabatic correction :math:`[1+(\omega\tau)^2]^{-\gamma}` from near unity to
  :math:`\sim 10^{-5}`, so that its exponent is tested and not merely its presence. Bound masses below the basic mass exercise
  the minimum taken in the satellite's half mass, a tiny bound mass exercises the fallback to the virial orbital frequency, and a
  path integral anti-aligned with the tidal field exercises the clamping of a negative rate to zero.

  The satellite's basic mass is not varied: with a fixed virial density contrast the internal orbital frequency, and so the rate,
  is independent of it, and such configurations would repeat another exactly.
  !!}
  use :: Display                   , only : displayMessage                , displayVerbositySet                 , verbosityLevelStandard
  use :: Error                     , only : Error_Handler_Register
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize    , nodeClassHierarchyInitialize        , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                    nodeComponentSatellite        , nodeComponentSpheroid               , treeNode
  use :: Input_Parameters          , only : inputParameters
  use :: ISO_Varying_String        , only : var_str
  use :: Node_Components           , only : Node_Components_Initialize    , Node_Components_Thread_Initialize   , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Tidal_Heating   , only : satelliteTidalHeatingRateClass
  use :: Satellites_Tidal_Fields   , only : satelliteTidalFieldStandard   , satelliteTidalFieldSphericalSymmetry
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric, assignment(=)
  use :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group              , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  ! The host halo, shared by every configuration.
  double precision, parameter :: massHost           =1.0d12, concentrationHost=10.0d0
  double precision, parameter :: timeNode           =13.8d0
  ! The reference configurations and rates, emitted by `satelliteTidalHeatingRate.py --fortran`. Positions are in Mpc, speeds in
  ! km/s, path-integrated tidal tensors in (km/s/Mpc)² Gyr, and heating rates in (km/s/Mpc)² Gyr⁻¹.
  integer         , parameter :: countConfigurations=16
  double precision, dimension(countConfigurations), parameter :: positionX                             =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01, 2.826842234878456d-02, 2.826842234878456d-02, 9.328103463593093d-02, 2.826842234878456d-02, 2.826842234878456d-02, 2.826842234878456d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 9.422807449594853d-03, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: positionY                             =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01,-4.711403724797426d-02,-4.711403724797426d-02, 0.000000000000000d+00,-4.711403724797426d-02,-4.711403724797426d-02,-4.711403724797426d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02,-1.570467908265809d-02, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: positionZ                             =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01, 7.538245959675884d-02, 7.538245959675884d-02, 0.000000000000000d+00, 7.538245959675884d-02, 7.538245959675884d-02, 7.538245959675884d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 2.512748653225294d-02, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: speed                                 =[ 2.352273917617680d+02, 1.764205438213260d+02, 5.880684794044201d+01, 7.056821752853041d+01, 4.704547835235361d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 2.352273917617680d+02, 1.764205438213260d+02]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXX                =[ 1.430775831020361d+05, 1.430775831020361d+05, 1.430775831020361d+05,-4.292327493061082d+04,-4.292327493061082d+04, 1.430775831020361d+05,-1.430775831020361d+05,-1.430775831020361d+05,-4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04,-1.430775831020361d+05, 1.430775831020361d+05,-4.292327493061082d+04, 1.430775831020361d+05]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXY                =[ 5.723103324081444d+04, 5.723103324081444d+04, 5.723103324081444d+04, 1.144620664816289d+05, 1.144620664816289d+05, 5.723103324081444d+04,-5.723103324081444d+04,-5.723103324081444d+04, 1.144620664816289d+05,-1.144620664816289d+05,-1.144620664816289d+05,-1.144620664816289d+05,-5.723103324081444d+04, 5.723103324081444d+04, 1.144620664816289d+05, 5.723103324081444d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXZ                =[-4.292327493061082d+04,-4.292327493061082d+04,-4.292327493061082d+04,-7.153879155101805d+04,-7.153879155101805d+04,-4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 4.292327493061082d+04,-4.292327493061082d+04,-7.153879155101805d+04,-4.292327493061082d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedYY                =[-7.153879155101805d+04,-7.153879155101805d+04,-7.153879155101805d+04,-8.584654986122165d+04,-8.584654986122165d+04,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04,-8.584654986122165d+04, 8.584654986122165d+04, 8.584654986122165d+04, 8.584654986122165d+04, 7.153879155101805d+04,-7.153879155101805d+04,-8.584654986122165d+04,-7.153879155101805d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedYZ                =[ 2.861551662040722d+04, 2.861551662040722d+04, 2.861551662040722d+04,-1.001543081714253d+05,-1.001543081714253d+05, 2.861551662040722d+04,-2.861551662040722d+04,-2.861551662040722d+04,-1.001543081714253d+05, 1.001543081714253d+05, 1.001543081714253d+05, 1.001543081714253d+05,-2.861551662040722d+04, 2.861551662040722d+04,-1.001543081714253d+05, 2.861551662040722d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedZZ                =[-7.153879155101805d+04,-7.153879155101805d+04,-7.153879155101805d+04, 1.287698247918325d+05, 1.287698247918325d+05,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 1.287698247918325d+05,-1.287698247918325d+05,-1.287698247918325d+05,-1.287698247918325d+05, 7.153879155101805d+04,-7.153879155101805d+04, 1.287698247918325d+05,-7.153879155101805d+04]
  double precision, dimension(countConfigurations), parameter :: massSatellite                         =[ 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10]
  double precision, dimension(countConfigurations), parameter :: massBound                             =[ 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 6.000000000000000d+09, 3.000000000000000d+09, 5.000000000000000d+08, 1.000000000000000d+03, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 3.000000000000000d+09]
  double precision, dimension(countConfigurations), parameter :: concentrationSatellite                =[ 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 5.000000000000000d+00, 3.000000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01]
  double precision, dimension(countConfigurations), parameter :: massSpheroid                          =[ 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 1.000000000000000d+09, 3.000000000000000d+09, 1.000000000000000d+09]
  double precision, dimension(countConfigurations), parameter :: radiusSpheroid                        =[ 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 0.000000000000000d+00, 2.000000000000000d-03, 5.000000000000000d-03, 2.000000000000000d-03]
  double precision, dimension(countConfigurations), parameter :: rateHeatingReference                  =[ 7.352871526342838d+05, 4.081543134638479d+04, 3.857328174759763d-01, 7.268422955010923d+03, 2.459168536178509d+05, 2.040771567319241d+05, 3.089602274814853d+05, 6.794783392315524d+04, 1.838192993311568d+04, 8.259565157905137d+02, 7.302336898509945d+00, 1.462768047066803d+05, 0.000000000000000d+00, 3.448814464945776d+04, 1.730093921002922d+06, 2.725535821311682d+02]
  double precision, dimension(countConfigurations), parameter :: tidalTensorXXReference                =[-6.756322663766456d+06,-5.630268886472050d+05,-2.233495095459988d+04,-1.539381380491868d+06,-1.539381380491868d+06, 2.132261567712710d+06,-1.539381380491868d+06,-1.539381380491868d+06,-1.539381380491868d+06,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-1.320663891540206d+07,-5.630268886472050d+05]
  double precision, dimension(countConfigurations), parameter :: tidalTensorXYReference                =[ 8.903253417750537d+06, 1.347644228179957d+06, 9.837116101857356d+04,-6.188162272254905d+05,-6.188162272254905d+05, 0.000000000000000d+00,-6.188162272254905d+05,-6.188162272254905d+05,-6.188162272254905d+05, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06,-4.088228610191579d+06, 1.347644228179957d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorXZReference                =[ 8.903253417750537d+06, 1.347644228179957d+06, 9.837116101857356d+04, 9.901059635607851d+05, 9.901059635607851d+05, 0.000000000000000d+00, 9.901059635607851d+05, 9.901059635607851d+05, 9.901059635607851d+05, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 6.541165776306526d+06, 1.347644228179957d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorYYReference                =[-6.756322663766456d+06,-5.630268886472050d+05,-2.233495095459988d+04,-8.793107381180120d+05,-8.793107381180120d+05,-1.910671116827163d+06,-8.793107381180120d+05,-8.793107381180120d+05,-8.793107381180120d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-8.845861731197711d+06,-5.630268886472050d+05]
  double precision, dimension(countConfigurations), parameter :: tidalTensorYZReference                =[ 8.903253417750537d+06, 1.347644228179957d+06, 9.837116101857356d+04,-1.650176605934642d+06,-1.650176605934642d+06, 0.000000000000000d+00,-1.650176605934642d+06,-1.650176605934642d+06,-1.650176605934642d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06,-1.090194296051088d+07, 1.347644228179957d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorZZReference                =[-6.756322663766456d+06,-5.630268886472050d+05,-2.233495095459988d+04, 7.296114526682650d+05, 7.296114526682650d+05,-1.910671116827163d+06, 7.296114526682650d+05, 7.296114526682650d+05, 7.296114526682650d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05,-5.630268886472050d+05, 1.783532655300398d+06,-5.630268886472050d+05]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalXXReference     =[ 1.232068841650502d+07, 6.292863038697622d+05,-1.041181902943021d+04,-1.486822268331937d+06, 7.965791599495369d+05, 5.709201145263611d+06,-1.210886929492296d+06,-1.210886929492296d+06,-1.210886929492296d+06, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05,-7.950727699408896d+06, 6.292863038697622d+05]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalXYReference     =[ 8.903253417750537d+06, 1.347644228179957d+06, 9.837116101857356d+04,-6.188162272254905d+05,-6.188162272254905d+05, 0.000000000000000d+00,-6.188162272254905d+05,-6.188162272254905d+05,-6.188162272254905d+05, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06,-4.088228610191579d+06, 1.347644228179957d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalXZReference     =[ 2.798026449802201d+07, 2.539957420696924d+06, 1.102942929437432d+05, 1.130263595987269d+06, 7.219334071404534d+06, 0.000000000000000d+00, 1.866091166226312d+06, 1.866091166226312d+06, 1.866091166226312d+06, 2.539957420696924d+06, 2.539957420696924d+06, 2.539957420696924d+06, 2.539957420696924d+06, 2.539957420696924d+06, 2.055692901895496d+07, 2.539957420696924d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalYYReference     =[ 3.139769949677650d+07, 1.821599496386729d+06, 1.511312895739462d+03,-4.529979394874554d+05, 1.806792475657339d+07, 1.666268460723739d+06, 1.785144253322966d+06, 1.785144253322966d+06, 1.785144253322966d+06, 1.821599496386729d+06, 1.821599496386729d+06, 1.821599496386729d+06, 1.821599496386729d+06, 1.821599496386729d+06, 3.378541813185795d+07, 1.821599496386729d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalYZReference     =[ 8.903253417750537d+06, 1.347644228179957d+06, 9.837116101857356d+04,-1.650176605934642d+06,-1.650176605934642d+06, 0.000000000000000d+00,-1.650176605934642d+06,-1.650176605934642d+06,-1.650176605934642d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06, 1.347644228179957d+06,-1.090194296051088d+07, 1.347644228179957d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorCentrifugalZZReference     =[ 1.232068841650502d+07, 6.292863038697622d+05,-1.041181902943021d+04, 1.103365139138890d+06, 1.734088640691826d+07,-1.910671116827163d+06, 3.065571993109671d+06, 3.065571993109671d+06, 3.065571993109671d+06, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05, 6.292863038697622d+05, 3.915890130236290d+07, 6.292863038697622d+05]
  double precision, dimension(countConfigurations), parameter :: tidalTensorRadialReference            =[ 1.105018417173462d+07, 2.132261567712708d+06, 1.744073710825472d+05, 2.132261567712710d+06, 2.132261567712710d+06, 2.132261567712710d+06, 2.132261567712710d+06, 2.132261567712710d+06, 2.132261567712710d+06, 2.132261567712708d+06, 2.132261567712708d+06, 2.132261567712708d+06, 2.132261567712708d+06, 2.132261567712708d+06, 1.105018417173464d+07, 2.132261567712708d+06]
  double precision, dimension(countConfigurations), parameter :: tidalTensorRadialCentrifugalReference =[ 4.920420633227757d+07, 4.516887952746642d+06, 1.982536349328866d+05, 2.558574366343266d+06, 2.107949706240411d+07, 5.709201145263611d+06, 4.796716559153688d+06, 4.796716559153688d+06, 4.796716559153688d+06, 4.516887952746642d+06, 4.516887952746642d+06, 4.516887952746642d+06, 4.516887952746642d+06, 4.516887952746642d+06, 5.368146403479030d+07, 4.516887952746642d+06]
  ! The assertion tolerance, justified rather than tuned.
  !
  ! Every ingredient of the rate that depends only on the geometry - the tidal tensor's orientation and contraction, and the
  ! half-mass radius as a fraction of the virial radius - was measured to agree with the reference to rounding. What remains
  ! is a difference of 7.5 x 10⁻⁶ in the virial radii themselves (a mean density differing by 2.3 x 10⁻⁵), not traced further.
  ! The adiabatic correction amplifies a fractional error in the internal orbital frequency by up to 2γ x²/(1+x²), approaching
  ! 5 at the largest x = ωτ in the grid, so that the rate inherits up to about 6 x 10⁻⁵. Measured: 4.7 x 10⁻⁵. 2 x 10⁻⁴ is
  ! allowed. Any error in the form of the rate - a missing factor of 2π in x, off-diagonal elements counted once, a wrong
  ! exponent - changes it by far more.
  double precision                                , parameter :: toleranceRate                        =2.0d-4
  ! The tidal tensor is compared element by element, at the same tolerance. Its off-diagonal elements pass through zero as the
  ! satellite's position approaches an axis, so each element is also given an absolute tolerance scaled to the xx element of
  ! its own configuration - without which a vanishing element would be compared relatively against zero.
  double precision                                , parameter :: toleranceTensor                      =2.0d-4, absoluteToleranceTensor  =1.0d-12
  class           (darkMatterHaloScaleClass      ), pointer   :: darkMatterHaloScale_
  class           (satelliteTidalHeatingRateClass), pointer   :: satelliteTidalHeatingRate_
  type            (treeNode                      ), pointer   :: nodeHost                                    , nodeSatellite
  class           (nodeComponentBasic            ), pointer   :: basicHost                                   , basicSatellite
  class           (nodeComponentDarkMatterProfile), pointer   :: profileHost                                 , profileSatellite
  class           (nodeComponentSatellite        ), pointer   :: satellite
  class           (nodeComponentSpheroid         ), pointer   :: spheroid
  type            (satelliteTidalFieldStandard         )      :: satelliteTidalFieldStandard_
  type            (satelliteTidalFieldSphericalSymmetry)      :: satelliteTidalFieldSphericalSymmetry_
  type            (tensorRank2Dimension3Symmetric      )      :: tidalTensorStandard                         , tidalTensorSpherical
  type            (inputParameters               )            :: parameters
  character       (len=128                       )            :: message
  integer                                                     :: iConfiguration
  double precision                                            :: radiusVirialHost                            , radiusVirialSatellite           , &
       &                                                         differenceMaximum                           , differenceRate                  , &
       &                                                         rateHeating                                 , tidalTensorRadialStandard       , &
       &                                                         tidalTensorRadialSpherical

  call displayVerbositySet              (verbosityLevelStandard)
  call Error_Handler_Register           (                      )
  parameters=inputParameters(var_str('testSuite/parameters/satelliteTidalHeatingRate.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
  <objectBuilder class="satelliteTidalHeatingRate" name="satelliteTidalHeatingRate_" source="parameters"/>
  !!]
  ! Both tidal field classes are constructed here rather than declared in the parameter file, so that the two can be compared
  ! with each other: they must agree, since the host is spherically symmetric.
  satelliteTidalFieldStandard_         =satelliteTidalFieldStandard         (      darkMatterHaloScale_)
  satelliteTidalFieldSphericalSymmetry_=satelliteTidalFieldSphericalSymmetry(1.0d0,darkMatterHaloScale_)
  call Unit_Tests_Begin_Group("Satellite tidal heating rate")
  differenceMaximum=0.0d0
  do iConfiguration=1,countConfigurations
     ! Build a fresh host and satellite for every configuration. A node memoizes the mass distributions built from its
     ! components, so reusing one and resetting its properties would leave the previous configuration's distribution in place.
     nodeHost         => treeNode                       (                 )
     nodeSatellite    => treeNode                       (                 )
     basicHost        => nodeHost     %basic            (autoCreate=.true.)
     profileHost      => nodeHost     %darkMatterProfile(autoCreate=.true.)
     basicSatellite   => nodeSatellite%basic            (autoCreate=.true.)
     profileSatellite => nodeSatellite%darkMatterProfile(autoCreate=.true.)
     satellite        => nodeSatellite%satellite        (autoCreate=.true.)
     ! Attach the satellite to its host, so that `mergesWith` and `isSatellite` resolve.
     nodeHost     %firstSatellite => nodeSatellite
     nodeSatellite%parent         => nodeHost
     call basicHost     %massSet(massHost                     )
     call basicHost     %timeSet(timeNode                     )
     call basicSatellite%massSet(massSatellite(iConfiguration))
     call basicSatellite%timeSet(timeNode                     )
     radiusVirialHost     =darkMatterHaloScale_%radiusVirial(nodeHost     )
     radiusVirialSatellite=darkMatterHaloScale_%radiusVirial(nodeSatellite)
     call profileHost     %scaleSet(radiusVirialHost     /concentrationHost                     )
     call profileSatellite%scaleSet(radiusVirialSatellite/concentrationSatellite(iConfiguration))
     ! Place the satellite. Only the speed enters the rate, so the velocity is set along the y-axis throughout.
     ! Give the satellite a stellar spheroid where the configuration asks for one, so that the total mass distribution - from
     ! which the circular velocity at the half-mass radius is taken - differs from the dark matter distribution.
     if (massSpheroid(iConfiguration) > 0.0d0) then
        spheroid => nodeSatellite%spheroid(autoCreate=.true.)
        call spheroid%massStellarSet(massSpheroid  (iConfiguration))
        call spheroid%radiusSet     (radiusSpheroid(iConfiguration))
     end if
     call satellite%boundMassSet                (massBound(iConfiguration))
     call satellite%positionSet                 ([positionX(iConfiguration),positionY(iConfiguration),positionZ(iConfiguration)])
     call satellite%velocitySet                 ([0.0d0                    ,speed    (iConfiguration),0.0d0                    ])
     call satellite%tidalTensorPathIntegratedSet(                                                                       &
          &                                      tensorRank2Dimension3Symmetric(                                        &
          &                                                                     tensorPathIntegratedXX(iConfiguration), &
          &                                                                     tensorPathIntegratedXY(iConfiguration), &
          &                                                                     tensorPathIntegratedXZ(iConfiguration), &
          &                                                                     tensorPathIntegratedYY(iConfiguration), &
          &                                                                     tensorPathIntegratedYZ(iConfiguration), &
          &                                                                     tensorPathIntegratedZZ(iConfiguration)  &
          &                                                                    )                                        &
          &                                     )
     ! Evaluate the rate.
     rateHeating=satelliteTidalHeatingRate_%heatingRate(nodeSatellite)
     ! Report the fractional difference of every comparison, so that a reader can see how much of the tolerance is in use. A
     ! test passing only because its tolerance is loose otherwise looks identical to one passing sharply. The rate clamped to zero
     ! has no fractional difference, and is compared absolutely below.
     if (rateHeatingReference(iConfiguration) > 0.0d0) then
        differenceRate   =abs(rateHeating/rateHeatingReference(iConfiguration)-1.0d0)
        differenceMaximum=max(differenceMaximum,differenceRate)
        write (message,'(a,i0,a,e12.5,a,e12.5)') 'configuration ',iConfiguration,': dQ/dt ',rateHeating,', fractional difference ',differenceRate
     else
        write (message,'(a,i0,a,e12.5)'        ) 'configuration ',iConfiguration,': dQ/dt ',rateHeating
     end if
     call displayMessage(trim(message))
     write (message,'(a,i0)') 'tidal heating rate, configuration ',iConfiguration
     call Assert(trim(message),rateHeating,rateHeatingReference(iConfiguration),absTol=1.0d-30,relTol=toleranceRate)
     ! Assert on the host's tidal tensor itself, from both tidal field classes. This is the chief input to the heating rate,
     ! and the rate uses it only through a contraction, which a wrongly oriented tensor can survive. Both classes must return
     ! the tensor in the frame in which the satellite's position is expressed: `sphericalSymmetry` formerly evaluated it along
     ! the x-axis, which returns the tensor of a frame rotating with the satellite.
     tidalTensorStandard        =satelliteTidalFieldStandard_         %tidalTensor      (nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.false.)
     tidalTensorSpherical       =satelliteTidalFieldSphericalSymmetry_%tidalTensor      (nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.false.)
     write (message,'(a,i0)') 'tidal tensor, standard, configuration ',iConfiguration
     call Assert(trim(message),                                                                                                          &
          &      [tidalTensorStandard %element(0,0),tidalTensorStandard %element(0,1),tidalTensorStandard %element(0,2),                 &
          &       tidalTensorStandard %element(1,1),tidalTensorStandard %element(1,2),tidalTensorStandard %element(2,2)],                &
          &      [tidalTensorXXReference(iConfiguration),tidalTensorXYReference(iConfiguration),tidalTensorXZReference(iConfiguration),  &
          &       tidalTensorYYReference(iConfiguration),tidalTensorYZReference(iConfiguration),tidalTensorZZReference(iConfiguration)], &
          &      relTol=toleranceTensor,absTol=absoluteToleranceTensor*abs(tidalTensorXXReference(iConfiguration)))
     write (message,'(a,i0)') 'tidal tensor, spherical symmetry, configuration ',iConfiguration
     call Assert(trim(message),                                                                                                          &
          &      [tidalTensorSpherical%element(0,0),tidalTensorSpherical%element(0,1),tidalTensorSpherical%element(0,2),                 &
          &       tidalTensorSpherical%element(1,1),tidalTensorSpherical%element(1,2),tidalTensorSpherical%element(2,2)],                &
          &      [tidalTensorXXReference(iConfiguration),tidalTensorXYReference(iConfiguration),tidalTensorXZReference(iConfiguration),  &
          &       tidalTensorYYReference(iConfiguration),tidalTensorYZReference(iConfiguration),tidalTensorZZReference(iConfiguration)], &
          &      relTol=toleranceTensor,absTol=absoluteToleranceTensor*abs(tidalTensorXXReference(iConfiguration)))
     ! Repeat with the centrifugal term included. It is built from the angular velocity, ω=r×v/r², and not from the velocity
     ! itself; the two agree only for a circular orbit, and no configuration here is circular.
     tidalTensorStandard        =satelliteTidalFieldStandard_         %tidalTensor      (nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.true. )
     tidalTensorSpherical       =satelliteTidalFieldSphericalSymmetry_%tidalTensor      (nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.true. )
     write (message,'(a,i0)') 'tidal tensor with centrifugal term, standard, configuration ',iConfiguration
     call Assert(trim(message),                                                                                                                                           &
          &      [tidalTensorStandard %element(0,0),tidalTensorStandard %element(0,1),tidalTensorStandard %element(0,2),                                                  &
          &       tidalTensorStandard %element(1,1),tidalTensorStandard %element(1,2),tidalTensorStandard %element(2,2)],                                                 &
          &      [tidalTensorCentrifugalXXReference(iConfiguration),tidalTensorCentrifugalXYReference(iConfiguration),tidalTensorCentrifugalXZReference(iConfiguration),  &
          &       tidalTensorCentrifugalYYReference(iConfiguration),tidalTensorCentrifugalYZReference(iConfiguration),tidalTensorCentrifugalZZReference(iConfiguration)], &
          &      relTol=toleranceTensor,absTol=absoluteToleranceTensor*abs(tidalTensorCentrifugalXXReference(iConfiguration)))
     write (message,'(a,i0)') 'tidal tensor with centrifugal term, spherical symmetry, configuration ',iConfiguration
     call Assert(trim(message),                                                                                                                                           &
          &      [tidalTensorSpherical%element(0,0),tidalTensorSpherical%element(0,1),tidalTensorSpherical%element(0,2),                                                  &
          &       tidalTensorSpherical%element(1,1),tidalTensorSpherical%element(1,2),tidalTensorSpherical%element(2,2)],                                                 &
          &      [tidalTensorCentrifugalXXReference(iConfiguration),tidalTensorCentrifugalXYReference(iConfiguration),tidalTensorCentrifugalXZReference(iConfiguration),  &
          &       tidalTensorCentrifugalYYReference(iConfiguration),tidalTensorCentrifugalYZReference(iConfiguration),tidalTensorCentrifugalZZReference(iConfiguration)], &
          &      relTol=toleranceTensor,absTol=absoluteToleranceTensor*abs(tidalTensorCentrifugalXXReference(iConfiguration)))
     ! The radial component, which is what every caller other than the tidal heating rate uses. With the centrifugal term it
     ! is ω², the *tangential* speed over the radius squared, not the total speed.
     tidalTensorRadialStandard  =satelliteTidalFieldStandard_         %tidalTensorRadial(nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.false.)
     tidalTensorRadialSpherical =satelliteTidalFieldSphericalSymmetry_%tidalTensorRadial(nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.false.)
     write (message,'(a,i0)') 'radial tidal field, configuration ',iConfiguration
     call Assert(trim(message),[tidalTensorRadialStandard,tidalTensorRadialSpherical],spread(tidalTensorRadialReference(iConfiguration),1,2),relTol=toleranceTensor)
     tidalTensorRadialStandard  =satelliteTidalFieldStandard_         %tidalTensorRadial(nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.true. )
     tidalTensorRadialSpherical =satelliteTidalFieldSphericalSymmetry_%tidalTensorRadial(nodeSatellite,nodeHost=nodeHost,includeCentrifugalAcceleration=.true. )
     write (message,'(a,i0)') 'radial tidal field with centrifugal term, configuration ',iConfiguration
     call Assert(trim(message),[tidalTensorRadialStandard,tidalTensorRadialSpherical],spread(tidalTensorRadialCentrifugalReference(iConfiguration),1,2),relTol=toleranceTensor)
     ! Detach before destroying, so that the host does not attempt to destroy the satellite a second time.
     nodeHost     %firstSatellite => null()
     nodeSatellite%parent         => null()
     call nodeSatellite%destroy()
     call nodeHost     %destroy()
     deallocate(nodeSatellite)
     deallocate(nodeHost     )
  end do
  write (message,'(a,e12.5)') 'largest fractional difference over all configurations: ',differenceMaximum
  call displayMessage(trim(message))
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  !![
  <objectDestructor name="darkMatterHaloScale_"      />
  <objectDestructor name="satelliteTidalHeatingRate_"/>
  !!]
end program Test_Satellite_Tidal_Heating_Rate
