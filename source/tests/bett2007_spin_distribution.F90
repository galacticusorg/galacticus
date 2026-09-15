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
Contains a program which tests the halo spin distribution of Bett et al. (2007) against independently computed reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Bett2007_Spin_Distribution
  !!{RST
  Tests the :galacticus-class:`haloSpinDistributionBett2007` halo spin distribution against reference values computed
  independently by the ``bett2007SpinDistribution.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository. The distribution is

  .. math::

     {\mathrm{d}p \over \mathrm{d}\lambda} = N \left({\lambda \over \lambda_0}\right)^3 \exp\left[ -\alpha \left({\lambda \over
     \lambda_0}\right)^{3/\alpha} \right] {1 \over \lambda},

  with :math:`\lambda_0` and :math:`\alpha` given by ``[lambda0]`` and ``[alpha]``.

  Note the trailing :math:`1/\lambda`: this is a density in :math:`\lambda` rather than in :math:`\ln \lambda`, so the leading
  power is effectively :math:`\lambda^2`. Substituting :math:`x = \alpha (\lambda/\lambda_0)^{3/\alpha}` turns the distribution
  into a gamma distribution of shape :math:`\alpha`, from which the normalization and moments follow in closed form:

  .. math::

     N = {3 \alpha^{\alpha-1} \over \Gamma(\alpha)}, \qquad \langle \lambda^n \rangle = \lambda_0^n \alpha^{-n\alpha/3}
     {\Gamma(\alpha + n\alpha/3) \over \Gamma(\alpha)}.

  The reference script confirms that the distribution as written integrates to unity, and that these closed forms agree with
  direct numerical quadrature to a relative accuracy of :math:`4\times 10^{-12}`.

  The distribution is asserted at a set of spins spanning the range occupied by physical haloes, which tests the normalization
  and the shape together, and then at a second pair of parameters, :math:`\lambda_0 = 0.035` and :math:`\alpha = 3`, for which
  the closed forms take the exact values :math:`N = 27/2` and :math:`\langle \lambda \rangle = \lambda_0`. That second case
  matters: it would detect a normalization which had been absorbed into a tabulation, and so is not merely a repeat of the
  first at different numbers.

  Evaluating the distribution requires a node, because the class divides the node's angular momentum by the halo's
  characteristic angular momentum scale to form the dimensionless spin. On the default (Peebles) path that scale integrates the
  energy of the halo's dark matter profile, so the node is given a dark matter profile with a physically sensible scale radius;
  the test then sets the angular momentum to place the spin at each required value.

  Tolerances are :math:`10^{-6}` relative throughout: both sides evaluate closed-form expressions involving only the gamma
  function and elementary operations, with no tabulation, interpolation or physical constants involved, so the only
  difference expected is in the evaluation of :math:`\Gamma(\alpha)` itself.
  !!}
  use            :: Cosmology_Functions        , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use            :: Dark_Matter_Halo_Scales    , only : darkMatterHaloScale                , darkMatterHaloScaleClass
  use            :: Dark_Matter_Halo_Spins     , only : Dark_Matter_Halo_Angular_Momentum_Scale
  use            :: Halo_Spin_Distributions    , only : haloSpinDistributionBett2007
  use            :: Display                    , only : displayVerbositySet                , verbosityLevelStandard
  use            :: Error                      , only : Error_Handler_Register
  use            :: Events_Hooks               , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities , only : Functions_Global_Set
  use            :: Galacticus_Nodes           , only : nodeClassHierarchyInitialize       , nodeComponentBasic               , &
  &                                                     nodeComponentDarkMatterProfile     , nodeComponentSpin                , &
  &                                                     treeNode
  use            :: ISO_Varying_String         , only : assignment(=)                      , varying_string
  use            :: Input_Parameters           , only : inputParameters
  use            :: Node_Components            , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, &
  &                                                     Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Unit_Tests                 , only : Assert                             , Unit_Tests_Begin_Group           , &
  &                                                     Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  class           (darkMatterHaloScaleClass    ), pointer                 :: darkMatterHaloScale_
  type            (haloSpinDistributionBett2007)                          :: distributionDefault  , distributionAlternate
  type            (varying_string              )                          :: parameterFile
  type            (inputParameters             )                          :: parameters
  ! The halo for which the distribution is evaluated. Its mass and epoch are arbitrary: the distribution depends on the node
  ! only through the dimensionless spin, which the test sets directly.
  double precision                              , parameter               :: massHalo             =1.0000000000d+12
  double precision                              , parameter               :: timeHalo             =1.3800000000d+01
  ! Parameters of the distribution. The first pair are the class defaults, which are those fit by Bett et al. (2007). The
  ! second pair are chosen so that the closed forms are exact: alpha=3 gives a normalization of 27/2 and a mean of lambda_0.
  double precision                              , parameter               :: lambda0Default       =4.3260000000d-02, alphaDefault  =2.5090000000d+00
  double precision                              , parameter               :: lambda0Alternate     =3.5000000000d-02, alphaAlternate=3.0000000000d+00
  ! Spins at which the distribution is evaluated, and the reference values.
  integer                                       , parameter               :: countSpinsDefault    =6
  double precision                              , dimension(countSpinsDefault) :: spinDefault      =[5.0000000000d-03,1.0000000000d-02,2.0000000000d-02,4.3260000000d-02,8.0000000000d-02,1.5000000000d-01]
  double precision                              , dimension(countSpinsDefault) :: distributionDefaultReference=[2.2945273727d+00,7.1812981274d+00,1.6375854294d+01,1.6898306259d+01,3.7914642223d+00,3.7878648339d-02]
  double precision                              , dimension(countSpinsDefault) :: distributionDefaultComputed
  integer                                       , parameter               :: countSpinsAlternate  =3
  double precision                              , dimension(countSpinsAlternate) :: spinAlternate  =[1.0000000000d-02,3.5000000000d-02,1.0000000000d-01]
  double precision                              , dimension(countSpinsAlternate) :: distributionAlternateReference=[1.3362177065d+01,1.9203583513d+01,5.9649321065d-01]
  double precision                              , dimension(countSpinsAlternate) :: distributionAlternateComputed
  integer                                                                 :: i

  ! Establish error handlers, so that errors reported by the GSL are trapped rather than aborting.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Bett et al. (2007) halo spin distribution")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/bett2007SpinDistribution.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
  !!]
  distributionDefault  =haloSpinDistributionBett2007(lambda0Default  ,alphaDefault  ,darkMatterHaloScale_)
  distributionAlternate=haloSpinDistributionBett2007(lambda0Alternate,alphaAlternate,darkMatterHaloScale_)

  ! Evaluate the distribution at each spin, for both parameter sets.
  do i=1,countSpinsDefault
     distributionDefaultComputed  (i)=distributionAtSpin(distributionDefault  ,spinDefault  (i))
  end do
  do i=1,countSpinsAlternate
     distributionAlternateComputed(i)=distributionAtSpin(distributionAlternate,spinAlternate(i))
  end do

  call Unit_Tests_Begin_Group("Distribution")
  call Assert("Bett et al. (2007) parameters"      ,distributionDefaultComputed  ,distributionDefaultReference  ,relTol=1.0d-6)
  call Assert("alpha = 3, where N = 27/2 is exact" ,distributionAlternateComputed,distributionAlternateReference,relTol=1.0d-6)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  !![
  <objectDestructor name="darkMatterHaloScale_"/>
  !!]

contains

  double precision function distributionAtSpin(distribution_,spin_)
    !!{RST
    Return the spin distribution evaluated at the given dimensionless spin, building the node for which it is evaluated. The
    class forms the dimensionless spin by dividing the node's angular momentum by the halo's characteristic angular momentum
    scale, so the angular momentum is set to the product of the required spin and that scale.
    !!}
    implicit none
    type            (haloSpinDistributionBett2007  ), intent(inout) :: distribution_
    double precision                                , intent(in   ) :: spin_
    type            (treeNode                      ), pointer       :: node
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentSpin             ), pointer       :: spin
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: angularMomentumScale

    node              => treeNode                  (                 )
    basic             => node    %basic            (autoCreate=.true.)
    spin              => node    %spin             (autoCreate=.true.)
    darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
    call basic%massSet            (massHalo)
    call basic%timeSet            (timeHalo)
    call basic%timeLastIsolatedSet(timeHalo)
    ! Give the halo a physically sensible scale radius. The angular momentum scale formed below integrates the energy of the
    ! dark matter profile, which is degenerate if the scale radius is left unset.
    call darkMatterProfile%scaleSet(2.0d-2)
    angularMomentumScale=Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterHaloScale_)
    call spin%angularMomentumSet(spin_*angularMomentumScale)
    distributionAtSpin=distribution_%distribution(node)
    call node%destroy()
    deallocate(node)
    return
  end function distributionAtSpin

end program Test_Bett2007_Spin_Distribution
