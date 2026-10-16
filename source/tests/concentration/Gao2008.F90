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
  Contains a program to test the :cite:t:`gao_redshift_2008` dark matter halo concentration algorithm.
  !!}

program Test_Concentration_Gao2008
  !!{RST
  Test the :galacticus-class:`darkMatterProfileConcentrationGao2008` dark matter halo concentration class.

  That class does not evaluate an expression taken from :cite:t:`gao_redshift_2008`. Their Table 1 tabulates the coefficients of

  .. math::

     \log_{10} c_{200} = A \log_{10} \left( M_{200} / h^{-1} \mathrm{M}_\odot \right) + B

  at five redshifts, for their relaxed halo sample, and the class interpolates those coefficients with smooth functions of
  expansion factor of its own. Two separate things therefore want checking, and are checked in two groups below.

  First, that the class evaluates the interpolating functions its own description gives. Those are reimplemented here and compared
  to machine precision. Because the coefficients apply to a mass in :math:`h^{-1} \mathrm{M}_\odot` with the :math:`h` of the
  simulation :cite:p:`gao_redshift_2008` used - not the :math:`h` of whatever cosmology is in use - and the cosmology here is
  given a different :math:`h`, this group also pins that convention. The coefficients are recovered from the class through its
  public interface, by evaluating the concentration at two masses: :math:`A` is the slope of :math:`\log_{10} c` in
  :math:`\log_{10} M` and :math:`B` follows by back-substitution. That the same pair is recovered from a third mass confirms that
  :math:`\log_{10} c` is exactly linear in :math:`\log_{10} M`, as the fitting form requires.

  Second, that those interpolating functions do reproduce Table 1 of :cite:t:`gao_redshift_2008`. This comparison is loose by
  construction - the interpolation is an approximation, agreeing with the tabulated concentrations to 2% at :math:`z=0` but only
  to 18% by :math:`z=2` - so it is made against a generous tolerance. It is still worth making: it would catch a coefficient
  transcribed wrongly, a sign error, or the wrong :math:`h`, none of which the first group would notice, since the first group
  compares the class only with its own description.

  Finally the two monotonicities the fit encodes: concentration falls with halo mass at fixed redshift, and falls with redshift at
  fixed mass.
  !!}
  use :: Cosmology_Functions                , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters               , only : cosmologyParametersSimple
  use :: Dark_Matter_Profiles_Concentration , only : darkMatterProfileConcentrationGao2008
  use :: Display                            , only : displayVerbositySet                  , verbosityLevelStandard
  use :: Events_Hooks                       , only : eventsHooksInitialize
  use :: Functions_Global_Utilities         , only : Functions_Global_Set
  use :: Galacticus_Nodes                   , only : mergerTree                           , nodeClassHierarchyFinalize       , nodeClassHierarchyInitialize      , nodeComponentBasic, &
       &                                             treeNode
  use :: Input_Parameters                   , only : inputParameters
  use :: Node_Components                    , only : Node_Components_Initialize            , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                         , only : Assert                               , Unit_Tests_Begin_Group           , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  type            (inputParameters                      )                :: parameters
  type            (cosmologyParametersSimple            ), pointer       :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda       ), pointer       :: cosmologyFunctions_
  type            (darkMatterProfileConcentrationGao2008), pointer       :: concentrationGao2008_
  type            (mergerTree                           )                :: tree
  type            (treeNode                             ), pointer       :: node
  class           (nodeComponentBasic                   ), pointer       :: basic
  ! The Hubble parameter of the cosmology used here. It is deliberately *not* the h = 0.73 of the simulation to which
  ! Gao et al. (2008) fit, so that the fit's use of that fixed value is tested rather than masked.
  double precision                                       , parameter     :: HubbleConstant       =67.36d0
  double precision                                       , parameter     :: littleHubbleGao2008  = 0.73d0
  ! Masses [M☉] at which the coefficients are recovered. The first two define the line; the third tests its linearity.
  double precision                                       , parameter     :: massLow              =1.0d12 , massHigh   =1.0d14, &
       &                                                                    massMiddle           =1.0d13
  ! Redshifts of Table 1 of Gao et al. (2008), and the coefficients tabulated there for their relaxed halo sample.
  integer                                                , parameter     :: countRedshifts       =5
  double precision                    , dimension(5)                     :: redshifts            =[0.0d0,0.5d0,1.0d0,2.0d0,3.0d0]
  double precision                    , dimension(5)                     :: coefficientAPaper    =[-0.138d0,-0.125d0,-0.092d0,-0.031d0,-0.004d0]
  double precision                    , dimension(5)                     :: coefficientBPaper    =[+2.646d0,+2.372d0,+1.891d0,+0.985d0,+0.577d0]
  ! Tolerances.
  !! The class is compared with its own description, so only round-off enters.
  double precision                                       , parameter     :: toleranceDescription =1.0d-12
  !! The interpolating functions are an approximation to Table 1, worst at high redshift: the concentration they give at 10¹²
  !! M☉ differs from the tabulated one by 2% at z=0, 12% at z=1, 18% at z=2 and 10% at z=3.
  double precision                                       , parameter     :: tolerancePaper       =2.5d-1
  double precision                    , dimension(5)                     :: coefficientA                  , coefficientB      , &
       &                                                                    coefficientAExpected          , coefficientBExpected, &
       &                                                                    concentrationClass            , concentrationPaper
  double precision                                                       :: expansionFactor               , logarithmExpansionFactor, &
       &                                                                    time                          , concentrationLow  , &
       &                                                                    concentrationHigh             , concentrationMiddle
  integer                                                                :: i

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Gao et al. (2008) concentrations")
  parameters=inputParameters('testSuite/parameters/concentration/Gao2008.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  allocate(cosmologyParameters_ )
  allocate(cosmologyFunctions_  )
  allocate(concentrationGao2008_)
  cosmologyParameters_ =cosmologyParametersSimple            (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=HubbleConstant)
  cosmologyFunctions_  =cosmologyFunctionsMatterLambda       (cosmologyParameters_=cosmologyParameters_)
  ! Scatter is set to zero so that the mean relation is returned, and no random number generator is required.
  concentrationGao2008_=darkMatterProfileConcentrationGao2008(scatter=0.0d0,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  node  => treeNode(hostTree=tree)
  basic => node    %basic   (autoCreate=.true.)
  ! Recover the fit coefficients at each of the redshifts tabulated by Gao et al. (2008).
  do i=1,countRedshifts
     expansionFactor         =+1.0d0             &
          &                   /(                 &
          &                     +1.0d0           &
          &                     +redshifts   (i) &
          &                    )
     time                    = cosmologyFunctions_%cosmicTime(expansionFactor)
     call basic%timeSet            (time)
     call basic%timeLastIsolatedSet(time)
     call basic%massSet            (massLow   )
     concentrationLow        = concentrationGao2008_%concentration(node)
     call basic%massSet            (massHigh  )
     concentrationHigh       = concentrationGao2008_%concentration(node)
     call basic%massSet            (massMiddle)
     concentrationMiddle     = concentrationGao2008_%concentration(node)
     ! A is the slope of log₁₀c in log₁₀M, and B follows from either point.
     coefficientA        (i) =+log10(concentrationHigh/concentrationLow)          &
          &                   /log10(massHigh         /massLow         )
     coefficientB        (i) =+log10(concentrationLow)                            &
          &                   -coefficientA(i)*log10(littleHubbleGao2008*massLow)
     ! The concentration at the middle mass must lie on that same line, since the fit is a power law in mass.
     call Assert('log₁₀c is linear in log₁₀M',concentrationMiddle,10.0d0**(coefficientA(i)*log10(littleHubbleGao2008*massMiddle)+coefficientB(i)),relTol=toleranceDescription)
     ! The interpolating functions of the class description.
     logarithmExpansionFactor=log10(cosmologyFunctions_%expansionFactor(time))
     coefficientAExpected(i) =-0.140d0*exp(-((logarithmExpansionFactor+0.05d0)/0.35d0)**2)
     coefficientBExpected(i) =+2.646d0*exp(-((logarithmExpansionFactor+0.00d0)/0.50d0)**2)
     ! The concentration at the low mass, as given by the class and by Table 1 of Gao et al. (2008).
     concentrationClass  (i) = concentrationLow
     concentrationPaper  (i) =10.0d0**(coefficientAPaper(i)*log10(littleHubbleGao2008*massLow)+coefficientBPaper(i))
  end do
  call Unit_Tests_Begin_Group("interpolating functions of the class description")
  call Assert('coefficient A',coefficientA,coefficientAExpected,relTol=toleranceDescription)
  call Assert('coefficient B',coefficientB,coefficientBExpected,relTol=toleranceDescription)
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group("Table 1 of Gao et al. (2008)")
  call Assert('concentration at 10¹² M☉',concentrationClass,concentrationPaper,relTol=tolerancePaper)
  call Unit_Tests_End_Group()
  ! The two monotonicities the fit encodes. A is negative at every redshift, so concentration falls with mass; and the
  ! concentration of a halo of fixed mass falls as redshift rises.
  call Unit_Tests_Begin_Group("monotonicity")
  call Assert('concentration falls with halo mass',all(coefficientA        <      0.0d0                        ),.true.)
  call Assert('concentration falls with redshift' ,all(concentrationClass(2:countRedshifts) < concentrationClass(1:countRedshifts-1)),.true.)
  call Unit_Tests_End_Group()
  ! Clean up.
  call node%destroy()
  deallocate(concentrationGao2008_)
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
end program Test_Concentration_Gao2008
