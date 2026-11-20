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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a program which tests the outer radius of the integral in the :galacticus-class:`stellarFeedbackOutflowsCreasey2013`
stellar feedback outflow class.
!!}

program Test_Stellar_Feedback_Outflows_Creasey2013
  !!{RST
  Tests the outer radius to which the :galacticus-class:`stellarFeedbackOutflowsCreasey2013` class integrates the outflow rate.
  :cite:t:`creasey_how_2013` integrate to infinite radius. For an untruncated Kennicutt-Schmidt law with exponent :math:`N=1.4`
  and an exponential disk whose gas fraction is independent of radius (as here, where gas and stars share a scale length), the
  integrand is proportional to :math:`x \exp[-(N-\mu)x]`, with :math:`x` the radius in units of the scale length and
  :math:`\mu=1.15`, so the fraction of the outflow rate beyond :math:`x_0` is :math:`(1+0.25x_0)\exp(-0.25x_0)`. The test
  checks that:

  * with all parameters at their defaults, the outflow rate is within 1% of that integrated to (effectively) infinite radius,
    taken to be 200 scale lengths;
  * integrated to 10 scale lengths (the outer radius used previously), the outflow rate is :math:`1-3.5\exp(-2.5)=0.7127` of that
    integrated to infinite radius.
  !!}
  use :: Display                                  , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize            , nodeComponentBasic               , nodeComponentDisk                  , treeNode
  use :: ISO_Varying_String                       , only : assignment(=)                           , varying_string
  use :: Input_Parameters                         , only : inputParameters
  use :: Node_Components                          , only : Node_Components_Initialize              , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Star_Formation_Rate_Surface_Density_Disks, only : starFormationRateSurfaceDensityDisksClass
  use :: Stellar_Feedback_Outflows                , only : stellarFeedbackOutflowsClass            , stellarFeedbackOutflowsCreasey2013
  use :: Unit_Tests                               , only : Assert                                  , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (varying_string                           )          :: parameterFile
  type            (inputParameters                          )          :: parameters
  class           (stellarFeedbackOutflowsClass             ), pointer :: stellarFeedbackOutflowsDefault_
  class           (starFormationRateSurfaceDensityDisksClass), pointer :: starFormationRateSurfaceDensityDisks_
  type            (stellarFeedbackOutflowsCreasey2013       )          :: stellarFeedbackOutflowsInfinite_     , stellarFeedbackOutflowsTen_
  type            (treeNode                                 ), pointer :: node
  class           (nodeComponentBasic                       ), pointer :: basic
  class           (nodeComponentDisk                        ), pointer :: disk
  ! Star formation and energy input rates [M☉/Gyr and arbitrary units] - only their ratio matters.
  double precision                                           , parameter :: rateStarFormation      =1.0d0, rateEnergyInput=1.0d0
  double precision                                                       :: rateOutflowDefault           , rateOutflowInfinite  , &
       &                                                                    rateOutflowTen               , rateOutflowExpulsive

  call displayVerbositySet(verbosityLevelStandard)
  parameterFile='testSuite/parameters/stellarFeedbackOutflowsCreasey2013.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="stellarFeedbackOutflows"              name="stellarFeedbackOutflowsDefault_"       source="parameters"/>
  <objectBuilder class="starFormationRateSurfaceDensityDisks" name="starFormationRateSurfaceDensityDisks_" source="parameters"/>
  !!]
  stellarFeedbackOutflowsInfinite_=stellarFeedbackOutflowsCreasey2013(mu=1.15d0,nu=0.16d0,beta0=13.0d0,radiusOuterDimensionless=200.0d0,starFormationRateSurfaceDensityDisks_=starFormationRateSurfaceDensityDisks_)
  stellarFeedbackOutflowsTen_     =stellarFeedbackOutflowsCreasey2013(mu=1.15d0,nu=0.16d0,beta0=13.0d0,radiusOuterDimensionless= 10.0d0,starFormationRateSurfaceDensityDisks_=starFormationRateSurfaceDensityDisks_)
  ! Build a node with a disk of gas fraction 0.3.
  node  => treeNode  (                 )
  basic => node%basic(autoCreate=.true.)
  disk  => node%disk (autoCreate=.true.)
  call disk%massGasSet    (3.0d9 )
  call disk%massStellarSet(7.0d9 )
  call disk%radiusSet     (3.0d-3)
  call Unit_Tests_Begin_Group("Creasey et al. (2013) outflow rate outer radius")
  call stellarFeedbackOutflowsDefault_ %outflowRate(disk,rateStarFormation,rateEnergyInput,rateOutflowDefault ,rateOutflowExpulsive)
  call stellarFeedbackOutflowsInfinite_%outflowRate(disk,rateStarFormation,rateEnergyInput,rateOutflowInfinite,rateOutflowExpulsive)
  call stellarFeedbackOutflowsTen_     %outflowRate(disk,rateStarFormation,rateEnergyInput,rateOutflowTen     ,rateOutflowExpulsive)
  call Assert('outflow rate is positive'                                   ,rateOutflowInfinite > 0.0d0         ,.true.                         )
  call Assert('default outer radius gives the infinite-radius outflow rate',rateOutflowDefault/rateOutflowInfinite,1.0d0                ,relTol=1.0d-2)
  call Assert('outer radius of 10 scale lengths omits 29% of the outflow'  ,rateOutflowTen    /rateOutflowInfinite,1.0d0-3.5d0*exp(-2.5d0),relTol=1.0d-2)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  call node%destroy()
  deallocate(node)
  !![
  <objectDestructor name="stellarFeedbackOutflowsDefault_"      />
  <objectDestructor name="starFormationRateSurfaceDensityDisks_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Stellar_Feedback_Outflows_Creasey2013
