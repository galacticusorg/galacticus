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
Contains a program which tests the fraction of outflowing gas stripped from the hot halos of satellites.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Hot_Halo_Outflow_Stripping
  !!{RST
  Tests the fraction of outflowing gas stripped from the hot halos of satellites by the
  :galacticus-class:`hotHaloOutflowStrippingStandard` class, against reference values. Following :cite:t:`font_colours_2008`
  this fraction is :math:`\epsilon_\mathrm{strip}` times the fraction of the *unstripped* hot halo which lies beyond the outer
  radius of the stripped halo. For a :math:`\beta=2/3` profile with core radius :math:`r_\mathrm{c}` this is

  .. math::

     f = \epsilon_\mathrm{strip} \left[ 1 - {x_\mathrm{outer} - \tan^{-1} x_\mathrm{outer} \over x_\mathrm{virial} - \tan^{-1}
     x_\mathrm{virial}} \right],

  where :math:`x=r/r_\mathrm{c}`. The reference values were computed by numerically integrating the :math:`\beta=2/3` density
  profile (using ``scipy.integrate.quad``), for :math:`r_\mathrm{c}=0.3 r_\mathrm{virial}`, rather than from this closed form.

  The cases tested are a satellite whose hot halo outer radius is at the virial radius (for which no stripping occurs), satellites
  whose outer radius is stripped to several fractions of the virial radius, a satellite with an empty hot halo (for which the
  fraction depends only on the shape of the profile, and so must equal that of the corresponding non-empty halo), and an
  isolated halo (for which no stripping occurs). The truncated hot halo mass distribution contains no mass beyond the outer
  radius, so these cases fail if the fraction is computed from the truncated, rather than the unstripped, distribution.
  !!}
  use :: Abundances_Structure       , only : zeroAbundances
  use :: Cosmology_Functions        , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales    , only : darkMatterHaloScale                , darkMatterHaloScaleClass
  use :: Display                    , only : displayVerbositySet                , verbosityLevelStandard
  use :: Error                      , only : Error_Handler_Register
  use :: Events_Hooks               , only : eventsHooksInitialize
  use :: Functions_Global_Utilities , only : Functions_Global_Set
  use :: Galacticus_Nodes           , only : nodeClassHierarchyInitialize       , nodeComponentBasic               , &
       &                                     nodeComponentHotHalo               , treeNode
  use :: Hot_Halo_Outflows_Stripping, only : hotHaloOutflowStripping            , hotHaloOutflowStrippingClass
  use :: IO_HDF5                    , only : ioHDF5AccessInitialize
  use :: ISO_Varying_String         , only : assignment(=)                      , varying_string
  use :: Input_Parameters           , only : inputParameters
  use :: Node_Components            , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, &
       &                                     Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                 , only : Assert                             , Unit_Tests_Begin_Group           , &
       &                                     Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  class           (cosmologyFunctionsClass     ), pointer                 :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass    ), pointer                 :: darkMatterHaloScale_
  class           (hotHaloOutflowStrippingClass), pointer                 :: hotHaloOutflowStripping_
  type            (treeNode                    ), pointer                 :: nodeHost                        , nodeSatellite
  class           (nodeComponentBasic          ), pointer                 :: basic
  class           (nodeComponentHotHalo        ), pointer                 :: hotHalo
  type            (varying_string              )                          :: parameterFile
  type            (inputParameters             )                          :: parameters
  ! The stripping efficiency - must match that in the parameter file.
  double precision                              , parameter               :: efficiency               =0.25d0
  ! Properties of the test halo.
  double precision                              , parameter               :: massHalo                 =1.0d12, massHotFraction=0.1d0
  ! Test cases: ratios of the hot halo outer radius to the virial radius. The outer radius is limited to be no less than one
  ! tenth of the virial radius, so all cases lie above that.
  integer                                       , parameter               :: countCases               =4
  double precision                              , dimension(countCases)   :: radiusOuterFraction      =[1.0000000000d+0,7.0000000000d-1,5.0000000000d-1,2.0000000000d-1]
  ! Reference values of the stripped fraction divided by the efficiency.
  double precision                              , dimension(countCases)   :: fractionStrippedReference=[0.0000000000d+0,4.3162983635d-1,6.9021822762d-1,9.6170189848d-1]
  double precision                              , dimension(countCases)   :: fractionStrippedComputed
  double precision                                                        :: time                            , radiusVirial            , &
       &                                                                     fractionStrippedEmpty           , fractionStrippedIsolated
  integer                                                                 :: i

  ! Establish error handlers, so that errors reported by the GSL are trapped rather than aborting.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize the lock used to serialize access to the HDF5 library.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Hot halo outflow stripping")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/hotHaloOutflowStripping.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="cosmologyFunctions"      name="cosmologyFunctions_"      source="parameters"/>
  <objectBuilder class="darkMatterHaloScale"     name="darkMatterHaloScale_"     source="parameters"/>
  <objectBuilder class="hotHaloOutflowStripping" name="hotHaloOutflowStripping_" source="parameters"/>
  !!]
  time=cosmologyFunctions_%cosmicTime(1.0d0)
  ! Build a host halo, and a satellite within it. A node is a satellite if it has a parent of which it is not a child.
  nodeHost      => treeNode()
  nodeSatellite => treeNode()
  nodeSatellite%parent => nodeHost
  basic => nodeHost     %basic(autoCreate=.true.)
  call basic%massSet            (1.0d1*massHalo)
  call basic%timeSet            (time          )
  call basic%timeLastIsolatedSet(time          )
  basic => nodeSatellite%basic(autoCreate=.true.)
  call basic%massSet            (      massHalo)
  call basic%timeSet            (time          )
  call basic%timeLastIsolatedSet(time          )
  hotHalo      => nodeSatellite                     %hotHalo     (autoCreate=.true.)
  radiusVirial =  darkMatterHaloScale_%radiusVirial(nodeSatellite                  )
  call hotHalo%abundancesSet(zeroAbundances)
  ! Stripped fractions for a range of outer radii.
  call hotHalo%massSet(massHotFraction*massHalo)
  do i=1,countCases
     call hotHalo%outerRadiusSet(radiusOuterFraction(i)*radiusVirial)
     fractionStrippedComputed(i)=hotHaloOutflowStripping_%fractionStripped(nodeSatellite)/efficiency
  end do
  ! An empty hot halo.
  call hotHalo%massSet       (                       0.0d0             )
  call hotHalo%outerRadiusSet(radiusOuterFraction(3)*radiusVirial      )
  fractionStrippedEmpty   =hotHaloOutflowStripping_%fractionStripped(nodeSatellite)/efficiency
  ! An isolated halo.
  call hotHalo%massSet       (massHotFraction       *massHalo          )
  nodeSatellite%parent => null()
  fractionStrippedIsolated=hotHaloOutflowStripping_%fractionStripped(nodeSatellite)/efficiency
  call Assert("stripped fraction"                ,fractionStrippedComputed,fractionStrippedReference   ,relTol=1.0d-6,absTol=1.0d-9)
  call Assert("stripped fraction (empty halo)"   ,fractionStrippedEmpty   ,fractionStrippedReference(3),relTol=1.0d-6              )
  call Assert("stripped fraction (isolated halo)",fractionStrippedIsolated,0.0d0                                     ,absTol=1.0d-9)
  call nodeSatellite%destroy()
  call nodeHost     %destroy()
  deallocate(nodeSatellite)
  deallocate(nodeHost     )
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  !![
  <objectDestructor name="cosmologyFunctions_"     />
  <objectDestructor name="darkMatterHaloScale_"    />
  <objectDestructor name="hotHaloOutflowStripping_"/>
  !!]
end program Test_Hot_Halo_Outflow_Stripping
