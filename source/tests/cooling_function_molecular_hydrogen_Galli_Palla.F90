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
Contains a program which tests the H--H\ :math:`_2` cooling function of the
:galacticus-class:`coolingFunctionMolecularHydrogenGalliPalla` class.
!!}

program Test_Cooling_Function_Molecular_Hydrogen_Galli_Palla
  !!{RST
  Tests the H--H\ :math:`_2` cooling function of the :galacticus-class:`coolingFunctionMolecularHydrogenGalliPalla` class, and
  its logarithmic gradients with respect to density and temperature, against reference values computed independently of
  Galacticus by the ``referenceGalliPallaH2Cooling.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository, which writes the cooling function directly from
  :cite:t:`galli_chemistry_1998` and :cite:t:`hollenbach_molecule_1979`.

  The gas is half atomic hydrogen, with a molecular hydrogen fraction of :math:`10^{-3}`, at hydrogen densities of :math:`1`,
  :math:`10^4` and :math:`10^{10}` cm\ :math:`^{-3}` (below, near, and far above the critical density) and temperatures of 300, 1000
  and 3000 K. The critical density of :cite:t:`galli_chemistry_1998` is defined relative to the density of *atomic* hydrogen,
  the collision partner of H\ :math:`_2`, so that the cooling rate tends to the LTE emission per molecule times the molecular
  hydrogen density at high density. Normalizing by the total hydrogen density instead would reduce the cooling at the two higher
  densities by up to a factor of two for this gas.

  The temperature gradient is tested only at the highest density, where the cooling is in local thermodynamic equilibrium, so
  that the gradient is that of the LTE cooling function. There :math:`n^\mathrm{cr}/n(\mathrm{H}) \approx 10^{-6}`, so that the
  gradient of the low-density limit (which this test does not cover) contributes less than :math:`10^{-4}` of it.
  !!}
  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances                        , Chemicals_Index
  use :: Cooling_Functions            , only : coolingFunctionMolecularHydrogenGalliPalla
  use :: Display                      , only : displayVerbositySet                       , verbosityLevelStandard
  use :: Events_Hooks                 , only : eventsHooksInitialize
  use :: Functions_Global_Utilities   , only : Functions_Global_Set
  use :: Galacticus_Nodes             , only : nodeClassHierarchyInitialize              , nodeComponentBasic               , treeNode                           , nodeClassHierarchyFinalize
  use :: Input_Paths                  , only : inputPath                                 , pathTypeExec
  use :: Input_Parameters             , only : inputParameters
  use :: ISO_Varying_String           , only : assignment(=)                             , operator(//)                     , varying_string
  use :: Node_Components              , only : Node_Components_Initialize                , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Radiation_Fields             , only : radiationFieldNull
  use :: Unit_Tests                   , only : Assert                                    , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                  ), pointer      :: node
  class           (nodeComponentBasic                        ), pointer      :: basic
  type            (coolingFunctionMolecularHydrogenGalliPalla)               :: coolingFunction_
  type            (radiationFieldNull                        )               :: radiation
  type            (abundances                                )               :: gasAbundances
  type            (chemicalAbundances                        )               :: chemicalDensities
  type            (inputParameters                           )               :: parameters
  type            (varying_string                            )               :: parameterFile
  ! Composition of the gas, as fractions of the hydrogen density.
  double precision                                            , parameter    :: fractionAtomic           =0.5d0               , fractionMolecular=1.0d-3
  ! Densities [cm⁻³] and temperatures [K] at which the cooling function is tested.
  integer                                                     , parameter    :: countDensities           =3                   , countTemperatures=3
  double precision                                            , dimension(3) :: numberDensitiesHydrogen  =[1.0d0,1.0d4,1.0d10], temperatures     =[3.0d2,1.0d3,3.0d3]
  ! Reference values computed by referenceGalliPallaH2Cooling.py, ordered with density in the outer loop and temperature in the
  ! inner loop. Cooling rates per unit volume [erg cm⁻³ s⁻¹], and logarithmic gradients with respect to density (at fixed
  ! composition) and temperature.
  double precision                                            , dimension(9) :: coolingReference         =[5.6482120634463d-30,2.8822048323793d-28,1.2311492985422d-26,8.4464620615498d-23,1.4130232750083d-20,6.8460358568289d-19,9.9306705067642d-17,2.7717550040403d-14,1.5419824254126d-12], &
       &                                                                        slopeDensityReference    =[1.9994312356886d+00,1.9998960152634d+00,1.9999201581200d+00,1.1494571721116d+00,1.4902067772127d+00,1.5560242989210d+00,1.0000001757149d+00,1.0000009615752d+00,1.0000012523825d+00], &
       &                                                                        slopeTemperatureReference=[3.4351591418869d+00,3.3445905157237d+00,3.2856499161582d+00,3.7978148241535d+00,4.2784223115433d+00,2.9265484666752d+00,3.8615831570610d+00,5.1765573626206d+00,2.4767387263225d+00]
  ! Tolerances.
  !! The LTE and low-density limit cooling functions are tabulated in log₁₀(T) at 100 points per decade and interpolated
  !! linearly in their logarithms.
  double precision                                            , parameter    :: toleranceCooling         =1.0d-3
  !! The density gradient is computed analytically by the class, and by finite differences in the reference.
  double precision                                            , parameter    :: toleranceSlopeDensity    =1.0d-4
  !! The temperature gradient is computed analytically by the class, and by finite differences in the reference.
  double precision                                            , parameter    :: toleranceSlopeTemperature=1.0d-3
  double precision                                            , dimension(9) :: cooling                                     , slopeDensity          , &
       &                                                                        slopeTemperature
  integer                                                                    :: i                                           , j                     , &
       &                                                                        k
  integer                                                                    :: indexAtomicHydrogen                         , indexMolecularHydrogen

  ! Initialize.
  parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/coolingFunctionMolecularHydrogenGalliPalla.xml'
  parameters   =inputParameters(parameterFile)
  call displayVerbositySet              (verbosityLevelStandard)
  call eventsHooksInitialize            (                      )
  call Functions_Global_Set             (                      )
  call nodeClassHierarchyInitialize     (parameters            )
  call Node_Components_Initialize       (parameters            )
  call Node_Components_Thread_Initialize(parameters            )
  node  => treeNode      (                 )
  basic => node    %basic(autoCreate=.true.)
  call basic%timeSet            (13.8d0)
  call basic%timeLastIsolatedSet(13.8d0)
  coolingFunction_      =coolingFunctionMolecularHydrogenGalliPalla()
  radiation             =radiationFieldNull                        ()
  indexAtomicHydrogen   =Chemicals_Index("AtomicHydrogen"   )
  indexMolecularHydrogen=Chemicals_Index("MolecularHydrogen")
  ! Evaluate the cooling function and its gradients.
  k=0
  do i=1,countDensities
     call chemicalDensities%abundanceSet(indexAtomicHydrogen   ,fractionAtomic   *numberDensitiesHydrogen(i))
     call chemicalDensities%abundanceSet(indexMolecularHydrogen,fractionMolecular*numberDensitiesHydrogen(i))
     do j=1,countTemperatures
        k=k+1
        cooling         (k)=coolingFunction_%coolingFunction                   (node,numberDensitiesHydrogen(i),temperatures(j),gasAbundances,chemicalDensities,radiation)
        slopeDensity    (k)=coolingFunction_%coolingFunctionDensityLogSlope    (node,numberDensitiesHydrogen(i),temperatures(j),gasAbundances,chemicalDensities,radiation)
        slopeTemperature(k)=coolingFunction_%coolingFunctionTemperatureLogSlope(node,numberDensitiesHydrogen(i),temperatures(j),gasAbundances,chemicalDensities,radiation)
     end do
  end do
  call Unit_Tests_Begin_Group("Galli & Palla (1998) H-H₂ cooling function")
  call Assert('cooling function'                           ,cooling              ,coolingReference              ,relTol=toleranceCooling         )
  call Assert('density log slope'                          ,slopeDensity         ,slopeDensityReference         ,relTol=toleranceSlopeDensity    )
  call Assert('temperature log slope (LTE limit)'          ,slopeTemperature(7:9),slopeTemperatureReference(7:9),relTol=toleranceSlopeTemperature)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call node%destroy()
  deallocate(node)
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Cooling_Function_Molecular_Hydrogen_Galli_Palla
