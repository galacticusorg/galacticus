!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a galactic high-pass filter for ISM gas fraction (i.e. ISM gas mass to stellar mass ratio).
!!}

  !![
  <galacticFilter name="galacticFilterGasFractionISM">
   <description>
   A galactic high-pass filter for ISM gas fraction (i.e. ISM gas mass to stellar mass ratio).
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterGasFractionISM
     !!{
     A galactic high-pass filter class for ISM gas fraction (i.e. ISM gas mass to stellar mass ratio).
     !!}
     private
     double precision :: fractionGasThreshold
   contains
     procedure :: passes => gasFractionISMPasses
  end type galacticFilterGasFractionISM

  interface galacticFilterGasFractionISM
     !!{
     Constructors for the ``gasFractionISM'' galactic filter class.
     !!}
     module procedure gasFractionISMConstructorParameters
     module procedure gasFractionISMConstructorInternal
  end interface galacticFilterGasFractionISM

contains

  function gasFractionISMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``gasFractionISM'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterGasFractionISM)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: fractionGasThreshold
    
    ! Check and read parameters.
    !![
    <inputParameter>
      <name>fractionGasThreshold</name>
      <source>parameters</source>
      <description>The ISM gas fraction above which to pass.</description>
    </inputParameter>
    !!]
    self=galacticFilterGasFractionISM(fractionGasThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gasFractionISMConstructorParameters

  function gasFractionISMConstructorInternal(fractionGasThreshold) result(self)
    !!{
    Internal constructor for the ``gasFractionISM'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterGasFractionISM)                :: self
    double precision                              , intent(in   ) :: fractionGasThreshold
    !![
    <constructorAssign variables="fractionGasThreshold"/>
    !!]
    return
  end function gasFractionISMConstructorInternal

  logical function gasFractionISMPasses(self,node)
    !!{
    Implement a ISM gas fraction high-pass filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticFilterGasFractionISM), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    class           (nodeComponentDisk           ), pointer               :: disk
    class           (nodeComponentSpheroid       ), pointer               :: spheroid
    double precision                                                      :: massStellar, massGas

    disk        => node    %disk       ()
    spheroid    => node    %spheroid   ()
    massStellar = +disk    %massStellar() &
         &        +spheroid%massStellar()
    massGas     = +disk    %massGas    () &
         &        +spheroid%massGas    ()
    gasFractionISMPasses=massGas > self%fractionGasThreshold*massStellar
    return
  end function gasFractionISMPasses
