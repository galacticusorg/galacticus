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

  !!{
  Implementation of an active mass for star formation class in which the entire ISM is active.
  !!}

  !![
  <starFormationActiveMass name="starFormationActiveMassTotalISM">
   <description>An active mass for star formation class in which the entire ISM is active.</description>
  </starFormationActiveMass>
  !!]
  type, extends(starFormationActiveMassClass) :: starFormationActiveMassTotalISM
     !!{
     Implementation of n active mass for star formation class in which the entire ISM is active.
     !!}
     private
   contains
     procedure :: massActive => totalISMMassActive
  end type starFormationActiveMassTotalISM

  interface starFormationActiveMassTotalISM
     !!{
     Constructors for the \refClass{starFormationActiveMassTotalISM} active mass for star formation class.
     !!}
     module procedure totalISMConstructorParameters
  end interface starFormationActiveMassTotalISM

contains

  function totalISMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationActiveMassTotalISM} active mass for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(starFormationActiveMassTotalISM)                :: self
    type(inputParameters                ), intent(inout) :: parameters
    
    self=starFormationActiveMassTotalISM()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function totalISMConstructorParameters

  double precision function totalISMMassActive(self,component)
    !!{
    Returns the mass (in $\mathrm{M}_\odot$) of gas actively undergoing star formation in the given {\normalfont \ttfamily
    component}, assuming that the entire ISM is active.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    class(starFormationActiveMassTotalISM), intent(inout) :: self
    class(nodeComponent                  ), intent(inout) :: component

    select type (component)
    class is (nodeComponentDisk    )
       totalISMMassActive=component%massGas()
    class is (nodeComponentSpheroid)
       totalISMMassActive=component%massGas()
    class is (nodeComponentNSC     )
       totalISMMassActive=component%massGas()
    class default
       totalISMMassActive=0.0d0
       call Error_Report('unsupported class'//{introspection:location})
    end select
    return
  end function totalISMMassActive
