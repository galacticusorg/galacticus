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

  !![
  <outputTimes name="outputTimesSimulationSnapshots">
   <description>An output times class which matches output times to snapshot times of a simulation.</description>
   <runTimeFileDependencies paths="fileName"/>
  </outputTimes>
  !!]
  type, extends(outputTimesList) :: outputTimesSimulationSnapshots
     !!{
     Implementation of an output times class which matches output times to snapshot times of a simulation.
     !!}
     private
     type(varying_string) :: fileName
  end type outputTimesSimulationSnapshots

  interface outputTimesSimulationSnapshots
     !!{
     Constructors for the \refClass{outputTimesSimulationSnapshots} output times class.
     !!}
     module procedure simulationSnapshotsConstructorParameters
     module procedure simulationSnapshotsConstructorInternal
  end interface outputTimesSimulationSnapshots

contains

  function simulationSnapshotsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesSimulationSnapshots} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputTimesSimulationSnapshots)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    type (varying_string                )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file from which to read simulation snapshots.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=outputTimesSimulationSnapshots(fileName,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function simulationSnapshotsConstructorParameters

  function simulationSnapshotsConstructorInternal(fileName,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{outputTimesSimulationSnapshots} output times class.
    !!}
    use :: Error         , only : Error_Report
    use :: FoX_DOM       , only : destroy       , node
    use :: IO_XML        , only : XML_Array_Read, XML_Get_First_Element_By_Tag_Name, XML_Parse
    implicit none
    type   (outputTimesSimulationSnapshots)                        :: self
    type   (varying_string                ), intent(in   )         :: fileName
    class  (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    type   (node                          ), pointer               :: doc                , snapshots
    integer(c_size_t                      )                        :: i
    integer                                                        :: ioStatus
    !![
    <constructorAssign variables="fileName, *cosmologyFunctions_"/>
    !!]
    
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(self%fileName,iostat=ioStatus)
    if (ioStatus /= 0) call Error_Report('unable to find or parse the simulation definition file'//{introspection:location})
    snapshots => XML_Get_First_Element_By_Tag_Name(doc,'snapshots')
    call XML_Array_Read(snapshots,'snapshot',self%redshifts)
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    allocate(self%times(size(self%redshifts)))
    do i=1,size(self%redshifts)
       self%times(i)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshifts(i)))
    end do
    return
  end function simulationSnapshotsConstructorInternal
