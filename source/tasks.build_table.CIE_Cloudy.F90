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

  !![
  <task name="taskBuildTableCIECloudy">
   <description>A task which builds collisional ionization equilibrium tables using Cloudy.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildTableCIECloudy
     !!{
     Implementation of a task which builds collisional ionization equilibrium tables using Cloudy.
     !!}
     private
     type            (varying_string) :: fileNameCoolingFunction      , fileNameChemicalState
     double precision                 :: metallicityLogarithmicMaximum
     logical                          :: includeContinuum
   contains
     procedure :: perform            => buildTableCIECloudyPerform
     procedure :: requiresOutputFile => buildTableCIECloudyRequiresOutputFile
  end type taskBuildTableCIECloudy

  interface taskBuildTableCIECloudy
     !!{
     Constructors for the \refClass{taskBuildTableCIECloudy} task.
     !!}
     module procedure buildTableCIECloudyParameters
     module procedure buildTableCIECloudyInternal
  end interface taskBuildTableCIECloudy

contains

  function buildTableCIECloudyParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildTableCIECloudy} task class which takes a parameter set as input.
    !!}
    use :: Input_Paths     , only : inputPath      , pathTypeDataDynamic
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (taskBuildTableCIECloudy)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    type            (varying_string         )                :: fileNameCoolingFunction      , fileNameChemicalState
    double precision                                         :: metallicityLogarithmicMaximum
    logical                                                  :: includeContinuum

    !![
    <inputParameter>
      <name>fileNameCoolingFunction</name>
      <description>The file name to which the cooling function table should be stored.</description>
      <defaultValue>inputPath(pathTypeDataDynamic)//'cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fileNameChemicalState</name>
      <description>The file name to which the chemical state table should be stored.</description>
      <defaultValue>inputPath(pathTypeDataDynamic)//'chemicalState/chemical_state_Atomic_CIE_Cloudy.hdf5'</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>metallicityLogarithmicMaximum</name>
      <description>The maximum metallicity to tabulated, expressed as log-10 relative to Solar.</description>
      <defaultValue>1.5d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeContinuum</name>
      <description>If true include the cumulative fraction of total power emitted in the continuum.</description>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=taskBuildTableCIECloudy(fileNameCoolingFunction,fileNameChemicalState,metallicityLogarithmicMaximum,includeContinuum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function buildTableCIECloudyParameters

  function buildTableCIECloudyInternal(fileNameCoolingFunction,fileNameChemicalState,metallicityLogarithmicMaximum,includeContinuum) result(self)
    !!{
    Internal constructor for the \refClass{taskBuildTableCIECloudy} task class.
    !!}
    implicit none
    type            (taskBuildTableCIECloudy)                :: self
    type            (varying_string         ), intent(in   ) :: fileNameCoolingFunction      , fileNameChemicalState
    double precision                         , intent(in   ) :: metallicityLogarithmicMaximum
    logical                                  , intent(in   ) :: includeContinuum
    !![
    <constructorAssign variables="fileNameCoolingFunction, fileNameChemicalState, metallicityLogarithmicMaximum, includeContinuum"/>
    !!]

    return
  end function buildTableCIECloudyInternal

  subroutine buildTableCIECloudyPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display              , only : displayIndent                , displayMessage, displayUnindent
    use :: Error     , only : errorStatusSuccess
    use :: Interfaces_Cloudy_CIE, only : Interface_Cloudy_CIE_Tabulate
    implicit none
    class  (taskBuildTableCIECloudy), intent(inout), target   :: self
    integer                         , intent(  out), optional :: status

    call displayIndent  ('Begin task: tabulate collisional ionization equilibrium using Cloudy')
    call Interface_Cloudy_CIE_Tabulate(                                                                  &
         &                             metallicityMaximumLogarithmic=self%metallicityLogarithmicMaximum, &
         &                             fileNameCoolingFunction      =self%fileNameCoolingFunction      , &
         &                             fileNameChemicalState        =self%fileNameChemicalState        , &
         &                             versionFileFormat            =     1                            , &
         &                             includeContinuum             =self%includeContinuum               &
         &                            )
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: tabulate collisional ionization equilibrium using Cloudy')
    return
  end subroutine buildTableCIECloudyPerform

  logical function buildTableCIECloudyRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildTableCIECloudy), intent(inout) :: self
    !$GLC attributes unused :: self

    buildTableCIECloudyRequiresOutputFile=.false.
    return
  end function buildTableCIECloudyRequiresOutputFile
