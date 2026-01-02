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
  Implements calculations of attenuation of stellar spectra using the model of \cite{witt_multiple_2000}.
  !!}

  !![
  <enumeration>
   <name>wittGordon2000Model</name>
   <description>Enumerates the models available in the {\normalfont \ttfamily wittGordon2000} dust attenuation class.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="milkyWayShellTau3"/>
   <entry label="SMCShellTau3"     />
  </enumeration>
  !!]

  !![
  <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationWittGordon2000">
   <description>Returns the dust attenuation of stellar spectra according to the model of \cite{witt_multiple_2000}.</description>
  </stellarSpectraDustAttenuation>
  !!]
  type, extends(stellarSpectraDustAttenuationTabulated) :: stellarSpectraDustAttenuationWittGordon2000
     !!{
     A class implementing calculations of attenuation of stellar spectra using the model of \cite{witt_multiple_2000}.
     !!}
     private
     type(enumerationWittGordon2000ModelType) :: model
   contains
  end type stellarSpectraDustAttenuationWittGordon2000

  interface stellarSpectraDustAttenuationWittGordon2000
     !!{
     Constructors for the \refClass{stellarSpectraDustAttenuationWittGordon2000} stellar spectra dust attenuation class.
     !!}
     module procedure wittGordon2003ConstructorParameters
     module procedure wittGordon2003ConstructorInternal
  end interface stellarSpectraDustAttenuationWittGordon2000

contains

  function wittGordon2003ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily wittGordon2003} stellar spectra dust attenuation class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(stellarSpectraDustAttenuationWittGordon2000)                :: self
    type(inputParameters                            ), intent(inout) :: parameters
    type(varying_string                             )                :: model

    !![
    <inputParameter>
      <name>model</name>
      <defaultValue>var_str('MilkyWayShellTau3.0')</defaultValue>
      <description>The name of the model from \cite{witt_multiple_2000} to use in dust attenuation calculations.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarSpectraDustAttenuationWittGordon2000(enumerationWittGordon2000ModelEncode(char(model),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function wittGordon2003ConstructorParameters

  function wittGordon2003ConstructorInternal(model) result(self)
    !!{
    Constructor for the \refClass{stellarSpectraDustAttenuationWittGordon2000} stellar spectra dust attenuation class.
    !!}
    use :: Array_Utilities                 , only : Array_Reverse
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : opticalDepthToMagnitudes
    use :: Numerical_Constants_Units       , only : micronsToAngstroms
    use :: Table_Labels                    , only : extrapolationTypeExtrapolate
    implicit none
    type(stellarSpectraDustAttenuationWittGordon2000)                :: self
    type(enumerationWittGordon2000ModelType         ), intent(in   ) :: model

    ! Initialize fitting function parameters for the chosen model.
    if (.not.enumerationWittGordon2000ModelIsValid(model)) call Error_Report('invalid model'//{introspection:location})
    call self%attenuationTable%create(micronsToAngstroms/Array_Reverse([ 1000.0d0,   1142.0d0,   1285.0d0,   1428.0d0,   1571.0d0,  1714.0d0,  1857.0d0,  2000.0d0,   2142.0d0,  2285.0d0,   2428.0d0,   2571.0d0,   2714.0d0,   2857.0d0,   3000.0d0,   3776.0d0,   4754.0d0,   5985.0d0,   7535.0d0,   9487.0d0,  11943.0d0,  15036.0d0,  18929.0d0,  23830.0d0,  30001.0d0]),tableCount=1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
     self%model=model
     select case (model%ID)
     case (wittGordon2000ModelMilkyWayShellTau3%ID)
        call self%attenuationTable%populate(opticalDepthToMagnitudes*Array_Reverse([ 15.714d0,  11.754d0,   9.546d0,   8.340d0,   7.752d0,   7.527d0,   7.683d0,   8.529d0,   9.570d0,   8.730d0,   7.416d0,   6.582d0,   6.066d0,   5.715d0,   5.454d0,   4.581d0,   3.597d0,   2.727d0,   2.001d0,   1.320d0,   0.912d0,   0.630d0,   0.435d0,   0.300d0,   0.207d0]))
     case (wittGordon2000ModelSMCShellTau3     %ID)
        call self%attenuationTable%populate(opticalDepthToMagnitudes*Array_Reverse([ 29.025d0,  22.320d0,  18.204d0,  15.501d0,  13.608d0,  12.222d0,  11.100d0,  10.137d0,   9.303d0,   8.571d0,   7.926d0,   7.356d0,   6.846d0,   6.399d0,   6.093d0,   4.830d0,   3.663d0,   2.640d0,   1.890d0,   1.290d0,   0.816d0,   0.498d0,   0.333d0,   0.225d0,   0.099d0]))
     end select
     return
  end function wittGordon2003ConstructorInternal
