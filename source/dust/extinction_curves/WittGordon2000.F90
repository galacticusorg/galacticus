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
  Implements the dust attenuation curves of :cite:t:`witt_multiple_2000`.
  !!}

  !![
  <enumeration docformat="rst">
   <name>wittGordon2000Model</name>
   <description>
   Enumerates the models available in the ``wittGordon2000`` dust extinction curve class.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="milkyWayShellTau3"/>
   <entry label="SMCShellTau3"     />
  </enumeration>
  !!]

  !![
  <dustExtinctionCurve name="dustExtinctionCurveWittGordon2000" docformat="rst">
   <description>
   The dust curves of :cite:t:`witt_multiple_2000`, for a shell geometry of :math:`V`-band optical depth 3 with either
   Milky Way or Small Magellanic Cloud dust.

   Note that these are *attenuation* curves rather than extinction curves: they are the outcome of a radiative
   transfer calculation through a specified geometry, including scattering back into the line of sight, and so already
   embed an assumption about how dust and stars are distributed. Combining one with a ``dustAttenuation`` object that
   itself models a geometry therefore double counts that geometry. They are provided here for continuity with the
   models this class replaces.
   </description>
  </dustExtinctionCurve>
  !!]
  type, extends(dustExtinctionCurveTabulated) :: dustExtinctionCurveWittGordon2000
     !!{RST
     The dust attenuation curves of :cite:t:`witt_multiple_2000`.
     !!}
     private
     type(enumerationWittGordon2000ModelType) :: model
  end type dustExtinctionCurveWittGordon2000

  interface dustExtinctionCurveWittGordon2000
     !!{RST
     Constructors for the :galacticus-class:`dustExtinctionCurveWittGordon2000` dust extinction curve class.
     !!}
     module procedure wittGordon2000ConstructorParameters
     module procedure wittGordon2000ConstructorInternal
  end interface dustExtinctionCurveWittGordon2000

contains

  function wittGordon2000ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustExtinctionCurveWittGordon2000` dust extinction curve class which takes
    a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(dustExtinctionCurveWittGordon2000)                :: self
    type(inputParameters                  ), intent(inout) :: parameters
    type(varying_string                   )                :: model

    !![
    <inputParameter docformat="rst">
      <name>model</name>
      <defaultValue>var_str('MilkyWayShellTau3.0')</defaultValue>
      <description>
      The name of the model from :cite:t:`witt_multiple_2000` to use.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustExtinctionCurveWittGordon2000(enumerationWittGordon2000ModelEncode(char(model),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function wittGordon2000ConstructorParameters

  function wittGordon2000ConstructorInternal(model) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustExtinctionCurveWittGordon2000` dust extinction curve class.
    !!}
    use :: Array_Utilities                 , only : Array_Reverse
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : opticalDepthToMagnitudes
    use :: Numerical_Constants_Units       , only : micronsToAngstroms
    use :: Table_Labels                    , only : extrapolationTypeExtrapolate
    implicit none
    type(dustExtinctionCurveWittGordon2000  )                :: self
    type(enumerationWittGordon2000ModelType ), intent(in   ) :: model
    !![
    <constructorAssign variables="model"/>
    !!]

    if (.not.enumerationWittGordon2000ModelIsValid(model)) call Error_Report('invalid model'//{introspection:location})
    call self%attenuationTable%create(micronsToAngstroms/Array_Reverse([ 1000.0d0,   1142.0d0,   1285.0d0,   1428.0d0,   1571.0d0,  1714.0d0,  1857.0d0,  2000.0d0,   2142.0d0,  2285.0d0,   2428.0d0,   2571.0d0,   2714.0d0,   2857.0d0,   3000.0d0,   3776.0d0,   4754.0d0,   5985.0d0,   7535.0d0,   9487.0d0,  11943.0d0,  15036.0d0,  18929.0d0,  23830.0d0,  30001.0d0]),tableCount=1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
    select case (model%ID)
    case (wittGordon2000ModelMilkyWayShellTau3%ID)
       call self%attenuationTable%populate(opticalDepthToMagnitudes*Array_Reverse([ 15.714d0,  11.754d0,   9.546d0,   8.340d0,   7.752d0,   7.527d0,   7.683d0,   8.529d0,   9.570d0,   8.730d0,   7.416d0,   6.582d0,   6.066d0,   5.715d0,   5.454d0,   4.581d0,   3.597d0,   2.727d0,   2.001d0,   1.320d0,   0.912d0,   0.630d0,   0.435d0,   0.300d0,   0.207d0]))
    case (wittGordon2000ModelSMCShellTau3     %ID)
       call self%attenuationTable%populate(opticalDepthToMagnitudes*Array_Reverse([ 29.025d0,  22.320d0,  18.204d0,  15.501d0,  13.608d0,  12.222d0,  11.100d0,  10.137d0,   9.303d0,   8.571d0,   7.926d0,   7.356d0,   6.846d0,   6.399d0,   6.093d0,   4.830d0,   3.663d0,   2.640d0,   1.890d0,   1.290d0,   0.816d0,   0.498d0,   0.333d0,   0.225d0,   0.099d0]))
    end select
    return
  end function wittGordon2000ConstructorInternal
