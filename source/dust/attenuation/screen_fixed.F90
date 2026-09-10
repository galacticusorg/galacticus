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
  Implements a dust screen of fixed optical depth.
  !!}

  !![
  <dustAttenuation name="dustAttenuationScreenFixed" docformat="rst">
   <description>
   A uniform dust screen whose :math:`V`-band optical depth is a fixed parameter, independent of the properties of the
   galaxy. This is the appropriate choice when the optical depth is to be treated as a free parameter---when fitting
   an observed spectral energy distribution, for example---rather than predicted from a model of the interstellar
   medium.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationScreen) :: dustAttenuationScreenFixed
     !!{RST
     A dust screen of fixed optical depth.
     !!}
     private
     double precision :: depthOpticalV_
   contains
     procedure :: depthOpticalV => screenFixedDepthOpticalV
  end type dustAttenuationScreenFixed

  interface dustAttenuationScreenFixed
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationScreenFixed` dust attenuation class.
     !!}
     module procedure screenFixedConstructorParameters
     module procedure screenFixedConstructorInternal
  end interface dustAttenuationScreenFixed

contains

  function screenFixedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationScreenFixed` dust attenuation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationScreenFixed)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (dustExtinctionCurveClass  ), pointer       :: dustExtinctionCurve_
    double precision                                            :: depthOpticalV

    !![
    <inputParameter docformat="rst">
      <name>depthOpticalV</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The :math:`V`-band optical depth of the dust screen.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
    !!]
    self=dustAttenuationScreenFixed(depthOpticalV,dustExtinctionCurve_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"/>
    !!]
    return
  end function screenFixedConstructorParameters

  function screenFixedConstructorInternal(depthOpticalV_,dustExtinctionCurve_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationScreenFixed` dust attenuation class.
    !!}
    implicit none
    type            (dustAttenuationScreenFixed)                        :: self
    double precision                            , intent(in   )         :: depthOpticalV_
    class           (dustExtinctionCurveClass  ), intent(in   ), target :: dustExtinctionCurve_
    !![
    <constructorAssign variables="depthOpticalV_, *dustExtinctionCurve_"/>
    !!]

    return
  end function screenFixedConstructorInternal

  double precision function screenFixedDepthOpticalV(self,node,componentType) result(depthOpticalV)
    !!{RST
    Return the fixed :math:`V`-band optical depth of the screen.
    !!}
    implicit none
    class(dustAttenuationScreenFixed  ), intent(inout)         :: self
    type (treeNode                    ), intent(inout), target :: node
    type (enumerationComponentTypeType), intent(in   )         :: componentType
    !$GLC attributes unused :: node, componentType

    depthOpticalV=self%depthOpticalV_
    return
  end function screenFixedDepthOpticalV
