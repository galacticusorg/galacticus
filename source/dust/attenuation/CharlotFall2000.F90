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
  Implements the two-component dust attenuation model of :cite:t:`charlot_simple_2000`.
  !!}

  use :: Dust_Extinction_Curves, only : dustExtinctionCurvePowerLaw, wavelengthVBand

  !![
  <dustAttenuation name="dustAttenuationCharlotFall2000" docformat="rst">
   <description>
   The two-component dust attenuation model of :cite:t:`charlot_simple_2000`: young stars are attenuated both by the
   dust of the birth clouds in which they remain embedded and by the diffuse dust of the interstellar medium, while
   older stars, having escaped their birth clouds, are attenuated only by the diffuse component.

   This class carries no physics of its own. It is exactly

   .. code-block:: xml

      &lt;dustAttenuation value="sequence"&gt;
        &lt;dustAttenuation value="birthCloud"&gt;
          &lt;dustExtinctionCurve value="powerLaw"&gt;&lt;exponent value="0.7"/&gt;&lt;/dustExtinctionCurve&gt;
        &lt;/dustAttenuation&gt;
        &lt;dustAttenuation value="screenSurfaceDensityMetals"&gt;
          &lt;dustExtinctionCurve value="powerLaw"&gt;&lt;exponent value="0.7"/&gt;&lt;/dustExtinctionCurve&gt;
        &lt;/dustAttenuation&gt;
      &lt;/dustAttenuation&gt;

   and is provided because that is the combination users most often want and because it pins the canonical parameter
   values of the model in one place. Anything expressible here is expressible with
   :galacticus-class:`dustAttenuationSequence`; use that directly to vary the extinction curve of either component
   independently, or to build a model with more than two components.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationCharlotFall2000
     !!{RST
     The two-component dust attenuation model of :cite:t:`charlot_simple_2000`.
     !!}
     private
     type            (dustAttenuationBirthCloud                ) :: birthCloud_
     type            (dustAttenuationScreenSurfaceDensityMetals) :: screenISM_
     ! Retained so that the object can describe itself back into a parameter file; the physics is carried entirely by
     ! the two component attenuators above.
     double precision                                            :: coefficientBirthCloud, coefficientISM, &
          &                                                         timescale            , exponent_     , &
          &                                                         wavelengthReference
   contains
     procedure :: transmission => charlotFall2000Transmission
     procedure :: request      => charlotFall2000Request
  end type dustAttenuationCharlotFall2000

  interface dustAttenuationCharlotFall2000
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationCharlotFall2000` dust attenuation class.
     !!}
     module procedure charlotFall2000ConstructorParameters
     module procedure charlotFall2000ConstructorInternal
  end interface dustAttenuationCharlotFall2000

contains

  function charlotFall2000ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationCharlotFall2000` dust attenuation class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationCharlotFall2000)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: coefficientBirthCloud, coefficientISM, &
         &                                                             timescale            , exponent_     , &
         &                                                             wavelengthReference

    !![
    <inputParameter docformat="rst">
      <name>coefficientBirthCloud</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The :math:`V`-band optical depth of a birth cloud of local interstellar medium metallicity.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>coefficientISM</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A dimensionless multiplicative coefficient applied to the :math:`V`-band optical depth of the diffuse
      interstellar medium.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>timescale</name>
      <defaultValue>1.0d-2</defaultValue>
      <defaultSource>:cite:t:`charlot_simple_2000`</defaultSource>
      <description>
      The lifetime of a stellar birth cloud, in Gyr.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponent</name>
      <variable>exponent_</variable>
      <defaultValue>0.7d0</defaultValue>
      <defaultSource>:cite:t:`charlot_simple_2000`</defaultSource>
      <description>
      The exponent of the power-law extinction curve applied to both components.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthReference</name>
      <defaultValue>wavelengthVBand</defaultValue>
      <description>
      The wavelength, in Å, at which the power-law extinction curve is normalized to unity. Set it to
      :math:`5500\,`Å to reproduce the ``lmnstyStllrCF2000`` property extractor exactly.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustAttenuationCharlotFall2000(coefficientBirthCloud,coefficientISM,timescale,exponent_,wavelengthReference)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function charlotFall2000ConstructorParameters

  function charlotFall2000ConstructorInternal(coefficientBirthCloud,coefficientISM,timescale,exponent_,wavelengthReference) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationCharlotFall2000` dust attenuation class. Both
    components are given their own power-law extinction curve of the same exponent.
    !!}
    implicit none
    type            (dustAttenuationCharlotFall2000)                :: self
    double precision                                , intent(in   ) :: coefficientBirthCloud, coefficientISM, &
         &                                                             timescale            , exponent_     , &
         &                                                             wavelengthReference
    class(dustExtinctionCurveClass), pointer :: curvePowerLaw
    !![
    <constructorAssign variables="coefficientBirthCloud, coefficientISM, timescale, exponent_, wavelengthReference"/>
    !!]

    ! The curve is shared by both components, so it must be a reference-counted heap object rather than a local: each
    ! component takes a reference to it, and a stack object would be destroyed while those references were still held.
    ! A freshly allocated object carries no references, so one is established here for the reference this constructor
    ! itself holds, and released once both components have taken theirs.
    allocate(dustExtinctionCurvePowerLaw :: curvePowerLaw)
    select type (curvePowerLaw)
    type is (dustExtinctionCurvePowerLaw)
       curvePowerLaw=dustExtinctionCurvePowerLaw(exponent_,wavelengthReference)
    end select
    call curvePowerLaw%referenceCountReset()
    self%birthCloud_=dustAttenuationBirthCloud                (coefficientBirthCloud,timescale,curvePowerLaw)
    self%screenISM_ =dustAttenuationScreenSurfaceDensityMetals(coefficientISM                 ,curvePowerLaw)
    !![
    <objectDestructor name="curvePowerLaw"/>
    !!]
    return
  end function charlotFall2000ConstructorInternal

  function charlotFall2000Transmission(self,node,descriptors) result(transmission)
    !!{RST
    Return the transmission through both the birth cloud and diffuse interstellar medium components.
    !!}
    implicit none
    class           (dustAttenuationCharlotFall2000), intent(inout)                               :: self
    type            (treeNode                      ), intent(inout), target                       :: node
    type            (emissionDescriptor            ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                               , dimension(size(descriptors)) :: transmission

    transmission=+self%birthCloud_%transmission(node,descriptors) &
         &       *self%screenISM_ %transmission(node,descriptors)
    return
  end function charlotFall2000Transmission

  function charlotFall2000Request(self) result(request)
    !!{RST
    Return a decomposition request splitting populations at the birth cloud lifetime, and by component.
    !!}
    implicit none
    type (decompositionRequest          )                :: request
    class(dustAttenuationCharlotFall2000), intent(inout) :: self

    request=self%birthCloud_%request()
    return
  end function charlotFall2000Request
