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
   values of the model in one place. Both components take the mass of their dust from the same
   :galacticus-class:`dustPropertiesClass` object, so that changing the dust-to-metals ratio changes the optical depth
   of birth clouds and diffuse medium alike. Anything expressible here is expressible with
   :galacticus-class:`dustAttenuationSequence`; use that directly to vary the extinction curve of either component
   independently, or to build a model with more than two components.

   Its two phases of dust are the birth clouds and the diffuse interstellar medium, in that order, labeled
   ``birthCloud`` and ``screenSurfaceDensityMetals``, so that the energy each absorbs can be re-emitted separately.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationCharlotFall2000
     !!{RST
     The two-component dust attenuation model of :cite:t:`charlot_simple_2000`.
     !!}
     private
     type            (dustAttenuationBirthCloud                )          :: birthCloud_
     type            (dustAttenuationScreenSurfaceDensityMetals)          :: screenISM_
     class           (dustPropertiesClass                      ), pointer :: dustProperties_       => null()
     ! Retained so that the object can describe itself back into a parameter file; the physics is carried entirely by
     ! the two component attenuators above.
     double precision                                                     :: coefficientBirthCloud          , coefficientISM, &
          &                                                                  timescale                      , exponent_     , &
          &                                                                  wavelengthReference
   contains
     final     ::                       charlotFall2000Destructor
     procedure :: transmission       => charlotFall2000Transmission
     procedure :: request            => charlotFall2000Request
     procedure :: countPhases        => charlotFall2000CountPhases
     procedure :: labelPhase         => charlotFall2000LabelPhase
     procedure :: transmissionPhases => charlotFall2000TransmissionPhases
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
    class           (dustPropertiesClass           ), pointer       :: dustProperties_
    double precision                                                :: coefficientBirthCloud, coefficientISM, &
         &                                                             timescale            , exponent_     , &
         &                                                             wavelengthReference

    !![
    <inputParameter docformat="rst">
      <name>coefficientBirthCloud</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A dimensionless multiplicative coefficient applied to the gas column of birth clouds, in units of the column for
      which a cloud of local interstellar medium metallicity has unit :math:`V`-band optical depth with the default dust
      properties. With those defaults it is therefore the :math:`V`-band optical depth of such a cloud.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>coefficientISM</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A dimensionless multiplicative coefficient applied to the :math:`V`-band optical depth of the diffuse
      interstellar medium, representing the effects of geometry. It does not change the mass of dust.
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
      :math:`5500\,\text{Å}` to reproduce the ``lmnstyStllrCF2000`` property extractor exactly.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustProperties" name="dustProperties_" source="parameters"/>
    !!]
    self=dustAttenuationCharlotFall2000(coefficientBirthCloud,coefficientISM,timescale,exponent_,wavelengthReference,dustProperties_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustProperties_"/>
    !!]
    return
  end function charlotFall2000ConstructorParameters

  function charlotFall2000ConstructorInternal(coefficientBirthCloud,coefficientISM,timescale,exponent_,wavelengthReference,dustProperties_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationCharlotFall2000` dust attenuation class. Both
    components are given their own power-law extinction curve of the same exponent.
    !!}
    implicit none
    type            (dustAttenuationCharlotFall2000)                        :: self
    double precision                                , intent(in   )         :: coefficientBirthCloud, coefficientISM, &
         &                                                                     timescale            , exponent_     , &
         &                                                                     wavelengthReference
    class           (dustPropertiesClass           ), intent(in   ), target :: dustProperties_
    class           (dustExtinctionCurveClass      ), pointer               :: curvePowerLaw
    !![
    <constructorAssign variables="coefficientBirthCloud, coefficientISM, timescale, exponent_, wavelengthReference, *dustProperties_"/>
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
    self%birthCloud_=dustAttenuationBirthCloud                (coefficientBirthCloud*densitySurfaceGasDepthOpticalVUnitMilkyWay,timescale,curvePowerLaw,dustProperties_)
    self%screenISM_ =dustAttenuationScreenSurfaceDensityMetals(coefficientISM                                                   ,curvePowerLaw,dustProperties_)
    !![
    <objectDestructor name="curvePowerLaw"/>
    !!]
    return
  end function charlotFall2000ConstructorInternal

  subroutine charlotFall2000Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationCharlotFall2000` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationCharlotFall2000), intent(inout) :: self

    !![
    <objectDestructor name="self%dustProperties_"/>
    !!]
    return
  end subroutine charlotFall2000Destructor

  function charlotFall2000Transmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through both the birth cloud and diffuse interstellar medium components.
    !!}
    implicit none
    class           (dustAttenuationCharlotFall2000), intent(inout)                               :: self
    type            (treeNode                      ), intent(inout), target                       :: node
    type            (emissionDescriptor            ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                               , dimension(size(descriptors)) :: transmission
    double precision                                , intent(in   ), optional                     :: inclination
    !$GLC attributes unused :: inclination

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

  integer function charlotFall2000CountPhases(self) result(countPhases)
    !!{RST
    Return the number of phases of dust: birth clouds and the diffuse interstellar medium.
    !!}
    implicit none
    class(dustAttenuationCharlotFall2000), intent(inout) :: self
    !$GLC attributes unused :: self

    countPhases=2
    return
  end function charlotFall2000CountPhases

  function charlotFall2000LabelPhase(self,indexPhase) result(label)
    !!{RST
    Return the label of a phase of dust: that of the birth cloud component for the first, and of the diffuse
    interstellar medium for the second.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (varying_string                )                :: label
    class  (dustAttenuationCharlotFall2000), intent(inout) :: self
    integer                                , intent(in   ) :: indexPhase

    select case (indexPhase)
    case (1)
       label=self%birthCloud_%labelPhase(1)
    case (2)
       label=self%screenISM_ %labelPhase(1)
    case default
       label=''
       call Error_Report('phase index out of range'//{introspection:location})
    end select
    return
  end function charlotFall2000LabelPhase

  function charlotFall2000TransmissionPhases(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through the birth clouds and through the diffuse interstellar medium, in that order.
    !!}
    implicit none
    double precision                                , allocatable  , dimension(:,:) :: transmission
    class           (dustAttenuationCharlotFall2000), intent(inout)                 :: self
    type            (treeNode                      ), intent(inout), target         :: node
    type            (emissionDescriptor            ), intent(in   ), dimension(:  ) :: descriptors
    double precision                                , intent(in   ), optional       :: inclination
    !$GLC attributes unused :: inclination

    allocate(transmission(size(descriptors),2))
    transmission(:,1)=self%birthCloud_%transmission(node,descriptors)
    transmission(:,2)=self%screenISM_ %transmission(node,descriptors)
    return
  end function charlotFall2000TransmissionPhases
