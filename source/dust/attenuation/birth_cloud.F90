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
  Implements attenuation of young stellar populations by the dust of their birth clouds.
  !!}

  !![
  <dustAttenuation name="dustAttenuationBirthCloud" docformat="rst">
   <description>
   Attenuation of stellar populations which remain embedded in the dense clouds from which they formed. Stars younger
   than ``timescale`` are attenuated by a screen of :math:`V`-band optical depth

   .. math::

      \tau_\mathrm{V} = C \frac{Z}{Z_\mathrm{ISM}},

   proportional to the metallicity :math:`Z` of the gas of the component in which they formed, relative to that of the
   local interstellar medium; stars older than ``timescale`` are not attenuated at all by this class. This is the
   birth cloud term of the :cite:t:`charlot_simple_2000` model, in which the coefficient :math:`C` (``coefficient``)
   is the optical depth of a birth cloud of local interstellar medium metallicity.

   Note that the optical depth depends only on metallicity, not on the surface density of the gas: a birth cloud is a
   local structure, and its column is not set by the global structure of the galaxy.

   A parcel of emission whose age is unresolved spans all ages, and is therefore treated as old and left unattenuated
   rather than being attenuated as though it were entirely young. Combine this class with a
   :galacticus-class:`dustAttenuationScreenClass` through :galacticus-class:`dustAttenuationSequence` to obtain the
   full two-component model.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationBirthCloud
     !!{RST
     Attenuation of young stellar populations by their birth clouds.
     !!}
     private
     class           (dustExtinctionCurveClass), pointer :: dustExtinctionCurve_ => null()
     double precision                                    :: coefficient                   , timescale
   contains
     final     ::                 birthCloudDestructor
     procedure :: transmission => birthCloudTransmission
     procedure :: request      => birthCloudRequest
  end type dustAttenuationBirthCloud

  interface dustAttenuationBirthCloud
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationBirthCloud` dust attenuation class.
     !!}
     module procedure birthCloudConstructorParameters
     module procedure birthCloudConstructorInternal
  end interface dustAttenuationBirthCloud

contains

  function birthCloudConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationBirthCloud` dust attenuation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationBirthCloud)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (dustExtinctionCurveClass ), pointer       :: dustExtinctionCurve_
    double precision                                           :: coefficient         , timescale

    !![
    <inputParameter docformat="rst">
      <name>coefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The :math:`V`-band optical depth of a birth cloud of local interstellar medium metallicity.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>timescale</name>
      <defaultValue>1.0d-2</defaultValue>
      <defaultSource>:cite:t:`charlot_simple_2000`</defaultSource>
      <description>
      The lifetime of a stellar birth cloud, in Gyr. Stellar populations younger than this are attenuated by their
      birth cloud; older populations are not.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
    !!]
    self=dustAttenuationBirthCloud(coefficient,timescale,dustExtinctionCurve_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"/>
    !!]
    return
  end function birthCloudConstructorParameters

  function birthCloudConstructorInternal(coefficient,timescale,dustExtinctionCurve_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationBirthCloud` dust attenuation class.
    !!}
    implicit none
    type            (dustAttenuationBirthCloud)                        :: self
    double precision                           , intent(in   )         :: coefficient         , timescale
    class           (dustExtinctionCurveClass ), intent(in   ), target :: dustExtinctionCurve_
    !![
    <constructorAssign variables="coefficient, timescale, *dustExtinctionCurve_"/>
    !!]

    return
  end function birthCloudConstructorInternal

  subroutine birthCloudDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationBirthCloud` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationBirthCloud), intent(inout) :: self

    !![
    <objectDestructor name="self%dustExtinctionCurve_"/>
    !!]
    return
  end subroutine birthCloudDestructor

  function birthCloudTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through the dust of a stellar birth cloud.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeMax, componentTypeMin
    implicit none
    class           (dustAttenuationBirthCloud), intent(inout)                                               :: self
    type            (treeNode                 ), intent(inout), target                                       :: node
    type            (emissionDescriptor       ), intent(in   ), dimension(:                                ) :: descriptors
    double precision                                    , intent(in   ), optional                            :: inclination
    double precision                                          , dimension(size(descriptors)                ) :: transmission
    double precision                                          , dimension(componentTypeMin:componentTypeMax) :: depthOpticalV
    logical                                                   , dimension(componentTypeMin:componentTypeMax) :: computed
    integer                                                                                                  :: i            , indexComponent
    double precision                                                                                         :: massGas      , radius        , &
         &                                                                                                      metallicity
    !$GLC attributes unused :: inclination

    computed=.false.
    do i=1,size(descriptors)
       ! Populations older than the birth cloud lifetime have escaped their birth clouds. A parcel whose age is
       ! unresolved spans all ages and so falls here, and is left unattenuated.
       if (descriptors(i)%ageMaximum > self%timescale) then
          transmission(i)=1.0d0
          cycle
       end if
       indexComponent=descriptors(i)%componentType%ID
       if (.not.computed(indexComponent)) then
          call componentGasProperties(node,descriptors(i)%componentType,massGas,radius,metallicity)
          depthOpticalV(indexComponent)=+self%coefficient    &
               &                        *metallicity         &
               &                        /metallicityISMLocal
          computed     (indexComponent)=.true.
       end if
       transmission(i)=exp(                                                                          &
            &              -depthOpticalV(indexComponent)                                            &
            &              *self%dustExtinctionCurve_%attenuationRelative(descriptors(i)%wavelength) &
            &             )
    end do
    return
  end function birthCloudTransmission

  function birthCloudRequest(self) result(request)
    !!{RST
    Return a decomposition request splitting populations at the birth cloud lifetime.
    !!}
    implicit none
    type (decompositionRequest      )                :: request
    class(dustAttenuationBirthCloud ), intent(inout) :: self

    request%ageBoundaries    =[self%timescale]
    request%resolveComponents=.true.
    return
  end function birthCloudRequest
