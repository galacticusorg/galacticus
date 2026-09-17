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
  Implements a property extractor for the spectral energy distribution of thermal emission from dust.
  !!}

  use :: Cosmology_Functions           , only : cosmologyFunctionsClass
  use :: Dust_Attenuations             , only : dustAttenuationClass     , gaussLegendreRule
  use :: Dust_Emission_Spectra         , only : dustEmissionSpectrumClass, dustEmissionSpectrumList
  use :: Dust_Properties               , only : dustPropertiesClass
  use :: Stellar_Luminosities_Structure, only : enumerationFrameType

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSEDDustEmission" docformat="rst">
   <description>
   A property extractor which returns the spectral energy distribution of thermal emission from the dust of a galaxy,
   :math:`L_\nu` in :math:`L_\odot\,\hbox{Hz}^{-1}`, on a grid of wavelengths from ``[wavelengthMinimum]`` to
   ``[wavelengthMaximum]`` at resolution :math:`\lambda/\Delta\lambda` of ``[resolution]``, in the ``[frame]`` (``rest``
   or ``observed``) frame. The grid and its conventions are those of :galacticus-class:`nodePropertyExtractorSED`, so
   that the two may be summed directly.

   The dust is heated by the light of the ``[nodePropertyExtractor]``---normally a
   :galacticus-class:`nodePropertyExtractorMulti` holding :galacticus-class:`nodePropertyExtractorSED` extractors for
   continuum and emission line extractors for lines, although a single extractor may be given directly. Each child is
   decomposed into parcels of emission and attenuated by the ``[dustAttenuation]`` object exactly as by
   :galacticus-class:`nodePropertyExtractorDustAttenuation`, which
   gives the luminosity absorbed by each phase of dust, averaged over orientation where the attenuator depends on it
   (with Gauss-Legendre quadrature of order ``[orderInclinationAverage]``). Absorbed spectra are integrated over
   frequency on each child's own wavelengths, so each spectral child must be computed in the rest frame, and must span
   at least ``[wavelengthHeatingMinimum]`` to ``[wavelengthHeatingMaximum]``---by default the Lyman limit to
   :math:`3\,\mu\hbox{m}`---both of which are checked at construction. The accuracy of the absorbed luminosity follows
   that child's resolution. Line luminosities are converted to :math:`L_\odot` using their units. Only the lines a child
   is asked for heat the dust, and nebular continuum is not included.

   Absorbed luminosity is summed over components into one luminosity for each phase of dust for the whole galaxy. That
   heats dust whose mass is the fraction ``[fractionsMassPhase]`` of the total dust mass of the galaxy, which is the sum
   over its components given by the ``[dustProperties]`` object. By default all of the dust mass is given to the last
   phase---the diffuse interstellar medium of :galacticus-class:`dustAttenuationCharlotFall2000`, for example---so
   earlier phases, such as birth clouds, should be given spectra whose temperature is fixed rather than set by energy
   balance. One ``[dustEmissionSpectrum]`` must be given for each phase, in the order of the phases, and each re-emits
   the luminosity its phase absorbs.

   Within each wavelength bin the emission is averaged over the resolution element as by
   :galacticus-class:`nodePropertyExtractorSED`. In the ``observed`` frame the wavelengths are observed-frame
   wavelengths, the emission is evaluated at the corresponding rest-frame wavelengths, and :math:`L_\nu` is not
   rescaled, again as by :galacticus-class:`nodePropertyExtractorSED`.

   The property is named ``dustEmissionSED:`` followed by the name of the attenuator and any ``[appendSuffix]``. Setting
   ``[outputPhases]`` also emits the emission of each phase separately, named with the label of the phase appended.
   Infrared emission is assumed not to be absorbed again.
   </description>
   <linkedList type="dustEmissionSpectrumList" variable="dustEmissionSpectra" next="next" object="dustEmissionSpectrum_" objectType="dustEmissionSpectrumClass"/>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSEDDustEmission
     !!{RST
     A property extractor for the spectral energy distribution of thermal emission from dust.
     !!}
     private
     class           (nodePropertyExtractorClass), pointer                   :: nodePropertyExtractor_   => null()
     class           (dustAttenuationClass      ), pointer                   :: dustAttenuation_         => null()
     class           (dustPropertiesClass       ), pointer                   :: dustProperties_          => null()
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_      => null()
     type            (dustEmissionSpectrumList  ), pointer                   :: dustEmissionSpectra      => null()
     type            (enumerationFrameType      )                            :: frame
     type            (varying_string            )                            :: appendSuffix
     logical                                                                 :: outputPhases
     integer                                                                 :: orderInclinationAverage
     double precision                                                        :: wavelengthMinimum                 , wavelengthMaximum       , &
          &                                                                     resolution                        , factorWavelength        , &
          &                                                                     wavelengthHeatingMinimum          , wavelengthHeatingMaximum
     ! The fraction of the dust mass heated by each phase; the quadrature rule used to average absorption over orientation;
     ! and the quadrature rule, on the unit interval in the logarithm of wavelength, used to average over a resolution
     ! element.
     double precision                            , allocatable, dimension(:) :: fractionsMassPhase                , cosineInclination       , &
          &                                                                     weight                            , abscissaeBin            , &
          &                                                                     weightsBin
   contains
     !![
     <methods docformat="rst">
       <method method="wavelengths" description="Return the wavelengths, in Å, of the centers of the bins of the spectrum, in the frame in which it is computed."/>
       <method method="countBins"   description="Return the number of bins in the spectrum."                                                                         />
       <method method="countElements" description="Return the number of properties emitted."                                                                         />
     </methods>
     !!]
     final     ::                       sedDustEmissionDestructor
     procedure :: wavelengths        => sedDustEmissionWavelengths
     procedure :: countBins          => sedDustEmissionCountBins
     procedure :: countElements      => sedDustEmissionCountElements
     procedure :: columnDescriptions => sedDustEmissionColumnDescriptions
     procedure :: size               => sedDustEmissionSize
     procedure :: elementCount       => sedDustEmissionElementCount
     procedure :: extract            => sedDustEmissionExtract
     procedure :: names              => sedDustEmissionNames
     procedure :: descriptions       => sedDustEmissionDescriptions
     procedure :: unitsInSI          => sedDustEmissionUnitsInSI
     procedure :: units              => sedDustEmissionUnits
  end type nodePropertyExtractorSEDDustEmission

  interface nodePropertyExtractorSEDDustEmission
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorSEDDustEmission` property extractor class.
     !!}
     module procedure sedDustEmissionConstructorParameters
     module procedure sedDustEmissionConstructorInternal
  end interface nodePropertyExtractorSEDDustEmission

  ! Order of the Gauss-Legendre rule used to average the emission over a resolution element. The emission of dust varies
  ! smoothly on the scale of a resolution element, so a low order suffices.
  integer, parameter :: sedDustEmissionOrderBin=8

contains

  function sedDustEmissionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorSEDDustEmission` property extractor class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters              , only : inputParameter        , inputParameters
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode
    implicit none
    type            (nodePropertyExtractorSEDDustEmission)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass          ), pointer                     :: nodePropertyExtractor_
    class           (dustAttenuationClass                ), pointer                     :: dustAttenuation_
    class           (dustPropertiesClass                 ), pointer                     :: dustProperties_
    class           (cosmologyFunctionsClass             ), pointer                     :: cosmologyFunctions_
    type            (dustEmissionSpectrumList            ), pointer                     :: dustEmissionSpectrum_
    double precision                                      , allocatable  , dimension(:) :: fractionsMassPhase
    type            (varying_string                      )                              :: frame                   , appendSuffix
    double precision                                                                    :: wavelengthMinimum       , wavelengthMaximum       , &
         &                                                                                 resolution              , wavelengthHeatingMinimum, &
         &                                                                                 wavelengthHeatingMaximum
    integer                                                                             :: orderInclinationAverage , i
    logical                                                                             :: outputPhases

    !![
    <inputParameter docformat="rst">
      <name>frame</name>
      <defaultValue>var_str('rest')</defaultValue>
      <description>
      The frame (``rest`` or ``observed``) for which to compute the spectrum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthMinimum</name>
      <defaultValue>1.0d4</defaultValue>
      <description>
      The minimum wavelength, in Å, of the spectrum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthMaximum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The maximum wavelength, in Å, of the spectrum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>resolution</name>
      <defaultValue>100.0d0</defaultValue>
      <description>
      The resolution, :math:`\lambda/\Delta\lambda`, of the spectrum. It must be positive.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthHeatingMinimum</name>
      <defaultValue>912.0d0</defaultValue>
      <description>
      The shortest wavelength, in Å, which each spectral child must reach, so that all of the light heating the dust is
      counted. The default is the Lyman limit: shortward of it, light is absorbed by gas.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthHeatingMaximum</name>
      <defaultValue>3.0d4</defaultValue>
      <description>
      The longest wavelength, in Å, which each spectral child must reach, so that all of the light heating the dust is
      counted.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>orderInclinationAverage</name>
      <defaultValue>8</defaultValue>
      <description>
      The order of the Gauss-Legendre quadrature used to average absorbed luminosity over orientation, where the
      attenuator depends on orientation.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputPhases</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, the emission of each phase of dust is emitted separately, as well as the total.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>appendSuffix</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      An extra label appended to the names of the properties emitted, after the name of the attenuator. It must be set
      where two of these extractors use the same class of attenuator. ``none`` appends nothing.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    <objectBuilder class="dustAttenuation"       name="dustAttenuation_"       source="parameters"/>
    <objectBuilder class="dustProperties"        name="dustProperties_"        source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !!]
    if (parameters%isPresent('fractionsMassPhase')) then
       allocate(fractionsMassPhase(parameters%count('fractionsMassPhase')))
       !![
       <inputParameter docformat="rst">
         <name>fractionsMassPhase</name>
         <description>
         The fraction of the total dust mass of the galaxy heated by each phase of dust, in the order of the phases.
         These must sum to one. If absent, all of the dust mass is given to the last phase.
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       allocate(fractionsMassPhase(dustAttenuation_%countPhases()))
       fractionsMassPhase=0.0d0
       if (size(fractionsMassPhase) > 0) fractionsMassPhase(size(fractionsMassPhase))=1.0d0
    end if
    ! Build the emission spectra, one for each phase of dust.
    self                 %dustEmissionSpectra => null()
    dustEmissionSpectrum_                     => null()
    do i=1,parameters%copiesCount('dustEmissionSpectrum',zeroIfNotPresent=.true.)
       if (associated(dustEmissionSpectrum_)) then
          allocate(dustEmissionSpectrum_%next)
          dustEmissionSpectrum_ => dustEmissionSpectrum_%next
       else
          allocate(self%dustEmissionSpectra)
          dustEmissionSpectrum_ => self%dustEmissionSpectra
       end if
       !![
       <objectBuilder class="dustEmissionSpectrum" name="dustEmissionSpectrum_%dustEmissionSpectrum_" source="parameters" copy="i" />
       !!]
    end do
    self%nodePropertyExtractor_ => nodePropertyExtractor_
    self%dustAttenuation_       => dustAttenuation_
    self%dustProperties_        => dustProperties_
    self%cosmologyFunctions_    => cosmologyFunctions_
    !![
    <referenceCountIncrement owner="self" object="nodePropertyExtractor_"/>
    <referenceCountIncrement owner="self" object="dustAttenuation_"      />
    <referenceCountIncrement owner="self" object="dustProperties_"       />
    <referenceCountIncrement owner="self" object="cosmologyFunctions_"   />
    !!]
    self%frame                   =enumerationFrameEncode(char(frame),includesPrefix=.false.)
    self%wavelengthMinimum       =wavelengthMinimum
    self%wavelengthMaximum       =wavelengthMaximum
    self%resolution              =resolution
    self%wavelengthHeatingMinimum=wavelengthHeatingMinimum
    self%wavelengthHeatingMaximum=wavelengthHeatingMaximum
    self%orderInclinationAverage =orderInclinationAverage
    self%outputPhases            =outputPhases
    self%appendSuffix            =appendSuffix
    self%fractionsMassPhase      =fractionsMassPhase
    call sedDustEmissionInitialize(self)
    !![
    <inputParametersValidate source="parameters" multiParameters="dustEmissionSpectrum"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    <objectDestructor name="dustAttenuation_"   />
    <objectDestructor name="dustProperties_"    />
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function sedDustEmissionConstructorParameters

  function sedDustEmissionConstructorInternal(frame,wavelengthMinimum,wavelengthMaximum,resolution,wavelengthHeatingMinimum,wavelengthHeatingMaximum,orderInclinationAverage,outputPhases,appendSuffix,fractionsMassPhase,nodePropertyExtractor_,dustAttenuation_,dustEmissionSpectra,dustProperties_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorSEDDustEmission` property extractor class.
    ``nodePropertyExtractor_`` supplies the light heating the dust: either a :galacticus-class:`nodePropertyExtractorMulti`
    holding the child extractors, or a single child extractor.
    !!}
    implicit none
    type            (nodePropertyExtractorSEDDustEmission)                              :: self
    type            (enumerationFrameType                ), intent(in   )               :: frame
    double precision                                      , intent(in   )               :: wavelengthMinimum       , wavelengthMaximum       , &
         &                                                                                 resolution              , wavelengthHeatingMinimum, &
         &                                                                                 wavelengthHeatingMaximum
    integer                                               , intent(in   )               :: orderInclinationAverage
    logical                                               , intent(in   )               :: outputPhases
    type            (varying_string                      ), intent(in   )               :: appendSuffix
    double precision                                      , intent(in   ), dimension(:) :: fractionsMassPhase
    class           (nodePropertyExtractorClass          ), intent(in   ), target       :: nodePropertyExtractor_
    class           (dustAttenuationClass                ), intent(in   ), target       :: dustAttenuation_
    type            (dustEmissionSpectrumList            ), intent(in   ), target       :: dustEmissionSpectra
    class           (dustPropertiesClass                 ), intent(in   ), target       :: dustProperties_
    class           (cosmologyFunctionsClass             ), intent(in   ), target       :: cosmologyFunctions_
    type            (dustEmissionSpectrumList            ), pointer                     :: dustEmissionSpectrum_
    !![
    <constructorAssign variables="frame, wavelengthMinimum, wavelengthMaximum, resolution, wavelengthHeatingMinimum, wavelengthHeatingMaximum, orderInclinationAverage, outputPhases, appendSuffix, fractionsMassPhase, *nodePropertyExtractor_, *dustAttenuation_, *dustProperties_, *cosmologyFunctions_"/>
    !!]

    self                 %dustEmissionSpectra => dustEmissionSpectra
    dustEmissionSpectrum_                     => dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       !![
       <referenceCountIncrement owner="dustEmissionSpectrum_" object="dustEmissionSpectrum_"/>
       !!]
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    call sedDustEmissionInitialize(self)
    return
  end function sedDustEmissionConstructorInternal

  subroutine sedDustEmissionInitialize(self)
    !!{RST
    Validate the configuration of a :galacticus-class:`nodePropertyExtractorSEDDustEmission` object, reporting any
    problem at construction rather than part way through a run, and build its quadrature rules.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nodePropertyExtractorSEDDustEmission), intent(inout) :: self
    type            (dustEmissionSpectrumList            ), pointer       :: dustEmissionSpectrum_
    type            (multiExtractorList                  ), pointer       :: extractor_                  , extractorSingle
    ! Tolerance on the mass fractions summing to one, loose enough to admit values written to a few significant figures.
    double precision                                      , parameter     :: toleranceSum         =1.0d-6
    integer                                                               :: countSpectra
    double precision                                                      :: wavelengthMinimum           , wavelengthMaximum

    if (self%resolution               <= 0.0d0                                                                    ) &
         & call Error_Report('`resolution` must be positive'                  //{introspection:location})
    if (self%wavelengthMinimum        <= 0.0d0 .or. self%wavelengthMaximum        <= self%wavelengthMinimum       ) &
         & call Error_Report('the wavelength range of the spectrum is invalid'//{introspection:location})
    if (self%wavelengthHeatingMinimum <= 0.0d0 .or. self%wavelengthHeatingMaximum <= self%wavelengthHeatingMinimum) &
         & call Error_Report('the wavelength range of the heating light is invalid'//{introspection:location})
    if (self%orderInclinationAverage  <  1                            ) call Error_Report('`orderInclinationAverage` must be positive'//{introspection:location})
    ! The factor by which the extreme wavelengths of a resolution element differ from its central wavelength, as for
    ! `nodePropertyExtractorSED`.
    self%factorWavelength=(1.0d0+sqrt(1.0d0+4.0d0*self%resolution**2))/2.0d0/self%resolution
    call gaussLegendreRule(self%orderInclinationAverage,self%cosineInclination,self%weight    )
    call gaussLegendreRule(sedDustEmissionOrderBin     ,self%abscissaeBin     ,self%weightsBin)
    ! One emission spectrum is needed for each phase of dust, and one mass fraction.
    countSpectra          =  0
    dustEmissionSpectrum_ => self%dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       countSpectra          =  countSpectra+1
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    if (countSpectra                 /= self%dustAttenuation_%countPhases())                                                                &
         & call Error_Report('one `dustEmissionSpectrum` must be given for each phase of dust of the attenuator'//{introspection:location})
    if (size(self%fractionsMassPhase) /= self%dustAttenuation_%countPhases())                                                               &
         & call Error_Report('one mass fraction must be given for each phase of dust of the attenuator'         //{introspection:location})
    if (any(self%fractionsMassPhase < 0.0d0) .or. abs(sum(self%fractionsMassPhase)-1.0d0) > toleranceSum)                                   &
         & call Error_Report('`fractionsMassPhase` must be non-negative and sum to one'                         //{introspection:location})
    ! Check the children: each must be able to decompose its luminosity, and a spectrum must be in the rest frame and span
    ! the heating light.
    call sedDustEmissionChildren(self,extractor_,extractorSingle)
    do while (associated(extractor_))
       if (.not.extractor_%extractor_%supportsAttenuation())                                                               &
            & call Error_Report(                                                                                           &
            &                   'property extractor "'//extractor_%extractor_%objectType()//'" does not support dust'   // &
            &                   ' attenuation, so can not supply light to heat dust'                                    // &
            &                   {introspection:location}                                                                   &
            &                  )
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorSED   )
          if (.not.child%isRestFrame())                                                                                    &
               & call Error_Report('spectra heating dust must be computed in the rest frame'//{introspection:location})
          call child%wavelengthRange(wavelengthMinimum,wavelengthMaximum)
          if (wavelengthMinimum > self%wavelengthHeatingMinimum .or. wavelengthMaximum < self%wavelengthHeatingMaximum)    &
               & call Error_Report(                                                                                        &
               &                   'spectra heating dust must span at least [wavelengthHeatingMinimum] to'              // &
               &                   ' [wavelengthHeatingMaximum], so that all of the light heating the dust is counted'  // &
               &                   {introspection:location}                                                                &
               &                  )
       class is (nodePropertyExtractorScalar)
          ! Supported.
       class is (nodePropertyExtractorTuple )
          ! Supported.
       class default
          call Error_Report(                                                                                               &
               &            'property extractor "'//extractor_%extractor_%objectType()//'" can not supply light to heat'// &
               &            ' dust - only spectra, and scalar or tuple luminosities, are supported'                     // &
               &            {introspection:location}                                                                       &
               &           )
       end select
       extractor_ => extractor_%next
    end do
    if (associated(extractorSingle)) deallocate(extractorSingle)
    return
  end subroutine sedDustEmissionInitialize

  subroutine sedDustEmissionChildren(self,extractors,extractorSingle)
    !!{RST
    Return a list of the child extractors supplying the light which heats the dust. If the heating-light extractor is a
    :galacticus-class:`nodePropertyExtractorMulti`, its own list is returned. Otherwise it is the only child, and a list
    of one element is allocated to hold it and returned as ``extractorSingle`` also, which the caller must deallocate
    when done---doing so does not destroy the extractor it points to.
    !!}
    implicit none
    class(nodePropertyExtractorSEDDustEmission), intent(inout), target  :: self
    type (multiExtractorList                  ), intent(  out), pointer :: extractors, extractorSingle

    extractorSingle => null()
    ! Match the multi extractor exactly: an extension of it, such as a dust attenuation extractor, is a single child.
    select type (multi_ => self%nodePropertyExtractor_)
    type is (nodePropertyExtractorMulti)
       extractors => multi_%extractors
    class default
       allocate(extractorSingle)
       extractorSingle%extractor_ => self%nodePropertyExtractor_
       extractorSingle%next       => null()
       extractors                 => extractorSingle
    end select
    return
  end subroutine sedDustEmissionChildren

  subroutine sedDustEmissionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorSEDDustEmission` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSEDDustEmission), intent(inout) :: self
    type(dustEmissionSpectrumList            ), pointer       :: dustEmissionSpectrum_, dustEmissionSpectrumNext

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%dustAttenuation_"      />
    <objectDestructor name="self%dustProperties_"       />
    <objectDestructor name="self%cosmologyFunctions_"   />
    !!]
    if (associated(self%dustEmissionSpectra)) then
       dustEmissionSpectrum_ => self%dustEmissionSpectra
       do while (associated(dustEmissionSpectrum_))
          dustEmissionSpectrumNext => dustEmissionSpectrum_%next
          !![
          <objectDestructor name="dustEmissionSpectrum_%dustEmissionSpectrum_"/>
          !!]
          deallocate(dustEmissionSpectrum_)
          dustEmissionSpectrum_ => dustEmissionSpectrumNext
       end do
    end if
    return
  end subroutine sedDustEmissionDestructor

  integer function sedDustEmissionCountBins(self) result(countBins)
    !!{RST
    Return the number of bins in the spectrum, computed as by :galacticus-class:`nodePropertyExtractorSED`.
    !!}
    implicit none
    class(nodePropertyExtractorSEDDustEmission), intent(inout) :: self

    countBins=int(log(self%wavelengthMaximum/self%wavelengthMinimum)/log(self%factorWavelength)/2.0d0)+1
    return
  end function sedDustEmissionCountBins

  integer function sedDustEmissionCountElements(self) result(countElements)
    !!{RST
    Return the number of properties emitted: the total, and optionally the emission of each phase.
    !!}
    implicit none
    class(nodePropertyExtractorSEDDustEmission), intent(inout) :: self

    countElements=1
    if (self%outputPhases) countElements=countElements+self%dustAttenuation_%countPhases()
    return
  end function sedDustEmissionCountElements

  function sedDustEmissionWavelengths(self) result(wavelengths)
    !!{RST
    Return the wavelengths, in Å, of the centers of the bins of the spectrum, in the frame in which it is computed. As
    for :galacticus-class:`nodePropertyExtractorSED`, adjacent bins are separated by a factor of the square of the factor
    by which the extremes of a resolution element differ from its center.
    !!}
    implicit none
    double precision                                      , allocatable  , dimension(:) :: wavelengths
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)               :: self
    integer                                                                             :: i

    allocate(wavelengths(self%countBins()))
    do i=1,size(wavelengths)
       wavelengths(i)=self%wavelengthMinimum*self%factorWavelength**(2*i-1)
    end do
    return
  end function sedDustEmissionWavelengths

  integer function sedDustEmissionElementCount(self,time) result(elementCount)
    !!{RST
    Return the number of properties emitted.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDDustEmission), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: time

    elementCount=self%countElements()
    return
  end function sedDustEmissionElementCount

  function sedDustEmissionSize(self,time) result(size)
    !!{RST
    Return the number of bins in the spectrum.
    !!}
    implicit none
    integer         (c_size_t                            )                :: size
    class           (nodePropertyExtractorSEDDustEmission), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: time

    size=self%countBins()
    return
  end function sedDustEmissionSize

  function sedDustEmissionExtract(self,node,time,instance) result(sed)
    !!{RST
    Return the spectrum of thermal emission from the dust of the galaxy.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galactic_Structure_Options      , only : componentTypeDisk, componentTypeNuclearStarCluster, componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    use :: Stellar_Luminosities_Structure  , only : frameObserved    , frameRest
    implicit none
    double precision                                      , dimension(:,:), allocatable :: sed
    class           (nodePropertyExtractorSEDDustEmission), intent(inout) , target      :: self
    type            (treeNode                            ), intent(inout) , target      :: node
    double precision                                      , intent(in   )               :: time
    type            (multiCounter                        ), intent(inout) , optional    :: instance
    type            (multiExtractorList                  ), pointer                     :: extractor_           , extractorSingle
    type            (dustEmissionSpectrumList            ), pointer                     :: dustEmissionSpectrum_
    double precision                                      , dimension(:,:), allocatable :: absorbed
    double precision                                      , dimension(:  ), allocatable :: luminosityAbsorbed   , wavelengthsChild , &
         &                                                                                 frequenciesChild     , unitsChild       , &
         &                                                                                 wavelengthsOutput    , wavelengthsSample, &
         &                                                                                 luminositySample     , luminosityBin
    integer                                                                             :: countPhases          , countBins        , &
         &                                                                                 countSamples         , i                , &
         &                                                                                 j                    , k
    double precision                                                                    :: massDust             , expansionFactor  , &
         &                                                                                 normalizationBin
    !$GLC attributes unused :: instance

    countPhases =self%dustAttenuation_%countPhases()
    countBins   =self%countBins                   ()
    countSamples=size(self%abscissaeBin)
    allocate(sed               (countBins,self%countElements()))
    allocate(luminosityAbsorbed(countPhases                   ))
    sed               =0.0d0
    luminosityAbsorbed=0.0d0
    ! Find the luminosity absorbed by each phase of dust, summed over all children, in L☉.
    call sedDustEmissionChildren(self,extractor_,extractorSingle)
       do while (associated(extractor_))
          call dustAbsorbedLuminosities(self%dustAttenuation_,extractor_%extractor_,node,time,absorbed=absorbed,cosineInclination=self%cosineInclination,weight=self%weight)
          select type (child => extractor_%extractor_)
          class is (nodePropertyExtractorSED   )
             ! Integrate the absorbed spectrum, L_ν in L☉ Hz⁻¹, over frequency by the trapezoidal rule.
             wavelengthsChild=child%wavelengths(time)
             frequenciesChild=speedLight*metersToAngstroms/wavelengthsChild
             do k=1,countPhases
                do j=2,size(frequenciesChild)
                   luminosityAbsorbed(k)=+luminosityAbsorbed(k)                          &
                        &                +0.5d0                                          &
                        &                *(absorbed(j,k)+absorbed(j-1,k))                &
                        &                *abs(frequenciesChild(j-1)-frequenciesChild(j))
                end do
             end do
          class is (nodePropertyExtractorScalar)
             ! Convert the absorbed luminosity to L☉.
             do k=1,countPhases
                luminosityAbsorbed(k)=+luminosityAbsorbed(  k) &
                     &                +absorbed          (1,k) &
                     &                *child%unitsInSI   (   ) &
                     &                /luminositySolar
             end do
          class is (nodePropertyExtractorTuple )
             ! Convert each absorbed luminosity to L☉.
             unitsChild=child%unitsInSI(time)
             do k=1,countPhases
                luminosityAbsorbed(k)=+luminosityAbsorbed(k)         &
                     &                +sum(absorbed(:,k)*unitsChild) &
                     &                /luminositySolar
             end do
          class default
             call Error_Report('unsupported child extractor'//{introspection:location})
          end select
          extractor_ => extractor_%next
       end do
    if (associated(extractorSingle)) deallocate(extractorSingle)
    ! Nothing more to do if no light is absorbed.
    if (all(luminosityAbsorbed <= 0.0d0)) return
    ! The dust mass of the galaxy, summed over its components.
    massDust=+self%dustProperties_%massDust(node,componentTypeDisk              ) &
         &   +self%dustProperties_%massDust(node,componentTypeSpheroid          ) &
         &   +self%dustProperties_%massDust(node,componentTypeNuclearStarCluster)
    ! Find the rest-frame wavelengths at which to evaluate the emission: the abscissae of the quadrature rule within each
    ! resolution element, uniformly spaced in the logarithm of wavelength between its extremes.
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=1.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select
    wavelengthsOutput=self%wavelengths()
    allocate(wavelengthsSample(countBins*countSamples))
    do i=1,countBins
       do j=1,countSamples
          wavelengthsSample((i-1)*countSamples+j)=+expansionFactor                                           &
               &                                  *wavelengthsOutput(i)                                      &
               &                                  *self%factorWavelength**(2.0d0*self%abscissaeBin(j)-1.0d0)
       end do
    end do
    ! The mean of L_ν over a resolution element, as defined by `nodePropertyExtractorSED`, is
    ! (f-1/f)⁻¹ ∫ L_ν dλ/λ, with the integral taken between the extremes of the element, a range of 2 ln f in ln λ.
    normalizationBin=+2.0d0                         &
         &           *log(self%factorWavelength)    &
         &           /(                             &
         &             +      self%factorWavelength &
         &             -1.0d0/self%factorWavelength &
         &            )
    allocate(luminosityBin(countBins))
    k                     =  0
    dustEmissionSpectrum_ => self%dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       k               =k+1
       luminositySample=dustEmissionSpectrum_%dustEmissionSpectrum_%luminosity(                                     &
            &                                                                  wavelengthsSample                  , &
            &                                                                  luminosityAbsorbed(k)              , &
            &                                                                  self%fractionsMassPhase(k)*massDust, &
            &                                                                  time                                 &
            &                                                                 )
       do i=1,countBins
          luminosityBin(i)=+normalizationBin                                                           &
               &           *sum(self%weightsBin*luminositySample((i-1)*countSamples+1:i*countSamples))
       end do
       sed(:,1)=sed(:,1)+luminosityBin
       if (self%outputPhases) sed(:,1+k)=luminosityBin
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    return
  end function sedDustEmissionExtract

  subroutine sedDustEmissionNames(self,names,time)
    !!{RST
    Return the names of the properties emitted.
    !!}
    use :: ISO_Varying_String, only : operator(//), operator(/=)
    implicit none
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)                             :: self
    double precision                                      , intent(in   ), optional                   :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: names
    type            (varying_string                      )               , dimension(:) , allocatable :: labels
    type            (varying_string                      )                                            :: name
    integer                                                                                           :: k
    !$GLC attributes unused :: time

    name="dustEmissionSED:"//self%dustAttenuation_%objectType(short=.true.)
    if (self%appendSuffix /= 'none') name=name//":"//self%appendSuffix
    allocate(names(self%countElements()))
    names(1)=name
    if (self%outputPhases) then
       ! Allocate explicitly: varying_string has a defined assignment, which does not allocate on assignment.
       allocate(labels(self%dustAttenuation_%countPhases()))
       labels=dustPhaseLabels(self%dustAttenuation_)
       do k=1,size(labels)
          names(1+k)=name//":"//labels(k)
       end do
    end if
    return
  end subroutine sedDustEmissionNames

  subroutine sedDustEmissionDescriptions(self,descriptions,time)
    !!{RST
    Return descriptions of the properties emitted.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)                             :: self
    double precision                                      , intent(in   ), optional                   :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: descriptions
    type            (varying_string                      )               , dimension(:) , allocatable :: labels
    integer                                                                                           :: k
    !$GLC attributes unused :: time

    allocate(descriptions(self%countElements()))
    descriptions(1)="Spectral energy density (SED), dL/dν, of thermal emission from dust, heated as absorbed by the '"//self%dustAttenuation_%objectType(short=.true.)//"' model [L☉ Hz⁻¹]."
    if (self%outputPhases) then
       allocate(labels(self%dustAttenuation_%countPhases()))
       labels=dustPhaseLabels(self%dustAttenuation_)
       do k=1,size(labels)
          descriptions(1+k)="Spectral energy density (SED), dL/dν, of thermal emission from the '"//labels(k)//"' phase of dust [L☉ Hz⁻¹]."
       end do
    end if
    return
  end subroutine sedDustEmissionDescriptions

  subroutine sedDustEmissionColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnits,time)
    !!{RST
    Return column descriptions of the spectrum: the wavelengths of its bins.
    !!}
    use            :: Numerical_Constants_Units, only : metersToAngstroms
    use            :: Units_MetaData           , only : unitType
    implicit none
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)                            :: self
    double precision                                      , intent(in   ), optional                  :: time
    type            (varying_string                      ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                      , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                      ), intent(  out)                            :: valuesDescription
    type            (unitType                            ), intent(  out)                            :: valuesUnits
    integer                                                                                          :: i
    character       (len=18                              )                                           :: label
    !$GLC attributes unused :: time

    values=self%wavelengths()
    allocate(descriptions(size(values)))
    do i=1,size(descriptions)
       write (label,'(a2,1x,e12.6,1x,a1)') "λ=",values(i),"Å"
       descriptions(i)=trim(label)
    end do
    valuesDescription=var_str('Wavelengths at which the SED is tabulated [in units of Å].')
    valuesUnits      =unitType(1.0d0/metersToAngstroms,"Angstroms","angstrom")
    return
  end subroutine sedDustEmissionColumnDescriptions

  function sedDustEmissionUnitsInSI(self,time) result(unitsInSI)
    !!{RST
    Return the units of the properties emitted in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    double precision                                      , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)               :: self
    double precision                                      , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(unitsInSI(self%countElements()))
    unitsInSI=luminositySolar
    return
  end function sedDustEmissionUnitsInSI

  function sedDustEmissionUnits(self,time) result(units)
    !!{RST
    Return the units of the properties emitted.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    use :: Units_MetaData                  , only : unitType
    implicit none
    type            (unitType                            ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorSEDDustEmission), intent(inout)             :: self
    double precision                                      , intent(in   ), optional   :: time
    !$GLC attributes unused :: time

    allocate(units(self%countElements()))
    units=unitType(luminositySolar,description='L☉',quantity='solLum')
    return
  end function sedDustEmissionUnits
