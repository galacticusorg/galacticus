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
Implements a property extractor class for broadband stellar luminosities computed from an SED.
!!}

  use :: Cosmology_Functions           , only : cosmologyFunctionsClass
  use :: ISO_Varying_String            , only : varying_string
  use :: Output_Times                  , only : outputTimesClass
  use :: Numerical_Interpolation       , only : interpolator
  use :: Stellar_Luminosities_Structure, only : enumerationFrameType
  
  type :: filterResponse
     !!{RST
     Type used to hold pointers to filter response functions.
     !!}
     type(interpolator), pointer :: interpolator_ => null()
  end type filterResponse
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityStellarFromSED" docformat="rst">
   <description>
   A property extractor that computes broadband stellar luminosities in AB units by integrating an SED (provided by a :galacticus-class:`nodePropertyExtractorSED` property extractor object) through a set of broadband filter response functions specified by ``filterNames``. For each filter, the filter transmission curve is convolved with the galaxy SED to compute :math:`L_\nu` in the AB zero-point system. This approach allows luminosities to be derived from an already-computed SED without re-running the full stellar population synthesis, making it efficient when many filters are required from a single SED.

   Luminosities are computed in the frame given by ``frame`` (``rest`` or ``observed``). If the SED is in the same frame the filter response is integrated over it directly. Observed-frame luminosities can also be computed from a rest-frame SED, by blueshifting the filter response into the rest frame of the galaxy:

   .. math::

      L = \frac{\int R(\lambda/a) M_z(\lambda,z) S(\lambda) \mathrm{d}\lambda/\lambda}{\int R(\lambda) L_\mathrm{AB} \mathrm{d}\lambda/\lambda},

   where :math:`R` is the filter response, :math:`S` the rest-frame SED, :math:`a` the expansion factor, and :math:`L_\mathrm{AB}` the AB zero point. The measure :math:`\mathrm{d}\lambda/\lambda` is unchanged by the blueshift, so this is identical to integrating the filter over the observed-frame SED.

   Here :math:`M_z(\lambda,z)` is the redshift-dependent factor of the postprocessing of the SED (e.g. attenuation by the :term:`IGM`). Where that postprocessing is separable (see the ``isSeparable`` method of :galacticus-class:`stellarPopulationSpectraPostprocessorClass`), and the star formation history has fixed ages (see :galacticus-class:`starFormationHistoryFixedAges`), the SED is computed with only the redshift-independent factor of the postprocessing applied, and :math:`M_z(\lambda,z)` is applied to the filter instead. The SED is then independent of redshift, so a single template serves galaxies at all redshifts---which makes this the efficient route to observed-frame luminosities for lightcone output. Otherwise, the SED is computed with its full postprocessing at the redshift of the galaxy, and :math:`M_z(\lambda,z)=1`.

   The SED must span the range of wavelengths covered by each filter (blueshifted into the rest frame if appropriate), as otherwise the luminosity would be silently truncated. This is reported as an error.
   </description>
   <stateStorable>
    <functionClass variables="nodePropertyExtractor_"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLuminosityStellarFromSED
     !!{RST
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     type            (varying_string          ), allocatable, dimension(:  ) :: filterNames
     integer                                   , allocatable, dimension(:  ) :: luminosityIndex
     type            (filterResponse          ), allocatable, dimension(:  ) :: filterResponses
     double precision                          , allocatable, dimension(:,:) :: wavelengthRanges
     class           (nodePropertyExtractorSED), pointer                     :: nodePropertyExtractor_ => null()
     class           (cosmologyFunctionsClass ), pointer                     :: cosmologyFunctions_    => null()
     type            (enumerationFrameType    )                              :: frame
   contains
     final     ::                luminosityStellarFromSEDDestructor
     procedure :: elementCount => luminosityStellarFromSEDElementCount
     procedure :: extract      => luminosityStellarFromSEDExtract
     procedure :: quantity     => luminosityStellarFromSEDQuantity
     procedure :: names        => luminosityStellarFromSEDNames
     procedure :: descriptions => luminosityStellarFromSEDDescriptions
     procedure :: unitsInSI    => luminosityStellarFromSEDUnitsInSI
     procedure :: units        => luminosityStellarFromSEDUnits
  end type nodePropertyExtractorLuminosityStellarFromSED

  interface nodePropertyExtractorLuminosityStellarFromSED
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED` property extractor class.
     !!}
     module procedure luminosityStellarFromSEDConstructorParameters
     module procedure luminosityStellarFromSEDConstructorInternal
  end interface nodePropertyExtractorLuminosityStellarFromSED

contains

  function luminosityStellarFromSEDConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED` property extractor class which takes a parameter set as input.
    !!}
    use :: Error                         , only : Error_Report
    use :: Input_Parameters              , only : inputParameter        , inputParameters
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode
    implicit none
    type (nodePropertyExtractorLuminosityStellarFromSED)                             :: self
    type (inputParameters                              ), intent(inout)              :: parameters
    class(nodePropertyExtractorClass                   ), pointer                    :: nodePropertyExtractor_
    class(cosmologyFunctionsClass                      ), pointer                    :: cosmologyFunctions_
    type (varying_string                               ), dimension(:) , allocatable :: filterNames
    type (varying_string                               )                             :: frame

    allocate(filterNames(parameters%count('filterNames')))
    !![
    <inputParameter docformat="rst">
      <name>filterNames</name>
      <source>parameters</source>
      <description>
      The filters to select.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>frame</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>
      The frame (``rest`` or ``observed``) for which to compute luminosities.
      </description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorSED)
       self=nodePropertyExtractorLuminosityStellarFromSED(filterNames,enumerationFrameEncode(char(frame),includesPrefix=.false.),nodePropertyExtractor_,cosmologyFunctions_)
    class default
       call Error_Report('an "SED" nodePropertyExtractor is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    <objectDestructor name="cosmologyFunctions_"   />
    !!]
    return
  end function luminosityStellarFromSEDConstructorParameters

  function luminosityStellarFromSEDConstructorInternal(filterNames,frame,nodePropertyExtractor_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED` property extractor class.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Error                         , only : Error_Report
    use            :: Instruments_Filters           , only : Filter_Get_Index, Filter_Response_Function, Filter_Extent
    use            :: Stellar_Luminosities_Structure, only : frameRest
    implicit none
    type   (nodePropertyExtractorLuminosityStellarFromSED)                              :: self
    type   (varying_string                               ), intent(in   ), dimension(:) :: filterNames
    type   (enumerationFrameType                         ), intent(in   )               :: frame
    class  (nodePropertyExtractorSED                     ), intent(in   ), target       :: nodePropertyExtractor_
    class  (cosmologyFunctionsClass                      ), intent(in   ), target       :: cosmologyFunctions_
    integer                                                                             :: i                     , indexFilter
    !![
    <constructorAssign variables="filterNames, frame, *nodePropertyExtractor_, *cosmologyFunctions_"/>
    !!]

    ! A rest-frame luminosity can not be computed from an observed-frame SED.
    if (self%frame == frameRest .and. .not.self%nodePropertyExtractor_%isRestFrame()) &
         & call Error_Report('rest-frame luminosities can not be computed from an observed-frame SED'//{introspection:location})

    ! Get filter functions for all requested filters.
    allocate(self%filterResponses (size(filterNames)  ))
    allocate(self%wavelengthRanges(size(filterNames),2))
    do i=1,size(filterNames)
       indexFilter                              =  Filter_Get_Index        (filterNames(i          ))
       self%filterResponses (i  )%interpolator_ => Filter_Response_Function(            indexFilter )
       self%wavelengthRanges(i,:)               =  Filter_Extent           (            indexFilter )
    end do
    return
  end function luminosityStellarFromSEDConstructorInternal
  
  subroutine luminosityStellarFromSEDDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityStellarFromSED), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%cosmologyFunctions_"   />
    !!]
    return
  end subroutine luminosityStellarFromSEDDestructor

  integer function luminosityStellarFromSEDElementCount(self,time)
    !!{RST
    Return the number of elements in the ``luminosityStellarFromSED`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    !$GLC attributes unused :: time

    luminosityStellarFromSEDElementCount=size(self%filterNames)
    return
  end function luminosityStellarFromSEDElementCount

  function luminosityStellarFromSEDExtract(self,node,time,instance) result(luminosities)
    !!{RST
    Implement a stellar luminosity output analysis property extractor.
    !!}
    use :: Stellar_Luminosities_Structure, only : frameObserved
    use :: Table_Labels                  , only : extrapolationTypeZero
    use :: Numerical_Integration         , only : GSL_Integ_Gauss15    , integrator
    implicit none
    double precision                                               , allocatable   , dimension(:) :: luminosities
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout) , target       :: self
    type            (treeNode                                     ), intent(inout) , target       :: node
    double precision                                               , intent(in   )                :: time
    type            (multiCounter                                 ), intent(inout) , optional     :: instance
    double precision                                               , dimension(:  ), allocatable  :: wavelengths
    double precision                                               , dimension(:,:), allocatable  :: sed
    type            (integrator                                   ), allocatable                  :: integrator_      , integratorAB_
    type            (interpolator                                 ), allocatable                  :: interpolatorSED   
    integer                                                                                       :: i                , errorStatus
    double precision                                                                              :: expansionFactor  , redshift
    logical                                                                                       :: separateRedshift

    allocate(luminosities(size(self%filterNames)))
    ! Observed-frame luminosities from a rest-frame SED are computed by blueshifting the filter response into the rest frame of
    ! the galaxy. Otherwise the SED is in the frame of the luminosity, and the filter is integrated over it directly.
    if (self%frame == frameObserved .and. self%nodePropertyExtractor_%isRestFrame()) then
       expansionFactor =self%cosmologyFunctions_%expansionFactor            (time           )
       redshift        =self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
       ! Where possible, compute the SED with only the redshift-independent factor of its postprocessing applied, and apply the
       ! redshift-dependent factor to the filter instead. The SED can then be computed from a single template.
       separateRedshift=self%nodePropertyExtractor_%isRedshiftSeparable()
    else
       expansionFactor =1.0d0
       redshift        =0.0d0
       separateRedshift=.false.
    end if
    wavelengths=self%nodePropertyExtractor_%wavelengths(time)
    ! Check that the SED spans the range of each filter. Beyond its tabulated range the SED is extrapolated as zero, so a
    ! filter which extends beyond it would give a silently truncated luminosity.
    do i=1,size(self%filterNames)
       if     (                                                                              &
            &   size(wavelengths)              == 0                                          &
            &  .or.                                                                          &
            &   wavelengths(1                ) >  self%wavelengthRanges(i,1)*expansionFactor &
            &  .or.                                                                          &
            &   wavelengths(size(wavelengths)) <  self%wavelengthRanges(i,2)*expansionFactor &
            & ) call reportCoverageError(i)
    end do
    if (separateRedshift) then
       sed=self%nodePropertyExtractor_%extractRedshiftIndependent(node,time         )
    else
       sed=self%nodePropertyExtractor_%extract                   (node,time,instance)
    end if
    ! A galaxy with no SED has zero luminosity.
    if (all(sed(:,1) == 0.0d0)) then
       luminosities=0.0d0
       return
    end if
    allocate(interpolatorSED)
    allocate(integrator_    )
    allocate(integratorAB_  )
    interpolatorSED=interpolator(wavelengths                  ,                  sed(:,1),extrapolationType=extrapolationTypeZero                                )
    integrator_    =integrator  (integrandFilteredLuminosity  ,toleranceRelative=4.0d-3  ,integrationRule  =GSL_Integ_Gauss15    ,intervalsMaximum=10000_c_size_t)
    integratorAB_  =integrator  (integrandFilteredLuminosityAB,toleranceRelative=4.0d-3                                                                          )
    do i=1,size(self%filterNames)
       luminosities(i)=+integrator_  %integrate(self%wavelengthRanges(i,1)*expansionFactor,self%wavelengthRanges(i,2)*expansionFactor,status=errorStatus) &
            &          /integratorAB_%integrate(self%wavelengthRanges(i,1)                ,self%wavelengthRanges(i,2)                                   )
    end do
    deallocate(interpolatorSED)
    deallocate(integrator_    )
    deallocate(integratorAB_  )
    return

  contains

    double precision function integrandFilteredLuminosity(wavelength)
      !!{RST
      Integrand for the luminosity through a given filter.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      ! The factor of 1/λ appears since we want to integrate F_ν (dν / ν) and dν = -c/λ² dλ. Note that we follow the convention of
      ! Hogg et al. (2002) and assume that the filter response gives the fraction of incident photons received by the detector at
      ! a given wavelength, multiplied by the relative photon response (which will be 1 for a photon-counting detector such as a
      ! CCD, or proportional to the photon energy for a bolometer/calorimeter type detector).
      !
      ! When computing an observed-frame luminosity from a rest-frame SED, the wavelength here is in the rest frame, and the
      ! filter response is evaluated at the corresponding observed wavelength, λ/a. The measure dλ/λ is unchanged by this
      ! transformation.
      integrandFilteredLuminosity=+self%filterResponses(i)%interpolator_  %interpolate(wavelength/expansionFactor) &
           &                      *                        interpolatorSED%interpolate(wavelength                ) &
           &                      /                                                    wavelength
      if (separateRedshift)                                                                                            &
           & integrandFilteredLuminosity=+                            integrandFilteredLuminosity                      &
           &                             *self%nodePropertyExtractor_%multiplierRedshift         (wavelength,redshift)
      return
    end function integrandFilteredLuminosity

    subroutine reportCoverageError(i)
      !!{RST
      Report that the SED does not span the range of a filter.
      !!}
      use :: Display           , only : displayGreen, displayReset
      use :: Error             , only : Error_Report
      use :: ISO_Varying_String, only : operator(//), var_str
      implicit none
      integer                  , intent(in   ) :: i
      character(len=12        )                :: label
      type     (varying_string)                :: message

      message=var_str('the SED does not cover the wavelength range of filter "')//self%filterNames(i)//'"'
      !$omp critical(gfortranInternalIO)
      write (label,'(f12.3)') self%wavelengthRanges(i,1)*expansionFactor
      message=message//char(10)//'  filter range: '//trim(adjustl(label))//'Å - '
      write (label,'(f12.3)') self%wavelengthRanges(i,2)*expansionFactor
      message=message//trim(adjustl(label))//'Å'
      if (expansionFactor < 1.0d0) then
         write (label,'(f12.6)') redshift
         message=message//' (blueshifted to the rest frame of a galaxy at z='//trim(adjustl(label))//')'
      end if
      if (size(wavelengths) > 0) then
         write (label,'(f12.3)') wavelengths(1)
         message=message//char(10)//'     SED range: '//trim(adjustl(label))//'Å - '
         write (label,'(f12.3)') wavelengths(size(wavelengths))
         message=message//trim(adjustl(label))//'Å'
      else
         message=message//char(10)//'     SED range: empty'
      end if
      !$omp end critical(gfortranInternalIO)
      message=message//char(10)//displayGreen()//'HELP:'//displayReset()//' widen the range of wavelengths (`wavelengthMinimum` and `wavelengthMaximum`) of the SED'
      if (self%frame == frameObserved .and. self%nodePropertyExtractor_%isRestFrame()) &
           & message=message//char(10)//'      to cover the filter blueshifted to the rest frame of galaxies at all redshifts of interest'
      call Error_Report(message//{introspection:location})
      return
    end subroutine reportCoverageError

    double precision function integrandFilteredLuminosityAB(wavelength)
      !!{RST
      Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
      !!}
      use :: Numerical_Constants_Astronomical, only : luminositySolar, luminosityZeroPointAB
      implicit none
      double precision, intent(in   ) :: wavelength
      ! Luminosity of a zeroth magnitude (AB) source in Solar luminosities per Hz.
      double precision, parameter     :: luminosityZeroPointABSolar=luminosityZeroPointAB/luminositySolar

      integrandFilteredLuminosityAB=+self%filterResponses(i)%interpolator_%interpolate(wavelength) &
           &                        *luminosityZeroPointABSolar                                    &
           &                        /                                                  wavelength
      return
    end function integrandFilteredLuminosityAB
    
  end function luminosityStellarFromSEDExtract

  function luminosityStellarFromSEDQuantity(self)
    !!{RST
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: luminosityStellarFromSEDQuantity
    class(nodePropertyExtractorLuminosityStellarFromSED), intent(inout) :: self
    !$GLC attributes unused :: self

    luminosityStellarFromSEDQuantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function luminosityStellarFromSEDQuantity

  subroutine luminosityStellarFromSEDNames(self,time,names)
    !!{RST
    Return the name of the luminosityStellarFromSED property.
    !!}
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)                            :: self
    double precision                                               , intent(in   )                            :: time
    type            (varying_string                               ), intent(inout), dimension(:), allocatable :: names
    type            (varying_string                               )               , dimension(:), allocatable :: name
    integer                                                                                                   :: i
    
    allocate(names(size(self%filterNames)))
    call self%nodePropertyExtractor_%names(name,time)
    do i=1,size(names)
       names(i)=name(1)//":"//self%filterNames(i)//":"//enumerationFrameDecode(self%frame,includePrefix=.false.)
    end do
    return
  end subroutine luminosityStellarFromSEDNames

  subroutine luminosityStellarFromSEDDescriptions(self,time,descriptions)
    !!{RST
    Return a description of the luminosityStellarFromSED property.
    !!}
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)                            :: self
    double precision                                               , intent(in   )                            :: time
    type            (varying_string                               ), intent(inout), dimension(:), allocatable :: descriptions
    type            (varying_string                               )               , dimension(:), allocatable :: description
    integer                                                                                                   :: i

    allocate(descriptions(size(self%filterNames)))
    call self%nodePropertyExtractor_%descriptions(description,time)
    do i=1,size(descriptions)
       descriptions(i)="Luminosity in "//self%filterNames(i)//" filter in the "//enumerationFrameDecode(self%frame,includePrefix=.false.)//" frame in units of the AB-magnitude system zero-point; derived from: "//description(1)
    end do
    return
  end subroutine luminosityStellarFromSEDDescriptions

  function luminosityStellarFromSEDUnitsInSI(self,time) result(unitsInSI)
    !!{RST
    Return the units of the luminosityStellarFromSED property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    implicit none
    double precision                                               , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)              :: self
    double precision                                               , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(unitsInSI(size(self%filterNames)))
    unitsInSI=luminosityZeroPointAB
    return
  end function luminosityStellarFromSEDUnitsInSI

  function luminosityStellarFromSEDUnits(self,time) result(units)
    !!{RST
    Return the units of the luminosityStellarFromSED properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                                     ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)             :: self
    double precision                                               , intent(in   )             :: time
    double precision                                               , dimension(:), allocatable :: siValues
    integer                                                                                    :: i

    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i),description='AB-magnitude zero point',quantity='4.465920e17 W/Hz')
    end do
    return
  end function luminosityStellarFromSEDUnits
