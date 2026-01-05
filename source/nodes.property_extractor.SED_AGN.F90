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
  Implements a property extractor class for the SED of the AGN.
  !!}
  use :: Cosmology_Functions   , only : cosmologyFunctionsClass
  use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSEDAGN">
    <description>
      A property extractor class for the SED of the AGN. The spectrum is computed using the provided
      \refClass{accretionDiskSpectraClass} object, and will be output between wavelengths {\normalfont \ttfamily
      [wavelengthMinimum]} and {\normalfont \ttfamily [wavelengthMaximum]}. If {\normalfont \ttfamily [resolution]} is set to a
      positive value then this specifies the resolution, $\lambda/\Delta\lambda$, at which to compute the SED. If {\normalfont
      \ttfamily [resolution]} is non-positive then the SED will be output at the full native resolution provided by the
      \refClass{accretionDiskSpectraClass} object. The frame for the SED, {\normalfont \ttfamily rest} or {\normalfont \ttfamily
      observed}, is specified by {\normalfont \ttfamily [frame]}. Note that using {\normalfont \ttfamily observed} merely means
      that the extracted spectrum is evaluated at wavelength $\lambda_\mathrm{r} = \lambda_\mathrm{o} / (1+z)$ where
      $\lambda_\mathrm{o}$ is the observed wavelength (and is the wavelength returned by the `columnDescriptions` method),
      $\lambda_\mathrm{r}$ is the rest-frame wavelength, and $z$ is redshift---\emph{no adjustment} is made for the boost in observed
      flux due to the $\mathrm{d}\lambda_\mathrm{o}/\mathrm{d}\lambda_\mathrm{r}$ term which appears when computing observed-frame
      fluxes.
    </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSEDAGN
     !!{
     A property extractor class for the SED of the AGN.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer                   :: cosmologyFunctions_   => null()
     class           (accretionDiskSpectraClass), pointer                   :: accretionDiskSpectra_ => null()
     integer                                                                :: countWavelengths
     double precision                           , allocatable, dimension(:) :: wavelengths_
     double precision                                                       :: wavelengthMinimum              , wavelengthMaximum, &
          &                                                                    resolution                     , factorWavelength
     type            (enumerationFrameType     )                            :: frame
   contains
     !![
     <methods>
       <method description="Return an array of the wavelengths at which the SED is computed." method="wavelengths"/>
     </methods>
     !!]
     final     ::                            sedAGNDestructor
     procedure :: wavelengths             => sedAGNWavelengths
     procedure :: columnDescriptions      => sedAGNColumnDescriptions
     procedure :: size                    => sedAGNSize
     procedure :: elementCount            => sedAGNElementCount
     procedure :: extract                 => sedAGNExtract
     procedure :: names                   => sedAGNNames
     procedure :: descriptions            => sedAGNDescriptions
     procedure :: unitsInSI               => sedAGNUnitsInSI
  end type nodePropertyExtractorSEDAGN
  
  interface nodePropertyExtractorSEDAGN
     !!{
     Constructors for the \refClass{nodePropertyExtractorSEDAGN} output analysis class.
     !!}
     module procedure sedAGNConstructorParameters
     module procedure sedAGNConstructorInternal
  end interface nodePropertyExtractorSEDAGN
      
contains

  function sedAGNConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSEDAGN} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters              , only : inputParameters
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode
    implicit none
    type            (nodePropertyExtractorSEDAGN)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    class           (accretionDiskSpectraClass  ), pointer       :: accretionDiskSpectra_
    type            (varying_string             )                :: frame
    double precision                                             :: wavelengthMinimum    , wavelengthMaximum, &
         &                                                          resolution
    
    !![
    <inputParameter>
      <name>frame</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>The frame ({\normalfont \ttfamily rest} or {\normalfont \ttfamily observed}) for which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>resolution</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The resolution, $\lambda/\Delta\lambda$, at which to compute the SED. If a negative value is given the SED will be computed at the full resolution provided by the \refClass{accretionDiskSpectraClass} object.</description>
    </inputParameter>
    <objectBuilder class="accretionDiskSpectra" name="accretionDiskSpectra_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !!]
    self=nodePropertyExtractorSEDAGN(enumerationFrameEncode(char(frame),includesPrefix=.false.),wavelengthMinimum,wavelengthMaximum,resolution,cosmologyFunctions_,accretionDiskSpectra_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="accretionDiskSpectra_"/>
    !!]
    return
  end function sedAGNConstructorParameters

  function sedAGNConstructorInternal(frame,wavelengthMinimum,wavelengthMaximum,resolution,cosmologyFunctions_,accretionDiskSpectra_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSEDAGN} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nodePropertyExtractorSEDAGN)                        :: self
    type            (enumerationFrameType       ), intent(in   )         :: frame
    class           (accretionDiskSpectraClass  ), intent(in   ), target :: accretionDiskSpectra_
    class           (cosmologyFunctionsClass    ), intent(in   ), target :: cosmologyFunctions_
    double precision                             , intent(in   )         :: wavelengthMinimum    , wavelengthMaximum, &
         &                                                                  resolution
    !![
    <constructorAssign variables="frame, wavelengthMinimum, wavelengthMaximum, resolution, *cosmologyFunctions_, *accretionDiskSpectra_"/>
    !!]
    
    if (resolution > 0.0d0) then
       ! Compute the factor by which the minimum/maximum wavelength in a resolution element differ from the central wavelength.
       self%factorWavelength=(1.0d0+sqrt(1.0d0+4.0d0*resolution**2))/2.0d0/resolution
    else
       ! Get lists of all wavelengths available for the spectrum.
       call self%accretionDiskSpectra_%wavelengths(self%countWavelengths,self%wavelengths_)
    end if
    return
  end function sedAGNConstructorInternal

  subroutine sedAGNDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorSEDAGN} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSEDAGN), intent(inout) :: self
    
    !![
    <objectDestructor name="self%accretionDiskSpectra_"/>
    <objectDestructor name="self%cosmologyFunctions_"  />
    !!]
    return
  end subroutine sedAGNDestructor

  integer function sedAGNElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily sedAGN} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDAGN), intent(inout) :: self
    double precision                             , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    sedAGNElementCount=1
    return
  end function sedAGNElementCount

  function sedAGNSize(self,time) result(size)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily sedAGN} property extractors.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved 
    implicit none
    integer         (c_size_t                   )                              :: size
    class           (nodePropertyExtractorSEDAGN), intent(inout)               :: self
    double precision                             , intent(in   )               :: time
    logical                                      , allocatable  , dimension(:) :: selection
    double precision                                                           :: expansionFactor

    if (self%resolution < 0.0d0) then
       ! Full resolution, so the number of wavelengths is simply the total number available within the wavelength range.
       allocate(selection(self%countWavelengths))
       select case (self%frame%ID)
       case (frameRest    %ID)
          expansionFactor=1.0d0
          selection      =                                                             &
               &           self%wavelengths_                 >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_                 <= self%wavelengthMaximum
       case (frameObserved%ID)
          expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
          selection      =                                                             &
               &           self%wavelengths_/expansionFactor >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_/expansionFactor <= self%wavelengthMaximum
       case default
          expansionFactor=1.0d0
          size           =0_c_size_t
          call Error_Report('unknown frame'//{introspection:location})
       end select
       size=count(selection)
    else
       ! The number of wavelengths must be computed.
       size=int(log(self%wavelengthMaximum/self%wavelengthMinimum)/log(self%factorWavelength)/2.0d0)+1
    end if
    return
  end function sedAGNSize

  function sedAGNExtract(self,node,time,instance) result(sed)
    !!{
    Implement a {\normalfont \ttfamily sedAGN} property extractor.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved
    implicit none
    double precision                             , dimension(:,:), allocatable :: sed
    class           (nodePropertyExtractorSEDAGN), intent(inout) , target      :: self
    type            (treeNode                   ), intent(inout) , target      :: node
    double precision                             , intent(in   )               :: time
    type            (multiCounter               ), intent(inout) , optional    :: instance
    double precision                             , dimension(:  ), allocatable :: wavelengths
    integer         (c_size_t                   )                              :: countWavelengths, i
    double precision                                                           :: expansionFactor
    !$GLC attributes unused :: instance

    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select
    countWavelengths=self%size(time)
    allocate(sed        (countWavelengths,1))
    allocate(wavelengths(countWavelengths  ))
    ! Get wavelengths and shift them to the rest frame of the AGN if necessary.
    wavelengths=+self%wavelengths    (time) &
         &      *     expansionFactor
    ! Evaluate the SEDAGN at each wavelength.
    do i=1,countWavelengths
       sed(i,1)=self%accretionDiskSpectra_%spectrum(node,wavelengths(i))
    end do
    return
  end function sedAGNExtract

  subroutine sedAGNNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily sedAGN} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDAGN), intent(inout)                             :: self
    double precision                             , intent(in   ), optional                   :: time
    type            (varying_string             ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(1))
    names(1)="agnSED"
    return
  end subroutine sedAGNNames

  subroutine sedAGNDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily agnSED} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDAGN), intent(inout)                             :: self
    double precision                             , intent(in   ), optional                   :: time
    type            (varying_string             ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(1))
    descriptions(1)="Spectral energy density (SED), dL/dν for the AGN [L☉ Hz⁻¹]."
    return
  end subroutine sedAGNDescriptions

  function sedAGNWavelengths(self,time) result(wavelengths)
    !!{
    Return wavelengths at which the SED is tabulated.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved
    implicit none
    double precision                             , dimension(:) , allocatable :: wavelengths
    class           (nodePropertyExtractorSEDAGN), intent(inout)              :: self
    double precision                             , intent(in   )              :: time
    integer         (c_size_t                   )                             :: i          , j
    double precision                                                          :: wavelength , expansionFactor

    allocate(wavelengths(self%size(time)))
    j=0
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select
    do i=1,size(wavelengths)
       if (self%resolution < 0.0d0) then
          ! Full resolution SED.
          j=j+1
          do while (self%wavelengths_(j)/expansionFactor < self%wavelengthMinimum)
             j=j+1
          end do
          wavelength=self%wavelengths_(j)/expansionFactor
       else
          ! Finite resolution SED.
          if (i == 1) then
             wavelength=+self%wavelengthMinimum    &
                  &     *self%factorWavelength
          else
             wavelength=+     wavelength           &
                  &     *self%factorWavelength **2
          end if
       end if
       wavelengths(i)=wavelength
    end do
    return
  end function sedAGNWavelengths

  subroutine sedAGNColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily sedAGN} property.
    !!}
    use :: Numerical_Constants_Units, only : metersToAngstroms
    implicit none
    class           (nodePropertyExtractorSEDAGN), intent(inout)                            :: self
    double precision                             , intent(in   ), optional                  :: time
    type            (varying_string             ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                             , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string             ), intent(  out)                            :: valuesDescription
    double precision                             , intent(  out)                            :: valuesUnitsInSI
    integer         (c_size_t                   )                                           :: i
    character       (len=18                     )                                           :: label
    !$GLC attributes unused :: self, time

    allocate(descriptions(self%size(time)))
    allocate(values      (self%size(time)))
    values=self%wavelengths(time)
    do i=1,size(descriptions)      
       write (label,'(a2,1x,e12.6,1x,a1)') "λ=",values(i),"Å"
       descriptions(i)=trim(label)
    end do
    valuesDescription=var_str('Wavelengths at which the SED is tabulated [in units of Å].')
    valuesUnitsInSI  =1.0d0/metersToAngstroms
    return
  end subroutine sedAGNColumnDescriptions

  function sedAGNUnitsInSI(self,time) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily sedAGN} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    double precision                             , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorSEDAGN), intent(inout)               :: self
    double precision                             , intent(in   ), optional     :: time
    !$GLC attributes unused :: self, time

    allocate(unitsInSI(1))
    unitsInSI(1)=luminositySolar
    return
  end function sedAGNUnitsInSI
