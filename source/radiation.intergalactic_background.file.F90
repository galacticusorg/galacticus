!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a class for intergalactic background light.
  
  use, intrinsic :: ISO_C_Binding
  use            :: Cosmology_Functions
  use            :: FGSL               , only : fgsl_interp_accel

  !# <radiationField name="radiationFieldIntergalacticBackgroundFile">
  !#  <description>A radiation field class for intergalactic background light with properties read from file.</description>
  !# </radiationField>
  type, extends(radiationFieldIntergalacticBackground) :: radiationFieldIntergalacticBackgroundFile
     !% A radiation field class for intergalactic background light with properties read from file.
     private
     class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_
     type            (varying_string         )                              :: fileName
     double precision                                                       :: time
     integer                                                                :: spectraTimesCount       , spectraWavelengthsCount
     double precision                         , allocatable, dimension(:  ) :: spectraTimes            , spectraWavelengths
     double precision                         , allocatable, dimension(:,:) :: spectra
     logical                                                                :: interpolationReset      , interpolationResetTimes
     type            (fgsl_interp_accel      )                              :: interpolationAccelerator, interpolationAcceleratorTimes
     integer         (c_size_t               )                              :: iTime
     double precision                         , dimension(0:1)              :: hTime
   contains
     !@ <objectMethods>
     !@  <object>radiationFieldIntergalacticBackground</object>
     !@  <objectMethod>
     !@   <method>timeSet</method>
     !@   <type>\void</type>
     !@   <arguments>\doublezero\ time\argin</arguments>
     !@   <description>Set the time for the radiation field.</description>
     !@  </objectMethod>
     !@ </objectMethods>
     final     ::            intergalacticBackgroundFileDestructor
     procedure :: flux    => intergalacticBackgroundFileFlux
     procedure :: timeSet => intergalacticBackgroundFileTimeSet
  end type radiationFieldIntergalacticBackgroundFile

  interface radiationFieldIntergalacticBackgroundFile
     !% Constructors for the {\normalfont \ttfamily intergalacticBackgroundFile} radiation field class.
     module procedure intergalacticBackgroundFileConstructorParameters
     module procedure intergalacticBackgroundFileConstructorInternal
  end interface radiationFieldIntergalacticBackgroundFile
  
  ! Current file format version for intergalactic background radiation files.
  integer, parameter :: intergalacticBackgroundFileFormatVersionCurrent=1

contains

  function intergalacticBackgroundFileConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily intergalacticBackgroundFile} radiation field class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type (radiationFieldIntergalacticBackgroundFile)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    type (varying_string                           )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file from which to read intergalactic background light properties.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=radiationFieldIntergalacticBackgroundFile(fileName,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function intergalacticBackgroundFileConstructorParameters

  function intergalacticBackgroundFileConstructorInternal(fileName,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily intergalacticBackgroundFile} radiation field class.
    use FoX_DOM
    use Galacticus_Error
    use Memory_Management
    use Array_Utilities
    use Galacticus_Paths
    use IO_XML
    implicit none
    type   (radiationFieldIntergalacticBackgroundFile)                        :: self
    type   (varying_string                           ), intent(in   )         :: fileName
    class  (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    type   (node                                     ), pointer               :: doc                , datum     , &
         &                                                                       spectrum           , wavelength
    type   (nodeList                                 ), pointer               :: spectraList
    integer                                                                   :: fileFormatVersion  , iSpectrum , &
         &                                                                       status             , jSpectrum
    logical                                                                   :: timesIncreasing
    !# <constructorAssign variables="fileName, *cosmologyFunctions_"/>

    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(fileName),iostat=status)
    if (status /= 0) call Galacticus_Error_Report('Unable to find or parse data file'//{introspection:location})
    ! Check the file format version of the file.
    datum => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(datum,fileFormatVersion)
    if (fileFormatVersion /= intergalacticBackgroundFileFormatVersionCurrent) call Galacticus_Error_Report('file format version is out of date'//{introspection:location})
    ! Get a list of all spectra.
    spectraList => getElementsByTagname(doc,"spectra")
    ! Get the wavelengths.
    wavelength => XML_Get_First_Element_By_Tag_Name(doc,"wavelengths")
    call XML_Array_Read(wavelength,"datum",self%spectraWavelengths)
    self%spectraWavelengthsCount=size(self%spectraWavelengths)
    ! Allocate array for spectra.
    self%spectraTimesCount=getLength(spectraList)
    call allocateArray(self%spectra     ,[self%spectraWavelengthsCount,self%spectraTimesCount])
    call allocateArray(self%spectraTimes,[                             self%spectraTimesCount])
    ! Read times.
    do iSpectrum=1,self%spectraTimesCount
       ! Get the data.
       spectrum => item(spectraList,iSpectrum-1)
       ! Extract the redshift.
       call extractDataContent(XML_Get_First_Element_By_Tag_Name(spectrum,"redshift"),self%spectraTimes(iSpectrum))
       ! Convert redshift to a time.
       self%spectraTimes(iSpectrum)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%spectraTimes(iSpectrum)))
    end do
    ! Check if the times are monotonically ordered.
    if (.not.Array_Is_Monotonic(self%spectraTimes)) call Galacticus_Error_Report('spectra must be monotonically ordered in time'//{introspection:location})
    timesIncreasing=Array_Is_Monotonic(self%spectraTimes,direction=directionIncreasing)
    ! Reverse times if necessary.
    if (.not.timesIncreasing) self%spectraTimes=Array_Reverse(self%spectraTimes)
    ! Read spectra into arrays.
    do iSpectrum=1,self%spectraTimesCount
       ! Determine where to store this spectrum, depending on whether the times were stored in increasing or decreasing order.
       if (timesIncreasing) then
          jSpectrum=                        +iSpectrum
       else
          jSpectrum=self%spectraTimesCount+1-iSpectrum
       end if
       ! Get the data.
       spectrum => item(spectraList,iSpectrum-1)
       ! Check that we have the correct number of data.
       if (XML_Array_Length(spectrum,"datum") /= self%spectraWavelengthsCount) call Galacticus_Error_Report('all spectra must contain the same number of wavelengths'//{introspection:location})
       ! Extract the data.
       call XML_Array_Read_Static(spectrum,"datum",self%spectra(:,jSpectrum))
    end do
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    self%interpolationReset     =.true.
    self%interpolationResetTimes=.true.
    return
  end function intergalacticBackgroundFileConstructorInternal

  subroutine intergalacticBackgroundFileDestructor(self)
    !% Destructor for the {\normalfont \ttfamily intergalacticBackgroundFile} radiation field class.
    implicit none
    type(radiationFieldIntergalacticBackgroundFile), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine intergalacticBackgroundFileDestructor
  
  double precision function intergalacticBackgroundFileFlux(self,wavelength,node)
    !% Return the flux in the intergalactic background radiation field.
    use Numerical_Interpolation
    implicit none
    class           (radiationFieldIntergalacticBackgroundFile), intent(inout)  :: self
    double precision                                           , intent(in   )  :: wavelength
    type            (treeNode                                 ), intent(inout)  :: node
    double precision                                           , dimension(0:1) :: hWavelength
    integer         (c_size_t                                 )                 :: iWavelength, jTime, &
         &                                                                         jWavelength
    !GCC$ attributes unused :: node

    intergalacticBackgroundFileFlux=0.0d0
    ! Return if out of range.
    if     (                                                                    &
         &   self%time  < self%spectraTimes      (1                           ) &
         &  .or.                                                                &
         &   self%time  > self%spectraTimes      (self%spectraTimesCount      ) &
         &  .or.                                                                &
         &   wavelength < self%spectraWavelengths(1                           ) &
         &  .or.                                                                &
         &   wavelength > self%spectraWavelengths(self%spectraWavelengthsCount) &
         & ) return
    ! Find interpolation in the array of wavelengths.
    iWavelength=Interpolate_Locate                 (self%spectraWavelengths,self%interpolationAccelerator,wavelength,reset=self%interpolationReset)
    hWavelength=Interpolate_Linear_Generate_Factors(self%spectraWavelengths,iWavelength                  ,wavelength                              )
    do jTime=0,1
       do jWavelength=0,1
          if (self%iTime+jTime > 0)                                                                          &
               & intergalacticBackgroundFileFlux=+intergalacticBackgroundFileFlux                            &
               &                                 +self%hTime      (                        jTime           ) &
               &                                 *     hWavelength(jWavelength                             ) &
               &                                 *self%spectra    (jWavelength+iWavelength,jTime+self%iTime)
       end do
    end do
    return
  end function intergalacticBackgroundFileFlux

  subroutine intergalacticBackgroundFileTimeSet(self,time)
    !% Set the time of the intergalactic backgroun radiation field.
    use Numerical_Interpolation
    implicit none
    class           (radiationFieldIntergalacticBackgroundFile), intent(inout) :: self
    double precision                                           , intent(in   ) :: time

    self%time=time
    ! Find interpolating factors in the array of times.
    self%iTime=Interpolate_Locate                 (self%spectraTimes,self%interpolationAcceleratorTimes,time,reset=self%interpolationResetTimes)
    self%hTime=Interpolate_Linear_Generate_Factors(self%spectraTimes,self%iTime                        ,time                                   )
    return
  end subroutine intergalacticBackgroundFileTimeSet
