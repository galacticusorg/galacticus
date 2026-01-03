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
Implements a stellar mass output analysis property extractor class.
!!}

  use :: ISO_Varying_String     , only : varying_string
  use :: Output_Times           , only : outputTimesClass
  use :: Numerical_Interpolation, only : interpolator
  
  type :: filterResponse
     !!{
     Type used to hold pointers to filter response functions.
     !!}
     type(interpolator), pointer :: interpolator_ => null()
  end type filterResponse
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityStellarFromSED">
   <description>A stellar luminosity output analysis property extractor class.</description>
   <deepCopy>
    <functionClass variables="nodePropertyExtractor_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="nodePropertyExtractor_"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLuminosityStellarFromSED
     !!{
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     type            (varying_string          ), allocatable, dimension(:  ) :: filterNames
     integer                                   , allocatable, dimension(:  ) :: luminosityIndex
     type            (filterResponse          ), allocatable, dimension(:  ) :: filterResponses
     double precision                          , allocatable, dimension(:,:) :: wavelengthRanges
     class           (nodePropertyExtractorSED), pointer                     :: nodePropertyExtractor_ => null()
   contains
     final     ::                luminosityStellarFromSEDDestructor
     procedure :: elementCount => luminosityStellarFromSEDElementCount
     procedure :: extract      => luminosityStellarFromSEDExtract
     procedure :: quantity     => luminosityStellarFromSEDQuantity
     procedure :: names        => luminosityStellarFromSEDNames
     procedure :: descriptions => luminosityStellarFromSEDDescriptions
     procedure :: unitsInSI    => luminosityStellarFromSEDUnitsInSI
  end type nodePropertyExtractorLuminosityStellarFromSED

  interface nodePropertyExtractorLuminosityStellarFromSED
     !!{
     Constructors for the \refClass{nodePropertyExtractorLuminosityStellarFromSED} output analysis class.
     !!}
     module procedure luminosityStellarFromSEDConstructorParameters
     module procedure luminosityStellarFromSEDConstructorInternal
  end interface nodePropertyExtractorLuminosityStellarFromSED

contains

  function luminosityStellarFromSEDConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorLuminosityStellarFromSED} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorLuminosityStellarFromSED)                             :: self
    type (inputParameters                              ), intent(inout)              :: parameters
    class(nodePropertyExtractorClass                   ), pointer                    :: nodePropertyExtractor_
    type (varying_string                               ), dimension(:) , allocatable :: filterNames

    allocate(filterNames(parameters%count('filterNames')))
    !![
    <inputParameter>
      <name>filterNames</name>
      <source>parameters</source>
      <description>The filters to select.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorSED)
       self=nodePropertyExtractorLuminosityStellarFromSED(filterNames,nodePropertyExtractor_)
    class default
       call Error_Report('an "SED" nodePropertyExtractor is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function luminosityStellarFromSEDConstructorParameters

  function luminosityStellarFromSEDConstructorInternal(filterNames,nodePropertyExtractor_) result(self)
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: Instruments_Filters, only : Filter_Get_Index, Filter_Response_Function, Filter_Extent
    implicit none
    type   (nodePropertyExtractorLuminosityStellarFromSED)                              :: self
    type   (varying_string                               ), intent(in   ), dimension(:) :: filterNames
    class  (nodePropertyExtractorSED                     ), intent(in   ), target       :: nodePropertyExtractor_
    integer                                                                             :: i                     , indexFilter
    !![
    <constructorAssign variables="filterNames, *nodePropertyExtractor_"/>
    !!]

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
    !!{
    Destructor for the \refClass{nodePropertyExtractorLuminosityStellarFromSED} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityStellarFromSED), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine luminosityStellarFromSEDDestructor

  integer function luminosityStellarFromSEDElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily luminosityStellarFromSED} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    !$GLC attributes unused :: time

    luminosityStellarFromSEDElementCount=size(self%filterNames)
    return
  end function luminosityStellarFromSEDElementCount

  function luminosityStellarFromSEDExtract(self,node,time,instance) result(luminosities)
    !!{
    Implement a stellar luminosity output analysis property extractor.
    !!}
    use :: Table_Labels         , only : extrapolationTypeZero
    use :: Numerical_Integration, only : GSL_Integ_Gauss15    , integrator
    implicit none
    double precision                                               , allocatable   , dimension(:) :: luminosities
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout) , target       :: self
    type            (treeNode                                     ), intent(inout) , target       :: node
    double precision                                               , intent(in   )                :: time
    type            (multiCounter                                 ), intent(inout) , optional     :: instance
    double precision                                               , dimension(:  ), allocatable  :: wavelengths
    double precision                                               , dimension(:,:), allocatable  :: sed
    type            (integrator                                   ), allocatable                  :: integrator_    , integratorAB_
    type            (interpolator                                 ), allocatable                  :: interpolatorSED   
    integer                                                                                       :: i              , errorStatus

    allocate(luminosities   (size(self%filterNames)))
    allocate(interpolatorSED                        )
    allocate(integrator_                            )
    allocate(integratorAB_                          )
    wavelengths     =  self        %nodePropertyExtractor_%wavelengths(                              time                                                                                              )
    sed             =  self        %nodePropertyExtractor_%extract    (node                         ,time                      ,                  instance                                             )
    interpolatorSED =  interpolator                                   (wavelengths                  ,                  sed(:,1),extrapolationType=extrapolationTypeZero                                )
    integrator_     =  integrator                                     (integrandFilteredLuminosity  ,toleranceRelative=4.0d-3  ,integrationRule  =GSL_Integ_Gauss15    ,intervalsMaximum=10000_c_size_t)
    integratorAB_   =  integrator                                     (integrandFilteredLuminosityAB,toleranceRelative=4.0d-3                                                                          )
    do i=1,size(self%filterNames)
       luminosities(i)=+integrator_  %integrate(self%wavelengthRanges(i,1),self%wavelengthRanges(i,2),status=errorStatus) &
            &          /integratorAB_%integrate(self%wavelengthRanges(i,1),self%wavelengthRanges(i,2)                   )
    end do
    deallocate(interpolatorSED)
    deallocate(integrator_    )
    deallocate(integratorAB_  )
    return

  contains

    double precision function integrandFilteredLuminosity(wavelength)
      !!{
      Integrand for the luminosity through a given filter.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      ! The factor of 1/λ appears since we want to integrate F_ν (dν / ν) and dν = -c/λ² dλ. Note that we follow the convention of
      ! Hogg et al. (2002) and assume that the filter response gives the fraction of incident photons received by the detector at
      ! a given wavelength, multiplied by the relative photon response (which will be 1 for a photon-counting detector such as a
      ! CCD, or proportional to the photon energy for a bolometer/calorimeter type detector).
      integrandFilteredLuminosity=+self%filterResponses(i)%interpolator_  %interpolate(wavelength) &
           &                      *                        interpolatorSED%interpolate(wavelength) &
           &                      /                                                    wavelength
      return
    end function integrandFilteredLuminosity

    double precision function integrandFilteredLuminosityAB(wavelength)
      !!{
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
    !!{
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
    !!{
    Return the name of the luminosityStellarFromSED property.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)                            :: self
    double precision                                               , intent(in   )                            :: time
    type            (varying_string                               ), intent(inout), dimension(:), allocatable :: names
    type            (varying_string                               )               , dimension(:), allocatable :: name
    integer                                                                                                   :: i
    
    allocate(names(size(self%filterNames)))
    call self%nodePropertyExtractor_%names(name,time)
    do i=1,size(names)
       names(i)=name(1)//":"//self%filterNames(i)
    end do
    return
  end subroutine luminosityStellarFromSEDNames

  subroutine luminosityStellarFromSEDDescriptions(self,time,descriptions)
    !!{
    Return a description of the luminosityStellarFromSED property.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromSED), intent(inout)                            :: self
    double precision                                               , intent(in   )                            :: time
    type            (varying_string                               ), intent(inout), dimension(:), allocatable :: descriptions
    type            (varying_string                               )               , dimension(:), allocatable :: description
    integer                                                                                                   :: i

    allocate(descriptions(size(self%filterNames)))
    call self%nodePropertyExtractor_%descriptions(description,time)
    do i=1,size(descriptions)
       descriptions(i)="Luminosity in "//self%filterNames(i)//" filter in units of the AB-magnitude system zero-point; derived from: "//description(1)
    end do
    return
  end subroutine luminosityStellarFromSEDDescriptions

  function luminosityStellarFromSEDUnitsInSI(self,time) result(unitsInSI)
    !!{
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
