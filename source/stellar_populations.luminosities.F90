!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of stellar population luminosities in the AB magnitude system.

module Stellar_Population_Luminosities
  !% Implements calculations of stellar population luminosities in the AB magnitude system.
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  use Abundances_Structure
  private
  public :: Stellar_Population_Luminosity

  type luminosityTable
     !% Structure for holding tables of simple stellar population luminosities.
     integer                                         :: agesCount,metallicitiesCount
     logical,          allocatable, dimension(:)     :: isTabulated
     double precision, allocatable, dimension(:)     :: age,metallicity
     double precision, allocatable, dimension(:,:,:) :: luminosity
     ! Interpolation structures.
     logical                                         :: resetAge=.true., resetMetallicity=.true.
     type(fgsl_interp_accel)                         :: interpolationAcceleratorAge,interpolationAcceleratorMetallicity
    end type luminosityTable
  
  ! Array of simple stellar population luminosity tables.
  type(luminosityTable), allocatable, dimension(:) :: luminosityTables

  ! Module global variables used in integrations.
  double precision          :: ageTabulate,redshiftTabulate
  integer                   :: filterIndexTabulate,imfIndexTabulate
  type(abundancesStructure) :: abundancesTabulate
  !$omp threadprivate(ageTabulate,redshiftTabulate,abundancesTabulate,filterIndexTabulate,imfIndexTabulate)

  ! Flag indicating if this module has been initialized yet.
  logical                   :: moduleInitialized=.false.
  
  ! Tolerance used in integrations.
  double precision          :: stellarPopulationLuminosityIntegrationToleranceRelative

contains
  
  function Stellar_Population_Luminosity(luminosityIndex,filterIndex,imfIndex,abundances,age,redshift)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\tt abundances} and {\tt age} drawn from IMF
    !% specified by {\tt imfIndex} and observed through the filter specified by {\tt filterIndex}.
    use Memory_Management
    use Abundances_Structure
    use Stellar_Population_Spectra
    use Instruments_Filters
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    use Input_Parameters
    implicit none
    integer,                   intent(in)                                    :: luminosityIndex(:),filterIndex(:),imfIndex
    double precision,          intent(in)                                    :: age(:),redshift(:)
    type(abundancesStructure), intent(in)                                    :: abundances
    double precision,                       dimension(size(luminosityIndex)) :: Stellar_Population_Luminosity
    type(luminosityTable),     allocatable, dimension(:)                     :: luminosityTablesTemporary
    double precision,          allocatable, dimension(:,:,:)                 :: luminosityTemporary
    logical,                   allocatable, dimension(:)                     :: isTabulatedTemporary
    double precision,                       dimension(2)                     :: wavelengthRange
    double precision,                       dimension(0:1)                   :: hAge,hMetallicity
    integer                                                                  :: iLuminosity,iAge,iMetallicity,jAge,jMetallicity
    logical                                                                  :: computeTable
    double precision                                                         :: normalization,metallicity,ageLast
    type(c_ptr)                                                              :: parameterPointer                         
    type(fgsl_function)                                                      :: integrandFunction
    type(fgsl_integration_workspace)                                         :: integrationWorkspace

    ! Determine if we have created space for this IMF yet.
    !$omp critical (Luminosity_Tables_Initialize)
    if (.not.moduleInitialized) then
       ! Read the parameter controlling integration tolerance.
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityIntegrationToleranceRelative</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The relative tolerance used when integrating the flux of stellar populations through filters.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityIntegrationToleranceRelative',stellarPopulationLuminosityIntegrationToleranceRelative,defaultValue=1.0d-3)

       ! Flag that this module is now initialized.
       moduleInitialized=.true.
    end if

    if (allocated(luminosityTables)) then
       if (size(luminosityTables) < imfIndex) then
          call Move_Alloc(luminosityTables,luminosityTablesTemporary)
          allocate(luminosityTables(imfIndex))
          luminosityTables(1:size(luminosityTablesTemporary))=luminosityTablesTemporary
          deallocate(luminosityTablesTemporary)
          call Memory_Usage_Record(sizeof(luminosityTables(1)),blockCount=0)
       end if
    else
       allocate(luminosityTables(imfIndex))
       call Memory_Usage_Record(sizeof(luminosityTables))
    end if
    
    ! Determine if we have tabulated luminosities for this luminosityIndex in this IMF yet.
    do iLuminosity=1,size(luminosityIndex)
       if (allocated(luminosityTables(imfIndex)%isTabulated)) then
          if (size(luminosityTables(imfIndex)%isTabulated) >= luminosityIndex(iLuminosity)) then
             computeTable=.not.luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))
          else
             call Move_Alloc (luminosityTables(imfIndex)%isTabulated,isTabulatedTemporary)
             call Move_Alloc (luminosityTables(imfIndex)%luminosity ,luminosityTemporary )
             call Alloc_Array(luminosityTables(imfIndex)%isTabulated,[luminosityIndex(iLuminosity)])
             call Alloc_Array(luminosityTables(imfIndex)%luminosity ,[luminosityIndex(iLuminosity)&
                  &,luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%metallicitiesCount])
             luminosityTables(imfIndex)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
             luminosityTables(imfIndex)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndex(iLuminosity))=.false.
             luminosityTables(imfIndex)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
             call Dealloc_Array(isTabulatedTemporary)
             call Dealloc_Array(luminosityTemporary)
             computeTable=.true.
          end if
       else
          call Alloc_Array(luminosityTables(imfIndex)%isTabulated,[luminosityIndex(iLuminosity)])
          luminosityTables(imfIndex)%isTabulated=.false.
          ! Since we have not yet tabulated any luminosities yet for this IMF, we need to get a list of suitable metallicities and
          ! ages at which to tabulate.
          call Stellar_Population_Spectrum_Tabulation(imfIndex,luminosityTables(imfIndex)%agesCount &
               &,luminosityTables(imfIndex)%metallicitiesCount,luminosityTables(imfIndex)%age&
               &,luminosityTables(imfIndex)%metallicity)
          where (luminosityTables(imfIndex)%metallicity > 0.0d0)
             luminosityTables(imfIndex)%metallicity=dlog10(luminosityTables(imfIndex)%metallicity/metallicitySolar)
          elsewhere
             luminosityTables(imfIndex)%metallicity=logMetallicityZero
          end where
          call Alloc_Array(luminosityTables(imfIndex)%luminosity,[luminosityIndex(iLuminosity)&
               &,luminosityTables(imfIndex)%agesCount ,luminosityTables(imfIndex)%metallicitiesCount])
          computeTable=.true.
       end if
       
       ! If we haven't, do so now.
       if (computeTable) then
          
          ! Get wavelength extent of the filter.
          wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
          
          ! Integrate over the wavelength range.
          filterIndexTabulate=filterIndex(iLuminosity)
          imfIndexTabulate=imfIndex
          redshiftTabulate=redshift(iLuminosity)
          do iAge=1,luminosityTables(imfIndex)%agesCount
             ageTabulate=luminosityTables(imfIndex)%age(iAge)
             do iMetallicity=1,luminosityTables(imfIndex)%metallicitiesCount
                call abundancesTabulate%metallicitySet(luminosityTables(imfIndex)%metallicity(iMetallicity) &
                     &,metallicityType=logarithmicByMassSolar)
                luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),iAge,iMetallicity) &
                     &=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand,parameterPointer &
                     &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
                     &=stellarPopulationLuminosityIntegrationToleranceRelative,integrationRule =FGSL_Integ_Gauss15,maxIntervals&
                     &=10000)
                call Integrate_Done(integrandFunction,integrationWorkspace)
             end do
          end do
          
          ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
          normalization=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand_AB,parameterPointer &
               &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
               &=stellarPopulationLuminosityIntegrationToleranceRelative)
          call Integrate_Done(integrandFunction,integrationWorkspace)
          
          ! Normalize the luminosity.
          luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:) &
               &=luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:)/normalization
          
          ! Flag that calculations have been performed for this filter.
          luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))=.true.
       end if
    end do
 
    ! Get interpolation in metallicity.
    metallicity=Abundances_Get_Metallicity(abundances,metallicityType=logarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(imfIndex)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(imfIndex)%metallicity(luminosityTables(imfIndex)%metallicitiesCount)) then
       iMetallicity=luminosityTables(imfIndex)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(imfIndex)%metallicitiesCount,luminosityTables(imfIndex)%metallicity &
            &,luminosityTables(imfIndex)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(imfIndex)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%metallicitiesCount &
            &,luminosityTables(imfIndex)%metallicity ,iMetallicity,metallicity)
    end if
    
    ! Do the interpolation.
    Stellar_Population_Luminosity(:)=0.0d0
    do iLuminosity=1,size(luminosityIndex)
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosity) > 0.0d0) then
          ! Get interpolation in age if the age for this luminosity differs from the previous one.
          if (iLuminosity == 1 .or. age(iLuminosity) /= ageLast) then
             ! Check for out of range age.
             if (age(iLuminosity) > luminosityTables(imfIndex)%age(luminosityTables(imfIndex)%agesCount)) call&
                  & Galacticus_Error_Report('Stellar_Population_Luminosity','age exceeds the maximum tabulated')
             iAge=Interpolate_Locate(luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%age &
                  &,luminosityTables(imfIndex)%interpolationAcceleratorAge,age(iLuminosity),luminosityTables(imfIndex)%resetAge)
             hAge=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%age,iAge&
                  &,age(iLuminosity))
             ageLast=age(iLuminosity)
          end if
          do jAge=0,1  
             do jMetallicity=0,1
                Stellar_Population_Luminosity(iLuminosity)=Stellar_Population_Luminosity(iLuminosity)&
                     &+luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),iAge +jAge,iMetallicity+jMetallicity)&
                     &*hAge(jAge)*hMetallicity(jMetallicity)
             end do
          end do
       end if
    end do
    !$omp end critical (Luminosity_Tables_Initialize)

    ! Prevent interpolation from returning negative fluxes.
    Stellar_Population_Luminosity=max(Stellar_Population_Luminosity,0.0d0)
    
    return
  end function Stellar_Population_Luminosity
  
  function Filter_Luminosity_Integrand(wavelength,parameterPointer) bind(c)
    !% Integrand for the luminosity through a given filter.
    use Stellar_Population_Spectra
    use Stellar_Population_Spectra_Postprocess
    use Instruments_Filters
    implicit none
    real(c_double)         :: Filter_Luminosity_Integrand
    real(c_double),  value :: wavelength
    type(c_ptr),     value :: parameterPointer
    double precision       :: wavelengthRedshifted

    ! If this luminosity is for a redshifted spectrum, then we shift wavelength at which we sample the stellar population spectrum
    ! to be a factor of (1+z) smaller. We therefore integrate over the stellar SED at shorter wavelengths, since these will be
    ! shifted into the filter by z=0. Factor of 1/wavelength appears since we want to integrate F_nu (dnu / nu) and dnu =
    ! -c/lambda^2 dlambda. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the
    ! fraction of incident photons received by the detector at a given wavelength, multiplied by the relative photon response
    ! (which will be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a
    ! bolometer/calorimeter type detector).
    wavelengthRedshifted=wavelength/(1.0d0+redshiftTabulate)
    Filter_Luminosity_Integrand=Filter_Response(filterIndexTabulate,wavelength)*Stellar_Population_Spectrum(abundancesTabulate &
         &,ageTabulate,wavelengthRedshifted,imfIndexTabulate)*Stellar_Population_Spectrum_PostProcess(wavelengthRedshifted,redshiftTabulate)/wavelength
    return
  end function Filter_Luminosity_Integrand
  
  function Filter_Luminosity_Integrand_AB(wavelength,parameterPointer) bind(c)
    !% Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
    use Stellar_Population_Spectra
    use Instruments_Filters
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    implicit none
    real(c_double)              :: Filter_Luminosity_Integrand_AB
    real(c_double),   value     :: wavelength
    type(c_ptr),      value     :: parameterPointer
    ! Luminosity of a zeroth magintude (AB) source in Solar luminosities per Hz.
    double precision, parameter :: luminosityZeroPointABSolar=luminosityZeroPointAB/luminositySolar

    Filter_Luminosity_Integrand_AB=Filter_Response(filterIndexTabulate,wavelength)*luminosityZeroPointABSolar/wavelength
    return
  end function Filter_Luminosity_Integrand_AB
  
end module Stellar_Population_Luminosities
