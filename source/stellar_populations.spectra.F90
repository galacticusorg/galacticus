!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements calculations of stellar population spectra.

module Stellar_Population_Spectra
  !% Implements calculations of stellar population spectra.
  use Abundances_Structure
  use ISO_Varying_String
  !# <include directive="stellarPopulationSpectraMethod" type="moduleUse">
  include 'stellar_populations.spectra.modules.inc'
  !# </include>
  implicit none
  private
  public :: Stellar_Population_Spectrum, Stellar_Population_Spectrum_Tabulation

  ! Flag to indicate if this module has been initialized.
  logical                                                      :: stellarPopulationSpectraInitialized=.false.

  ! Name of cooling time available method used.
  type     (varying_string                          )          :: stellarPopulationSpectraMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Population_Spectrum_Get_Template), pointer :: Stellar_Population_Spectrum_Get    =>null()
  abstract interface
     double precision function Stellar_Population_Spectrum_Get_Template(abundancesStellar,age,wavelength,imfIndex)
       import abundances
       type            (abundances), intent(in   ) :: abundancesStellar
       double precision            , intent(in   ) :: age              , wavelength
       integer                     , intent(in   ) :: imfIndex
     end function Stellar_Population_Spectrum_Get_Template
  end interface
  procedure(Stellar_Population_Tabulation_Get_Template), pointer :: Stellar_Population_Spectrum_Tabulation_Get=>null()
  abstract interface
     subroutine Stellar_Population_Tabulation_Get_Template(imfIndex,agesCount,metallicitiesCount,ages,metallicity)
       integer                                    , intent(in   ) :: imfIndex
       integer                                    , intent(  out) :: agesCount, metallicitiesCount
       double precision, allocatable, dimension(:), intent(  out) :: ages     , metallicity
     end subroutine Stellar_Population_Tabulation_Get_Template
  end interface

contains

  subroutine Stellar_Population_Spectrum_Initialize
    !% Initialize the stellar population spectra module
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.stellarPopulationSpectraInitialized) then
       !$omp critical(Stellar_Population_Spectrum_Initialization)
       if (.not.stellarPopulationSpectraInitialized) then
          ! Get the cooling function method parameter.
          !@ <inputParameter>
          !@   <name>stellarPopulationSpectraMethod</name>
          !@   <defaultValue>Conroy-White-Gunn2009</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of stellar population spectra.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarPopulationSpectraMethod',stellarPopulationSpectraMethod,defaultValue='Conroy-White-Gunn2009')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarPopulationSpectraMethod" type="functionCall" functionType="void">
          !#  <functionArgs>stellarPopulationSpectraMethod,Stellar_Population_Spectrum_Get,Stellar_Population_Spectrum_Tabulation_Get</functionArgs>
          include 'stellar_populations.spectra.inc'
          !# </include>
          if (.not.(associated(Stellar_Population_Spectrum_Get).and.associated(Stellar_Population_Spectrum_Tabulation_Get))) call&
               & Galacticus_Error_Report('Stellar_Population_Spectrum','method ' //char(stellarPopulationSpectraMethod)//' is&
               & unrecognized')
          stellarPopulationSpectraInitialized=.true.
       end if
       !$omp end critical(Stellar_Population_Spectrum_Initialization)
    end if
    return
  end subroutine Stellar_Population_Spectrum_Initialize

  double precision function Stellar_Population_Spectrum(abundancesStellar,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms).
    implicit none
    type            (abundances), intent(in   ) :: abundancesStellar
    double precision            , intent(in   ) :: age              , wavelength
    integer                     , intent(in   ) :: imfIndex

    ! Initialize the module.
    call Stellar_Population_Spectrum_Initialize

    ! Get the spectrum using the selected method.
    Stellar_Population_Spectrum=Stellar_Population_Spectrum_Get(abundancesStellar,age,wavelength,imfIndex)

    return
  end function Stellar_Population_Spectrum

  subroutine Stellar_Population_Spectrum_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    implicit none
    integer                                    , intent(in   ) :: imfIndex
    integer                                    , intent(  out) :: agesCount, metallicitiesCount
    double precision, allocatable, dimension(:), intent(  out) :: age      , metallicity

    ! Initialize the module.
    call Stellar_Population_Spectrum_Initialize

    ! Get the tabulation using the selected method.
    call Stellar_Population_Spectrum_Tabulation_Get(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    return
  end subroutine Stellar_Population_Spectrum_Tabulation

end module Stellar_Population_Spectra
