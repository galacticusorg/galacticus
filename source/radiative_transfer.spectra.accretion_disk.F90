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
  
  use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass

  !# <radiativeTransferSpectrum name="radiativeTransferSpectrumAccretionDisk">
  !#  <description>A photon spectrum class for accretionDisk spectrums.</description>
  !# </radiativeTransferSpectrum>
  type, extends(radiativeTransferSpectrumClass) :: radiativeTransferSpectrumAccretionDisk
     !% Implementation of a black body spectrum for radiative transfer calculations.
     private
     class           (accretionDiskSpectraClass), pointer :: accretionDiskSpectra_
     double precision                                     :: massBlackHole        , accretionRateEddington, &
          &                                                  accretionRate
   contains
     final     ::               accretionDiskDestructor
     procedure :: luminosity => accretionDiskLuminosity
     procedure :: spectrum   => accretionDiskSpectrum
  end type radiativeTransferSpectrumAccretionDisk
  
  interface radiativeTransferSpectrumAccretionDisk
     !% Constructors for the {\normalfont \ttfamily accretionDisk} radiative transfer spectrum class.
     module procedure accretionDiskConstructorParameters
     module procedure accretionDiskConstructorInternal
  end interface radiativeTransferSpectrumAccretionDisk
      
contains
      
  function accretionDiskConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily accretionDisk} radiative transfer spectrum class which takes a parameter set as
    !% input.
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferSpectrumAccretionDisk)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (accretionDiskSpectraClass             ), pointer       :: accretionDiskSpectra_
    double precision                                                        :: massBlackHole        , accretionRateEddington
    
    !# <inputParameter>
    !#   <name>massBlackHole</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d6</defaultValue>
    !#   <description>The mass of the black hole at the center of the accretion disk.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>accretionRateEddington</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-1</defaultValue>
    !#   <description>Accretion rate onto the black hole in units of the Eddington rate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="accretionDiskSpectra" name="accretionDiskSpectra_" source="parameters"/>
    self=radiativeTransferSpectrumAccretionDisk(massBlackHole,accretionRateEddington,accretionDiskSpectra_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="accretionDiskSpectra_"/>
    return
  end function accretionDiskConstructorParameters

  function accretionDiskConstructorInternal(massBlackHole,accretionRateEddington,accretionDiskSpectra_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily accretionDisk} radiative transfer spectrum class.
    use :: Numerical_Constants_Astronomical, only : gigaYear
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : gravitationalConstant , speedLight, thomsonCrossSection
    implicit none
    type            (radiativeTransferSpectrumAccretionDisk)                        :: self
    class           (accretionDiskSpectraClass             ), intent(in   ), target :: accretionDiskSpectra_
    double precision                                        , intent(in   )         :: massBlackHole        , accretionRateEddington
    !# <constructorAssign variables="massBlackHole, accretionRateEddington, *accretionDiskSpectra_"/>

    self%accretionRate=+4.0d0                  &
         &             *Pi                     &
         &             *gravitationalConstant  &
         &             *massBlackHole          &
         &             *massHydrogenAtom       &
         &             *gigaYear               &
         &             /thomsonCrossSection    &
         &             /speedLight             &
         &             *accretionRateEddington
    return
  end function accretionDiskConstructorInternal

  subroutine accretionDiskDestructor(self)
    !% Destructor for the {\normalfont \ttfamily accretionDisk} radiative transfer spectrum class.
    implicit none
    type(radiativeTransferSpectrumAccretionDisk), intent(inout) :: self

    !# <objectDestructor name="self%accretionDiskSpectra_"/>
    return
  end subroutine accretionDiskDestructor

  double precision function accretionDiskLuminosity(self,wavelengthMinimum,wavelengthMaximum)
    !% Compute the luminosity in the given wavelength range for an accretion disk.
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    class           (radiativeTransferSpectrumAccretionDisk), intent(inout) :: self
    double precision                                        , intent(in   ) :: wavelengthMinimum    , wavelengthMaximum
    type            (fgsl_function                         )                :: integrandFunction
    type            (fgsl_integration_workspace            )                :: integrationWorkspace

    accretionDiskLuminosity=+Integrate(                                        &
         &                                               wavelengthMinimum   , &
         &                                               wavelengthMaximum   , &
         &                                               integrand           , &
         &                                               integrandFunction   , &
         &                                               integrationWorkspace, &
         &                             toleranceAbsolute=0.0d+0              , &
         &                             toleranceRelative=1.0d-2                &
         &                            )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
    
  contains

    double precision function integrand(wavelength)
      !% Integrand over black body spectrum.
      implicit none
      double precision, intent(in   ) :: wavelength

      integrand=self%spectrum(wavelength)
      return
    end function integrand

  end function accretionDiskLuminosity

  double precision function accretionDiskSpectrum(self,wavelength)
    !% Return the spectrum of the accretion disk.
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Units   , only : angstromsPerMeter
    implicit none
    class           (radiativeTransferSpectrumAccretionDisk), intent(inout) :: self
    double precision                                        , intent(in   ) :: wavelength

    accretionDiskSpectrum=+self%accretionDiskSpectra_%spectrum(self%accretionRate,1.0d0,wavelength) &
         &                *speedLight                                                               &
         &                *angstromsPerMeter                                                        &
         &                /wavelength**2
    return
  end function accretionDiskSpectrum
