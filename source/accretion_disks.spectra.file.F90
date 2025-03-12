!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of the accretion disk spectra class for tabulated spectra read from file.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass
  use :: Numerical_Interpolation   , only : interpolator

  !![
  <accretionDiskSpectra name="accretionDiskSpectraFile">
   <description>Accretion disk spectra are interpolated from tables read from file.</description>
   <runTimeFileDependencies paths="fileName"/>
  </accretionDiskSpectra>
  !!]

  type, extends(accretionDiskSpectraClass) :: accretionDiskSpectraFile
     !!{
     An accretion disk spectra class which interpolates in spectra read from file.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer                     :: blackHoleAccretionRate_ => null()
     class           (accretionDisksClass        ), pointer                     :: accretionDisks_         => null()
     type            (varying_string             )                              :: fileName
     double precision                             , allocatable, dimension(:  ) :: luminosity                       , wavelength
     double precision                             , allocatable, dimension(:,:) :: SED
     type            (interpolator               )                              :: interpolatorLuminosity           , interpolatorWavelength
   contains
     !![
     <methods>
       <method description="Load a file of AGN spectra." method="loadFile" />
     </methods>
     !!]
     final     ::                     fileDestructor
     procedure :: spectrumNode     => fileSpectrumNode
     procedure :: spectrumMassRate => fileSpectrumMassRate
     procedure :: loadFile         => fileLoadFile
  end type accretionDiskSpectraFile

  interface accretionDiskSpectraFile
     !!{
     Constructors for the {\normalfont \ttfamily file} accretion disk spectra class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface accretionDiskSpectraFile

  ! Supported file format.
  integer :: fileFormatCurrent=2

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily file} accretion disk spectra class which takes a parameter set as input.
    !!}
    implicit none
    type (accretionDiskSpectraFile   )                :: self
    type (inputParameters            ), intent(inout) :: parameters
    type (varying_string             )                :: fileName
    class(blackHoleAccretionRateClass), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass        ), pointer       :: accretionDisks_

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of a file from which to read tabulated spectra of accretion disks.</description>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=accretionDiskSpectraFile(char(fileName),blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily file} accretion disk spectra class.
    !!}
    implicit none
    type     (accretionDiskSpectraFile   )                        :: self
    character(len=*                      )        , intent(in   ) :: fileName
    class    (blackHoleAccretionRateClass), target, intent(in   ) :: blackHoleAccretionRate_
    class    (accretionDisksClass        ), target, intent(in   ) :: accretionDisks_
    !![
    <constructorAssign variables="fileName, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    call self%loadFile(fileName)
    return
  end function fileConstructorInternal

  subroutine fileDestructor(self)
    !!{
    Default destructor for the {\normalfont \ttfamily file} accretion disk spectra class.
    !!}
    implicit none
    type(accretionDiskSpectraFile), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    return
  end subroutine fileDestructor

  subroutine fileLoadFile(self,fileName)
    !!{
    Load a file of AGN spectra.
    !!}
    use :: Error       , only : Error_Report
    use :: HDF5_Access , only : hdf5Access
    use :: IO_HDF5     , only : hdf5Object
    use :: Table_Labels, only : extrapolationTypeZero, extrapolationTypeFix
    implicit none
    class    (accretionDiskSpectraFile), intent(inout) :: self
    character(len=*                   ), intent(in   ) :: fileName
    integer                                            :: fileFormatFile
    type     (hdf5Object              )                :: spectraFile

    ! Open the file.
    !$ call hdf5Access%set()
    spectraFile=hdf5Object(fileName,readOnly=.true.)
    ! Check file format.
    call spectraFile%readAttribute('fileFormat',fileFormatFile,allowPseudoScalar=.true.)
    if (fileFormatFile /= fileFormatCurrent) call Error_Report('file format mismatch'//{introspection:location})
    ! Read datasets.
    call spectraFile%readDataset('wavelength'          ,self%wavelength)
    call spectraFile%readDataset('bolometricLuminosity',self%luminosity)
    call spectraFile%readDataset('SED'                 ,self%SED       )
    !$ call hdf5Access%unset()
    ! Convert luminosities to logarithmic form for interpolation.
    self%luminosity=log(self%luminosity)
    ! Build interpolators.
    self%interpolatorLuminosity=interpolator(self%luminosity,extrapolationType=extrapolationTypeFix )
    self%interpolatorWavelength=interpolator(self%wavelength,extrapolationType=extrapolationTypeZero)
    return
  end subroutine fileLoadFile

  double precision function fileSpectrumNode(self,node,wavelength)
    !!{
    Return the accretion disk spectrum for tabulated spectra.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (accretionDiskSpectraFile), intent(inout)  :: self
    type            (treeNode                ), intent(inout)  :: node
    double precision                          , intent(in   )  :: wavelength
    class           (nodeComponentBlackHole  ), pointer        :: blackHole
    double precision                                           :: rateAccretionSpheroid, rateAccretionHotHalo           , &
         &                                                        rateAccretion        , rateAccretionNuclearStarCluster, &
         &                                                        efficiencyRadiative

    blackHole => node%blackHole()
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion      =+                    rateAccretionSpheroid                                                         &
         &              +                    rateAccretionHotHalo                                                          &
         &              +                    rateAccretionNuclearStarCluster
    efficiencyRadiative=self%accretionDisks_%efficiencyRadiative  (blackHole,rateAccretion                               )
    fileSpectrumNode   =self                %spectrum             (          rateAccretion,efficiencyRadiative,wavelength)
    return
  end function fileSpectrumNode

  double precision function fileSpectrumMassRate(self,accretionRate,efficiencyRadiative,wavelength)
    !!{
    Return the accretion disk spectrum for tabulated spectra.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : gigaYear  , luminositySolar   , massSolar
    use            :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    class           (accretionDiskSpectraFile), intent(inout)  :: self
    double precision                          , intent(in   )  :: accretionRate       , efficiencyRadiative, &
         &                                                        wavelength
    double precision                          , dimension(0:1) :: hLuminosity         , hWavelength
    integer         (c_size_t                )                 :: iLuminosity         , iWavelength        , &
         &                                                        jLuminosity         , jWavelength
    double precision                                           :: luminosityBolometric

    ! Initialize to zero spectrum.
    fileSpectrumMassRate=0.0d0
    ! Get the bolometric luminosity.
    luminosityBolometric=+efficiencyRadiative    &
         &               *accretionRate          &
         &               *massSolar              &
         &               *speedLight         **2 &
         &               /gigaYear               &
         &               /luminositySolar
    ! Return on non-positive luminosities.
    if (luminosityBolometric <= 0.0d0) return
    ! Assume zero flux outside of the tabulated wavelength range.
    if     (                                                     &
         &   wavelength < self%wavelength(                   1 ) &
         & .or.                                                  &
         &   wavelength > self%wavelength(size(self%wavelength)) &
         & ) return
    ! Get the interpolating factors.
    call self%interpolatorLuminosity%linearFactors(log(luminosityBolometric),iLuminosity,hLuminosity)
    call self%interpolatorWavelength%linearFactors(    wavelength           ,iWavelength,hWavelength)
    ! Do the interpolation.
    do jLuminosity=0,1
       do jWavelength=0,1
          fileSpectrumMassRate=+fileSpectrumMassRate                      &
               &               +self        %SED(                         &
               &                                 iWavelength+jWavelength, &
               &                                 iLuminosity+jLuminosity  &
               &                                )                         &
               &               *hLuminosity     (            jLuminosity) &
               &               *hWavelength     (            jWavelength)
       end do
    end do
    ! Prevent interpolation from returning negative fluxes.
    fileSpectrumMassRate=max(fileSpectrumMassRate,0.0d0)
    return
  end function fileSpectrumMassRate

