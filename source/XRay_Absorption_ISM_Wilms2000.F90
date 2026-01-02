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
Contains a program which wraps the {\normalfont \ttfamily dotbvabs} function (which implements the model of \citealt{wilms_absorption_2000}) from
\href{https://heasarc.gsfc.nasa.gov/xanadu/xspec/}{\normalfont \scshape XSpec} to produce a table of X-ray absorption cross-sections in the
\gls{ism}. This program assumes that various files from \href{https://heasarc.gsfc.nasa.gov/xanadu/xspec/}{\normalfont \scshape XSpec} have been
downloaded into the {\normalfont \ttfamily aux/XSpec} folder---usually this program will be run automatically as needed by the {\normalfont \ttfamily
Galacticus::ISMCrossSections} module.
!!}

! Add explicit dependencies on the XSpec files.
!: aux/XSpec/dotbvabs.o aux/XSpec/gphoto.o aux/XSpec/j4save.o aux/XSpec/phfit2.o
!/ exclude

program XRay_Absorption_ISM_Wilms2000
  !!{
  Wraps the {\normalfont \ttfamily dotbvabs} function (which implements the model of \citealt{wilms_absorption_2000}) from
  \href{https://heasarc.gsfc.nasa.gov/xanadu/xspec/}{\normalfont \scshape XSpec} to produce a table of X-ray absorption
  cross-sections in the \gls{ism}. This program assumes that various files from
  \href{https://heasarc.gsfc.nasa.gov/xanadu/xspec/}{\normalfont \scshape XSpec} have been downloaded into the {\normalfont
  \ttfamily aux/XSpec} folder---usually this program will be run automatically as needed by the {\normalfont \ttfamily
  Galacticus::ISMColumnDensity} module.
  !!}
  use :: Atomic_Cross_Sections_Compton, only : Atomic_Cross_Section_Compton
  use :: Dates_and_Times              , only : Formatted_Date_and_Time
  use :: IO_HDF5                      , only : hdf5Object
  use :: Numerical_Constants_Prefixes , only : kilo
  use :: Numerical_Constants_Units    , only : electronVolt
  use :: Numerical_Ranges             , only : Make_Range
  implicit none
  integer                                                                , parameter :: energyCount       =1000
  integer                                                                , parameter :: metallicityCount  =100
  double precision                                                       , parameter :: energyMaximum     =20.0d0       , energyMinimum     =0.01d0     !   keV
  double precision                                                       , parameter :: metallicityMaximum=0.1d0        , metallicityMinimum=1.0d-4
  ! The metallicity of the ISM computed consistently from Wilms et al. abundances and element weights.
  double precision                                                       , parameter :: metallicityIsm    =0.012000314d0
  ! Electrons per hydrogen atom contributed by helium and metals for ISM metallicity.
  double precision                                                       , parameter :: electronsHelium   =0.195447444d0
  double precision                                                       , parameter :: electronsMetals   =0.008383587d0
  double precision            , dimension(0:energyCount                 )            :: energy
  double precision            , dimension(              metallicityCount)            :: metallicity
  real                        , dimension(  energyCount                 )            :: photar                          , photer                   , &
       &                                                                                sigavg                          , siggas                   , &
       &                                                                                siggrains                       , sigmol
  double precision            , dimension(  energyCount,metallicityCount)            :: crossSection
  double precision            , dimension(                            42)            :: parameters
  type            (hdf5Object)                                                       :: myDataset                       , outputFile
  integer                                                                            :: iMetallicity
  double precision                                                                   :: electronNumber

  ! Create array of energies.
  energy(0            )=energyMinimum
  energy(1:energyCount)=Make_Range(energyMinimum     ,energyMaximum     ,energyCount     ,rangeType=rangeTypeLogarithmic)

  ! Create array of metallicities.
  metallicity          =Make_Range(metallicityMinimum,metallicityMaximum,metallicityCount,rangeType=rangeTypeLogarithmic)

  ! Specify values of the parameters to be passed to the dotbvabs function. The default values from Wilms et al. (2000) are used.
  parameters( 1   )=1.000d0 ! Hydrogen column density (units of 10²² cm⁻²).
  parameters( 2   )=1.000d0 ! Helium abundance relative to Milky Way ISM.
  parameters(19   )=0.200d0 ! Fraction of hydrogen in molecular form.
  parameters(20   )=1.000d0 ! Density of dust (units of g/cm³).
  parameters(21   )=0.025d0 ! Minimum thickness of dust.
  parameters(22   )=0.250d0 ! Maximum thickness of dust.
  parameters(23   )=3.500d0 ! Power-law index of dust distribution.
  parameters(42   )=0.000d0 ! Redshift.
  ! Element depletion factors.
  parameters(24:41)=[1.0d0,1.0d0,0.5d0,1.0d0,0.6d0,1.0d0,0.25d0,0.2d0,0.02d0,0.1d0,0.6d0,0.5d0,1.0d0,0.003d0,0.03d0,0.3d0,0.05d0&
       &,0.04d0]

  ! Loop over metallicities.
  do iMetallicity=1,metallicityCount
     ! Set the abundance scaling factors such that metal abundance is scaled by metallicity.
     parameters(3:18)=metallicity(iMetallicity)/metallicityIsm
     ! Call dotbvabs to compute the absorption cross-sections.
     call dotbvabs(real(energy),energyCount,real(parameters),0,photar,photer,siggas,siggrains,sigmol,sigavg)
     ! Store the cross sections.
     crossSection(:,iMetallicity)=dble(siggas+siggrains+sigmol)*1.0d-18
     ! Compute Compton scattering contribution.
     electronNumber=1.0d0+electronsHelium+electronsMetals*metallicity(iMetallicity)/metallicityIsm
     crossSection(:,iMetallicity)=crossSection(:,iMetallicity)+electronNumber*Atomic_Cross_Section_Compton(energy)
  end do

  ! Open the output file.
  call outputFile%openFile('data/atomic/Interstellar_Absorption_Wilms_2000.hdf5',overWrite=.true.,chunkSize=1024_hsize_t,compressionLevel=9)
  ! Write energy table.
  call outputFile%writeDataset(energy(1:energyCount),datasetName="energy",comment="Photon energy in keV",datasetReturned=myDataset)
  call myDataset %writeAttribute(kilo*electronVolt,"unitsInSI")
  call myDataset %close()
  ! Write metallicity table.
  call outputFile%writeDataset(metallicity,datasetName="metallicity",comment="Metallicity")
  ! Write crossSection table.
  call outputFile%writeDataset(crossSection,datasetName="crossSection",comment="Absorption cross section in cm²",datasetReturned=myDataset)
  call myDataset %writeAttribute(1.0d-4,"unitsInSI")
  call myDataset %close()
  ! Add meta-data.
  call outputFile%writeAttribute("Created by Galacticus using dotbvabs function from XSpec","description")
  call outputFile%writeAttribute(Formatted_Date_and_Time(),"timestamp")
  call outputFile%writeAttribute("Wilms, Allen & McCray (2000; ApJ, 542, 914; http://adsabs.harvard.edu/abs/2000ApJ...542..914W)","source")
  ! Close the output file.
  call outputFile%close()
end program XRay_Absorption_ISM_Wilms2000

subroutine xwrite(msg,i)
  !!{
  Message display function required by {\normalfont \ttfamily dotbvabs}.
  !!}
  implicit none
  character(len=*), intent(in   ) :: msg
  integer         , intent(in   ) :: i

  write (0,*) msg
  return
end subroutine xwrite

subroutine xermsg(a,b,c,i,j)
  !!{
  Error message function required by {\normalfont \ttfamily dotbvabs}.
  !!}
  use, intrinsic :: ISO_Fortran_Env, only : output_unit
  use :: Error, only : Error_Report
  implicit none
  character(len=*), intent(in   ) :: a, b, c
  integer         , intent(in   ) :: i, j

  write (output_unit,*) a
  write (output_unit,*) b
  write (output_unit,*) c
  call Error_Report('error thrown by XSpec functions'//{introspection:location})
  return
end subroutine xermsg

real function fgabnd(c)
  !!{
  Function to return the abundance (relative to hydrogen) of elements. Required by {\normalfont \ttfamily dotbvabs}.
  !!}
  use :: Error, only : Error_Report
  implicit none
  character(len=2), intent(in   ) :: c

  fgabnd=0.0
  if (c == 'H ') fgabnd=12.00 ! H
  if (c == 'He') fgabnd=10.99 ! He
  if (c == 'C ') fgabnd= 8.38 ! C
  if (c == 'N ') fgabnd= 7.88 ! N
  if (c == 'O ') fgabnd= 8.69 ! O
  if (c == 'Ne') fgabnd= 7.94 ! Ne
  if (c == 'Na') fgabnd= 6.16 ! Na
  if (c == 'Mg') fgabnd= 7.40 ! Mg
  if (c == 'Al') fgabnd= 6.33 ! Al
  if (c == 'Si') fgabnd= 7.27 ! Si
  if (c == 'S ') fgabnd= 7.09 ! S
  if (c == 'Cl') fgabnd= 5.12 ! Cl
  if (c == 'Ar') fgabnd= 6.41 ! Ar
  if (c == 'Ca') fgabnd= 6.20 ! Ca
  if (c == 'Cr') fgabnd= 5.51 ! Cr
  if (c == 'Fe') fgabnd= 7.43 ! Fe
  if (c == 'Co') fgabnd= 4.92 ! Co
  if (c == 'Ni') fgabnd= 6.05 ! Ni
  if (fgabnd == 0.0) call Error_Report('unknown element'//{introspection:location})
  fgabnd=10.0**(fgabnd-12.00)
  return
end function fgabnd
