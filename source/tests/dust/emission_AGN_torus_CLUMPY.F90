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
Contains a program which tests the CLUMPY AGN torus dust emission class.
!!}

program Test_Dust_Emission_AGN_Torus_CLUMPY
  !!{RST
  Tests the :galacticus-class:`dustEmissionSpectrumAGNTorusCLUMPY` dust emission spectrum class, using a small synthetic
  tabulation written here in the format of the CLUMPY tabulation, so that the expected results are known exactly. Each
  parameter is tabulated at three values, and the emission at each wavelength is a linear function of :math:`N_0`,
  :math:`\ln Y`, :math:`q`, :math:`\sigma` and :math:`\ln \tau_\mathrm{V}`, with a different dependence on each at each
  wavelength. Interpolation of the tabulation, linear in :math:`N_0`, :math:`q` and :math:`\sigma` and in the logarithms
  of :math:`Y` and :math:`\tau_\mathrm{V}`, then reproduces that function exactly, wherever the parameters fall, which
  would not be so were the parameters confused with each other, or interpolated on the wrong scale.

  The checks are made against quantities computed here, independently of the class:

  * the emission spectrum integrates to the luminosity absorbed;
  * the ratios of :math:`\nu L_\nu` at the tabulated wavelengths are those of the function, both for parameters on
    tabulated values and for parameters between them;
  * the emission is zero outside the tabulated wavelengths.
  !!}
  use :: Display                     , only : displayVerbositySet               , verbosityLevelStandard
  use :: Dust_Emission_Spectra       , only : dustEmissionSpectrumAGNTorusCLUMPY
  use :: Error                       , only : Error_Handler_Register
  use :: Events_Hooks                , only : eventsHooksInitialize
  use :: File_Utilities              , only : Directory_Make
  use :: HDF5_Access                 , only : hdf5Access
  use :: Input_Paths                 , only : inputPath                         , pathTypeDataDynamic
  use :: IO_HDF5                     , only : hdf5File                          , ioHDF5AccessInitialize
  use :: ISO_Varying_String          , only : assignment(=)                     , char                  , operator(//)        , var_str          , &
       &                                      varying_string
  use :: Numerical_Constants_Physical, only : speedLight
  use :: Numerical_Constants_Units   , only : metersToAngstroms  , micronsToAngstroms
  use :: Unit_Tests                  , only : Assert                            , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! The synthetic tabulation: wavelengths (in μm), and the tabulated values of each parameter.
  integer                                             , parameter                                                    :: countWavelength    =6     , countParameter     =3
  double precision                                    , parameter, dimension(countWavelength)                        :: wavelengthTable    =[0.01d0,0.1d0,1.0d0,10.0d0,100.0d0,1000.0d0]
  double precision                                    , parameter, dimension(countParameter )                        :: countCloudsTable   =[ 1.0d0, 5.0d0, 10.0d0], &
       &                                                                                                                ratioRadiiTable    =[ 5.0d0,20.0d0,100.0d0], &
       &                                                                                                                exponentTable      =[ 0.0d0, 1.0d0,  2.0d0], &
       &                                                                                                                widthTable         =[15.0d0,30.0d0, 45.0d0], &
       &                                                                                                                depthOpticalTable  =[10.0d0,40.0d0,160.0d0]
  ! Coefficients of the emission, λL_λ per unit luminosity, at each wavelength: a constant, and the dependence on each
  ! of N₀, ln Y, q, σ and ln τ_V.
  double precision                                    , parameter, dimension(countWavelength)                        :: emissionBase       =[ 0.020d0, 0.050d0, 0.300d0, 1.000d0, 0.400d0,0.010d0], &
       &                                                                                                                emissionCountClouds=[ 0.001d0, 0.002d0, 0.010d0, 0.030d0, 0.020d0,0.001d0], &
       &                                                                                                                emissionRatioRadii =[ 0.000d0, 0.001d0,-0.010d0, 0.050d0, 0.100d0,0.005d0], &
       &                                                                                                                emissionExponent   =[ 0.002d0,-0.003d0, 0.020d0,-0.050d0, 0.080d0,0.002d0], &
       &                                                                                                                emissionWidth      =[ 1.0d-4 , 2.0d-4 , 1.0d-3 , 4.0d-3 ,-2.0d-3 ,1.0d-4 ], &
       &                                                                                                                emissionDepth      =[-0.002d0,-0.005d0,-0.020d0, 0.030d0, 0.050d0,0.001d0]
  ! Number of points in the fine grid used for the independent integral.
  integer                                             , parameter                                                    :: countWavelengthFine=200001
  double precision                                    , parameter                                                    :: luminosityAbsorbed =1.0d10
  type            (dustEmissionSpectrumAGNTorusCLUMPY), pointer                                                      :: emission_
  type            (varying_string                    )                                                               :: fileName                  , pathDirectory
  double precision                                    , dimension(countWavelength,countParameter,countParameter, &
       &                                                          countParameter ,countParameter,countParameter)     :: emissionTable
  double precision                                    , dimension(countWavelengthFine)                               :: wavelengthsFine           , luminosity
  double precision                                    , dimension(countWavelength    )                               :: emissionExpected          , emissionComputed
  double precision                                    , dimension(5                  )                               :: parametersModel
  double precision                                                                                                   :: integral
  integer                                                                                                            :: i                         , j               , &
       &                                                                                                                k                         , l               , &
       &                                                                                                                m                         , n               , &
       &                                                                                                                iCase
  character       (len=128                           )                                                               :: label

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call eventsHooksInitialize ()
  call ioHDF5AccessInitialize()
  ! Build the synthetic tabulation, indexed (wavelength, optical depth, width, exponent, ratio of radii, number of
  ! clouds).
  do n=1,countParameter
     do m=1,countParameter
        do l=1,countParameter
           do k=1,countParameter
              do j=1,countParameter
                 emissionTable(:,j,k,l,m,n)=emissionFunction([countCloudsTable(n),ratioRadiiTable(m),exponentTable(l),widthTable(k),depthOpticalTable(j)])
              end do
           end do
        end do
     end do
  end do
  fileName     ='CLUMPY_unitTest.hdf5'
  pathDirectory=inputPath(pathTypeDataDynamic)//'dust/agnTorus'
  call Directory_Make(pathDirectory)
  !$ call hdf5Access%set()
  hdf5WriteScope : block
    type(hdf5File) :: file
    file=hdf5File(char(pathDirectory//'/'//fileName),overWrite=.true.)
    call file%writeDataset(wavelengthTable  ,'wavelength'           )
    call file%writeDataset(countCloudsTable ,'countCloudsEquatorial')
    call file%writeDataset(ratioRadiiTable  ,'ratioRadii'           )
    call file%writeDataset(exponentTable    ,'exponentRadial'       )
    call file%writeDataset(widthTable       ,'widthAngular'         )
    call file%writeDataset(depthOpticalTable,'depthOpticalCloudV'   )
    call file%writeDataset(emissionTable    ,'emissionDustMean'     )
  end block hdf5WriteScope
  !$ call hdf5Access%unset()
  ! Fine grid of wavelength, in Å, spanning the tabulated wavelengths.
  do i=1,countWavelengthFine
     wavelengthsFine(i)=wavelengthTable(1)*micronsToAngstroms*exp(log(wavelengthTable(countWavelength)/wavelengthTable(1))*dble(i-1)/dble(countWavelengthFine-1))
  end do
  ! Parameters on tabulated values, between them, and at the ends of the tabulated ranges.
  do iCase=1,3
     select case (iCase)
     case (1)
        parametersModel=[5.0d0, 20.0d0,1.0d0,30.0d0, 40.0d0]
     case (2)
        parametersModel=[3.0d0, 40.0d0,1.5d0,20.0d0, 30.0d0]
     case (3)
        parametersModel=[1.0d0,100.0d0,2.0d0,15.0d0,160.0d0]
     end select
     write (label,'(a,5(1x,f5.1))') "Dust emission: CLUMPY AGN torus, N0 Y q sigma tauV =",parametersModel
     call Unit_Tests_Begin_Group(trim(label))
     allocate(emission_)
     !![
     <referenceConstruct object="emission_" constructor="dustEmissionSpectrumAGNTorusCLUMPY(fileName,var_str('none'),parametersModel(1),parametersModel(4),parametersModel(5),parametersModel(2),parametersModel(3))"/>
     !!]
     ! Integrate L_ν over frequency on a fine grid, as ∫ νL_ν dln λ.
     luminosity=emission_%luminosity(wavelengthsFine,luminosityAbsorbed,0.0d0,13.8d0)
     integral  =trapezoid(log(wavelengthsFine),luminosity*speedLight*metersToAngstroms/wavelengthsFine)
     call Assert('emission integrates to the luminosity absorbed',integral,luminosityAbsorbed,relTol=1.0d-6)
     ! The shape: νL_ν at the tabulated wavelengths, relative to that at 10 μm, is in the ratio of the function.
     emissionExpected=emissionFunction(parametersModel)
     emissionComputed=emission_%luminosity(wavelengthTable*micronsToAngstroms,luminosityAbsorbed,0.0d0,13.8d0)*speedLight*metersToAngstroms/(wavelengthTable*micronsToAngstroms)
     call Assert('emission has the shape of the function',emissionComputed/emissionComputed(4),emissionExpected/emissionExpected(4),relTol=1.0d-9)
     luminosity(1:2)=emission_%luminosity([0.5d0*wavelengthTable(1),2.0d0*wavelengthTable(countWavelength)]*micronsToAngstroms,luminosityAbsorbed,0.0d0,13.8d0)
     call Assert('no emission beyond the tabulated wavelengths',luminosity(1:2),[0.0d0,0.0d0],absTol=0.0d0)
     !![
     <objectDestructor name="emission_"/>
     !!]
     call Unit_Tests_End_Group()
  end do
  call Unit_Tests_Finish()

contains

  function emissionFunction(parametersModel) result(emission)
    !!{RST
    The emission, :math:`\lambda L_\lambda` per unit luminosity, at each tabulated wavelength for the given parameters
    (:math:`N_0`, :math:`Y`, :math:`q`, :math:`\sigma`, :math:`\tau_\mathrm{V}`).
    !!}
    implicit none
    double precision, intent(in   ), dimension(5              ) :: parametersModel
    double precision               , dimension(countWavelength) :: emission

    emission=+emissionBase                                &
         &   +emissionCountClouds*    parametersModel(1)  &
         &   +emissionRatioRadii *log(parametersModel(2)) &
         &   +emissionExponent   *    parametersModel(3)  &
         &   +emissionWidth      *    parametersModel(4)  &
         &   +emissionDepth      *log(parametersModel(5))
    return
  end function emissionFunction

  double precision function trapezoid(x,y) result(integral)
    !!{RST
    Integrate ``y`` over ``x`` by the trapezoidal rule.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y

    integral=sum(0.5d0*(y(2:size(y))+y(1:size(y)-1))*(x(2:size(x))-x(1:size(x)-1)))
    return
  end function trapezoid

end program Test_Dust_Emission_AGN_Torus_CLUMPY
