!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Implements a triaxiality modifier for power spectra in the halo model of clustering based on the results of \cite{smith_triaxial_2005}.

  !# <haloModelPowerSpectrumModifier name="haloModelPowerSpectrumModifierTriaxiality">
  !#  <description>A triaxiality modifier for power spectra in the halo model of clustering based on the results of \cite{smith_triaxial_2005}.</description>
  !# </haloModelPowerSpectrumModifier>

  use Tables

  type, extends(haloModelPowerSpectrumModifierClass) :: haloModelPowerSpectrumModifierTriaxiality
   contains
     procedure :: modify => triaxialityModify
  end type haloModelPowerSpectrumModifierTriaxiality
  
  interface haloModelPowerSpectrumModifierTriaxiality
     !% Constructor for the triaxiality halo model power spectra modifier class.
     module procedure triaxialityConstructor
  end interface haloModelPowerSpectrumModifierTriaxiality

  ! Tabulated results read from figures in Smith et al. (2005).
  double precision, parameter                  :: triaxialityWavenumberMinimum=1.0d-2
  double precision, parameter                  :: triaxialityWavenumberMaximum=1.0d+2
  integer         , parameter                  :: triaxialityWavenumberCount  =20
  double precision, parameter, dimension(20  ) :: triaxialityTwoHalo      =                                             &
       &                                                   [                                                            &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +0.99994d0, &
       &                                                    +0.99953d0, +0.99909d0, +0.99862d0, +0.99634d0, +0.99150d0, &
       &                                                    +0.98190d0, +0.97345d0, +0.96753d0, +0.96269d0, +0.95766d0, &
       &                                                    +0.95275d0, +0.94701d0, +0.94095d0, +0.93382d0, +0.92410d0  &
       &                                                   ]
  double precision, parameter, dimension(20,4) :: triaxialityOneHalo      =                                             &
       &                                           reshape(                                                             &
       &                                                   [                                                            &
       ! 1e11 < Mhalo/Msun < 1e12
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.99920d0, +0.99744d0, &
       &                                                    +0.99372d0, +0.98678d0, +0.97315d0, +0.95042d0, +0.92762d0, &
       &                                                    +0.93043d0, +0.92965d0, +0.92329d0, +0.92329d0, +0.92328d0, &
       ! 1e12 < Mhalo/Msun < 1e13
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.99632d0, +0.99203d0, &
       &                                                    +0.98111d0, +0.96274d0, +0.93784d0, +0.92888d0, +0.94435d0, &
       &                                                    +0.93548d0, +0.93460d0, +0.93642d0, +0.94121d0, +0.94834d0, &
       ! 1e13 < Mhalo/Msun < 1e14
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.98701d0, +0.97392d0, &
       &                                                    +0.95379d0, +0.94002d0, +0.95533d0, +0.95836d0, +0.96092d0, &
       &                                                    +0.97041d0, +0.98241d0, +0.99732d0, +1.01862d0, +1.04115d0, &
       ! 1e14 < Mhalo/Msun
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.97246d0, +0.96142d0, &
       &                                                    +0.97073d0, +0.99659d0, +0.99625d0, +1.01229d0, +1.03127d0, &
       &                                                    +1.05791d0, +1.08701d0, +1.11753d0, +1.14804d0, +1.17855d0  &
       &                                                   ]                                                          , &
       &                                                   [20,4]                                                       &
       &                                                  )
  double precision, parameter, dimension(   4) :: triaxialityMass         =[1.0d12,1.0d13,1.0d14,huge(1.0d0)]
  logical                                      :: triaxialityInitialized  =.false.
  type(table1DLogarithmicLinear)               :: triaxialityTable
  !$omp threadprivate(triaxialityInitialized,triaxialityTable)

contains

  function triaxialityConstructor()
    !% Default constructor for the triaxiality halo model power spectra modifier class.
    use Input_Parameters
    implicit none
    type(haloModelPowerSpectrumModifierTriaxiality) :: triaxialityConstructor

    return
  end function triaxialityConstructor
  
  subroutine triaxialityModify(self,wavenumber,term,powerSpectrum,powerSpectrumCovariance,mass)
    !% Applies a triaxiality modification to a halo model power spectrum based on the results of \cite{smith_triaxial_2005}.
    use Cosmology_Parameters
    use Vectors
    use Galacticus_Error
    implicit none
    class           (haloModelPowerSpectrumModifierTriaxiality), intent(inout)                           :: self
    double precision                                           , intent(in   ), dimension(:  )           :: wavenumber
    integer                                                    , intent(in   )                           :: term
    double precision                                           , intent(inout), dimension(:  )           :: powerSpectrum
    double precision                                           , intent(inout), dimension(:,:), optional :: powerSpectrumCovariance
    double precision                                           , intent(in   )                , optional :: mass
    class           (cosmologyParametersClass                 ), pointer                                 :: cosmologyParameters_
    double precision                                           , parameter                               :: covarianceFraction     =0.4d0
    double precision                                           , allocatable  , dimension(:  )           :: covariance
    integer                                                                                              :: i                            , tableIndex
    double precision                                                                                     :: termPower

    ! Mass is required.
    if (.not.present(mass)) call Galacticus_Error_Report('triaxialityModify','mass is required')
    ! Initialize tables if necessary.
    if (.not.triaxialityInitialized) then
       call triaxialityTable%create(                              &
            &                       triaxialityWavenumberMinimum, &
            &                       triaxialityWavenumberMaximum, &
            &                       triaxialityWavenumberCount  , &
            &                       5                           , &
            &                       extrapolationTypeExtrapolate  &
            &                      )
       do i=1,4
          call triaxialityTable%populate(triaxialityOneHalo(:,i),i)
       end do
       call    triaxialityTable%populate(triaxialityTwoHalo     ,5)
       triaxialityInitialized=.true.
    end if
    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    ! Determine table to use.
    select case (term)
    case (termOneHalo)
       tableIndex=1
       do while (mass < triaxialityMass(tableIndex)/cosmologyParameters_%HubbleConstant(unitsLittleH))
          tableIndex=tableIndex+1
       end do
       termPower=1.0d0
    case (termTwoHalo)
       tableIndex=5
       termPower =0.5d0 ! Since the two-halo term is accumulated as the square-root of the power spectrum.
    end select
    ! Compute covariance if required.
    if (present(powerSpectrumCovariance)) then
       allocate(covariance(size(powerSpectrum)))
       do i=1,size(wavenumber)
          covariance(i)=+covarianceFraction                                                                 &
               &        *powerSpectrum(i)                                                                   &
               &        *(                                                                                  &
               &          +triaxialityTable%interpolate(                                                    &
               &                                        +wavenumber(i)                                      &
               &                                        /cosmologyParameters_%HubbleConstant(unitsLittleH), &
               &                                        table=tableIndex                                    &
               &                                       )**termPower                                         &
               &          -1.0d0                                                                            &
               &         )
       end do
       powerSpectrumCovariance=powerSpectrumCovariance*Vector_Outer_Product(covariance,covariance)
    end if
    ! Compute the modification.
    do i=1,size(wavenumber)
       powerSpectrum(i)=+powerSpectrum(i)                                                                 &
            &           *triaxialityTable%interpolate(                                                    &
            &                                         +wavenumber(i)                                      &
            &                                         /cosmologyParameters_%HubbleConstant(unitsLittleH), &
            &                                         table=tableIndex                                    &
            &                                        )**termPower
    end do
    return
  end subroutine triaxialityModify
