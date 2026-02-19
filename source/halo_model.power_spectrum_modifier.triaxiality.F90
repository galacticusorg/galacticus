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
  Implements a triaxiality modifier for power spectra in the halo model of clustering based on the results of \cite{smith_triaxial_2005}.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Tables              , only : table1DLogarithmicLinear

  !![
  <haloModelPowerSpectrumModifier name="haloModelPowerSpectrumModifierTriaxiality">
   <description>
    A halo model power spectrum modifier class which attempts to modify power spectra to approximately account for the effects
    of halo triaxiality using the results of \cite{smith_triaxial_2005}. Specifically, the one- and two-halo power spectra are
    multiplied by a correction factor, $\Delta^2_\mathrm{triax}/\Delta^2_\mathrm{sphere}$, derived from the lower panels of
    Figures 3 and 2 of \cite{smith_triaxial_2005} respectively for their ``JS02'' profile model. Given the uncertainty in this
    correction, the power spectrum covariance (if provided) is incremented by $\epsilon^2 (
    [\Delta^2_\mathrm{triax}/\Delta^2_\mathrm{sphere}-1] \otimes [\Delta^2_\mathrm{triax}/\Delta^2_\mathrm{sphere}-1)$ where
    $\epsilon=0.4$ is chosen to approximate the difference between ``continuity'' and ``JS02'' profiles in
    \cite{smith_triaxial_2005}.
   </description>
  </haloModelPowerSpectrumModifier>
  !!]
  type, extends(haloModelPowerSpectrumModifierClass) :: haloModelPowerSpectrumModifierTriaxiality
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     type (table1DLogarithmicLinear)          :: triaxialityTable
   contains
     final     ::           triaxialityDestructor
     procedure :: modify => triaxialityModify
  end type haloModelPowerSpectrumModifierTriaxiality

  interface haloModelPowerSpectrumModifierTriaxiality
     !!{
     Constructors for the \refClass{haloModelPowerSpectrumModifierTriaxiality} halo model power spectrum modifier class.
     !!}
     module procedure triaxialityConstructorParameters
     module procedure triaxialityConstructorInternal
  end interface haloModelPowerSpectrumModifierTriaxiality

  ! Tabulated results read from figures in Smith et al. (2005).
  double precision, parameter                  :: wavenumberMinimum=1.0d-2
  double precision, parameter                  :: wavenumberMaximum=1.0d+2
  integer         , parameter                  :: countWavenumbers =20
  double precision, parameter, dimension(20  ) :: factorTwoHalo    =                                                    &
       &                                                   [                                                            &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +0.99994d0, &
       &                                                    +0.99953d0, +0.99909d0, +0.99862d0, +0.99634d0, +0.99150d0, &
       &                                                    +0.98190d0, +0.97345d0, +0.96753d0, +0.96269d0, +0.95766d0, &
       &                                                    +0.95275d0, +0.94701d0, +0.94095d0, +0.93382d0, +0.92410d0  &
       &                                                   ]
  double precision, parameter, dimension(20,4) :: factorOneHalo    =                                                    &
       &                                           reshape(                                                             &
       &                                                   [                                                            &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, & ! 1e11 < Mhalo/Msun < 1e12
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.99920d0, +0.99744d0, &
       &                                                    +0.99372d0, +0.98678d0, +0.97315d0, +0.95042d0, +0.92762d0, &
       &                                                    +0.93043d0, +0.92965d0, +0.92329d0, +0.92329d0, +0.92328d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, & ! 1e12 < Mhalo/Msun < 1e13
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.99632d0, +0.99203d0, &
       &                                                    +0.98111d0, +0.96274d0, +0.93784d0, +0.92888d0, +0.94435d0, &
       &                                                    +0.93548d0, +0.93460d0, +0.93642d0, +0.94121d0, +0.94834d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, & ! 1e13 < Mhalo/Msun < 1e14
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.98701d0, +0.97392d0, &
       &                                                    +0.95379d0, +0.94002d0, +0.95533d0, +0.95836d0, +0.96092d0, &
       &                                                    +0.97041d0, +0.98241d0, +0.99732d0, +1.01862d0, +1.04115d0, &
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, +1.00000d0, & ! 1e14 < Mhalo/Msun
       &                                                    +1.00000d0, +1.00000d0, +1.00000d0, +0.97246d0, +0.96142d0, &
       &                                                    +0.97073d0, +0.99659d0, +0.99625d0, +1.01229d0, +1.03127d0, &
       &                                                    +1.05791d0, +1.08701d0, +1.11753d0, +1.14804d0, +1.17855d0  &
       &                                                   ]                                                          , &
       &                                                   [20,4]                                                       &
       &                                                  )
  double precision, parameter, dimension(   4) :: massTable         =[1.0d12,1.0d13,1.0d14,huge(1.0d0)]

contains

  function triaxialityConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily triaxiality} hot halo outflow reincorporation class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloModelPowerSpectrumModifierTriaxiality)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyParametersClass                 ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=haloModelPowerSpectrumModifierTriaxiality(cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function triaxialityConstructorParameters

  function triaxialityConstructorInternal(cosmologyParameters_) result(self)
    !!{
    Default constructor for the triaxiality hot halo outflow reincorporation class.
    !!}
    use :: Table_Labels, only : extrapolationTypeExtrapolate
    implicit none
    type   (haloModelPowerSpectrumModifierTriaxiality)                        :: self
    class  (cosmologyParametersClass                 ), intent(in   ), target :: cosmologyParameters_
    integer                                                                   :: i
    !![
    <constructorAssign variables="*cosmologyParameters_"/>
    !!]

    call self%triaxialityTable%create(                                          &
         &                            wavenumberMinimum                       , &
         &                            wavenumberMaximum                       , &
         &                            countWavenumbers                        , &
         &                            5                                       , &
         &                            spread(extrapolationTypeExtrapolate,1,2)  &
         &                           )
    do i=1,4
       call self%triaxialityTable%populate(factorOneHalo(:,i),i)
    end do
    call    self%triaxialityTable%populate(factorTwoHalo     ,5)
    return
  end function triaxialityConstructorInternal

  subroutine triaxialityDestructor(self)
    !!{
    Destructor for the \refClass{haloModelPowerSpectrumModifierTriaxiality} halo model power spectrum modifier class.
    !!}
    implicit none
    type(haloModelPowerSpectrumModifierTriaxiality), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine triaxialityDestructor

  subroutine triaxialityModify(self,wavenumber,term,powerSpectrum,powerSpectrumCovariance,mass)
    !!{
    Applies a triaxiality modification to a halo model power spectrum based on the results of \cite{smith_triaxial_2005}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    use :: Vectors             , only : Vector_Outer_Product
    implicit none
    class           (haloModelPowerSpectrumModifierTriaxiality), intent(inout)                           :: self
    double precision                                           , intent(in   ), dimension(:  )           :: wavenumber
    type            (enumerationHaloModelTermType             ), intent(in   )                           :: term
    double precision                                           , intent(inout), dimension(:  )           :: powerSpectrum
    double precision                                           , intent(inout), dimension(:,:), optional :: powerSpectrumCovariance
    double precision                                           , intent(in   )                , optional :: mass
    double precision                                           , parameter                               :: covarianceFraction     =0.4d0
    double precision                                           , allocatable  , dimension(:  )           :: covariance
    integer                                                                                              :: i                            , tableIndex

    ! Mass is required.
    if (.not.present(mass)) call Error_Report('mass is required'//{introspection:location})
    ! Determine table to use.
    select case (term%ID)
    case (haloModelTermOneHalo%ID)
       tableIndex=1
       do while (mass*self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH) < massTable(tableIndex))
          tableIndex=tableIndex+1
       end do
    case (haloModelTermTwoHalo%ID)
       tableIndex=5
    end select
    ! Compute covariance if required.
    if (present(powerSpectrumCovariance)) then
       allocate(covariance(size(powerSpectrum)))
       select case (term%ID)
       case (haloModelTermOneHalo%ID)
          do i=1,size(wavenumber)
             covariance(i)=+covarianceFraction                                                                                      &
                  &        *powerSpectrum(i)                                                                                        &
                  &        *(                                                                                                       &
                  &          +     self%triaxialityTable%interpolate(                                                               &
                  &                                                  +wavenumber(i)                                                 &
                  &                                                  /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH), &
                  &                                                  table=tableIndex                                               &
                  &                                                 )                                                               &
                  &          -1.0d0                                                                                                 &
                  &         )
          end do
       case (haloModelTermTwoHalo%ID)
          do i=1,size(wavenumber)
             covariance(i)=+covarianceFraction                                                                                      &
                  &        *powerSpectrum(i)                                                                                        &
                  &        *(                                                                                                       &
                  &          +sqrt(                                                                                                 &
                  &                self%triaxialityTable%interpolate(                                                               &
                  &                                                  +wavenumber(i)                                                 &
                  &                                                  /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH), &
                  &                                                  table=tableIndex                                               &
                  &                                                 )                                                               &
                  &               )                                                                                                 &
                  &          -1.0d0                                                                                                 &
                  &         )
          end do
       end select
       powerSpectrumCovariance=powerSpectrumCovariance*Vector_Outer_Product(covariance,covariance)
    end if
    ! Compute the modification.
    do i=1,size(wavenumber)
       select case (term%ID)
       case (haloModelTermOneHalo%ID)
          powerSpectrum(i)=+powerSpectrum(i)                                                                                      &
               &           *     self%triaxialityTable%interpolate(                                                               &
               &                                                   +wavenumber(i)                                                 &
               &                                                   /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH), &
               &                                                   table=tableIndex                                               &
               &                                                  )
       case (haloModelTermTwoHalo%ID)
          powerSpectrum(i)=+powerSpectrum(i)                                                                                      &
               &           *sqrt(                                                                                                 &
               &                 self%triaxialityTable%interpolate(                                                               &
               &                                                   +wavenumber(i)                                                 &
               &                                                   /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH), &
               &                                                   table=tableIndex                                               &
               &                                                  )                                                               &
               &                )
       end select
    end do
    return
  end subroutine triaxialityModify
