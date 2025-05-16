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

  !![
  <radiativeTransferConvergence name="radiativeTransferConvergenceLycEscape">
   <description>A task which performs radiative transfer.</description>
  </radiativeTransferConvergence>
  !!]
  type, extends(radiativeTransferConvergenceClass) :: radiativeTransferConvergenceLycEscape
     !!{
     Implementation of a radiative transfer convergence class based on the recombination rate of hydrogen.
     !!}
     private
     double precision :: toleranceRelative
     double precision :: escapeRateTotal  , escapeRateTotalPrevious
   contains
     procedure :: testConvergence     => lycEscapeTestConvergence
     procedure :: photonPacketEscapes => lycEscapePhotonPacketEscapes
  end type radiativeTransferConvergenceLycEscape
  
  interface radiativeTransferConvergenceLycEscape
     !!{
     Constructors for the \refClass{radiativeTransferConvergenceLycEscape} radiative transfer matter class.
     !!}
     module procedure lycEscapeConstructorParameters
     module procedure lycEscapeConstructorInternal
  end interface radiativeTransferConvergenceLycEscape
  
contains

  function lycEscapeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferConvergenceLycEscape} radiative transfer matter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (radiativeTransferConvergenceLycEscape)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: toleranceRelative
    
    !![
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The relative tolerance in hydrogen Lyc escape required to declare convergence.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=radiativeTransferConvergenceLycEscape(toleranceRelative)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function lycEscapeConstructorParameters

  function lycEscapeConstructorInternal(toleranceRelative) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferConvergenceLycEscape} radiative transfer matter class.
    !!}
    implicit none
    type            (radiativeTransferConvergenceLycEscape)                :: self
    double precision                                       , intent(in   ) :: toleranceRelative
    !![
    <constructorAssign variables="toleranceRelative"/>
    !!]

    self%escapeRateTotal        =      0.0d0
    self%escapeRateTotalPrevious=-huge(0.0d0)
    return
  end function lycEscapeConstructorInternal
  
  subroutine lycEscapeTestConvergence(self,radiativeTransferMatter_,properties,statusCell,converged)
    !!{
    Test convergence in the computational domain cell.
    !!}
    use :: Display                   , only : displayMessage               , verbosityLevelStandard
    use :: MPI_Utilities             , only : mpiSelf
    use :: Radiative_Transfer_Matters, only : radiativeTransferMatterAtomic, radiativeTransferPropertiesMatterAtomic
    implicit none
    class           (radiativeTransferConvergenceLycEscape), intent(inout) :: self
    class           (radiativeTransferMatterClass         ), intent(inout) :: radiativeTransferMatter_
    class           (radiativeTransferPropertiesMatter    ), intent(inout) :: properties
    type            (enumerationStatusCellType            ), intent(in   ) :: statusCell
    logical                                                , intent(  out) :: converged
    character       (len=128                              )                :: message
    double precision                                                       :: escapeRateTotal

    ! Test convergence for the last cell.
    if (statusCell == statusCellLast) then
       ! Sum escape rate over all processes.
       escapeRateTotal=mpiSelf%sum(self%escapeRateTotal)
       converged      =                              abs(escapeRateTotal-self%escapeRateTotalPrevious) &
            &          <                                                                               &
            &           self%toleranceRelative*0.5d0*   (escapeRateTotal+self%escapeRateTotalPrevious)
       if (mpiSelf%isMaster()) then
          write (message,'(a,e12.6,a)') 'hydrogen Lyc escape rate = ',escapeRateTotal,' s⁻¹ ('
          if (.not.converged) then
             message=trim(message)//'not converged)'
          else
             message=trim(message)//    'converged)'
          end if
          call displayMessage(trim(message),verbosityLevelStandard)
       end if
       ! Reset accumulated recombination rate.
       self%escapeRateTotalPrevious=escapeRateTotal
       self%escapeRateTotal        =0.0d0
    else
       converged=.true.
    end if
    return
  end subroutine lycEscapeTestConvergence

  subroutine lycEscapePhotonPacketEscapes(self,photonPacket)
    !!{
    Process an escaping photon packet.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    use :: Numerical_Constants_Atomic      , only : lymanSeriesLimitWavelengthHydrogen_atomic
    use :: Numerical_Constants_Physical    , only : plancksConstant                          , speedLight
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    implicit none
    class           (radiativeTransferConvergenceLycEscape), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass   ), intent(inout) :: photonPacket
    double precision                                                       :: energyPhoton

    if (photonPacket%wavelength() < lymanSeriesLimitWavelengthHydrogen_atomic) then
       energyPhoton        =+plancksConstant                           &
            &               *speedLight                                &
            &               *metersToAngstroms                         &
            &               /photonPacket%wavelength                ()
       self%escapeRateTotal=+self%escapeRateTotal                      &
            &               +photonPacket%luminosity                () &
            &               *luminositySolar                           &
            &               /energyPhoton       
    end if
    return
  end subroutine lycEscapePhotonPacketEscapes
