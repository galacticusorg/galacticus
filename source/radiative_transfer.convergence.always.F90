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
  <radiativeTransferConvergence name="radiativeTransferConvergenceAlways">
   <description>A convergence criterion for radiative transfer which always passes.</description>
  </radiativeTransferConvergence>
  !!]
  type, extends(radiativeTransferConvergenceClass) :: radiativeTransferConvergenceAlways
     !!{
     A convergence criterion for radiative transfer which always passes.
     !!}
     private
   contains
     procedure :: testConvergence     => alwaysTestConvergence
     procedure :: photonPacketEscapes => alwaysPhotonPacketEscapes
  end type radiativeTransferConvergenceAlways
  
  interface radiativeTransferConvergenceAlways
     !!{
     Constructors for the \refClass{radiativeTransferConvergenceAlways} radiative transfer matter class.
     !!}
     module procedure alwaysConstructorParameters
  end interface radiativeTransferConvergenceAlways
  
contains  

  function alwaysConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferConvergenceAlways} radiative transfer matter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(radiativeTransferConvergenceAlways)                :: self
    type(inputParameters                   ), intent(inout) :: parameters
    
    self=radiativeTransferConvergenceAlways()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function alwaysConstructorParameters

  subroutine alwaysTestConvergence(self,radiativeTransferMatter_,properties,statusCell,converged)
    !!{
    Test convergence in the computational domain cell.
    !!}
    implicit none
    class  (radiativeTransferConvergenceAlways), intent(inout) :: self
    class  (radiativeTransferMatterClass      ), intent(inout) :: radiativeTransferMatter_
    class  (radiativeTransferPropertiesMatter ), intent(inout) :: properties
    type   (enumerationStatusCellType         ), intent(in   ) :: statusCell
    logical                                    , intent(  out) :: converged
    !$GLC attributes unused :: self, radiativeTransferMatter_, properties, statusCell
    
    converged=.true.
    return
  end subroutine alwaysTestConvergence

  subroutine alwaysPhotonPacketEscapes(self,photonPacket)
    !!{
    Process an escaping photon packet.
    !!}
    implicit none
    class(radiativeTransferConvergenceAlways), intent(inout) :: self
    class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    !$GLC attributes unused :: self, photonPacket

    return
  end subroutine alwaysPhotonPacketEscapes
