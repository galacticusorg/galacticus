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

  !!{RST
  Implements a node operator class used solely for testing.
  !!}

  !![
  <nodeOperator name="nodeOperatorTestInactivePropertyRead" docformat="rst">
   <description>
   A node operator, used solely for testing, that deliberately reads the stellar luminosities of the disk
   component during evaluation of the derivatives of {\normalfont \ttfamily active} properties. When the disk
   stellar luminosities are being solved as {\normalfont \ttfamily inactive} properties (and a Jacobian-based
   ODE solver is in use, so that {\normalfont \ttfamily inactive} properties are marked as such) this is an
   erroneous read of a stale property value, and (in debugging builds) is trapped by the run-time check
   introduced for \href{https://github.com/galacticusorg/galacticus/issues/128}{issue \#128}. Used by the
   {\normalfont \ttfamily test-inactive-property-read.py} test.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTestInactivePropertyRead
     !!{RST
     A node operator class, used solely for testing, that deliberately reads an inactive property value during
     active-property differential evolution.
     !!}
     private
   contains
     procedure :: differentialEvolution => testInactivePropertyReadDifferentialEvolution
  end type nodeOperatorTestInactivePropertyRead

  interface nodeOperatorTestInactivePropertyRead
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorTestInactivePropertyRead` node operator class.
     !!}
     module procedure testInactivePropertyReadConstructorParameters
  end interface nodeOperatorTestInactivePropertyRead

contains

  function testInactivePropertyReadConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorTestInactivePropertyRead` node operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorTestInactivePropertyRead)                :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=nodeOperatorTestInactivePropertyRead()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function testInactivePropertyReadConstructorParameters

  subroutine testInactivePropertyReadDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Deliberately read the (inactive) stellar luminosities of the disk during active-property differential
    evolution, to exercise the run-time inactive-property-read check.
    !!}
    use :: Galacticus_Nodes              , only : propertyActive     , nodeComponentDisk
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class    (nodeOperatorTestInactivePropertyRead), intent(inout), target  :: self
    type     (treeNode                            ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure(interruptTask                       ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    class    (nodeComponentDisk                   )               , pointer :: disk
    type     (stellarLuminosities                 )                         :: luminosities
    !$GLC attributes unused :: self, interrupt, functionInterrupt

    ! Act only during evaluation of the derivatives of active properties - this is precisely the phase during which
    ! reading an inactive property's (stale) value is an error.
    if (.not.propertyActive(propertyType)) return
    disk => node%disk()
    ! Deliberately read the stellar luminosities. If they are being solved as inactive properties this read is
    ! erroneous, and (in debugging builds) is trapped by the run-time check.
    luminosities=disk%luminositiesStellar()
    return
  end subroutine testInactivePropertyReadDifferentialEvolution
