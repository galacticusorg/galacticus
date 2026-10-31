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

!+    Contributions to this file made by: Andrew Robertson, Codex, Claude.

!!{RST
Implements the standard angular-momentum-transfer prescription for spheroids formed by bar instabilities.
!!}

  !![
  <barInstabilitySpheroidAngularMomentum name="barInstabilitySpheroidAngularMomentumAngularMomentumTransfer" docformat="rst">
   <description>
   Assigns to the spheroid a fixed fraction of the angular momentum transfer rate associated with a bar instability. This reproduces the historical Galacticus prescription.
   </description>
  </barInstabilitySpheroidAngularMomentum>
  !!]
  type, extends(barInstabilitySpheroidAngularMomentumClass) :: barInstabilitySpheroidAngularMomentumAngularMomentumTransfer
     !!{RST
     Standard angular-momentum-transfer prescription for a bar-built spheroid.
     !!}
   contains
     procedure :: rate => angularMomentumTransferRate
  end type barInstabilitySpheroidAngularMomentumAngularMomentumTransfer

  interface barInstabilitySpheroidAngularMomentumAngularMomentumTransfer
     module procedure angularMomentumTransferConstructorParameters
  end interface barInstabilitySpheroidAngularMomentumAngularMomentumTransfer

contains

  function angularMomentumTransferConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the angular-momentum-transfer prescription from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(barInstabilitySpheroidAngularMomentumAngularMomentumTransfer)                :: self
    type(inputParameters                                             ), intent(inout) :: parameters

    self=barInstabilitySpheroidAngularMomentumAngularMomentumTransfer()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function angularMomentumTransferConstructorParameters

  double precision function angularMomentumTransferRate(self,node,rateMassGasTransfer,rateMassStellarTransfer,rateAngularMomentumDisk,rateAngularMomentumTransfer,fractionAngularMomentumRetainedSpheroid)
    !!{RST
    Return the spheroid pseudo-angular momentum rate.
    !!}
    implicit none
    class           (barInstabilitySpheroidAngularMomentumAngularMomentumTransfer), intent(inout)         :: self
    type            (treeNode                                                    ), intent(inout), target :: node
    double precision                                                              , intent(in   )         :: rateMassGasTransfer                    , &
         &                                                                                                   rateMassStellarTransfer                , &
         &                                                                                                   rateAngularMomentumDisk                , &
         &                                                                                                   rateAngularMomentumTransfer            , &
         &                                                                                                   fractionAngularMomentumRetainedSpheroid
    !$GLC attributes unused :: self, node, rateMassGasTransfer, rateMassStellarTransfer, rateAngularMomentumDisk

    angularMomentumTransferRate=fractionAngularMomentumRetainedSpheroid*rateAngularMomentumTransfer
    return
  end function angularMomentumTransferRate
