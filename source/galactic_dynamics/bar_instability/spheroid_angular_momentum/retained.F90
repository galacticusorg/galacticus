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

!+    Contributions to this file made by: Andrew Robertson, Codex.

!!{RST
Implements the standard retained-angular-momentum prescription for spheroids formed by bar instabilities.
!!}

  !![
  <barInstabilitySpheroidAngularMomentum name="barInstabilitySpheroidAngularMomentumRetained" docformat="rst">
   <description>
   Assigns to the spheroid a fixed fraction of the angular momentum transfer rate associated with a bar instability. This reproduces the historical Galacticus prescription.
   </description>
  </barInstabilitySpheroidAngularMomentum>
  !!]
  type, extends(barInstabilitySpheroidAngularMomentumClass) :: barInstabilitySpheroidAngularMomentumRetained
     !!{RST
     Standard retained-angular-momentum prescription for a bar-built spheroid.
     !!}
   contains
     procedure :: rate => retainedRate
  end type barInstabilitySpheroidAngularMomentumRetained

  interface barInstabilitySpheroidAngularMomentumRetained
     module procedure retainedConstructorParameters
  end interface barInstabilitySpheroidAngularMomentumRetained

contains

  function retainedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the retained-angular-momentum prescription from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(barInstabilitySpheroidAngularMomentumRetained)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=barInstabilitySpheroidAngularMomentumRetained()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function retainedConstructorParameters

  double precision function retainedRate(self,node,rateMassGasTransfer,rateMassStellarTransfer,rateAngularMomentumDisk,rateAngularMomentumTransfer,fractionAngularMomentumRetainedSpheroid)
    !!{RST
    Return the spheroid pseudo-angular momentum rate.
    !!}
    implicit none
    class           (barInstabilitySpheroidAngularMomentumRetained), intent(inout)         :: self
    type            (treeNode                                     ), intent(inout), target :: node
    double precision                                               , intent(in   )         :: rateMassGasTransfer                    , &
         &                                                                                    rateMassStellarTransfer                , &
         &                                                                                    rateAngularMomentumDisk                , &
         &                                                                                    rateAngularMomentumTransfer            , &
         &                                                                                    fractionAngularMomentumRetainedSpheroid
    !$GLC attributes unused :: self, node, rateMassGasTransfer, rateMassStellarTransfer, rateAngularMomentumDisk

    retainedRate=fractionAngularMomentumRetainedSpheroid*rateAngularMomentumTransfer
    return
  end function retainedRate
