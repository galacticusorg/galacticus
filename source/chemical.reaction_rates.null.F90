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

!% Contains a module which implements a null calculation of chemical reaction rates.

module Chemical_Reaction_Rates_Null
  !% Implements a null calculation of chemical reaction rates.
  implicit none
  private
  public :: Chemical_Reaction_Rates_Null_Initialize, Chemical_Reaction_Rates_Null_Compute

contains

  !# <chemicalReactionRates>
  !#  <unitName>Chemical_Reaction_Rates_Null_Initialize</unitName>
  !# </chemicalReactionRates>
  subroutine Chemical_Reaction_Rates_Null_Initialize(chemicalReactionRatesMethods)
    !% Initializes the null chemical reaction network module.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(in   ) :: chemicalReactionRatesMethods(:)

    return
  end subroutine Chemical_Reaction_Rates_Null_Initialize

  !# <chemicalRatesCompute>
  !#  <unitName>Chemical_Reaction_Rates_Null_Compute</unitName>
  !# </chemicalRatesCompute>
  subroutine Chemical_Reaction_Rates_Null_Compute(temperature,chemicalDensity,radiation,chemicalRates)
    !% Compute rates of change of chemical abundances due to reactions involving chemical hydrogen species.
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    type            (chemicalAbundances), intent(in   ) :: chemicalDensity
    double precision                    , intent(in   ) :: temperature
    type            (radiationStructure), intent(in   ) :: radiation
    type            (chemicalAbundances), intent(inout) :: chemicalRates

    return
  end subroutine Chemical_Reaction_Rates_Null_Compute

end module Chemical_Reaction_Rates_Null
