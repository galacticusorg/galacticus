!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements calculations of chemical reaction rates.

module Chemical_Reaction_Rates
  !% Implements calculations of chemical reaction rates.
  use Chemical_Abundances_Structure
  use ISO_Varying_String
  implicit none
  private
  public :: Chemical_Reaction_Rate

  ! Flag to indicate if this module has been initialized.
  logical                                            :: chemicalReactionRateInitialized=.false.

  ! Name of chemical reaction rates methods used.
  type   (varying_string), allocatable, dimension(:) :: chemicalReactionRateMethods

contains

  subroutine Chemical_Reaction_Rates_Initialize
    !% Initialize the chemical reaction rates module.
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="chemicalReactionRates" type="moduleUse">
    include 'chemical.reaction_rates.modules.inc'
    !# </include>
    implicit none
    integer :: chemicalReactionRatesCount

    ! Initialize if necessary.
    if (.not.chemicalReactionRateInitialized) then
       !$omp critical(Chemical_Reaction_Rates_Initialization)
       if (.not.chemicalReactionRateInitialized) then
          ! Get the chemical reaction rates method parameter.
          !@ <inputParameter>
          !@   <name>chemicalReactionRateMethods</name>
          !@   <defaultValue>hydrogenNetwork</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The names of the methods to be used for computing chemical reaction rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1..*</cardinality>
          !@ </inputParameter>
          chemicalReactionRatesCount=max(1,Get_Input_Parameter_Array_Size('chemicalReactionRatesMethods'))
          allocate(chemicalReactionRateMethods(chemicalReactionRatesCount))
          call Memory_Usage_Record(sizeof(chemicalReactionRateMethods))
          call Get_Input_Parameter('chemicalReactionRateMethods',chemicalReactionRateMethods,defaultValue=['null'])

          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="chemicalReactionRates" type="functionCall" functionType="void">
          !#  <functionArgs>chemicalReactionRateMethods</functionArgs>
          include 'chemical.reaction_rates.inc'
          !# </include>
          chemicalReactionRateInitialized=.true.
       end if
       !$omp end critical(Chemical_Reaction_Rates_Initialization)
    end if
    return
  end subroutine Chemical_Reaction_Rates_Initialize

  subroutine Chemical_Reaction_Rate(chemicalRates,temperature,chemicalDensity,radiation)
    !% Return chemical reaction rates at the given temperature for the specified set of chemical densities (in cm$^{-3}$) and radiation
    !% field. Units of the returned rates are cm$^-3$ s$^{-1}$.
    use Abundances_Structure
    use Radiation_Structure
    !# <include directive="chemicalRatesCompute" type="moduleUse">
    include 'chemical.reaction_rates.compute.modules.inc'
    !# </include>
    implicit none
    type            (chemicalAbundances), intent(inout) :: chemicalRates
    double precision                    , intent(in   ) :: temperature
    type            (chemicalAbundances), intent(in   ) :: chemicalDensity
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module.
    call Chemical_Reaction_Rates_Initialize

    call chemicalRates%reset()
    !# <include directive="chemicalRatesCompute" type="functionCall" functionType="void">
    !#  <functionArgs>temperature,chemicalDensity,radiation,chemicalRates</functionArgs>
    include 'chemical.reaction_rates.compute.inc'
    !# </include>

    return
  end subroutine Chemical_Reaction_Rate

end module Chemical_Reaction_Rates
