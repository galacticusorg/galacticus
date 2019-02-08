!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a cooling function class which implements cooling due to Compton scattering off of \gls{cmb} photons.

  use Chemical_States, only : chemicalStateClass, chemicalState

  !# <coolingFunction name="coolingFunctionCMBCompton">
  !#  <description>Class providing a cooling function due to Compton scattering off of \gls{cmb} photons.</description>
  !# </coolingFunction>
  type, extends(coolingFunctionClass) :: coolingFunctionCMBCompton
     !% A cooling function class which implements cooling due to Compton scattering off of \gls{cmb} photons.
     private
     class(chemicalStateClass), pointer :: chemicalState_ => null()
   contains
     final     ::                                       cmbComptonDestructor
     procedure :: coolingFunction                    => cmbComptonCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => cmbComptonCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => cmbComptonCoolingFunctionDensityLogSlope
  end type coolingFunctionCMBCompton

  interface coolingFunctionCMBCompton
     !% Constructors for the ``CMB Compton'' cooling function class.
     module procedure cmbComptonConstructorParameters
     module procedure cmbComptonConstructorInternal
  end interface coolingFunctionCMBCompton

contains

  function cmbComptonConstructorParameters(parameters) result(self)
    !% Constructor for the ``CMB Compton'' cooling function class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (coolingFunctionCMBCompton)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(chemicalStateClass       ), pointer       :: chemicalState_

    !# <objectBuilder class="chemicalState" name="chemicalState_" source="parameters"/>
    self=coolingFunctionCMBCompton(chemicalState_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="chemicalState_"/>
    return
  end function cmbComptonConstructorParameters

  function cmbComptonConstructorInternal(chemicalState_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily cmbCompton} cooling function class.
    use Input_Parameters
    implicit none
    type (coolingFunctionCMBCompton)                        :: self
    class(chemicalStateClass       ), intent(in   ), target :: chemicalState_
    !# <constructorAssign variables="*chemicalState_"/>

    return
  end function cmbComptonConstructorInternal

  subroutine cmbComptonDestructor(self)
    !% Destructor for the {\normalfont \ttfamily cmbCompton} cooling function class.
    implicit none
    type(coolingFunctionCMBCompton), intent(inout) :: self

    !# <objectDestructor name="self%chemicalState_"/>
    return
  end subroutine cmbComptonDestructor

  double precision function cmbComptonCoolingFunction(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function due to Compton scattering off of \gls{cmb} photons.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Fields
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    class           (coolingFunctionCMBCompton), intent(inout) :: self
    double precision                           , intent(in   ) :: numberDensityHydrogen                                    , temperature
    type            (abundances               ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances       ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass      ), intent(inout) :: radiation
    double precision                           , parameter     :: comptonRateNormalization            =+4.0d0                             &
         &                                                                                             *thomsonCrossSection               &
         &                                                                                             *radiationConstant                 &
         &                                                                                             *boltzmannsConstant                &
         &                                                                                             /electronMass                      &
         &                                                                                             /speedLight                        &
         &                                                                                             /ergs
    double precision                                           :: temperatureCosmicMicrowaveBackground
    !GCC$ attributes unused :: self, chemicalDensities

    ! Find the cosmic microwave background radiation field.
    temperatureCosmicMicrowaveBackground=cmbComptonTemperature(radiation)
    ! Compute the Compton cooling rate.
    cmbComptonCoolingFunction=+comptonRateNormalization                                     &
         &                    *  self%chemicalState_%electronDensity(                       &
         &                                                           numberDensityHydrogen, &
         &                                                           temperature          , &
         &                                                           gasAbundances        , &
         &                                                           radiation              &
         &                                                          )                       &
         &                    *  temperatureCosmicMicrowaveBackground**4                    &
         &                    *(                                                            &
         &                      +temperature                                                &
         &                      -temperatureCosmicMicrowaveBackground                       &
         &                     )
    return
  end function cmbComptonCoolingFunction

  double precision function cmbComptonCoolingFunctionDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic gradient with respect to density of the cooling function due to Compton scattering off of \gls{cmb}
    !% photons.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Fields
    implicit none
    class           (coolingFunctionCMBCompton), intent(inout) :: self
    double precision                           , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances               ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances       ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass      ), intent(inout) :: radiation
    !GCC$ attributes unused :: self, chemicalDensities
    
    ! Slope depends only on the behavior of electron density with density.
    cmbComptonCoolingFunctionDensityLogSlope=+self%chemicalState_%electronDensityDensityLogSlope(                       &
         &                                                                                       numberDensityHydrogen, &
         &                                                                                       temperature          , &
         &                                                                                       gasAbundances        , &
         &                                                                                       radiation              &
         &                                                                                      )
    return
  end function cmbComptonCoolingFunctionDensityLogSlope
  
  double precision function cmbComptonCoolingFunctionTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic gradient with respect to temperature of the cooling function due to Compton scattering off of
    !% \gls{cmb} photons.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Fields
    implicit none
    class           (coolingFunctionCMBCompton), intent(inout) :: self
    double precision                           , intent(in   ) :: numberDensityHydrogen               , temperature
    type            (abundances               ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances       ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass      ), intent(inout) :: radiation
    double precision                                           :: temperatureCosmicMicrowaveBackground
    !GCC$ attributes unused :: self, chemicalDensities
    
    ! Find the cosmic microwave background radiation field.
    temperatureCosmicMicrowaveBackground=cmbComptonTemperature(radiation)
    ! Compute the logarithmic slope.
    cmbComptonCoolingFunctionTemperatureLogSlope=                                                &
         & +  self     %chemicalState_%electronDensityTemperatureLogSlope(                       &
         &                                                                numberDensityHydrogen, &
         &                                                                temperature          , &
         &                                                                gasAbundances        , &
         &                                                                radiation              &
         &                                                               )                       &
         & +  temperature                                                                        &
         & /(                                                                                    &
         &   +temperature                                                                        &
         &   -temperatureCosmicMicrowaveBackground                                               &
         &  )
    return
  end function cmbComptonCoolingFunctionTemperatureLogSlope

  double precision function cmbComptonTemperature(radiation)
    !% Return the temperature of the cosmic microwave background.
    use Radiation_Fields
    use Galacticus_Error
    implicit none
    class(radiationFieldClass), intent(inout) :: radiation
    type (radiationFieldList ), pointer       :: radiationField_
    
    select type (radiation)
    class is (radiationFieldCosmicMicrowaveBackground)
       cmbComptonTemperature=radiation%temperature()
    class is (radiationFieldSummation                )
       cmbComptonTemperature =  0.0d0
       radiationField_       => radiation%list()
       do while (associated(radiationField_))
          select type (radiation_ => radiationField_%radiationField_)
          class is (radiationFieldCosmicMicrowaveBackground)
             cmbComptonTemperature=radiation_%temperature()
          end select
          radiationField_ => radiationField_%next
       end do
    class default
       cmbComptonTemperature=0.0d0
       call Galacticus_Error_Report('unable to find cosmic microwave background radiation field'//{introspection:location})
    end select
    return
  end function cmbComptonTemperature
