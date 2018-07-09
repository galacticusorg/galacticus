!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implements a cooling function class which sums over other cooling functions.

  !# <coolingFunction name="coolingFunctionSummation" defaultThreadPrivate="yes">
  !#  <description>Class providing a cooling function which sums over other cooling functions.</description>
  !# </coolingFunction>

  type, public :: coolantList
     class(coolingFunctionClass), pointer :: coolingFunction
     type (coolantList         ), pointer :: next            => null()
  end type coolantList

  type, extends(coolingFunctionClass) :: coolingFunctionSummation
     !% A cooling function class which sums over other cooling functions.
     private
     type(coolantList), pointer :: coolants
   contains
     final     ::                                       summationDestructor
     procedure :: coolingFunction                    => summationCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => summationCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => summationCoolingFunctionDensityLogSlope
     procedure :: descriptor                         => summationDescriptor
     procedure :: deepCopy                           => summationDeepCopy
  end type coolingFunctionSummation

  interface coolingFunctionSummation
     !% Constructors for the ``summation'' cooling function class.
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface coolingFunctionSummation

contains

  function summationConstructorParameters(parameters)
    !% Constructor for the ``summation'' cooling function class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (coolingFunctionSummation)                :: summationConstructorParameters
    type   (inputParameters         ), intent(inout) :: parameters
    type   (coolantList             ), pointer       :: coolant
    integer                                          :: i

    !$omp critical(coolingFunctionSummationInitialize)
    coolant => null()
    do i=1,parameters%copiesCount('coolingFunctionMethod',zeroIfNotPresent=.true.)
       if (associated(coolant)) then
          allocate(coolant                       %next    )
          coolant => coolant                       %next
       else
          allocate(summationConstructorParameters%coolants)
          coolant => summationConstructorParameters%coolants
       end if
       coolant%coolingFunction => coolingFunction(parameters,i)
    end do
    !$omp end critical(coolingFunctionSummationInitialize)
    return
  end function summationConstructorParameters
  
  function summationConstructorInternal(coolants)
    !% Internal constructor for the ``summation'' cooling function class.
    implicit none
    type(coolingFunctionSummation)                        :: summationConstructorInternal
    type(coolantList             ), target, intent(in   ) :: coolants

    summationConstructorInternal%coolants => coolants
    return
  end function summationConstructorInternal
  
  subroutine summationDestructor(self)
    !% Destructor for the ``summation'' cooling function class.
    implicit none
    type(coolingFunctionSummation), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ! Nothing to do.
    return
  end subroutine summationDestructor

  double precision function summationCoolingFunction(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function due to Compton scattering off of \gls{cmb} photons.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    double precision                          , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    type            (radiationStructure      ), intent(in   ) :: radiation
    type            (coolantList             ), pointer       :: coolant
    
    summationCoolingFunction =  0.0d0
    coolant                  => self%coolants
    do while (associated(coolant))
       summationCoolingFunction=                                                                &
            &                   +summationCoolingFunction                                       &
            &                   +coolant%coolingFunction%coolingFunction(                       &
            &                                                            numberDensityHydrogen, &
            &                                                            temperature          , &
            &                                                            gasAbundances        , &
            &                                                            chemicalDensities    , &
            &                                                            radiation              &
            &                                                           )
       coolant => coolant%next
    end do
    return
  end function summationCoolingFunction

  double precision function summationCoolingFunctionDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic gradient with respect to density of the cooling function due to Compton scattering off of \gls{cmb}
    !% photons.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Galacticus_Error
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    double precision                          , intent(in   ) :: numberDensityHydrogen  , temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    type            (radiationStructure      ), intent(in   ) :: radiation
    type            (coolantList             ), pointer       :: coolant
    double precision                                          :: coolingFunction        , coolingFunctionCumulative, &
         &                                                       coolingFunctionGradient

    coolingFunctionCumulative =  0.0d0
    coolingFunctionGradient   =  0.0d0
    coolant                   => self%coolants
    do while (associated(coolant))
       coolingFunction          =+coolant%coolingFunction%coolingFunction               (                       &
            &                                                                            numberDensityHydrogen, &
            &                                                                            temperature          , &
            &                                                                            gasAbundances        , &
            &                                                                            chemicalDensities    , &
            &                                                                            radiation              &
            &                                                                           )
       coolingFunctionCumulative=+coolingFunctionCumulative                                                     &
            &                    +coolingFunction
       coolingFunctionGradient  =+coolingFunctionGradient                                                       &
            &                    +coolingFunction                                                               &
            &                    /numberDensityHydrogen                                                         &
            &                    *coolant%coolingFunction%coolingFunctionDensityLogSlope(                       &
            &                                                                            numberDensityHydrogen, &
            &                                                                            temperature          , &
            &                                                                            gasAbundances        , &
            &                                                                            chemicalDensities    , &
            &                                                                            radiation              &
            &                                                                           )
       coolant => coolant%next
    end do
    if (coolingFunctionCumulative > 0.0d0) then
       summationCoolingFunctionDensityLogSlope=+coolingFunctionGradient   &
            &                                  *numberDensityHydrogen     &
            &                                  /coolingFunctionCumulative
    else
       summationCoolingFunctionDensityLogSlope=0.0d0
       if (coolingFunctionGradient /= 0.0d0) call Galacticus_Error_Report('cooling function is zero but has non-zero gradient with density - logarithmic slope is undefined'//{introspection:location})
    end if
    return
  end function summationCoolingFunctionDensityLogSlope
  
  double precision function summationCoolingFunctionTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic gradient with respect to temperature of the cooling function due to Compton scattering off of
    !% \gls{cmb} photons.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Galacticus_Error
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    double precision                          , intent(in   ) :: numberDensityHydrogen                       , temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    type            (radiationStructure      ), intent(in   ) :: radiation
    type            (coolantList             ), pointer       :: coolant
    double precision                                          :: coolingFunction        , coolingFunctionCumulative, &
         &                                                       coolingFunctionGradient

    coolingFunctionCumulative =  0.0d0
    coolingFunctionGradient   =  0.0d0
    coolant                   => self%coolants
    do while (associated(coolant))
       coolingFunction          =+coolant%coolingFunction%coolingFunction                   (                       &
            &                                                                                numberDensityHydrogen, &
            &                                                                                temperature          , &
            &                                                                                gasAbundances        , &
            &                                                                                chemicalDensities    , &
            &                                                                                radiation              &
            &                                                                               )
       coolingFunctionCumulative=+coolingFunctionCumulative                                                         &
            &                    +coolingFunction
       coolingFunctionGradient  =+coolingFunctionGradient                                                           &
            &                    +coolingFunction                                                                   &
            &                    /temperature                                                                       &
            &                    *coolant%coolingFunction%coolingFunctionTemperatureLogSlope(                       &
            &                                                                                numberDensityHydrogen, &
            &                                                                                temperature          , &
            &                                                                                gasAbundances        , &
            &                                                                                chemicalDensities    , &
            &                                                                                radiation              &
            &                                                                               )
       coolant => coolant%next
    end do
    if (coolingFunctionCumulative > 0.0d0) then
       summationCoolingFunctionTemperatureLogSlope=+coolingFunctionGradient   &
            &                                      *temperature               &
            &                                      /coolingFunctionCumulative
    else
       summationCoolingFunctionTemperatureLogSlope=0.0d0
       if (coolingFunctionGradient /= 0.0d0) call Galacticus_Error_Report('cooling function is zero but has non-zero gradient with temperature - logarithmic slope is undefined'//{introspection:location})
    end if
    return
  end function summationCoolingFunctionTemperatureLogSlope

  subroutine summationDescriptor(self,descriptor,includeMethod)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters
    use FoX_DOM
    implicit none
    class  (coolingFunctionSummation), intent(inout)           :: self
    type   (inputParameters         ), intent(inout)           :: descriptor
    logical                          , intent(in   ), optional :: includeMethod
    type   (coolantList             ), pointer                 :: coolant
    type   (inputParameters         )                          :: subParameters

    if (.not.present(includeMethod).or.includeMethod) call descriptor%addParameter("coolingFunctionMethod","summation")
    subParameters=descriptor%subparameters("coolingFunctionMethod")
    coolant       => self%coolants
    do while (associated(coolant))
       call coolant%coolingFunction%descriptor(subParameters)
       coolant => coolant%next
    end do
    return
  end subroutine summationDescriptor

  subroutine summationDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily summation} cooling function class.
    use Galacticus_Error
    implicit none
    class(coolingFunctionSummation), intent(inout) :: self
    class(coolingFunctionClass    ), intent(  out) :: destination
    type (coolantList             ), pointer       :: coolant_   , coolantDestination_, &
         &                                            coolantNew_

    call self%coolingFunctionClass%deepCopy(destination)
    select type (destination)
    type is (coolingFunctionSummation)
       destination%coolants => null          ()
       coolantDestination_  => null          ()
       coolant_             => self%coolants
       do while (associated(coolant_))
          allocate(coolantNew_)
          if (associated(coolantDestination_)) then
             coolantDestination_%next       => coolantNew_
             coolantDestination_            => coolantNew_             
          else
             destination          %coolants => coolantNew_
             coolantDestination_            => coolantNew_
          end if
          allocate(coolantNew_%coolingFunction,mold=coolant_%coolingFunction)
          call coolant_%coolingFunction%deepCopy(coolantNew_%coolingFunction)
          coolant_ => coolant_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine summationDeepCopy
