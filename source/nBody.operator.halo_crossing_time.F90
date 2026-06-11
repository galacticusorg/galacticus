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

!!{
Implements an N-body data operator which computes the crossing time of halos.
!!}

  !![
  <nbodyOperator name="nbodyOperatorHaloCrossingTime">
   <description>An N-body data operator which computes and stores the crossing time of halos defined as $t_\mathrm{cross} = 2 r_\mathrm{vir}/V_\mathrm{vir}$.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorHaloCrossingTime
     !!{
     An N-body data operator which computes the crossing time of halos.
     !!}
     private
   contains
     procedure :: operate => haloCrossingTimeOperate
  end type nbodyOperatorHaloCrossingTime

  interface nbodyOperatorHaloCrossingTime
     !!{
     Constructors for the \refClass{nbodyOperatorHaloCrossingTime} N-body operator class.
     !!}
     module procedure haloCrossingTimeConstructorParameters
  end interface nbodyOperatorHaloCrossingTime

contains

  function haloCrossingTimeConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorHaloCrossingTime} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorHaloCrossingTime)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=nbodyOperatorHaloCrossingTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloCrossingTimeConstructorParameters

  subroutine haloCrossingTimeOperate(self,simulations)
    !!{
    Compute the crossing times of halos.
    !!}
    use :: Display                         , only : displayIndent                 , displayUnindent  , verbosityLevelStandard
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, MpcPerKmPerSToGyr
    implicit none
    class           (nbodyOperatorHaloCrossingTime), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(:  ) :: simulations
    double precision                               , pointer      , dimension(:  ) :: radiusVirial, massVirial, &
         &                                                                            timeCrossing
    integer         (c_size_t                     )                                :: iSimulation

    call displayIndent('compute halo crossing times',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Retrieve required properties.
       radiusVirial => simulations(iSimulation)%propertiesReal%value('radiusVirial')
       massVirial   => simulations(iSimulation)%propertiesReal%value('massVirial'  )
       ! Allocate workspace.
       allocate(timeCrossing(size(massVirial)))
       ! Compute crossing times.
       !$omp workshare
       timeCrossing=+2.0d0                                   &
            &       *sqrt(                                   &
            &             +radiusVirial                  **3 &
            &             /gravitationalConstant_internal    &
            &             /massVirial                        &
            &            )                                   &
            &       *MpcPerKmPerSToGyr
       !$omp end workshare
       ! Store results.
       call simulations(iSimulation)%propertiesReal%set('timeCrossing',timeCrossing)
       nullify(timeCrossing)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine haloCrossingTimeOperate
