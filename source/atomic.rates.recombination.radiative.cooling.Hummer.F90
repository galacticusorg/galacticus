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
  An implementation of atomic recombination cooling rates based on the tabulated results of \cite{hummer_total_1994} and \cite{hummer_recombination_1998}.
  !!}

  use :: Atomic_Rates_Recombination_Radiative, only : atomicRecombinationRateRadiativeClass
  use :: Tables                              , only : table1DLogarithmicLinear
  
  !![
  <atomicRecombinationRateRadiativeCooling name="atomicRecombinationRateRadiativeCoolingHummer">
   <description>Atomic radiative cooling rates based on the tabulated results of \cite{hummer_total_1994} and \cite{hummer_recombination_1998}. For non-hydrogenic or helium-like ions an optional default of $\beta = \gamma \alpha$ can be returned where $\alpha$ is the corresponding radiative recombination coefficient and $\gamma$ is a parameter.</description>
  </atomicRecombinationRateRadiativeCooling>
  !!]
  type, extends(atomicRecombinationRateRadiativeCoolingClass) :: atomicRecombinationRateRadiativeCoolingHummer
     !!{
     A recombination cooling rate class assuming a thermal electron distribution.
     !!}
     private
     class           (atomicRecombinationRateRadiativeClass), pointer      :: atomicRecombinationRateRadiative_ => null()
     double precision                                                      :: gamma
     type            (table1DLogarithmicLinear             ), dimension(2) :: coefficientTable
   contains
     final     ::         hummerDestructor
     procedure :: rate => hummerRate
  end type atomicRecombinationRateRadiativeCoolingHummer

  interface atomicRecombinationRateRadiativeCoolingHummer
     !!{
     Constructors for the \refClass{atomicRecombinationRateRadiativeCoolingHummer} atomic radiative recombination class.
     !!}
     module procedure hummerConstructorParameters
     module procedure hummerConstructorInternal
  end interface atomicRecombinationRateRadiativeCoolingHummer

  !![
  <enumeration>
   <name>sequence</name>
   <description>Enumeration of isoelectronic sequences.</description>
   <indexing>0</indexing>
   <entry label="none"    />
   <entry label="hydrogen"/>
   <entry label="helium"  />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>level</name>
   <description>Enumeration of levels.</description>
   <indexing>1</indexing>
   <entry label="1"/>
   <entry label="B"/>
  </enumeration>
  !!]

contains

  function hummerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicRecombinationRateRadiativeCoolingHummer} atomic radiative recombination class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (atomicRecombinationRateRadiativeCoolingHummer)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (atomicRecombinationRateRadiativeClass        ), pointer       :: atomicRecombinationRateRadiative_
    double precision                                                               :: gamma

    !![
    <inputParameter>
      <name>gamma</name>
      <description>The multiplicative factor, $\gamma$, used to compute the cooling coefficient in cases of non-hydrogenic or helium-like ions.</description>
      <source>parameters</source>
      <defaultValue>0.67d0</defaultValue>
    </inputParameter>
    <objectBuilder class="atomicRecombinationRateRadiative" name="atomicRecombinationRateRadiative_" source="parameters"/>
    !!]
    self=atomicRecombinationRateRadiativeCoolingHummer(gamma,atomicRecombinationRateRadiative_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicRecombinationRateRadiative_"/>
    !!]
    return
  end function hummerConstructorParameters
  
  function hummerConstructorInternal(gamma,atomicRecombinationRateRadiative_) result(self)
    !!{
    Internal constructor for the \refClass{atomicRecombinationRateRadiativeCoolingHummer} atomic radiative recombination class.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Table_Labels    , only : extrapolationTypeExtrapolate
    implicit none
    type            (atomicRecombinationRateRadiativeCoolingHummer)                        :: self
    class           (atomicRecombinationRateRadiativeClass        ), intent(in   ), target :: atomicRecombinationRateRadiative_
    double precision                                               , intent(in   )         :: gamma
    !![
    <constructorAssign variables="gamma, *atomicRecombinationRateRadiative_"/>
    !!]

    ! Build the tables of cooling coefficients.
    call self%coefficientTable(sequenceHydrogen%ID)%create(1.0d1,1.00000d7,31,2,extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
    call self%coefficientTable(sequenceHelium  %ID)%create(1.0d1,2.51186d4,18,2,extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
    call self%coefficientTable(sequenceHydrogen%ID)%populate([1.646d-11,1.646d-11,1.646d-11,1.646d-11,1.646d-11,1.645d-11,1.644d-11,1.643d-11,1.641d-11,1.638d-11,1.633d-11,1.625d-11,1.613d-11,1.594d-11,1.565d-11,1.522d-11,1.460d-11,1.374d-11,1.260d-11,1.119d-11,9.571d-12,7.844d-12,6.146d-12,4.601d-12,3.295d-12,2.262d-12,1.494d-12,9.520d-13,5.878d-13,3.528d-13,2.066d-13],table=level1%ID)
    call self%coefficientTable(sequenceHydrogen%ID)%populate([8.287d-11,7.821d-11,7.356d-11,6.892d-11,6.430d-11,5.971d-11,5.515d-11,5.062d-11,4.614d-11,4.170d-11,3.734d-11,3.306d-11,2.888d-11,2.484d-11,2.098d-11,1.736d-11,1.402d-11,1.103d-11,8.442d-12,6.279d-12,4.539d-12,3.192d-12,2.185d-12,1.458d-12,9.484d-13,6.023d-13,3.738d-13,2.268d-13,1.348d-13,7.859d-14,4.499d-13],table=levelB%ID)
    call self%coefficientTable(sequenceHelium  %ID)%populate([1.569d-11,1.569d-11,1.569d-11,1.569d-11,1.569d-11,1.569d-11,1.570d-11,1.570d-11,1.571d-11,1.572d-11,1.574d-11,1.578d-11,1.583d-11,1.591d-11,1.602d-11,1.619d-11,1.641d-11,1.670d-11],table=level1%ID)
    call self%coefficientTable(sequenceHelium  %ID)%populate([8.347d-11,7.889d-11,7.430d-11,6.971d-11,6.512d-11,6.056d-11,5.603d-11,5.154d-11,4.710d-11,4.274d-11,3.847d-11,3.431d-11,3.031d-11,2.650d-11,2.291d-11,1.960d-11,1.660d-11,1.394d-11],table=levelB%ID)
    return
  end function hummerConstructorInternal

  subroutine hummerDestructor(self)
    !!{
    Destructor for the \refClass{atomicRecombinationRateRadiativeCoolingHummer} recombination cooling class.
    !!}
    implicit none
    type(atomicRecombinationRateRadiativeCoolingHummer), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicRecombinationRateRadiative_"/>
    !!]
    return
  end subroutine hummerDestructor

  double precision function hummerRate(self,atomicNumber,ionizationState,temperature,level)
    !!{
    Returns the cooling rate coefficient.
    !!}
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseA, recombinationCaseB
    use :: Error                               , only : Error_Report
    implicit none
    class           (atomicRecombinationRateRadiativeCoolingHummer), intent(inout)           :: self
    integer                                                        , intent(in   )           :: atomicNumber        , ionizationState
    double precision                                               , intent(in   )           :: temperature
    type            (enumerationRecombinationCaseType             ), intent(in   ), optional :: level
    type            (enumerationSequenceType                      )                          :: sequence_
    double precision                                                                         :: temperatureEffective, scaleFactor
    !![
    <optionalArgument name="level" defaultsTo="recombinationCaseA" />
    !!]

    ! Handle non-hydrogenic, non-helium-like ions separately.
    if (atomicNumber-ionizationState > 1) then
       ! Assume a fixed multiple of the corresponding recombination rate coefficient.
       hummerRate=+self%gamma                                                                                  &
            &     *self%atomicRecombinationRateRadiative_%rate(atomicNumber,ionizationState,temperature,level)
    else
       ! Determine which sequence to use.
       if      (atomicNumber-ionizationState == 0) then
          ! Hydrogen-like.
          sequence_           = sequenceHydrogen
          temperatureEffective=+temperature             &
               &               /dble(atomicNumber  )**2
          scaleFactor         =+dble(atomicNumber  )
       else if (atomicNumber-ionizationState == 1) then
          ! Helium-like. Use the hydrogenic approximation (e.g. Bautista & Kallman, 2000, ApJ,544, 581) to estimate the
          ! recombination cooling coefficient for helium-like ions.
          sequence_           = sequenceHelium
          temperatureEffective=+temperature             &
               &               /dble(atomicNumber-1)**2
          scaleFactor         =+dble(atomicNumber-1)
       else
          sequence_           = sequenceNone
          temperatureEffective=+temperature
          scaleFactor         =+1.0d0
          call Error_Report('expected hydrogen- or helium-like ion'//{introspection:location})
       end if
       ! Handle zero temperature.
       if (temperatureEffective <= 0.0d0) then
          hummerRate=0.0d0
          return
       end if
       ! Use tabulated solutions for hydrogenic and helium-like sequences.
       select case (level%ID)
       case (recombinationCaseA%ID)
          hummerRate=+self%coefficientTable(sequence_%ID)%interpolate(temperatureEffective,table=level1%ID) &
          &          +self%coefficientTable(sequence_%ID)%interpolate(temperatureEffective,table=levelB%ID)
       case (recombinationCaseB%ID)
          hummerRate=+self%coefficientTable(sequence_%ID)%interpolate(temperatureEffective,table=levelB%ID)
       case (1                 )
          hummerRate=+self%coefficientTable(sequence_%ID)%interpolate(temperatureEffective,table=level1%ID)
       case default
          hummerRate=0.0d0
          call Error_Report('only level 1, case A, and case B recombination is supported'//{introspection:location})
       end select
       ! Scale and divide out the thermal energy.
       hummerRate=+hummerRate                 &
            &     *scaleFactor                &
            &     /sqrt(temperatureEffective)
    end if
    return
  end function hummerRate
