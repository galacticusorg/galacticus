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
  Implements a transfer function accelerator class which tabulates a transfer function for rapid interpolation.
  !!}

  use :: Tables, only : table1DLinearLinear

  !![
  <transferFunction name="transferFunctionAccelerator">
   <description>A transfer function class which accelerates calculations of another transfer function class by tabulation for rapid interpolation.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionAccelerator
     !!{
     A transfer function class which accelerates calculations of another transfer function class by tabulation for rapid interpolation.
     !!}
     private
     type            (table1DLinearLinear  )          :: transferTable
     class           (transferFunctionClass), pointer :: transferFunction_            => null()
     double precision                                 :: wavenumberLogarithmicMinimum           , wavenumberLogarithmicMaximum
     integer                                          :: tablePointsPerDecade
     logical                                          :: tableInitialized             =  .false.
   contains
     final     ::                          acceleratorDestructor
     procedure :: value                 => acceleratorValue
     procedure :: logarithmicDerivative => acceleratorLogarithmicDerivative
     procedure :: halfModeMass          => acceleratorHalfModeMass
     procedure :: quarterModeMass       => acceleratorQuarterModeMass
     procedure :: fractionModeMass      => acceleratorFractionModeMass
     procedure :: epochTime             => acceleratorEpochTime
  end type transferFunctionAccelerator

  interface transferFunctionAccelerator
     !!{
     Constructors for the accelerator transfer function class.
     !!}
     module procedure acceleratorConstructorParameters
     module procedure acceleratorConstructorInternal
  end interface transferFunctionAccelerator

contains

  function acceleratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the accelerator transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    type  (transferFunctionAccelerator)                :: self
    type  (inputParameters            ), intent(inout) :: parameters
    class (transferFunctionClass      ), pointer       :: transferFunction_
    class (cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    integer                                            :: tablePointsPerDecade

    !![
    <inputParameter>
      <name>tablePointsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of points per decade of wavenumber at which to tabulate the transfer function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunction_"    source="parameters"/>
    !!]
    self=transferFunctionAccelerator(transferFunction_,cosmologyParameters_,tablePointsPerDecade)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="transferFunction_"/>
    !!]
    return
  end function acceleratorConstructorParameters

  function acceleratorConstructorInternal(transferFunction_,cosmologyParameters_,tablePointsPerDecade) result(self)
    !!{
    Internal constructor for the accelerator transfer function class.
    !!}
    implicit none
    type   (transferFunctionAccelerator)                        :: self
    class  (cosmologyParametersClass   ), intent(in   ), target :: cosmologyParameters_
    class  (transferFunctionClass      ), intent(in   ), target :: transferFunction_
    integer                             , intent(in   )         :: tablePointsPerDecade
    !![
    <constructorAssign variables="*cosmologyParameters_, *transferFunction_, tablePointsPerDecade"/>
    !!]

    self%tableInitialized            =.false.
    self%wavenumberLogarithmicMinimum=log(1.0d-6)
    self%wavenumberLogarithmicMaximum=log(1.0d+6)
    return
  end function acceleratorConstructorInternal

  subroutine acceleratorDestructor(self)
    !!{
    Destructor for the accelerator transfer function class.
    !!}
    implicit none
    type(transferFunctionAccelerator), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%transferFunction_"   />
    !!]
    if (self%tableInitialized) call self%transferTable%destroy()
    return
  end subroutine acceleratorDestructor

  double precision function acceleratorValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber
    double precision                                             :: wavenumberLogarithmic

    wavenumberLogarithmic=log(wavenumber)
    call acceleratorTabulate(self,wavenumberLogarithmic)
    acceleratorValue=exp(self%transferTable%interpolate(wavenumberLogarithmic))
    return
  end function acceleratorValue

  double precision function acceleratorLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber
    double precision                                             :: wavenumberLogarithmic

    wavenumberLogarithmic=log(wavenumber)
    call acceleratorTabulate(self,wavenumberLogarithmic)
    acceleratorLogarithmicDerivative=+self%transferTable%interpolateGradient(wavenumberLogarithmic)
    return
  end function acceleratorLogarithmicDerivative

  double precision function acceleratorEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionAccelerator), intent(inout) :: self

    acceleratorEpochTime=self%transferFunction_%epochTime()
    return
  end function acceleratorEpochTime

  double precision function acceleratorHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of two relative to a \gls{cdm} transfer function
    !!}
    implicit none
    class  (transferFunctionAccelerator), intent(inout), target   :: self
    integer                             , intent(  out), optional :: status

    acceleratorHalfModeMass=self%transferFunction_%halfModeMass(status)
    return
  end function acceleratorHalfModeMass

  double precision function acceleratorQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of four relative to a \gls{cdm} transfer function
    !!}
    implicit none
    class  (transferFunctionAccelerator), intent(inout), target   :: self
    integer                             , intent(  out), optional :: status

    acceleratorQuarterModeMass=self%transferFunction_%quarterModeMass(status)
    return
  end function acceleratorQuarterModeMass

  double precision function acceleratorFractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of four relative to a \gls{cdm} transfer function
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout), target   :: self
    double precision                             , intent(in   )           :: fraction
    integer                                      , intent(  out), optional :: status

    acceleratorFractionModeMass=self%transferFunction_%fractionModeMass(fraction,status)
    return
  end function acceleratorFractionModeMass

  subroutine acceleratorTabulate(self,wavenumberLogarithmic)
    !!{
    Tabulate the transfer function for rapid interpolation.
    !!}
    implicit none
    class           (transferFunctionAccelerator), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumberLogarithmic
    logical                                                      :: makeTable
    integer                                                      :: pointCount           , i

    makeTable=.not.self%tableInitialized
    if (.not.makeTable)                                                         &
         & makeTable= wavenumberLogarithmic < self%wavenumberLogarithmicMinimum &
         &           .or.                                                       &
         &            wavenumberLogarithmic > self%wavenumberLogarithmicMaximum
    if (makeTable) then
       self%wavenumberLogarithmicMinimum=min(self%wavenumberLogarithmicMinimum,wavenumberLogarithmic-1.0d0)
       self%wavenumberLogarithmicMaximum=max(self%wavenumberLogarithmicMaximum,wavenumberLogarithmic+1.0d0)
       pointCount=int((self%wavenumberLogarithmicMaximum-self%wavenumberLogarithmicMinimum)*dble(self%tablePointsPerDecade)/log(10.0d0))+1
       call self%transferTable%create(self%wavenumberLogarithmicMinimum,self%wavenumberLogarithmicMaximum,pointCount)
       do i=1,pointCount
          call self%transferTable%populate(log(self%transferFunction_%value(exp(self%transferTable%x(i)))),i)
       end do
       self%tableInitialized=.true.
    end if
    return
  end subroutine acceleratorTabulate

