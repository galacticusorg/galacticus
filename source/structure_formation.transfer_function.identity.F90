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
Implements an identity transfer function class.
!!}

  !![
  <transferFunction name="transferFunctionIdentity">
   <description>Provides an identity transfer function, i.e. $T(k)=1$ for all $k$.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionIdentity
     !!{
     A identity transfer function class.
     !!}
     private
     double precision :: time
   contains
     final     ::                          identityDestructor
     procedure :: value                 => identityValue
     procedure :: logarithmicDerivative => identityLogarithmicDerivative
     procedure :: halfModeMass          => identityHalfModeMass
     procedure :: quarterModeMass       => identityQuarterModeMass
     procedure :: epochTime             => identityEpochTime
     procedure :: descriptor            => identityDescriptor
  end type transferFunctionIdentity

  interface transferFunctionIdentity
     !!{
     Constructors for the identity transfer function class.
     !!}
     module procedure identityConstructorParameters
     module procedure identityConstructorInternal
  end interface transferFunctionIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the identity transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    use :: Cosmology_Functions , only : cosmologyFunctionsClass
    use :: Input_Parameters    , only : inputParameter          , inputParameters
    implicit none
    type            (transferFunctionIdentity)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                                          :: redshift

    !![
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which the transfer function is defined.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=transferFunctionIdentity(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
  return
  end function identityConstructorParameters

  function identityConstructorInternal(time,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the identity transfer function class.
    !!}
    implicit none
    type            (transferFunctionIdentity)                        :: self
    double precision                          , intent(in   )         :: time
    class           (cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="time, *cosmologyParameters_"/>
    !!]

    return
  end function identityConstructorInternal

  subroutine identityDestructor(self)
    !!{
    Destructor for the identity transfer function class.
    !!}
    implicit none
    type(transferFunctionIdentity), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine identityDestructor

  double precision function identityValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionIdentity), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    !$GLC attributes unused :: self, wavenumber

    identityValue=1.0d0
    return
  end function identityValue

  double precision function identityLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionIdentity), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    !$GLC attributes unused :: self, wavenumber

    identityLogarithmicDerivative=0.0d0
    return
  end function identityLogarithmicDerivative

  double precision function identityHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionIdentity), intent(inout), target   :: self
    integer                          , intent(  out), optional :: status
    !$GLC attributes unused :: self

    identityHalfModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function identityHalfModeMass

  double precision function identityQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionIdentity), intent(inout), target   :: self
    integer                          , intent(  out), optional :: status
    !$GLC attributes unused :: self

    identityQuarterModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function identityQuarterModeMass

  double precision function identityEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionIdentity), intent(inout) :: self

    identityEpochTime=self%time
    return
  end function identityEpochTime

  subroutine identityDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (transferFunctionIdentity), intent(inout)           :: self
    type   (inputParameters         ), intent(inout)           :: descriptor
    logical                          , intent(in   ), optional :: includeClass, includeFileModificationTimes
    !$GLC attributes unused :: self, includeFileModificationTimes

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('transferFunction','identity')
    return
  end subroutine identityDescriptor
