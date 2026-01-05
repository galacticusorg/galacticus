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
  An implementation of the intergalactic medium state class for a simplistic model of instantaneous and full reionization.
  !!}

  !![
  <intergalacticMediumState name="intergalacticMediumStateSimple">
   <description>
    An intergalactic medium state class which implements a simple model of reionization in which the universe is assumed to be
    fully neutral prior to the redshift given by {\normalfont \ttfamily [reionizationRedshift]} and fully ionized
    thereafter. The temperature is given by {\normalfont \ttfamily [preReionizationTemperature]} before reionization, and
    {\normalfont \ttfamily [reionizationTemperature]} thereafter.
   </description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateSimple
     !!{
     An \gls{igm} state class for a simple model in which the \gls{igm} is assumed to be instantaneously and fully reionized at
     a fixed redshift, and heated to a fixed temperature.
     !!}
     private
     double precision :: reionizationTime, reionizationTemperature, preReionizationTemperature
   contains
     final     ::                                simpleDestructor
     procedure :: electronFraction            => simpleElectronFraction
     procedure :: temperature                 => simpleTemperature
     procedure :: neutralHydrogenFraction     => simpleNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => simpleNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => simpleSinglyIonizedHeliumFraction
     procedure :: descriptor                  => simpleDescriptor
  end type intergalacticMediumStateSimple

  interface intergalacticMediumStateSimple
     !!{
     Constructors for the simple intergalactic medium state class.
     !!}
     module procedure simpleIGMConstructorParameters
     module procedure simpleIGMConstructorInternal
  end interface intergalacticMediumStateSimple

contains

  function simpleIGMConstructorParameters(parameters) result (self)
    !!{
    Constructor for the simple \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (intergalacticMediumStateSimple)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    double precision                                                :: reionizationRedshift      , reionizationTemperature, &
         &                                                             preReionizationTemperature

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>reionizationRedshift</name>
      <source>parameters</source>
      <variable>reionizationRedshift</variable>
      <defaultValue>9.97d0</defaultValue>
      <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
      <description>The redshift of reionization in the simple \gls{igm} state model.</description>
    </inputParameter>
    <inputParameter>
      <name>reionizationTemperature</name>
      <source>parameters</source>
      <variable>reionizationTemperature</variable>
      <defaultValue>1.0d4</defaultValue>
      <description>The post-reionization temperature (in units of Kelvin) in the simple \gls{igm} state model.</description>
    </inputParameter>
    <inputParameter>
      <name>preReionizationTemperature</name>
      <source>parameters</source>
      <variable>preReionizationTemperature</variable>
      <defaultValue>10.0d0</defaultValue>
      <description>The pre-reionization temperature (in units of Kelvin) in the simple \gls{igm} state model.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    ! Construct the object.
    self=intergalacticMediumStateSimple(reionizationRedshift,reionizationTemperature,preReionizationTemperature,cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function simpleIGMConstructorParameters

  function simpleIGMConstructorInternal(reionizationRedshift,reionizationTemperature,preReionizationTemperature,cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Constructor for the simple \gls{igm} state class.
    !!}
    implicit none
    type            (intergalacticMediumStateSimple)                        :: self
    double precision                                , intent(in   )         :: reionizationRedshift      , reionizationTemperature, &
         &                                                                     preReionizationTemperature
    class           (cosmologyFunctionsClass       ), intent(inout), target :: cosmologyFunctions_
    class           (cosmologyParametersClass      ), intent(inout), target :: cosmologyParameters_
    !![
    <constructorAssign variables="reionizationTemperature, preReionizationTemperature, *cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    self%reionizationTime=cosmologyFunctions_%cosmicTime                 (                      &
         &                cosmologyFunctions_%expansionFactorFromRedshift (                     &
         &                                                                 reionizationRedshift &
         &                                                                )                     &
         &                                                               )
    return
  end function simpleIGMConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the simple \gls{igm} state class.
    !!}
    implicit none
    type(intergalacticMediumStateSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine simpleDestructor

  double precision function simpleElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm} in the simple model.
    !!}
    use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial, hydrogenByMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassHelium      , atomicMassHydrogen
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleElectronFraction=+(                                                    &
            &                   +      hydrogenByMassPrimordial/atomicMassHydrogen  &
            &                   +2.0d0*  heliumByMassPrimordial/atomicMassHelium    &
            &                  )                                                    &
            &                 /(       hydrogenByMassPrimordial/atomicMassHydrogen)
    else
       simpleElectronFraction=0.0d0
    end if
    return
  end function simpleElectronFraction

  double precision function simpleNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm} in the simple model.
    !!}
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleNeutralHydrogenFraction=0.0d0
    else
       simpleNeutralHydrogenFraction=1.0d0
    end if
    return
  end function simpleNeutralHydrogenFraction

  double precision function simpleNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm} in the simple model.
    !!}
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleNeutralHeliumFraction=0.0d0
    else
       simpleNeutralHeliumFraction=1.0d0
    end if
    return
  end function simpleNeutralHeliumFraction

  double precision function simpleSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly-ionized helium fraction of the \gls{igm} in the simple model.
    !!}
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleSinglyIonizedHeliumFraction=0.0d0
    else
       simpleSinglyIonizedHeliumFraction=1.0d0
    end if
    return
  end function simpleSinglyIonizedHeliumFraction

  double precision function simpleTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm} in the simple model.
    !!}
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleTemperature=self%   reionizationTemperature
    else
       simpleTemperature=self%preReionizationTemperature
    end if
    return
  end function simpleTemperature

  subroutine simpleDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (intergalacticMediumStateSimple), intent(inout)           :: self
    type     (inputParameters               ), intent(inout)           :: descriptor
    logical                                  , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                        )                          :: parameterLabel
    type     (inputParameters               )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('intergalacticMediumState','simple')
    parameters=descriptor%subparameters('intergalacticMediumState')
    write (parameterLabel,'(e17.10)') self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%reionizationTime          ))
    call parameters%addParameter('reionizationRedshift'      ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)')                                                                                               self%reionizationTemperature
    call parameters%addParameter('reionizationTemperature'   ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)')                                                                                               self%preReionizationTemperature
    call parameters%addParameter('preReionizationTemperature',trim(adjustl(parameterLabel)))
    call self%cosmologyFunctions_ %descriptor(parameters,includeClass,includeFileModificationTimes)
    call self%cosmologyParameters_%descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine simpleDescriptor
