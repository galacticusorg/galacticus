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
  An implementation of fixed dark matter halo virial density contrasts.
  !!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  ! Enumeration for different reference densities.
  !![
  <enumeration>
   <name>fixedDensityType</name>
   <description>Specifies reference density type for fixed virial density contrast.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="critical" />
   <entry label="mean"     />
  </enumeration>
  !!]

  !![
  <virialDensityContrast name="virialDensityContrastFixed">
   <description>
    A dark matter halo virial density contrast class which uses a fixed virial density contrast of {\normalfont \ttfamily
    [densityContrastValue]}, defined relative to {\normalfont \ttfamily critical} or {\normalfont \ttfamily mean} density as
    specified by {\normalfont \ttfamily [densityType]}.
   </description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastFixed
     !!{
     A dark matter halo virial density contrast class assuming fixed contrast.
     !!}
     private
     class           (cosmologyParametersClass       ), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_  => null()
     double precision                                           :: densityContrastValue          , turnAroundOverVirialRadius
     type            (enumerationFixedDensityTypeType)          :: densityType
   contains
     final     ::                                fixedDestructor
     procedure :: densityContrast             => fixedDensityContrast
     procedure :: densityContrastRateOfChange => fixedDensityContrastRateOfChange
     procedure :: turnAroundOverVirialRadii   => fixedTurnAroundOverVirialRadii
  end type virialDensityContrastFixed

  interface virialDensityContrastFixed
     !!{
     Constructors for the \refClass{virialDensityContrastFixed} dark matter halo virial density contrast class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface virialDensityContrastFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastFixed} dark matter halo virial density contrast class that takes a parameter set as input.
    !!}
    use :: ISO_Varying_String, only : var_str       , varying_string
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (virialDensityContrastFixed)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    double precision                                            :: densityContrastValue, turnAroundOverVirialRadius
    type            (varying_string            )                :: densityType

    !![
    <inputParameter>
     <name>densityContrastValue</name>
     <source>parameters</source>
     <defaultValue>200.0d0</defaultValue>
     <description>The virial density contrast to use in the fixed value model.</description>
    </inputParameter>
    <inputParameter>
     <name>densityType</name>
     <source>parameters</source>
     <defaultValue>var_str('critical')</defaultValue>
     <description>The reference density to use in the fixed value virial density contrast model. Either of {\normalfont \ttfamily critical} and {\normalfont \ttfamily mean} are allowed.</description>
    </inputParameter>
    <inputParameter>
     <name>turnAroundOverVirialRadius</name>
     <source>parameters</source>
     <defaultValue>2.0d0</defaultValue>
     <description>The ratio of the turnaround to virial radii in the fixed value model.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=virialDensityContrastFixed(densityContrastValue,enumerationFixedDensityTypeEncode(char(densityType),includesPrefix=.false.),turnAroundOverVirialRadius,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(densityContrastValue,densityType,turnAroundOverVirialRadius,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastFixed} dark matter halo virial density contrast class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (virialDensityContrastFixed     )                        :: self
    double precision                                 , intent(in   )         :: densityContrastValue, turnAroundOverVirialRadius
    type            (enumerationFixedDensityTypeType), intent(in   )         :: densityType
    class           (cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="densityContrastValue, densityType, turnAroundOverVirialRadius, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    if (.not.enumerationFixedDensityTypeIsValid(densityType)) call Error_Report('invalid densityType'//{introspection:location})
    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the \refClass{virialDensityContrastFixed} virial density contrast class.
    !!}
    implicit none
    type(virialDensityContrastFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine fixedDestructor

  double precision function fixedDensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming a fixed contrast.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsStaticUniverse
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   )           :: mass
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    !$GLC attributes unused :: mass

    ! Set the density contrast.
    fixedDensityContrast=self%densityContrastValue
    ! If density contrast is specified relative to critical density, convert to mean density.
    if (self%densityType == fixedDensityTypeCritical) then
       select type (cosmologyFunctions_ => self%cosmologyFunctions_)
       class is (cosmologyFunctionsStaticUniverse)
          ! For a static universe, ``omegaMatterEpochal" is not well defined. So here we use the
          ! value of ``OmegaMatter" specified in the parameter file when doing the conversion.
          fixedDensityContrast=+fixedDensityContrast                                           &
               &               /self%cosmologyParameters_%omegaMatter       ()
       class default
          fixedDensityContrast=+fixedDensityContrast                                           &
               &               /self%cosmologyFunctions_ %omegaMatterEpochal(                  &
               &                 self%cosmologyFunctions_%epochTime          (                 &
               &                                                              time           , &
               &                                                              expansionFactor, &
               &                                                              collapsing       &
               &                                                             )                 &
               &                                                            )
       end select
    end if
    return
  end function fixedDensityContrast

  double precision function fixedDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming a fixed contrast.
    !!}
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   )           :: mass
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    double precision                                                      :: epochTime
    !$GLC attributes unused :: mass

    ! Zero rate of change for fixed density contrast.
    fixedDensityContrastRateOfChange=0.0d0
    ! If density contrast is defined relative to critical density, include the rate of change of
    ! critical mean density.
    if (self%densityType == fixedDensityTypeCritical) then
       epochTime=self%cosmologyFunctions_%epochTime(                 &
            &                                       time           , &
            &                                       expansionFactor, &
            &                                       collapsing       &
            &                                      )
       fixedDensityContrastRateOfChange=-self%densityContrastValue                                       &
            &                           *self%cosmologyFunctions_ %omegaMatterRateOfChange(epochTime)    &
            &                           /self%cosmologyFunctions_ %omegaMatterEpochal     (epochTime)**2
    end if
    return
  end function fixedDensityContrastRateOfChange

  double precision function fixedTurnAroundOverVirialRadii(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming a fixed contrast.
    !!}
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   )           :: mass
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    !$GLC attributes unused :: mass, time, expansionFactor, collapsing

    fixedTurnAroundOverVirialRadii=self%turnAroundOverVirialRadius
    return
  end function fixedTurnAroundOverVirialRadii
