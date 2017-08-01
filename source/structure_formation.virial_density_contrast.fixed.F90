!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of fixed dark matter halo virial density contrasts.

  use Cosmology_Functions

  !# <virialDensityContrast name="virialDensityContrastFixed">
  !#  <description>Fixed dark matter halo virial density contrasts.</description>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastClass) :: virialDensityContrastFixed
     !% A dark matter halo virial density contrast class assuming fixed contrast.
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     double precision                                    :: densityContrastValue
     integer                                             :: densityType
   contains
     final     ::                                fixedDestructor
     procedure :: densityContrast             => fixedDensityContrast
     procedure :: densityContrastRateOfChange => fixedDensityContrastRateOfChange
  end type virialDensityContrastFixed

  interface virialDensityContrastFixed
     !% Constructors for the {\normalfont \ttfamily fixed} dark matter halo virial density contrast class.
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface virialDensityContrastFixed

  ! Enumeration for different reference densities.
  !# <enumeration>
  !#  <name>fixedDensityType</name>
  !#  <description>Specifies reference density type for fixed virial density contrast.</description>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="critical" />
  !#  <entry label="mean"     />
  !# </enumeration>
  
contains

  function fixedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fixed} dark matter halo virial density contrast class that takes a parameter set as input.
    use Input_Parameters2
    use ISO_Varying_String
    implicit none
    type            (virialDensityContrastFixed)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    double precision                                            :: densityContrastValue
    type            (varying_string            )                :: densityType
    
    !# <inputParameter>
    !#  <name>densityContrastValue</name>
    !#  <source>parameters</source>
    !#  <defaultValue>200.0d0</defaultValue>
    !#  <description>The virial density contrast to use in the fixed value model.</description>
    !#  <type>float</type>
    !#  <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#  <name>densityType</name>
    !#  <source>parameters</source>
    !#  <defaultValue>var_str('critical')</defaultValue>
    !#  <description>The reference density to use in the fixed value virial density contrast model. Either of {\normalfont \ttfamily critical density} and {\normalfont \ttfamily mean density} are allowed.</description>
    !#  <type>string</type>
    !#  <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=virialDensityContrastFixed(densityContrastValue,enumerationFixedDensityTypeEncode(char(densityType),includesPrefix=.false.),cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(densityContrastValue,densityType,cosmologyFunctions_) result(self)
    !% Constructor for the {\normalfont \ttfamily fixed} dark matter halo virial density contrast class.
    use Galacticus_Error
    implicit none
    type            (virialDensityContrastFixed)                        :: self
    double precision                            , intent(in   )         :: densityContrastValue
    integer                                     , intent(in   )         :: densityType
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="densityContrastValue, densityType, *cosmologyFunctions_"/>

    if (.not.enumerationFixedDensityTypeIsValid(densityType)) call Galacticus_Error_Report('fixedConstructor','invalid densityType')
    return
  end function fixedConstructorInternal
  
  subroutine fixedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily fixed} virial density contrast class.
    implicit none
    type(virialDensityContrastFixed), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine fixedDestructor
  
  double precision function fixedDensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming a fixed contrast.
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   )           :: mass
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    !GCC$ attributes unused :: mass

    ! Set the density contrast.
    fixedDensityContrast=self%densityContrastValue
    ! If density contrast is specified relative to critical density, convert to mean density.
    if (self%densityType == fixedDensityTypeCritical)                                           &
         & fixedDensityContrast=+fixedDensityContrast                                           &
            &                   /self%cosmologyFunctions_ %omegaMatterEpochal(                  &
            &                     self%cosmologyFunctions_%epochTime          (                 &
            &                                                                  time           , &
            &                                                                  expansionFactor, &
            &                                                                  collapsing       &
            &                                                                 )                 &
            &                                                                )
    return
  end function fixedDensityContrast

  double precision function fixedDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming a fixed contrast.
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   )           :: mass
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    double precision                                                      :: epochTime
    !GCC$ attributes unused :: mass

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
