!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Contains a module which implements a transfer function class based on the fuzzy dark matter modifier of \cite{hu_fuzzy_2000}.

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !# <transferFunction name="transferFunctionHu2008Fuzzy">
  !#  <description>Provides a transfer function based on the fuzzy dark matter modifier of \cite{hu_fuzzy_2000}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionHu2008Fuzzy
     !% A transfer function class which modifies another transfer function using the fuzzy dark matter modifier of
     !% \cite{hu_fuzzy_2000}.
     private
     double precision                                    :: m22                           , time, &
          &                                                 redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
   contains
     final     ::                          hu2008FuzzyDestructor
     procedure :: value                 => hu2008FuzzyValue
     procedure :: logarithmicDerivative => hu2008FuzzyLogarithmicDerivative
     procedure :: halfModeMass          => hu2008FuzzyHalfModeMass
     procedure :: epochTime             => hu2008FuzzyEpochTime
  end type transferFunctionHu2008Fuzzy

  interface transferFunctionHu2008Fuzzy
     !% Constructors for the {\normalfont \ttfamily bode2001} transfer function class.
     module procedure hu2008FuzzyConstructorParameters
     module procedure hu2008FuzzyConstructorInternal
  end interface transferFunctionHu2008Fuzzy

contains

  function hu2008FuzzyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily bode2001} transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionHu2008Fuzzy)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (transferFunctionClass      ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    double precision                                             :: m22                 , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>m22</name>
    !#   <source>parameters</source>
    !#   <defaultValue>40.0d0</defaultValue>
    !#   <description>Dark matter particle mass in units of $10^{-22}$~eV.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
    !#   <description>The redshift of the epoch at which the transfer function is defined.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
     self=transferFunctionHu2008Fuzzy(transferFunctionCDM,m22,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function hu2008FuzzyConstructorParameters
  
  function hu2008FuzzyConstructorInternal(transferFunctionCDM,m22,time,cosmologyParameters_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily bode2001} transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionHu2008Fuzzy)                        :: self
    class           (transferFunctionClass      ), target, intent(in   ) :: transferFunctionCDM
    double precision                                     , intent(in   ) :: m22                 , time
    class           (cosmologyParametersClass   ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass    ), target, intent(in   ) :: cosmologyFunctions_
    !# <constructorAssign variables="*transferFunctionCDM, m22, time, *cosmologyParameters_, *cosmologyFunctions_"/>

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function hu2008FuzzyConstructorInternal

  subroutine hu2008FuzzyDestructor(self)
    !% Destructor for the {\normalfont \ttfamily bode2001} transfer function class.
    implicit none
    type(transferFunctionHu2008Fuzzy), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine hu2008FuzzyDestructor

  double precision function hu2008FuzzyValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionHu2008Fuzzy), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber
    double precision                                             :: x

    hu2008FuzzyValue=+self%transferFunctionCDM%value(wavenumber)
    x               =+1.61d0                     &
         &           /9.00d0                     &
         &           *self%m22  **(-4.0d0/9.0d0) &
         &           *wavenumber
    hu2008FuzzyValue =+hu2008FuzzyValue &
         &            *(                &
         &              +cos(x**3)      &
         &              /(              &
         &                +1.0d0        &
         &                +x**8         &
         &               )              &
         &             )   
    return
  end function hu2008FuzzyValue

  double precision function hu2008FuzzyLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionHu2008Fuzzy), intent(inout) :: self
    double precision                             , intent(in   ) :: wavenumber
    double precision                                             :: x

    hu2008FuzzyLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    x                               =+1.61d0                           &
         &                           /9.00d0                           &
         &                           *self%m22**(-4.0d0/9.0d0)         &
         &                           *wavenumber
    hu2008FuzzyLogarithmicDerivative=+hu2008FuzzyLogarithmicDerivative &
         &                           +(                                &
         &                             -8.0d0                          &
         &                             *    x** 7                      &
         &                             -3.0d0                          &
         &                             *(                              &
         &                               +  x**10                      &
         &                               +  x** 2                      &
         &                              )                              &
         &                             *tan(x** 3)                     &
         &                            )                                &
         &                           /(                                &
         &                             +1.0d0                          &
         &                             +    x** 8                      &
         &                            )
    return
  end function hu2008FuzzyLogarithmicDerivative

  double precision function hu2008FuzzyHalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionHu2008Fuzzy), intent(inout)           :: self
    integer                                      , intent(  out), optional :: status
    double precision                                                       :: matterDensity, wavenumberHalfMode

    matterDensity          =+self%cosmologyParameters_%OmegaMatter    () &
         &                  *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode     =+4.5d0                   &
         &                  *self%m22**(4.0d0/9.0d0)
    hu2008FuzzyHalfModeMass=+4.0d0                &
         &                  *Pi                   &
         &                  /3.0d0                &
         &                  *matterDensity        &
         &                  *(                    &
         &                    +Pi                 &
         &                    /wavenumberHalfMode &
         &                  )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function hu2008FuzzyHalfModeMass

  double precision function hu2008FuzzyEpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionHu2008Fuzzy), intent(inout) :: self

    hu2008FuzzyEpochTime=self%time
    return
  end function hu2008FuzzyEpochTime
