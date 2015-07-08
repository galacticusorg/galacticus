!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% A primordial power spectrum class which provides a power-law power spectrum.
  
  !# <powerSpectrumPrimordial name="powerSpectrumPrimordialPowerLaw">
  !#  <description>
  !#   Implements a power-law primordial power spectrum, possibly with a running index. It is defined by
  !#   \begin{equation}
  !#    P(k)\propto k^{n_{\mathrm s} + \ln(k/k_{\mathrm ref}) [\d n /\d\ln k]},
  !#   \end{equation}
  !#   where the parameters are specified by input parameters $n_{\mathrm s}\equiv${\normalfont \ttfamily [index]}, $k_{\mathrm ref}\equiv${\normalfont \ttfamily [wavenumberReference]} and $\d n / \d \ln k \equiv${\normalfont \ttfamily [running]}.
  !#  </description>
  !# </powerSpectrumPrimordial>
  type, extends(powerSpectrumPrimordialClass) :: powerSpectrumPrimordialPowerLaw
     !% A power-law primordial power spectrum class.
     private
     double precision :: index              , running, &
          &              wavenumberReference
   contains
     final     ::                          powerLawDestructor
     procedure :: power                 => powerLawPower
     procedure :: logarithmicDerivative => powerLawLogarithmicDerivative
     procedure :: descriptor            => powerLawDescriptor
  end type powerSpectrumPrimordialPowerLaw

  interface powerSpectrumPrimordialPowerLaw
     !% Constructors for the ``power-law'' primordial power spectrum class.
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface powerSpectrumPrimordialPowerLaw

contains

  function powerLawConstructorParameters(parameters)
    !% Constructor for the ``power-law'' primordial power spectrum class which takes a parameter set as input.
    use Galacticus_Display
    implicit none
    type(powerSpectrumPrimordialPowerLaw)                :: powerLawConstructorParameters
    type(inputParameters                ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>index</name>
    !#   <source>parameters</source>
    !#   <variable>powerLawConstructorParameters%index</variable>
    !#   <defaultValue>0.9646d0</defaultValue>
    !#   <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
    !#   <description>The index of the power-law primordial power spectrum.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>running</name>
    !#   <source>parameters</source>
    !#   <variable>powerLawConstructorParameters%running</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The running, $\d n_{\mathrm s} / \d \ln k$, of the power spectrum index.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>wavenumberReference</name>
    !#   <source>parameters</source>
    !#   <variable>powerLawConstructorParameters%wavenumberReference</variable>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>When a running power spectrum index is used, this is the wavenumber at which the index is equal to {\normalfont \ttfamily [index]}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(index,running,wavenumberReference)
    !% Internal constructor for the ``power-law'' primordial power spectrum class.
    implicit none
    type            (powerSpectrumPrimordialPowerLaw)                :: powerLawConstructorInternal
    double precision                                 , intent(in   ) :: index                      , running, &
         &                                                              wavenumberReference

    powerLawConstructorInternal%index              =index
    powerLawConstructorInternal%running            =running
    powerLawConstructorInternal%wavenumberReference=wavenumberReference
    return
  end function powerLawConstructorInternal

  elemental subroutine powerLawDestructor(self)
    !% Destructor for the ``power-law'' primordial power spectrum class.
    implicit none
    type(powerSpectrumPrimordialPowerLaw), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine powerLawDestructor

  double precision function powerLawPower(self,wavenumber)
    !% Return the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    implicit none
    class           (powerSpectrumPrimordialPowerLaw), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavenumber

    powerLawPower=+(                                  &
         &          +             wavenumber          &
         &          /        self%wavenumberReference &
         &         )**(                               &
         &             +self%index                    &
         &             +0.5d0                         &
         &             *self%running                  &
         &             *log(                          &
         &                  +     wavenumber          &
         &                  /self%wavenumberReference &
         &                 )                          &
         &            )
    return
  end function powerLawPower

  double precision function powerLawLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    implicit none
    class           (powerSpectrumPrimordialPowerLaw), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavenumber

    powerLawLogarithmicDerivative=+self%index                    &
         &                        +self%running                  &
         &                        *log(                          &
         &                             +     wavenumber          &
         &                             /self%wavenumberReference &
         &                            )
    return
  end function powerLawLogarithmicDerivative

  subroutine powerLawDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (powerSpectrumPrimordialPowerLaw), intent(inout) :: self
    type     (inputParameters                ), intent(inout) :: descriptor
    type     (node                           ), pointer       :: parameterNode
    type     (inputParameters                )                :: subParameters
    character(len=10                         )                :: parameterLabel

    call descriptor%addParameter("powerSpectrumPrimordialMethod","powerLaw")
    parameterNode => descriptor%node("powerSpectrumPrimordialMethod")
    subParameters=inputParameters(parameterNode)
    write (parameterLabel,'(f10.6)') self%index
    call subParameters%addParameter("index"              ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%running
    call subParameters%addParameter("running"            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%wavenumberReference
    call subParameters%addParameter("wavenumberReference",trim(adjustl(parameterLabel)))
    return
  end subroutine powerLawDescriptor
